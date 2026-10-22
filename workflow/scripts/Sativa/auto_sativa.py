#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SATIVA automation for Snakemake.
Detect taxonomic misidentifications on MSAs.
"""

import argparse
import json
import subprocess
import sys
import time
import logging
import shutil
import re
from pathlib import Path
from datetime import datetime
from typing import Optional, Tuple
from dataclasses import dataclass, asdict


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)
logger = logging.getLogger(__name__)


_SANITIZE_PATTERN = re.compile(r"[^A-Za-z0-9_.-]")
_TAXON_PLACEHOLDER = "NA"


def _sanitize_label(value: str) -> str:
    """Convert arbitrary taxon/barcode labels into a filesystem-safe token."""
    return _SANITIZE_PATTERN.sub("_", value)


def _clean_sativa_tmp(sativa_script: Path, output_name: Optional[str]) -> None:
    """Remove the scratch directory SATIVA leaves behind when it does not exit cleanly.

    sativa.py creates <sativa_home>/tmp/<name>_<timestamp>/ and deletes it itself when
    the analysis finishes or aborts on a user error, but not when it hits a fatal error,
    raises, or is killed (timeout, OOM). Those runs are the ones that fill the disk: a
    single leftover jplace reaches 900 MB.

    Matching is scoped to this task's own directory: the name it passed to sativa.py
    followed by sativa's millisecond timestamp. Anchoring on those digits matters, since
    a bare "<name>_*" prefix would also catch a sibling whose taxon merely extends this
    one ("Abies" would match "Abies alba") and delete the scratch of a task still running.
    """
    if not output_name:
        return
    tmp_root = sativa_script.resolve().parent / "tmp"
    if not tmp_root.is_dir():
        return
    own = re.compile(rf"^{re.escape(_sanitize_label(output_name))}_\d+$")
    for leftover in tmp_root.glob(f"{_sanitize_label(output_name)}_*"):
        if leftover.is_dir() and own.match(leftover.name):
            shutil.rmtree(leftover, ignore_errors=True)


def _normalize_tax_value(value: str) -> str:
    """Ensure taxonomy ranks are never empty to keep SATIVA Newick parsing happy.

    Spaces in taxon names were escaped as underscores in the MSA FASTA header
    (so MAFFT/trimAl do not truncate them). Restore them here: the .tax file is
    tab-separated and is not re-parsed as FASTA, so the full name (e.g.
    "Russula azurea") is safe and needed for CoL validation / integration.
    """
    value = (value or "").strip().replace("_", " ")
    return value if value else _TAXON_PLACEHOLDER


@dataclass
class SativaResult:
    """SATIVA execution result metadata."""
    barcode: str
    tax_level: str
    taxon: str
    n_seq: int
    aln_len: int
    status: str  # 'success', 'skipped', 'error'
    reason: Optional[str] = None
    runtime_sec: Optional[float] = None
    tool_version: str = "SATIVA"
    command: Optional[str] = None
    timestamp: str = ""

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()


def is_output_complete(output_dir: Path, taxon: str) -> bool:
    """
    Determine whether a previous SATIVA run finished successfully.

    Args:
        output_dir: SATIVA output directory.
        taxon: Taxon name (base filename to check).

    Returns:
        True if the run is considered complete, otherwise False.
    """
    if not output_dir.exists():
        return False

    # Check the summary file
    summary_path = output_dir / "sativa.summary.json"
    if not summary_path.exists():
        return False

    # Ensure the summary is valid JSON
    try:
        with summary_path.open('r') as f:
            result = json.load(f)

        status = result.get('status')

        # If skipped or error, we consider it complete (no rerun)
        if status in ['skipped', 'error']:
            return True

        # If success, ensure the .mis file exists.
        # SATIVA names it after the -n label, i.e. `{barcode}_{tax_level}_{taxon}.mis`
        # (see the run_sativa call), NOT `{taxon}.mis`. Looking for the latter never
        # matched, so this returned False for every successful task and the caller then
        # ran clean_incomplete_output() -> shutil.rmtree(output_dir), DELETING a finished
        # result before redoing it. Glob instead: the label is sanitised, and a clade
        # directory holds exactly one .mis.
        if status == 'success':
            mis_files = [f for f in output_dir.glob("*.mis") if f.name != "consensus.mis"]
            if mis_files:
                return True
            logger.warning(f"Summary indicates success but no .mis file in {output_dir}")
            return False

        return False

    except (json.JSONDecodeError, Exception) as e:
        logger.warning(f"Invalid or corrupted summary file: {e}")
        return False


def clean_incomplete_output(output_dir: Path):
    """
    Remove an incomplete or corrupted SATIVA output directory.

    Args:
        output_dir: Directory to delete.
    """
    if output_dir.exists():
        logger.warning(f"Cleaning incomplete/corrupted output directory: {output_dir}")
        try:
            shutil.rmtree(output_dir)
            logger.info(f"Successfully cleaned: {output_dir}")
        except Exception as e:
            logger.error(f"Failed to clean directory: {e}")
            raise


def quick_msa_stats(msa_path: Path) -> Tuple[int, int]:
    """
    Compute quick stats (sequence count and alignment length) from an MSA.

    Args:
        msa_path: FASTA alignment path.

    Returns:
        Tuple (num_sequences, alignment_length)
    """
    n_seq = 0
    max_len = 0
    current_seq_len = 0

    with msa_path.open('r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq_len > max_len:
                    max_len = current_seq_len
                current_seq_len = 0
                n_seq += 1
            else:
                current_seq_len += len(line)

        # Don't forget the last sequence
        if current_seq_len > max_len:
            max_len = current_seq_len

    return n_seq, max_len


def parse_msa_to_sativa_format(msa_path: Path, output_dir: Path, taxon: str) -> Tuple[Path, Path]:
    """
    Convert an MSA to SATIVA's .fasta and .tax inputs.

    Expected MSA header format:
    >seq_id|kingdom=X|phylum=Y|class=Z|order=A|family=B|genus=C|species=D
    SEQUENCE

    Args:
        msa_path: Input MSA file.
        output_dir: Destination directory.
        taxon: Taxon name used for generated filenames.

    Returns:
        Tuple (fasta_path, tax_path).
    """
    fasta_output = output_dir / f"{taxon}.fasta"
    tax_output = output_dir / f"{taxon}.tax"

    logger.info(f"Parsing MSA: {msa_path.name}")

    with msa_path.open('r', encoding='utf-8') as in_file, \
         fasta_output.open('w', encoding='utf-8') as fasta_out, \
         tax_output.open('w', encoding='utf-8') as tax_out:

        current_id = None
        current_sequence_parts = []

        for line in in_file:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Write previous sequence
                if current_id:
                    fasta_out.write(f">{current_id}\n")
                    fasta_out.write("".join(current_sequence_parts) + "\n")

                current_sequence_parts = []
                parts = line[1:].split('|')
                current_id = parts[0]

                # Parse taxonomic info
                tax_info = {}
                for item in parts[1:]:
                    if '=' in item:
                        key, value = item.split('=', 1)
                        tax_info[key] = value

                # Build taxonomy string in order
                tax_order = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                tax_string = ';'.join(_normalize_tax_value(tax_info.get(rank, "")) for rank in tax_order)

                tax_out.write(f"{current_id}\t{tax_string}\n")
            else:
                current_sequence_parts.append(line)

        # Write last sequence
        if current_id:
            fasta_out.write(f">{current_id}\n")
            fasta_out.write("".join(current_sequence_parts) + "\n")

    logger.info(f"Generated: {fasta_output.name}, {tax_output.name}")
    return fasta_output, tax_output


def run_sativa(
    fasta_path: Path,
    tax_path: Path,
    output_dir: Path,
    sativa_script: Path,
    mode: str = "ultrafast",
    confidence: float = 0.6,
    taxonomy_code: str = "BOT",
    threads: int = 4,
    timeout_seconds: Optional[int] = None,
    output_name: Optional[str] = None,
    rand_seed: Optional[int] = None
) -> bool:
    """
    Execute SATIVA with the provided inputs.

    Args:
        fasta_path: FASTA sequences.
        tax_path: Taxonomy file.
        output_dir: Working/output directory for SATIVA.
        sativa_script: Path to sativa.py.
        mode: SATIVA execution mode.
        confidence: Confidence threshold.
        taxonomy_code: Taxonomy code passed to SATIVA.
        threads: Number of threads for SATIVA.
        timeout_seconds: Optional timeout (seconds).
        output_name: Optional sanitized output stem.

    Returns:
        True on success, False otherwise.

    Raises:
        subprocess.TimeoutExpired: When the timeout is reached.
    """
    command = [
        "python3",
        str(sativa_script.resolve()),
        "-s", fasta_path.name,
        "-t", tax_path.name,
        "-m", mode,
        "-C", str(confidence),
        "-x", taxonomy_code,
        "-T", str(threads),
        "-S",  # Silent mode
    ]
    if output_name:
        command.extend(["-n", _sanitize_label(output_name)])
    if rand_seed is not None:
        command.extend(["-p", str(rand_seed)])

    logger.info(f"Running SATIVA: {' '.join(command)}")
    if timeout_seconds:
        logger.info(f"Timeout set to: {timeout_seconds}s ({timeout_seconds/3600:.1f}h)")

    try:
        process = subprocess.run(
            command,
            cwd=output_dir.resolve(),
            capture_output=True,
            text=True,
            check=True,
            encoding='utf-8',
            timeout=timeout_seconds
        )

        logger.info("SATIVA completed successfully")
        return True

    except subprocess.TimeoutExpired as e:
        logger.warning(f"SATIVA timeout after {timeout_seconds}s ({timeout_seconds/3600:.1f}h)")
        logger.warning(f"Process was killed due to timeout")
        raise

    except subprocess.CalledProcessError as e:
        logger.error(f"SATIVA failed with return code {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")

        # Save error log
        error_log = output_dir / f"{fasta_path.stem}_sativa_error.log"
        with error_log.open('w', encoding='utf-8') as f:
            f.write(f"Return code: {e.returncode}\n\n")
            f.write(f"STDOUT:\n{e.stdout}\n\n")
            f.write(f"STDERR:\n{e.stderr}\n")

        raise

    except Exception as e:
        logger.error(f"Unexpected error running SATIVA: {e}")
        raise

    finally:
        # Covers every exit sativa.py does not clean up after itself: fatal error,
        # uncaught exception, and the timeout above (where the child is killed and
        # never reaches its own clean_tempdir).
        _clean_sativa_tmp(sativa_script, output_name)


def main():
    """Main entry point (Snakemake or standalone)."""
    start_time = time.time()

    try:
        # Configuration Snakemake
        if 'snakemake' in globals():
            msa_path = Path(snakemake.input.msa)
            output_dir = Path(snakemake.params.outdir)
            taxon = snakemake.params.taxon
            barcode = snakemake.params.barcode
            tax_level = snakemake.params.tax_level
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            max_seqs = snakemake.params.max_seqs
            max_len = snakemake.params.max_len
            threads = snakemake.threads

            sativa_script = Path(snakemake.params.sativa_script)
            sativa_mode = snakemake.params.sativa_mode
            confidence = snakemake.params.confidence
            taxonomy_code = snakemake.params.taxonomy_code
            rand_seed = snakemake.params.rand_seed

            # Configure file logging
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s | %(levelname)s | %(name)s | %(message)s'
            ))
            logger.addHandler(file_handler)

        else:
            # Standalone mode for testing
            parser = argparse.ArgumentParser(
                description="SATIVA wrapper for Snakemake pipeline"
            )
            parser.add_argument("--msa", type=str, required=True, help="MSA input file")
            parser.add_argument("--outdir", type=str, required=True, help="Output directory")
            parser.add_argument("--taxon", type=str, required=True, help="Taxon name")
            parser.add_argument("--barcode", type=str, required=True, help="Barcode (ITS/ITS1/ITS2)")
            parser.add_argument("--tax-level", type=str, required=True, help="Taxonomic level")
            parser.add_argument("--log-file", type=str, required=True, help="Log file")
            parser.add_argument("--max-seqs", type=int, default=20000, help="Max sequences")
            parser.add_argument("--max-len", type=int, default=5000, help="Max alignment length")
            parser.add_argument("--threads", type=int, default=4, help="Number of threads")
            parser.add_argument("--sativa-script", type=str, required=True, help="Path to sativa.py")
            parser.add_argument("--sativa-mode", type=str, default="ultrafast", help="SATIVA mode")
            parser.add_argument("--confidence", type=float, default=0.6, help="Confidence threshold")
            parser.add_argument("--taxonomy-code", type=str, default="BOT", help="Taxonomy code")
            parser.add_argument("--rand-seed", type=int, default=None, help="Random seed forwarded to SATIVA/RAxML")
            parser.add_argument("--timeout", type=int, default=21600, help="SATIVA timeout (seconds)")

            args = parser.parse_args()

            msa_path = Path(args.msa)
            output_dir = Path(args.outdir)
            taxon = args.taxon
            barcode = args.barcode
            tax_level = args.tax_level
            log_file = args.log_file
            max_seqs = args.max_seqs
            max_len = args.max_len
            threads = args.threads
            sativa_script = Path(args.sativa_script)
            sativa_mode = args.sativa_mode
            confidence = args.confidence
            taxonomy_code = args.taxonomy_code
            rand_seed = args.rand_seed

            # Configure file logging
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s | %(levelname)s | %(name)s | %(message)s'
            ))
            logger.addHandler(file_handler)

        # Check if output is already complete or needs cleaning
        if is_output_complete(output_dir, taxon):
            logger.info(f"Output already complete for {taxon}, skipping execution")
            sys.exit(0)
        else:
            # Clean incomplete/corrupted output
            clean_incomplete_output(output_dir)

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info("=" * 60)
        logger.info(f"SATIVA Pipeline - {taxon}")
        logger.info(f"Barcode: {barcode} | Tax level: {tax_level}")
        logger.info(f"MSA input: {msa_path}")
        logger.info(f"Output dir: {output_dir}")
        logger.info("=" * 60)

        # Step 1: Quick stats check
        logger.info("Computing MSA statistics...")
        n_seq, aln_len = quick_msa_stats(msa_path)
        logger.info(f"MSA stats: {n_seq} sequences, {aln_len} bp alignment length")

        # Step 2: Check thresholds
        # Size-based skips (n_seq > max_seqs, aln_len > max_len) REMOVED on purpose
        # (2026-08-14): SATIVA must run on EVERY group regardless of sequence count or
        # alignment length. Only the genuine technical floors below are kept — SATIVA
        # cannot build a reference tree from 0 or 1 sequence. max_seqs / max_len are no
        # longer used to gate execution.
        skip_reasons = []
        if n_seq == 0:
            skip_reasons.append("MSA is empty (no sequences)")
        elif n_seq < 2:
            skip_reasons.append("fewer than 2 sequences (SATIVA requires at least 2)")

        if skip_reasons:
            logger.warning(f"Skipping: {'; '.join(skip_reasons)}")

            # Create skip marker
            skip_marker = output_dir / ".skipped"
            skip_marker.touch()

            # Create summary JSON
            result = SativaResult(
                barcode=barcode,
                tax_level=tax_level,
                taxon=taxon,
                n_seq=n_seq,
                aln_len=aln_len,
                status="skipped",
                reason="; ".join(skip_reasons),
                runtime_sec=time.time() - start_time
            )

            summary_path = output_dir / "sativa.summary.json"
            with summary_path.open('w', encoding='utf-8') as f:
                json.dump(asdict(result), f, indent=2)

            logger.info(f"Summary written to: {summary_path}")
            logger.info("Exiting with success (skip)")
            sys.exit(0)

        # Step 3: Parse MSA to SATIVA format
        logger.info("Parsing MSA to SATIVA format...")
        fasta_path, tax_path = parse_msa_to_sativa_format(msa_path, output_dir, taxon)

        # Step 4: Run SATIVA
        logger.info("Launching SATIVA...")
        sativa_start = time.time()

        # Get timeout from params (default: 6 hours = 21600 seconds)
        timeout_seconds = snakemake.params.get('timeout_seconds', 21600) if 'snakemake' in globals() else args.timeout

        run_sativa(
            fasta_path=fasta_path,
            tax_path=tax_path,
            output_dir=output_dir,
            sativa_script=sativa_script,
            mode=sativa_mode,
            confidence=confidence,
            taxonomy_code=taxonomy_code,
            threads=threads,
            timeout_seconds=timeout_seconds,
            output_name=f"{barcode}_{tax_level}_{taxon}",
            rand_seed=rand_seed
        )

        sativa_runtime = time.time() - sativa_start

        # Step 5: Create success summary
        command_str = f"sativa.py -s {fasta_path.name} -t {tax_path.name} -m {sativa_mode} -C {confidence} -x {taxonomy_code} -T {threads}"
        if rand_seed is not None:
            command_str += f" -p {rand_seed}"

        result = SativaResult(
            barcode=barcode,
            tax_level=tax_level,
            taxon=taxon,
            n_seq=n_seq,
            aln_len=aln_len,
            status="success",
            runtime_sec=sativa_runtime,
            command=command_str
        )

        summary_path = output_dir / "sativa.summary.json"
        with summary_path.open('w', encoding='utf-8') as f:
            json.dump(asdict(result), f, indent=2)

        total_time = time.time() - start_time
        logger.info("=" * 60)
        logger.info("SATIVA pipeline completed successfully")
        logger.info(f"SATIVA runtime: {sativa_runtime:.2f}s")
        logger.info(f"Total runtime: {total_time:.2f}s")
        logger.info(f"Summary: {summary_path}")
        logger.info("=" * 60)

    except subprocess.TimeoutExpired as e:
        logger.warning(f"SATIVA timeout reached")

        # Create timeout summary
        result = SativaResult(
            barcode=barcode if 'barcode' in locals() else "unknown",
            tax_level=tax_level if 'tax_level' in locals() else "unknown",
            taxon=taxon if 'taxon' in locals() else "unknown",
            n_seq=n_seq if 'n_seq' in locals() else 0,
            aln_len=aln_len if 'aln_len' in locals() else 0,
            status="timeout",
            reason=f"Execution exceeded {timeout_seconds}s ({timeout_seconds/3600:.1f}h) timeout",
            runtime_sec=time.time() - start_time
        )

        summary_path = output_dir / "sativa.summary.json"
        with summary_path.open('w', encoding='utf-8') as f:
            json.dump(asdict(result), f, indent=2)

        logger.info(f"Timeout summary written to: {summary_path}")
        logger.info("Exiting with success (timeout)")
        sys.exit(0)

    except subprocess.CalledProcessError as e:
        # SATIVA (sativa.py) itself failed on THIS MSA — typically a degenerate alignment
        # (too few informative sites, marginal genus-only/small clade). This is a per-MSA
        # outcome, not a pipeline failure: record an "error" summary and exit 0 so the run
        # continues (the sequences of this MSA are simply left unvalidated).
        logger.warning(f"SATIVA failed on this MSA (recorded as error, not fatal): {e}")
        result = SativaResult(
            barcode=barcode if 'barcode' in locals() else "unknown",
            tax_level=tax_level if 'tax_level' in locals() else "unknown",
            taxon=taxon if 'taxon' in locals() else "unknown",
            n_seq=n_seq if 'n_seq' in locals() else 0,
            aln_len=aln_len if 'aln_len' in locals() else 0,
            status="error",
            reason=f"SATIVA execution failed: {str(e)[:400]}",
            runtime_sec=time.time() - start_time
        )
        summary_path = output_dir / "sativa.summary.json"
        with summary_path.open('w', encoding='utf-8') as f:
            json.dump(asdict(result), f, indent=2)
        logger.info(f"Error summary written to: {summary_path}; exiting 0 (per-MSA SATIVA failure)")
        sys.exit(0)

    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)

        # Try to write error summary
        try:
            result = SativaResult(
                barcode=barcode if 'barcode' in locals() else "unknown",
                tax_level=tax_level if 'tax_level' in locals() else "unknown",
                taxon=taxon if 'taxon' in locals() else "unknown",
                n_seq=n_seq if 'n_seq' in locals() else 0,
                aln_len=aln_len if 'aln_len' in locals() else 0,
                status="error",
                reason=str(e)[:500],  # Truncate long errors
                runtime_sec=time.time() - start_time
            )

            summary_path = output_dir / "sativa.summary.json"
            with summary_path.open('w', encoding='utf-8') as f:
                json.dump(asdict(result), f, indent=2)
        except Exception as summary_err:
            # Best-effort error summary; don't mask the original failure.
            logger.warning(f"Could not write error summary: {summary_err}")

        sys.exit(1)


if __name__ == "__main__":
    main()
