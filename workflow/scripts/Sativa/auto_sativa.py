#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script SATIVA pour Snakemake
Détection des misidentifications taxonomiques sur des MSA.
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


# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)
logger = logging.getLogger(__name__)


_SANITIZE_PATTERN = re.compile(r"[^A-Za-z0-9_.-]")


def _sanitize_label(value: str) -> str:
    """Convert arbitrary taxon/barcode labels into a filesystem-safe token."""
    return _SANITIZE_PATTERN.sub("_", value)


@dataclass
class SativaResult:
    """Résultat d'exécution SATIVA."""
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
    Vérifie si une exécution SATIVA précédente s'est terminée correctement.

    Args:
        output_dir: Répertoire de sortie SATIVA
        taxon: Nom du taxon

    Returns:
        True si le résultat est complet, False sinon
    """
    if not output_dir.exists():
        return False

    # Vérifier le fichier summary
    summary_path = output_dir / "sativa.summary.json"
    if not summary_path.exists():
        return False

    # Vérifier que le summary est valide
    try:
        with summary_path.open('r') as f:
            result = json.load(f)

        status = result.get('status')

        # Si skipped ou error, c'est considéré comme complet (pas besoin de relancer)
        if status in ['skipped', 'error']:
            return True

        # Si success, vérifier que les fichiers .mis existent
        if status == 'success':
            mis_file = output_dir / f"{taxon}.mis"
            if mis_file.exists():
                return True
            else:
                logger.warning(f"Summary indicates success but .mis file missing: {mis_file}")
                return False

        return False

    except (json.JSONDecodeError, Exception) as e:
        logger.warning(f"Invalid or corrupted summary file: {e}")
        return False


def clean_incomplete_output(output_dir: Path):
    """
    Nettoie un dossier de sortie incomplet ou corrompu.

    Args:
        output_dir: Répertoire à nettoyer
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
    Calcule rapidement les stats d'un MSA sans tout charger en RAM.

    Args:
        msa_path: Chemin vers le fichier MSA FASTA

    Returns:
        Tuple (nombre_sequences, longueur_alignement)
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
    Parse un MSA en fichiers .fasta et .tax pour SATIVA.

    Format MSA attendu:
    >seq_id|kingdom=X|phylum=Y|class=Z|order=A|family=B|genus=C|species=D
    SEQUENCE

    Args:
        msa_path: Fichier MSA d'entrée
        output_dir: Répertoire de sortie
        taxon: Nom du taxon (pour nommer les fichiers)

    Returns:
        Tuple (fasta_path, tax_path)
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
                tax_string = ';'.join(tax_info.get(rank, '') for rank in tax_order)

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
    output_name: Optional[str] = None
) -> bool:
    """
    Exécute SATIVA.

    Args:
        fasta_path: Fichier FASTA des séquences
        tax_path: Fichier de taxonomie
        output_dir: Répertoire de sortie (cwd pour SATIVA)
        sativa_script: Chemin vers sativa.py
        mode: Mode d'exécution SATIVA
        confidence: Seuil de confiance
        taxonomy_code: Code taxonomique
        threads: Nombre de threads
        timeout_seconds: Timeout en secondes (None = pas de limite)

    Returns:
        True si succès, False sinon

    Raises:
        subprocess.TimeoutExpired: Si le timeout est atteint
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


def main():
    """Point d'entrée principal."""
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
        skip_reasons = []
        if n_seq > max_seqs:
            skip_reasons.append(f"n_seq ({n_seq}) > max_seqs ({max_seqs})")
        if aln_len > max_len:
            skip_reasons.append(f"aln_len ({aln_len}) > max_len ({max_len})")

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
        timeout_seconds = snakemake.params.get('timeout_seconds', 21600) if 'snakemake' in globals() else None

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
            output_name=f"{barcode}_{tax_level}_{taxon}"
        )

        sativa_runtime = time.time() - sativa_start

        # Step 5: Create success summary
        command_str = f"sativa.py -s {fasta_path.name} -t {tax_path.name} -m {sativa_mode} -C {confidence} -x {taxonomy_code} -T {threads}"

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
        except:
            pass

        sys.exit(1)


if __name__ == "__main__":
    main()
