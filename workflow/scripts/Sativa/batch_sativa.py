#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch SATIVA processor with parallel execution
Processes multiple MSA files in parallel using ThreadPoolExecutor
"""

import json
import logging
import time
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from datetime import datetime
from typing import List, Tuple, Dict, Any, Optional
import sys

# Import functions from auto_sativa
sys.path.insert(0, str(Path(__file__).parent))
from auto_sativa import (
    quick_msa_stats,
    parse_msa_to_sativa_format,
    run_sativa,
    SativaResult,
    is_output_complete,
    clean_incomplete_output
)
from dataclasses import asdict


# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)
logger = logging.getLogger(__name__)


def cleanup_old_temp_files(sativa_script_path: Path, retention_hours: int = 24, cleanup_enabled: bool = True):
    """
    Nettoie les fichiers temporaires SATIVA anciens.

    Args:
        sativa_script_path: Chemin vers le script SATIVA
        retention_hours: Durée de rétention des fichiers temp (en heures)
        cleanup_enabled: Active/désactive le nettoyage
    """
    if not cleanup_enabled:
        logger.info("Cleanup des fichiers temporaires désactivé")
        return

    tmp_dir = sativa_script_path.parent / "sativa" / "tmp"

    if not tmp_dir.exists():
        logger.info(f"Dossier temporaire inexistant : {tmp_dir}")
        return

    logger.info(f"Nettoyage des fichiers temporaires dans : {tmp_dir}")
    logger.info(f"Rétention : {retention_hours}h")

    now = time.time()
    cleaned = 0
    total_size = 0

    try:
        for folder in tmp_dir.iterdir():
            if not folder.is_dir():
                continue

            age_hours = (now - folder.stat().st_mtime) / 3600

            if age_hours > retention_hours:
                # Calculer la taille avant suppression
                folder_size = sum(f.stat().st_size for f in folder.rglob('*') if f.is_file())
                total_size += folder_size

                logger.debug(f"Suppression : {folder.name} (âge: {age_hours:.1f}h, taille: {folder_size / 1024 / 1024:.1f} MB)")

                shutil.rmtree(folder, ignore_errors=True)
                cleaned += 1

        if cleaned > 0:
            logger.info(f"✓ {cleaned} dossiers temporaires supprimés ({total_size / 1024 / 1024:.1f} MB libérés)")
        else:
            logger.info("Aucun fichier temporaire à nettoyer")

    except Exception as e:
        logger.warning(f"Erreur lors du nettoyage des fichiers temporaires : {e}")


def _infer_taxon_from_filename(msa_file: Path) -> str:
    """
    Derive taxon name from an MSA filename that ends with
    `_msa.fasta` or `_msa.trimmed.fasta`.
    """
    name = msa_file.name
    if name.endswith("_msa.trimmed.fasta"):
        taxon = name[: -len("_msa.trimmed.fasta")]
    elif name.endswith("_msa.fasta"):
        taxon = name[: -len("_msa.fasta")]
    else:
        taxon = msa_file.stem.replace("_msa", "")
    return taxon


def get_all_msa_files(msa_dir: Path, barcodes: List[str], tax_levels: List[str]) -> List[Tuple[str, str, str, Path]]:
    """
    Récupère tous les fichiers MSA disponibles.

    Returns:
        Liste de tuples (barcode, tax_level, taxon, msa_path)
    """
    msa_files = []

    if not msa_dir.exists():
        logger.warning(f"MSA directory does not exist: {msa_dir}")
        return []

    # Parcourir les barcodes
    for barcode in barcodes:
        barcode_dir = msa_dir / barcode
        if not barcode_dir.exists():
            continue

        # Parcourir les niveaux taxonomiques
        for tax_level in tax_levels:
            tax_dir = barcode_dir / f"alignments_{tax_level}"
            if not tax_dir.exists():
                continue

            # Trouver tous les fichiers MSA
            trimmed_files = sorted(tax_dir.glob("*_msa.trimmed.fasta"))
            candidate_files = trimmed_files if trimmed_files else sorted(tax_dir.glob("*_msa.fasta"))

            for msa_file in candidate_files:
                taxon = _infer_taxon_from_filename(msa_file)
                msa_files.append((barcode, tax_level, taxon, msa_file))

    return msa_files


def process_single_msa(
    msa_path: Path,
    barcode: str,
    tax_level: str,
    taxon: str,
    output_dir: Path,
    sativa_script: Path,
    sativa_mode: str,
    confidence: float,
    taxonomy_code: str,
    threads: int,
    max_seqs: int,
    max_len: int,
    log_dir: Path,
    timeout_seconds: Optional[int] = None
) -> Dict[str, Any]:
    """
    Traite un seul MSA avec SATIVA.

    Returns:
        Dictionnaire avec les résultats (summary JSON)
    """
    start_time = time.time()

    # Define output and log directories
    msa_output_dir = output_dir / barcode / f"alignments_{tax_level}" / taxon
    log_file = log_dir / barcode / f"{tax_level}_{taxon}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    # Configure file logging for this MSA
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)s | %(name)s | %(message)s'
    ))
    task_logger = logging.getLogger(f"sativa.{barcode}.{tax_level}.{taxon}")
    task_logger.addHandler(file_handler)
    task_logger.setLevel(logging.INFO)

    try:
        task_logger.info("=" * 60)
        task_logger.info(f"SATIVA Pipeline - {taxon}")
        task_logger.info(f"Barcode: {barcode} | Tax level: {tax_level}")
        task_logger.info(f"MSA input: {msa_path}")
        task_logger.info(f"Output dir: {msa_output_dir}")
        task_logger.info("=" * 60)

        # Step 0: Check if already processed or needs cleaning
        if is_output_complete(msa_output_dir, taxon):
            task_logger.info(f"Output already complete, loading existing results")
            summary_path = msa_output_dir / "sativa.summary.json"
            with summary_path.open('r') as f:
                existing_result = json.load(f)
            task_logger.info(f"Existing result status: {existing_result.get('status', 'unknown')}")
            return existing_result
        else:
            # Clean incomplete/corrupted output
            clean_incomplete_output(msa_output_dir)

        # Create output directory (after cleaning if needed)
        msa_output_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Quick stats check
        task_logger.info("Computing MSA statistics...")
        n_seq, aln_len = quick_msa_stats(msa_path)
        task_logger.info(f"MSA stats: {n_seq} sequences, {aln_len} bp alignment length")

        # Step 2: Check thresholds
        skip_reasons = []
        if n_seq > max_seqs:
            skip_reasons.append(f"n_seq ({n_seq}) > max_seqs ({max_seqs})")
        if aln_len > max_len:
            skip_reasons.append(f"aln_len ({aln_len}) > max_len ({max_len})")

        if skip_reasons:
            task_logger.warning(f"Skipping: {'; '.join(skip_reasons)}")

            # Create skip marker
            skip_marker = msa_output_dir / ".skipped"
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

            summary_path = msa_output_dir / "sativa.summary.json"
            with summary_path.open('w', encoding='utf-8') as f:
                json.dump(asdict(result), f, indent=2)

            task_logger.info(f"Summary written to: {summary_path}")
            return asdict(result)

        # Step 3: Parse MSA to SATIVA format
        task_logger.info("Parsing MSA to SATIVA format...")
        fasta_path, tax_path = parse_msa_to_sativa_format(msa_path, msa_output_dir, taxon)

        # Step 4: Run SATIVA
        task_logger.info("Launching SATIVA...")
        if timeout_seconds:
            task_logger.info(f"Timeout set to: {timeout_seconds}s ({timeout_seconds/3600:.1f}h)")
        sativa_start = time.time()

        run_sativa(
            fasta_path=fasta_path,
            tax_path=tax_path,
            output_dir=msa_output_dir,
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

        summary_path = msa_output_dir / "sativa.summary.json"
        with summary_path.open('w', encoding='utf-8') as f:
            json.dump(asdict(result), f, indent=2)

        total_time = time.time() - start_time
        task_logger.info("=" * 60)
        task_logger.info("SATIVA pipeline completed successfully")
        task_logger.info(f"SATIVA runtime: {sativa_runtime:.2f}s")
        task_logger.info(f"Total runtime: {total_time:.2f}s")
        task_logger.info(f"Summary: {summary_path}")
        task_logger.info("=" * 60)

        return asdict(result)

    except subprocess.TimeoutExpired as e:
        task_logger.warning(f"SATIVA timeout reached after {timeout_seconds}s ({timeout_seconds/3600:.1f}h)")

        # Create timeout summary
        result = SativaResult(
            barcode=barcode,
            tax_level=tax_level,
            taxon=taxon,
            n_seq=n_seq if 'n_seq' in locals() else 0,
            aln_len=aln_len if 'aln_len' in locals() else 0,
            status="timeout",
            reason=f"Execution exceeded {timeout_seconds}s ({timeout_seconds/3600:.1f}h) timeout",
            runtime_sec=time.time() - start_time
        )

        summary_path = msa_output_dir / "sativa.summary.json"
        with summary_path.open('w', encoding='utf-8') as f:
            json.dump(asdict(result), f, indent=2)

        task_logger.info(f"Timeout summary written to: {summary_path}")
        return asdict(result)

    except Exception as e:
        task_logger.error(f"Fatal error: {e}", exc_info=True)

        # Create error summary
        result = SativaResult(
            barcode=barcode,
            tax_level=tax_level,
            taxon=taxon,
            n_seq=n_seq if 'n_seq' in locals() else 0,
            aln_len=aln_len if 'aln_len' in locals() else 0,
            status="error",
            reason=str(e)[:500],
            runtime_sec=time.time() - start_time
        )

        summary_path = msa_output_dir / "sativa.summary.json"
        with summary_path.open('w', encoding='utf-8') as f:
            json.dump(asdict(result), f, indent=2)

        return asdict(result)

    finally:
        # Remove handler to avoid memory leaks
        task_logger.removeHandler(file_handler)
        file_handler.close()


def main():
    """Point d'entrée principal."""
    overall_start_time = time.time()

    try:
        # Configuration Snakemake
        if 'snakemake' in globals():
            input_msa_dir = Path(snakemake.params.input_msa_dir)
            output_dir = Path(snakemake.params.output_dir)
            log_dir = Path(snakemake.params.logs_dir)
            report_path = Path(snakemake.output.report)

            barcodes = snakemake.params.barcodes
            tax_levels = snakemake.params.tax_levels
            max_seqs = snakemake.params.max_seqs
            max_len = snakemake.params.max_len
            sativa_script = Path(snakemake.params.sativa_script)
            sativa_mode = snakemake.params.sativa_mode
            confidence = snakemake.params.confidence
            taxonomy_code = snakemake.params.taxonomy_code
            threads_per_job = snakemake.params.threads
            max_parallel_jobs = snakemake.params.max_parallel_jobs
            timeout_seconds = snakemake.params.get('timeout_seconds', 21600)  # Default: 6 hours

            main_log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            # Configure main logging
            file_handler = logging.FileHandler(main_log_file, mode="w")
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s | %(levelname)s | %(name)s | %(message)s'
            ))
            logger.addHandler(file_handler)
        else:
            raise RuntimeError("This script must be run from Snakemake")

        logger.info("=" * 80)
        logger.info("SATIVA BATCH PROCESSOR")
        logger.info("=" * 80)
        logger.info(f"Input MSA directory: {input_msa_dir}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Barcodes: {', '.join(barcodes)}")
        logger.info(f"Tax levels: {', '.join(tax_levels)}")
        logger.info(f"Max parallel jobs: {max_parallel_jobs}")
        logger.info(f"Threads per job: {threads_per_job}")
        logger.info(f"Timeout per task: {timeout_seconds}s ({timeout_seconds/3600:.1f}h)")
        logger.info("=" * 80)

        # Get all MSA files
        logger.info("Scanning for MSA files...")
        msa_files = get_all_msa_files(input_msa_dir, barcodes, tax_levels)
        logger.info(f"Found {len(msa_files)} MSA files to process")

        if not msa_files:
            logger.warning("No MSA files found!")
            # Create empty report
            report_data = {
                'total_msa_processed': 0,
                'summary': {'success': 0, 'skipped': 0, 'error': 0, 'timeout': 0},
                'by_barcode': {},
                'by_tax_level': {},
                'results': {'success': [], 'skipped': [], 'error': [], 'timeout': []}
            }
            report_path.parent.mkdir(parents=True, exist_ok=True)
            with report_path.open('w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Empty report written to {report_path}")
            return

        # Process MSA files in parallel
        logger.info(f"Starting parallel processing with {max_parallel_jobs} workers...")
        all_results = []

        with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
            futures = {}

            for barcode, tax_level, taxon, msa_path in msa_files:
                future = executor.submit(
                    process_single_msa,
                    msa_path=msa_path,
                    barcode=barcode,
                    tax_level=tax_level,
                    taxon=taxon,
                    output_dir=output_dir,
                    sativa_script=sativa_script,
                    sativa_mode=sativa_mode,
                    confidence=confidence,
                    taxonomy_code=taxonomy_code,
                    threads=threads_per_job,
                    max_seqs=max_seqs,
                    max_len=max_len,
                    log_dir=log_dir,
                    timeout_seconds=timeout_seconds
                )
                futures[future] = (barcode, tax_level, taxon)

            # Process results as they complete
            completed = 0
            for future in as_completed(futures):
                barcode, tax_level, taxon = futures[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    completed += 1
                    status = result.get('status', 'unknown')
                    logger.info(f"[{completed}/{len(msa_files)}] {barcode}/{tax_level}/{taxon}: {status}")
                except Exception as e:
                    logger.error(f"Error processing {barcode}/{tax_level}/{taxon}: {e}")
                    completed += 1

        # Cleanup temporary files
        logger.info("Cleaning up temporary files...")
        cleanup_enabled = snakemake.config.get('sativa', {}).get('cleanup_temp_files', True) if 'snakemake' in globals() else True
        retention_hours = snakemake.config.get('sativa', {}).get('temp_retention_hours', 24) if 'snakemake' in globals() else 24
        cleanup_old_temp_files(sativa_script, retention_hours, cleanup_enabled)

        # Generate report
        logger.info("Generating final report...")
        results_by_status = defaultdict(list)
        stats_by_barcode = defaultdict(lambda: {'total': 0, 'success': 0, 'skipped': 0, 'error': 0, 'timeout': 0})
        stats_by_tax_level = defaultdict(lambda: {'total': 0, 'success': 0, 'skipped': 0, 'error': 0, 'timeout': 0})

        for result in all_results:
            status = result.get('status', 'unknown')
            barcode = result.get('barcode', 'unknown')
            tax_level = result.get('tax_level', 'unknown')

            results_by_status[status].append(result)
            stats_by_barcode[barcode]['total'] += 1
            stats_by_barcode[barcode][status] += 1
            stats_by_tax_level[tax_level]['total'] += 1
            stats_by_tax_level[tax_level][status] += 1

        report_data = {
            'total_msa_processed': len(all_results),
            'summary': {
                'success': len(results_by_status['success']),
                'skipped': len(results_by_status['skipped']),
                'error': len(results_by_status['error']),
                'timeout': len(results_by_status['timeout'])
            },
            'by_barcode': dict(stats_by_barcode),
            'by_tax_level': dict(stats_by_tax_level),
            'results': dict(results_by_status)
        }

        report_path.parent.mkdir(parents=True, exist_ok=True)
        with report_path.open('w') as f:
            json.dump(report_data, f, indent=2)

        total_time = time.time() - overall_start_time
        logger.info("=" * 80)
        logger.info("BATCH PROCESSING COMPLETE")
        logger.info(f"Total MSA processed: {len(all_results)}")
        logger.info(f"Success: {len(results_by_status['success'])}")
        logger.info(f"Skipped: {len(results_by_status['skipped'])}")
        logger.info(f"Timeouts: {len(results_by_status['timeout'])}")
        logger.info(f"Errors: {len(results_by_status['error'])}")
        logger.info(f"Total runtime: {total_time:.2f}s")
        logger.info(f"Report written to: {report_path}")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Fatal error in batch processor: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
