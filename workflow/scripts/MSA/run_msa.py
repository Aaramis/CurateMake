#!/usr/bin/env python3
"""
Aggregate the Snakemake-driven MSA outputs into the legacy JSON/Markdown report.
"""

import json
import logging
import time
from pathlib import Path
from typing import Dict, Any, List

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def add_reassignment_section_to_report(report_path: str, stats: dict):
    if not stats:
        return

    with open(report_path, 'a') as f:
        f.write("\n\n## Statistiques de réaffectation inter-barcodes\n\n")
        f.write("Cette section montre comment les séquences initialement ")
        f.write("associées à un barcode ont été réassignées après classification.\n\n")

        for source in ['ITS', 'ITS1', 'ITS2']:
            f.write(f"\n### Séquences source : {source}\n\n")
            f.write("| Barcode cible | Nombre | Pourcentage |\n")
            f.write("|--------------|---------|-------------|\n")

            for target in ['ITS', 'ITS1', 'ITS2']:
                key = f'{source}_to_{target}'
                if key in stats:
                    count = stats[key]['count']
                    pct = stats[key]['percentage']
                    f.write(f"| {target} | {count:,} | {pct:.1f}% |\n")


def generate_markdown_report(report_data: Dict[str, Any], report_md_path: str):
    with open(report_md_path, 'w') as f:
        f.write("# Rapport MSA - Alignements Hiérarchiques\n\n")
        f.write(f"**Date** : {report_data['timestamp']}\n")
        f.write(f"**Status** : {report_data['status']}\n")
        f.write(f"**Temps de traitement** : {report_data['processing_time']}s\n\n")

        f.write("## Paramètres\n\n")
        for key, value in report_data['parameters'].items():
            f.write(f"- **{key}** : {value}\n")

        f.write("\n## Résultats par barcode\n\n")
        results = report_data.get('results', {})
        barcode_results = results.get('barcode_results', {})

        for barcode, data in barcode_results.items():
            f.write(f"\n### {barcode}\n\n")
            for key, value in data.items():
                f.write(f"- **{key}** : {value}\n")

        f.write(f"\n## Résumé global\n\n")
        f.write(f"- **Barcodes traités** : {results.get('total_barcodes_processed', 0)}\n")
        f.write(f"- **Phases complétées** : {results.get('phases_completed', 0)}\n")

    logger.info(f"Rapport markdown sauvegardé: {report_md_path}")


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def summarize_barcodes(task_manifest: Dict[str, Any], classification_manifest: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    barcode_summary: Dict[str, Dict[str, Any]] = {}
    barcode_info = task_manifest.get("barcodes", {})
    classification_info = classification_manifest.get("barcode_summary", {})

    for barcode, payload in barcode_info.items():
        tasks = payload.get("tasks", []) + payload.get("large_tasks", [])
        sequences_processed = sum(task.get("sequence_count", 0) for task in tasks)
        sequences_unprocessed = len(payload.get("unprocessed_sequence_ids", []))

        classification_counts = classification_info.get(barcode, {})
        total_classification = classification_counts.get("total_sequences", sequences_unprocessed)

        barcode_summary[barcode] = {
            "tasks_processed": len(tasks),
            "sequences_processed": sequences_processed,
            "sequences_unprocessed": total_classification,
            "classification_sequences": total_classification
        }

    return barcode_summary


def count_fasta_sequences(fasta_path: Path) -> int:
    """Count sequences in a FASTA without loading it entirely into memory."""
    count = 0
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def compute_final_msa_totals(msa_root: str, barcodes: List[str]) -> Dict[str, int]:
    """
    Aggregate the number of sequences present in trimmed MSA files for each barcode.
    Only *_msa.trimmed.fasta files are considered to reflect the actual alignments used downstream.
    """
    totals: Dict[str, int] = {}
    base_path = Path(msa_root)

    for barcode in barcodes:
        barcode_dir = base_path / barcode
        if not barcode_dir.exists():
            logger.warning(f"MSA directory not found for barcode {barcode}: {barcode_dir}")
            totals[barcode] = 0
            continue

        final_count = 0
        for alignment_dir in barcode_dir.glob("alignments_*"):
            if not alignment_dir.is_dir():
                continue
            for trimmed_file in alignment_dir.glob("*_msa.trimmed.fasta"):
                final_count += count_fasta_sequences(trimmed_file)

        totals[barcode] = final_count
        logger.info(f"Final trimmed sequences for {barcode}: {final_count}")

    return totals


def build_report_data(params: Dict[str, Any], manifest_path: str, classification_manifest_path: str,
                      classification_stats_path: str) -> Dict[str, Any]:
    task_manifest = load_json(manifest_path)
    classification_manifest = load_json(classification_manifest_path)
    classification_stats = load_json(classification_stats_path)

    barcode_summary = summarize_barcodes(task_manifest, classification_manifest)
    final_counts = compute_final_msa_totals(params["output_folder"], params["parameters"].get("barcodes", []))
    for barcode, final_count in final_counts.items():
        if barcode in barcode_summary:
            barcode_summary[barcode]["sequences_processed"] = final_count
        else:
            barcode_summary[barcode] = {
                "tasks_processed": 0,
                "sequences_processed": final_count,
                "sequences_unprocessed": 0,
                "classification_sequences": 0
            }
    results = {
        "barcode_results": barcode_summary,
        "global_assignments": classification_stats,
        "total_barcodes_processed": len(barcode_summary),
        "phases_completed": 2
    }

    processing_time = round(params["elapsed"], 2)
    report_data: Dict[str, Any] = {
        "source": "Snakemake_MSA_scatter",
        "input_csv": params["input_csv"],
        "output_folder": params["output_folder"],
        "parameters": params["parameters"],
        "results": results,
        "processing_time": processing_time,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "status": "success"
    }
    return report_data


def main():
    start_time = time.time()
    try:
        if 'snakemake' in globals():
            params = snakemake.params
            inputs = snakemake.input
            output_report = snakemake.output.report
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)

            manifest_path = inputs.task_manifest
            classification_manifest_path = inputs.classification_manifest
            classification_stats_path = inputs.classification_stats
        else:
            logger.warning("Running in standalone mode")
            params = {
                "csv_path": "data/full.csv",
                "db": "",
                "output_folder": "data/MSA",
                "min_seq_threshold": 20,
                "max_seq_threshold": 1000,
                "large_task_threshold": 2000,
                "min_domain": "species",
                "trimal_seq_overlap": 75,
                "trimal_res_overlap": 0.5,
                "domain": "family",
                "barcodes": ["ITS", "ITS1", "ITS2"],
                "hierarchy_levels": ["genus", "family", "order", "class", "phylum"],
                "tax_columns": ["kingdom", "phylum", "class", "order", "family", "genus", "species"],
                "id_col": "seq_id",
                "seq_col": "sequence",
                "barcode_col": "barcode",
                "max_thread_workers": 1,
                "classification_max_workers": 1,
                "classification_threads_per_worker": 1,
                "classification_batch_size": 1000,
                "logger_folder": "logs/MSA"
            }
            output_report = "data/processed/msa_run_report.json"
            manifest_path = "data/processed/msa_task_manifest.json"
            classification_manifest_path = "data/processed/msa_classification_manifest.json"
            classification_stats_path = "data/processed/msa_classification_stats.json"

        parameter_block = {
            "min_seq_threshold": params["min_seq_threshold"],
            "trimal_seq_overlap": params["trimal_seq_overlap"],
            "trimal_res_overlap": params["trimal_res_overlap"],
            "domain": params["domain"],
            "barcodes": params["barcodes"],
            "hierarchy_levels": params["hierarchy_levels"],
            "tax_columns": params["tax_columns"],
            "id_col": params["id_col"],
            "seq_col": params["seq_col"],
            "barcode_col": params["barcode_col"],
            "classification_batch_size": params["classification_batch_size"],
        }

        runtime_params = {
            "input_csv": params["csv_path"],
            "output_folder": params["output_folder"],
            "parameters": parameter_block,
            "elapsed": time.time() - start_time
        }

        report_data = build_report_data(
            runtime_params,
            manifest_path,
            classification_manifest_path,
            classification_stats_path
        )

        Path(output_report).parent.mkdir(parents=True, exist_ok=True)
        with open(output_report, "w", encoding="utf-8") as fh:
            json.dump(report_data, fh, indent=2)
        logger.info("Rapport JSON sauvegardé: %s", output_report)

        report_md_path = output_report.replace('.json', '.md')
        generate_markdown_report(report_data, report_md_path)

        global_assignments = report_data['results'].get('global_assignments', {})
        add_reassignment_section_to_report(
            report_md_path,
            global_assignments.get('reassignment_stats', {})
        )
        logger.info("Statistiques de réaffectation ajoutées au rapport")

    except Exception as exc:
        logger.critical("Erreur fatale: %s", exc, exc_info=True)
        raise


if __name__ == "__main__":
    main()
