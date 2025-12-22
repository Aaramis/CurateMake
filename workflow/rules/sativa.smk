"""
Snakemake rules for SATIVA: orchestrate per-MSA SATIVA runs, reporting, visualization,
and database integration.
"""

import json
import os
from pathlib import Path


SATIVA_MANIFEST = "data/processed/sativa_task_manifest.json"
SATIVA_REPORT = "data/processed/sativa_analysis_report.json"
SATIVA_OUTPUT_DIR = config["sativa"]["output_dir"]


def _load_sativa_manifest():
    manifest_path = checkpoints.sativa_task_manifest.get().output.manifest
    with open(manifest_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _sativa_tasks():
    manifest = _load_sativa_manifest()
    return manifest.get("tasks", [])


def _sativa_summary_path(task):
    return os.path.join(
        SATIVA_OUTPUT_DIR,
        task["barcode"],
        f"alignments_{task['tax_level']}",
        task["taxon"],
        "sativa.summary.json",
    )


def sativa_summary_inputs(_wildcards):
    return [_sativa_summary_path(task) for task in _sativa_tasks()]


def sativa_task_msa_path(wildcards):
    for task in _sativa_tasks():
        if (
            task["barcode"] == wildcards.barcode
            and task["tax_level"] == wildcards.level
            and task["taxon"] == wildcards.taxon
        ):
            return task["msa_path"]
    raise ValueError(f"Task not found in SATIVA manifest: {wildcards}")


# ---------------------------------------------------------------------------
# Manifest generation (checkpoint)
# ---------------------------------------------------------------------------
checkpoint sativa_task_manifest:
    input:
        msa_report="data/processed/msa_run_report.json"
    output:
        manifest=SATIVA_MANIFEST
    params:
        input_msa_dir=config["sativa"]["input_msa_dir"],
        barcodes=config["sativa"]["barcodes"],
        tax_levels=config["sativa"]["tax_levels"]
    log:
        "logs/sativa/manifest.log"
    benchmark:
        "benchmarks/sativa/manifest.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/prepare_sativa_manifest.py"


# ---------------------------------------------------------------------------
# Per-MSA SATIVA execution
# ---------------------------------------------------------------------------
rule sativa_task:
    input:
        manifest=SATIVA_MANIFEST,
        msa=lambda wildcards: sativa_task_msa_path(wildcards)
    output:
        summary=os.path.join(
            SATIVA_OUTPUT_DIR,
            "{barcode}",
            "alignments_{level}",
            "{taxon}",
            "sativa.summary.json",
        )
    params:
        outdir=lambda wildcards: os.path.join(SATIVA_OUTPUT_DIR, wildcards.barcode, f"alignments_{wildcards.level}", wildcards.taxon),
        taxon=lambda wildcards: wildcards.taxon,
        barcode=lambda wildcards: wildcards.barcode,
        tax_level=lambda wildcards: wildcards.level,
        max_seqs=config["sativa"]["max_seq_count"],
        max_len=config["sativa"]["max_alignment_length"],
        sativa_script=config["sativa"]["sativa_script_path"],
        sativa_mode=config["sativa"]["sativa_mode"],
        confidence=config["sativa"]["confidence_threshold"],
        taxonomy_code=config["sativa"]["taxonomy_code"],
        timeout_seconds=config["sativa"]["timeout_seconds"],
        rand_seed=config["sativa"].get("rand_seed")
    threads: config["sativa"]["threads"]
    resources:
        sativa_jobs=1
    log:
        "logs/sativa/tasks/{barcode}/{level}/{taxon}.log"
    benchmark:
        "benchmarks/sativa/tasks/{barcode}.{level}.{taxon}.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/auto_sativa.py"

# ---------------------------------------------------------------------------
# Aggregate report
# ---------------------------------------------------------------------------
rule sativa_analysis_report:
    input:
        manifest=SATIVA_MANIFEST,
        summaries=sativa_summary_inputs
    output:
        report=SATIVA_REPORT
    params:
        cleanup_temp_files=config["sativa"].get("cleanup_temp_files", True),
        temp_retention_hours=config["sativa"].get("temp_retention_hours", 24),
        sativa_script_path=config["sativa"]["sativa_script_path"]
    log:
        "logs/sativa/collect_report.log"
    benchmark:
        "benchmarks/sativa/collect_report.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/collect_sativa_results.py"


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
rule sativa_visualization:
    input:
        report=SATIVA_REPORT,
        integration_report=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json"
    output:
        report=config["sativa"]["visualization_report"]
    params:
        sativa_output_dir=SATIVA_OUTPUT_DIR
    log:
        "logs/sativa/sativa_visualization.log"
    benchmark:
        "benchmarks/sativa/sativa_visualization.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/visualisation_sativa.py"


# ---------------------------------------------------------------------------
# SATIVA integration into the database
# ---------------------------------------------------------------------------
rule sativa_integrate:
    input:
        sativa_report=SATIVA_REPORT,
        db=get_db_path()
    output:
        report=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json"
    params:
        db_path=get_db_path(),
        sativa_output_dir=SATIVA_OUTPUT_DIR,
        barcodes=config["sativa"]["barcodes"],
        tax_levels=config["sativa"]["tax_levels"]
    log:
        "logs/sativa/integrate_sativa.log"
    benchmark:
        "benchmarks/sativa/sativa_integrate.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/integrate_sativa.py"
