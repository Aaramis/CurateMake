"""
Snakemake rules for SATIVA: orchestrate per-MSA SATIVA runs, reporting, visualization,
and database integration.
"""

import json
import os
from pathlib import Path


SATIVA_MANIFEST = "results/processed/sativa_task_manifest.json"
SATIVA_REPORT = "results/processed/sativa_analysis_report.json"
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


def _sativa_integration_path(task):
    return os.path.join(
        SATIVA_OUTPUT_DIR,
        task["barcode"],
        f"alignments_{task['tax_level']}",
        task["taxon"],
        "integration.json",
    )


def sativa_integration_inputs(_wildcards):
    return [_sativa_integration_path(task) for task in _sativa_tasks()]


# NOTE on SATIVA memory (2026-08-14/15): a sativa_task peaks around ~20 GB RSS. The cost is
# driven by the reference tree epa-ng builds internally for the leave-one-out placement
# (ref.nwk + its conditional likelihood vectors), NOT by the size of the input MSA -- an
# estimate derived from the input file size underestimates it several-fold and is not safe
# to schedule on. Throttle by job count instead: `--resources sativa_jobs=N`, with
# N ~= memory_budget_GB / 20. LAUNCH.sh does this.
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
# SATIVA consumes the manifest produced by clustered_collect (already in the right
# format, {taxon}=msa_id filesystem-safe).
checkpoint sativa_task_manifest:
    input:
        manifest=CL_SATIVA_MANIFEST
    output:
        manifest=SATIVA_MANIFEST
    log:
        "logs/sativa/manifest.log"
    shell:
        "cp {input.manifest} {output.manifest}"


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
        rand_seed_arg=(f"--rand-seed {config['sativa']['rand_seed']}"
                       if config["sativa"].get("rand_seed") is not None else "")
    threads: config["sativa"]["threads"]
    resources:
        sativa_jobs=1
    log:
        "logs/sativa/tasks/{barcode}/{level}/{taxon}.log"
    benchmark:
        "benchmarks/sativa/tasks/{barcode}.{level}.{taxon}.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    # Run the real script file (shell) rather than Snakemake's `script:` temp-copy
    # mechanism, which races on the parallel FS under high job concurrency
    # (.snakemake/scripts/tmp*.auto_sativa.py missing / 0-byte).
    shell:
        """
        python workflow/scripts/Sativa/auto_sativa.py \
            --msa {input.msa} \
            --outdir {params.outdir} \
            --taxon {params.taxon} \
            --barcode {params.barcode} \
            --tax-level {params.tax_level} \
            --log-file {log} \
            --max-seqs {params.max_seqs} \
            --max-len {params.max_len} \
            --threads {threads} \
            --sativa-script {params.sativa_script} \
            --sativa-mode {params.sativa_mode} \
            --confidence {params.confidence} \
            --taxonomy-code {params.taxonomy_code} \
            --timeout {params.timeout_seconds} \
            {params.rand_seed_arg}
        """

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
# SATIVA integration into the database — scattered per clade
#
# Each clade is integrated as soon as its SATIVA job finishes, so the DB/CoL I/O
# (the bulk of integration wall-clock, ~0 CPU) overlaps with the remaining SATIVA
# compute of other clades instead of running as one job after the barrier.
# Concurrent writes are safe (add_taxonomy handles the insert race; bulk ops +
# writer_retry absorb contention). db_lock keeps the same lever as the format rules
# (no effect unless `--resources db_lock=N` is passed; use db_lock=1 to serialise
# SQLite). On PostgreSQL/YugabyteDB the jobs run concurrently under `-j`.
# ---------------------------------------------------------------------------
rule sativa_integrate_task:
    input:
        summary=os.path.join(
            SATIVA_OUTPUT_DIR, "{barcode}", "alignments_{level}", "{taxon}", "sativa.summary.json"
        ),
        db=get_db_path()
    output:
        report=os.path.join(
            SATIVA_OUTPUT_DIR, "{barcode}", "alignments_{level}", "{taxon}", "integration.json"
        )
    params:
        sativa_output_dir=SATIVA_OUTPUT_DIR
    resources:
        db_lock=GLOBAL_DB_LOCK
    retries: 3            # same transient temp-script race as sativa_task
    log:
        "logs/sativa/integrate_tasks/{barcode}/{level}/{taxon}.log"
    benchmark:
        "benchmarks/sativa/integrate_tasks/{barcode}.{level}.{taxon}.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/integrate_sativa.py"


# Barrier: aggregate the per-clade integration reports into the single report the
# downstream steps expect. Keeps the rule name/output stable for sativa_visualization.
rule sativa_integrate:
    input:
        clade_reports=sativa_integration_inputs
    output:
        report=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json"
    log:
        "logs/sativa/integrate_sativa.log"
    benchmark:
        "benchmarks/sativa/sativa_integrate.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/collect_sativa_integration.py"
