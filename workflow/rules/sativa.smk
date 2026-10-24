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

# --- consensus vote over independent splits --------------------------------------------
# See config/config.yaml (sativa.replicates) for the measurements. `replicates: 1` leaves
# the DAG exactly as it was: no reslice job, no vote job, and the merge below is a copy.
SATIVA_REPLICATES = int(config["sativa"].get("replicates", 1))
# kept as text: vote_sativa.py reads it as an exact Fraction ("2/3"), because a decimal
# would ask for 3 votes out of 3 where 2/3 means 2 (0.667 x 3 = 2.001).
SATIVA_VOTE_THRESHOLD = str(config["sativa"].get("vote_threshold", "2/3"))
SATIVA_BARCODES = config["sativa"]["barcodes"]


def _load_sativa_manifest():
    manifest_path = checkpoints.sativa_task_manifest.get().output.manifest
    with open(manifest_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _sativa_tasks():
    manifest = _load_sativa_manifest()
    return manifest.get("tasks", [])


def _native_tasks():
    """Replicate 0 only — the split the database and the reports are built on.

    Replicates exist to be voted on, not to be integrated N times: integrating each one
    would take sativa_warnings from 1.5 M rows to 15 M and repeat the integration work N
    times, while every single-split flag is ~63% artefact anyway (2026-08-26 measurement).
    """
    return [t for t in _sativa_tasks() if int(t.get("replicate", 0)) == 0]


def _sativa_summary_path(task):
    return os.path.join(
        SATIVA_OUTPUT_DIR,
        task["barcode"],
        f"alignments_{task['tax_level']}",
        task["taxon"],
        "sativa.summary.json",
    )


def sativa_summary_inputs(_wildcards):
    # Report on the native split so the figures stay comparable across runs; the replicate
    # jobs are pulled in by sativa_vote instead.
    return [_sativa_summary_path(task) for task in _native_tasks()]


def sativa_vote_inputs(wildcards):
    """Every replicate of one barcode — this is what forces the replicate SATIVA jobs."""
    return [_sativa_summary_path(t) for t in _sativa_tasks()
            if t["barcode"] == wildcards.barcode and t["tax_level"] == wildcards.level]


def _consensus_path(barcode, level):
    return os.path.join(SATIVA_OUTPUT_DIR, barcode, f"alignments_{level}", "consensus.mis")


def sativa_consensus_input(wildcards):
    """The clade integrates the consensus when there is one, its own .mis otherwise."""
    if SATIVA_REPLICATES <= 1:
        return []
    return [_consensus_path(wildcards.barcode, wildcards.level)]


def _sativa_integration_path(task):
    return os.path.join(
        SATIVA_OUTPUT_DIR,
        task["barcode"],
        f"alignments_{task['tax_level']}",
        task["taxon"],
        "integration.json",
    )


def sativa_integration_inputs(_wildcards):
    return [_sativa_integration_path(task) for task in _native_tasks()]


# NOTE on SATIVA memory (2026-08-14/15): a sativa_task peaks around ~20 GB RSS. The cost is
# driven by the reference tree epa-ng builds internally for the leave-one-out placement
# (ref.nwk + its conditional likelihood vectors), NOT by the size of the input MSA -- an
# estimate derived from the input file size underestimates it several-fold and is not safe
# to schedule on. Throttle by job count instead: `--resources sativa_jobs=N`, with
# N ~= memory_budget_GB / 20. Nothing applies it for you: without that flag the resource
# declared below is unlimited, and concurrency is floor(--cores / sativa.threads), which
# at --cores 120 is 15 tasks and about 300 GB. docs/usage.md lists the flag.
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
# SATIVA consumes the manifest produced by collect_groups (already in the right
# format, {taxon}=msa_id filesystem-safe).
def _replicate_manifests(_wildcards):
    if SATIVA_REPLICATES <= 1:
        return []
    return expand(
        os.path.join(EH_REPL_ROOT, "{barcode}", "seed{seed}", "sativa_manifest.json"),
        barcode=SATIVA_BARCODES, seed=range(1, SATIVA_REPLICATES),
    )


checkpoint sativa_task_manifest:
    input:
        manifest=CL_SATIVA_MANIFEST,
        replicates=_replicate_manifests
    output:
        manifest=SATIVA_MANIFEST
    params:
        replicate_args=lambda w, input: (
            "--replicate-manifest " + " ".join(input.replicates) if input.replicates else ""
        )
    log:
        "logs/sativa/manifest.log"
    shell:
        """
        python workflow/scripts/Sativa/merge_replicate_manifests.py \
            --native {input.manifest} {params.replicate_args} \
            --output {output.manifest} --log-file {log}
        """


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
        sativa_binary=config["sativa"].get("sativa_binary", "sativa-epang"),
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
        # The packaged tool, not the copy that used to live in this repository. Its own
        # environment: see envs/sativa_epang.yaml for why it cannot share ours.
        "../envs/sativa_epang.yaml"
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
            --sativa-binary {params.sativa_binary} \
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
        sativa_output_dir=SATIVA_OUTPUT_DIR
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
# Consensus vote across the independent splits
#
# One job per barcode: reads every replicate's .mis, votes rank by rank, and writes ONE
# consensus in the SATIVA .mis layout. integrate_sativa.py then consumes it through the
# very same SativaMisParser as a raw SATIVA output — there is no second path into the
# database. A sequence below the threshold is simply absent from the consensus, so the
# clade's .tax marks it VALIDATED, which is exactly the intended meaning: no reproducible
# objection to its current label.
#
# Not defined at all when replicates == 1, so the historical DAG is untouched.
# ---------------------------------------------------------------------------
if SATIVA_REPLICATES > 1:

    rule sativa_vote:
        input:
            manifest=SATIVA_MANIFEST,
            summaries=sativa_vote_inputs,
            # members of the NATIVE split: one clade per sequence, scaffold copies excluded
            members=os.path.join(EH_OUT, "{barcode}", "group_members.tsv")
        output:
            consensus=os.path.join(SATIVA_OUTPUT_DIR, "{barcode}", "alignments_{level}",
                                   "consensus.mis"),
            report=os.path.join(SATIVA_OUTPUT_DIR, "{barcode}", "alignments_{level}",
                                "vote_report.json")
        params:
            sativa_dir=SATIVA_OUTPUT_DIR,
            replicates=SATIVA_REPLICATES,
            threshold=SATIVA_VOTE_THRESHOLD
        log:
            "logs/sativa/vote/{barcode}.{level}.log"
        benchmark:
            "benchmarks/sativa/vote/{barcode}.{level}.benchmark.txt"
        conda:
            "../envs/sativa.yaml"
        shell:
            """
            python workflow/scripts/Sativa/vote_sativa.py \
                --manifest {input.manifest} \
                --sativa-dir {params.sativa_dir} \
                --barcode {wildcards.barcode} --level {wildcards.level} \
                --replicates {params.replicates} \
                --native-members {input.members} \
                --vote-threshold {params.threshold} \
                --output {output.consensus} --report {output.report} \
                --log-file {log}
            """


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
        consensus=sativa_consensus_input,
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
