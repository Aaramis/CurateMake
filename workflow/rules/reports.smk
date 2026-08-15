# workflow/rules/reports.smk
"""
Snakemake rules for generating HTML reports
"""

import os


_FLOW_CL = config["msa"].get("clustered", {})
_FLOW_COMMON_IN = dict(
    final_integration=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json",
    formatting_summary="results/processed/formatting_summary.json",
    itsx_deduplicate_pre="results/processed/itsx_deduplicate_pre.json",
    itsx_deduplicate_post="results/processed/itsx_deduplicate_post.json",
    itsx_process="results/processed/itsx_process.json",
    itsx_prepare="results/processed/itsx_prepare_report.json",
    msa_prepare="results/processed/msa_prepare_report.json",
)

rule sequence_flow_report:
    """Flow report (bottom-up: reads the bottom-up summary)."""
    input:
        bottomup_summary=_FLOW_CL.get("summary", "results/MSA/clustered/clustered_summary.json"),
        **_FLOW_COMMON_IN
    output:
        report="reports/html/sequence_flow_report.html"
    params:
        processed_dir="results/processed"
    log:
        "logs/reports/sequence_flow_report.log"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/generate_sequence_flow_report.py"


rule sativa_taxonomy_report:
    """
    Generate a detailed SATIVA taxonomy report with a 4-way comparison and KPIs.
    Runs after SATIVA integration.
    """
    input:
        final_integration=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json",
        db=get_db_path(),  # Uses the helper that auto-detects SQLite or PostgreSQL
        sativa_report=config["sativa"]["visualization_report"]
    output:
        report="reports/html/sativa_taxonomy_report.html"
    log:
        "logs/reports/sativa_taxonomy_report.log"
    benchmark:
        "benchmarks/reports/sativa_taxonomy_report.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/generate_sativa_taxonomy_report.py"


rule generate_reports_index:
    """
    Generate an HTML index to navigate all pipeline reports
    """
    input:
        msa_report=config["msa"]["output_report"],
        sequence_flow_report="reports/html/sequence_flow_report.html",
        benchmark_report="reports/html/benchmark_report.html",
        sativa_report=config["sativa"]["visualization_report"],
        sativa_taxonomy_report="reports/html/sativa_taxonomy_report.html"
    params:
        reports_dir="reports/html"
    output:
        index="reports/html/index.html"
    log:
        "logs/reports/generate_index.log"
    benchmark:
        "benchmarks/reports/generate_index.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/generate_index_report.py"
