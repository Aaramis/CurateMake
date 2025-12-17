# workflow/rules/reports.smk
"""
Règles Snakemake pour la génération de rapports HTML
"""

import os


rule sequence_flow_report:
    """
    Génère un rapport visuel interactif du flux de séquences à travers le pipeline
    """
    input:
        final_integration=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json",
        formatting_summary="data/processed/formatting_summary.json",
        itsx_deduplicate_pre="data/processed/itsx_deduplicate_pre.json",
        itsx_deduplicate_post="data/processed/itsx_deduplicate_post.json",
        itsx_process="data/processed/itsx_process.json",
        itsx_prepare="data/processed/itsx_prepare_report.json",
        msa_prepare="data/processed/msa_prepare_report.json",
        msa_run="data/processed/msa_run_report.json"
    output:
        report="reports/html/sequence_flow_report.html"
    params:
        processed_dir="data/processed"
    log:
        "logs/reports/sequence_flow_report.log"
    benchmark:
        "benchmarks/reports/sequence_flow_report.benchmark.txt"
    conda:
        "../envs/sativa.yaml"  # Utilise l'environnement SATIVA qui contient plotly
    script:
        "../scripts/generate_sequence_flow_report.py"


rule sativa_taxonomy_report:
    """
    Génère un rapport taxonomique SATIVA détaillé avec comparaison 4-way et KPIs.
    S'exécute après SATIVA integration.
    """
    input:
        final_integration=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json",
        db=get_db_path(),  # Utilise la fonction qui détecte automatiquement SQLite ou PostgreSQL
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
    Génère un index HTML pour naviguer entre tous les rapports du pipeline
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
    script:
        "../scripts/generate_index_report.py"
