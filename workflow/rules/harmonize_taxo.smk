# workflow/rules/harmonize_taxo.smk
"""
Rules Snakemake pour le matching taxonomique avec Catalogue of Life (COL).
"""

import os

rule col_taxonomic_matching:
    input:
        formatting_summary="data/processed/formatting_summary.json"
        # deduplication="data/processed/taxonomy_deduplication_report.json"
    params:
        db=get_db_path(),
        batch_size=config["CoL_api"]["batch_size"]
    output:
        report="data/processed/col_matching_report.json"
    log:
        "logs/harmonize_taxo/col_matching.log"
    benchmark:
        "benchmarks/harmonize_taxo/col_taxonomic_matching.benchmark.txt"
    conda:
        "../envs/harmonize_taxo.yaml"
    script:
        "../scripts/col_taxonomic_matcher.py"