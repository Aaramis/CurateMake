# workflow/rules/harmonize_taxo.smk
"""
Règles Snakemake pour le matching taxonomique avec Catalogue of Life (COL).
"""

import os

rule col_taxonomic_matching:
    input:
        formatting_summary="data/processed/formatting_summary.json"
        # deduplication="data/processed/taxonomy_deduplication_report.json"
    params:
        db=os.path.join(config["database"]["path"], config["database"]["name"]),
    output:
        report="data/processed/col_matching_report.json"
    log:
        "logs/harmonize_taxo/col_matching.log"
    conda:
        "../envs/harmonize_taxo.yaml"
    script:
        "../scripts/col_taxonomic_matcher.py"