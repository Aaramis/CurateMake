# workflow/rules/benchmarking.smk
"""
Règle Snakemake pour générer le rapport de benchmarking.
"""

rule generate_benchmark_report:
    """Génère un rapport HTML interactif à partir de tous les fichiers benchmark disponibles."""
    input:
        # Dépend de toutes les étapes qui génèrent des benchmarks
        final_integration=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json",
        itsx_deduplicate_post="data/processed/itsx_deduplicate_post.json",
        msa_run_report="data/processed/msa_run_report.json",
        sativa_visualization=config["sativa"]["visualization_report"],
        msa_report=config["msa"]["output_report"],
        sequence_flow="reports/html/sequence_flow_report.html",
        taxonomy_sunburst=config["taxonomy_sunburst"]["output"]
    params:
        benchmarks_dir="benchmarks"
    output:
        report="reports/html/benchmark_report.html"
    log:
        "logs/benchmarking/generate_report.log"
    conda:
        "../envs/sativa.yaml"  # Contient pandas, plotly
    script:
        "../scripts/generate_benchmark_report.py"
