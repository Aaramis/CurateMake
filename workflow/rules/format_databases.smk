# workflow/rules/format_databases.smk

# Règle générique : input depuis config.yaml, script spécifique par DB
rule format_database:
    input:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        source = lambda wc: config["databases"][wc.database]["path"]
    output:
        report = "data/processed/{database}_formatted_report.json",
    log:
        "logs/format_databases/{database}_formatting.log"
    conda:
        "../envs/format_databases.yaml"
    params:
        database_name = "{database}",
        database_config = lambda wc: config["databases"][wc.database],
        sequence_validation = lambda wc: config["sequence_validation"]
    script:
        "../scripts/format_{wildcards.database}.py"


# Agrégation : déclenche toutes les formatations puis produit un rapport
rule format_all_databases:
    input:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        formatted = expand("data/processed/{database}_formatted.fasta", database=DATABASES)
    output:
        summary = "data/processed/formatting_summary.json"
    log:
        "logs/format_databases/all_formatting.log"
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/format_all_databases.py"



