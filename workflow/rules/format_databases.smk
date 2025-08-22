# workflow/rules/format_databases.smk

# Règle de déduplication des taxonomies
rule deduplicate_taxonomies:
    input:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        formatted = expand("data/processed/{database}_formatted_report.json", database=DATABASES),
    output:
        report = "data/processed/taxonomy_deduplication_report.json"
    log:
        "logs/format_databases/taxonomy_deduplication.log"
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/deduplicate_taxonomies.py"


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

rule format_all_databases:
    input:
        deduplication = "data/processed/taxonomy_deduplication_report.json"
    output:
        summary = "data/processed/formatting_summary.json"
    log:
        "logs/format_databases/all_formatting.log"
    params:
        formatted = expand("data/processed/{database}_formatted_report.json", database=DATABASES)
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/format_all_databases.py"

