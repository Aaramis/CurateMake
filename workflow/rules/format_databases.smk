# workflow/rules/format_databases.smk

# Règle UNITE
rule format_unite_database:
    input:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        source = config["databases"]["unite"]["path"]
    output:
        report = "data/processed/unite_formatted_report.json"
    log:
        "logs/format_databases/unite_formatting.log"
    conda:
        "../envs/format_databases.yaml"
    params:
        database_name = "unite",
        database_config = config["databases"]["unite"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_unite.py"


# Règle PLANITS
rule format_planits_database:
    input:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        source = config["databases"]["planits"]["path"]
    output:
        report = "data/processed/planits_formatted_report.json"
    log:
        "logs/format_databases/planits_formatting.log"
    conda:
        "../envs/format_databases.yaml"
    params:
        database_name = "planits",
        database_config = config["databases"]["planits"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_planits.py"

# Règle Bold
rule format_bold_database:
    input:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        source = config["databases"]["bold"]["path"]
    output:
        report = "data/processed/bold_formatted_report.json"
    log:
        "logs/format_databases/bold_formatting.log"
    conda:
        "../envs/format_databases.yaml"
    params:
        database_name = "bold",
        database_config = config["databases"]["bold"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_bold.py"


# Règle de déduplication des taxonomies
rule deduplicate_taxonomies:
    input:
        formatted = [
            # "data/processed/unite_formatted_report.json",
            # "data/processed/planits_formatted_report.json",
            "data/processed/bold_formatted_report.json"
        ]
    output:
        report = "data/processed/taxonomy_deduplication_report.json"
    params:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
    log:
        "logs/format_databases/taxonomy_deduplication.log"
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/deduplicate_taxonomies.py"


# Règle finale
rule format_all_databases:
    input:
        deduplication = "data/processed/taxonomy_deduplication_report.json"
    output:
        summary = "data/processed/formatting_summary.json"
    log:
        "logs/format_databases/all_formatting.log"
    params:
        formatted = [
            "data/processed/unite_formatted_report.json",
            "data/processed/planits_formatted_report.json"
        ]
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/format_all_databases.py"
