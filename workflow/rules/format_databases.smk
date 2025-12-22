# workflow/rules/format_databases.smk
import os

def get_db_lock():
    """
    Pick an appropriate database lock level for the selected backend.

    SQLite: return 1 (strict serialization—single writer)
    PostgreSQL: return 999 (practically unlimited thanks to MVCC)
    """
    db_type = config.get("database", {}).get("type", "sqlite")

    if db_type == "sqlite":
        return 1  # Séquentiel obligatoire (SQLite = un seul writer)
    else:  # postgresql
        return 999  # Parallélisme autorisé (PostgreSQL MVCC)

GLOBAL_DB_LOCK = get_db_lock()

def get_caledna_files():
    """Retourne la liste des fichiers CALeDNA depuis la config."""
    return config["databases"]["caledna"]["files"]

def get_caledna_basenames():
    """Retourne les noms de base des fichiers CALeDNA pour les wildcards."""
    files = get_caledna_files()
    return [os.path.splitext(os.path.basename(f))[0] for f in files]


def formatted_database_reports():
    """Chemins des rapports formatés (toutes bases + CALeDNA)."""
    reports = [
        "data/processed/unite_formatted_report.json",
        "data/processed/planits_formatted_report.json",
        "data/processed/bold_formatted_report.json",
    ]
    reports.extend(
        expand(
            "data/processed/caledna_{caledna_file}_formatted_report.json",
            caledna_file=get_caledna_basenames(),
        )
    )
    return reports

# Rule UNITE
rule format_unite_database:
    input:
        source = config["databases"]["unite"]["path"],
        db_ready = rules.create_database.output.report
    output:
        report = "data/processed/unite_formatted_report.json"
    log:
        "logs/format_databases/unite_formatting.log"
    benchmark:
        "benchmarks/format_databases/format_unite_database.benchmark.txt"
    conda:
        "../envs/format_databases.yaml"
    resources: db_lock=GLOBAL_DB_LOCK
    params:
        database_name = "unite",
        database_config = config["databases"]["unite"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_unite.py"

# Rule PLANITS
rule format_planits_database:
    input:
        source = config["databases"]["planits"]["path"],
        db_ready = rules.create_database.output.report
    output:
        report = "data/processed/planits_formatted_report.json"
    log:
        "logs/format_databases/planits_formatting.log"
    benchmark:
        "benchmarks/format_databases/format_planits_database.benchmark.txt"
    conda:
        "../envs/format_databases.yaml"
    resources: db_lock=GLOBAL_DB_LOCK
    params:
        database_name = "planits",
        database_config = config["databases"]["planits"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_planits.py"

# Rule Bold
rule format_bold_database:
    input:
        source = config["databases"]["bold"]["path"],
        db_ready = rules.create_database.output.report
    output:
        report = "data/processed/bold_formatted_report.json"
    log:
        "logs/format_databases/bold_formatting.log"
    benchmark:
        "benchmarks/format_databases/format_bold_database.benchmark.txt"
    conda:
        "../envs/format_databases.yaml"
    resources: db_lock=GLOBAL_DB_LOCK
    params:
        database_name = "bold",
        database_config = config["databases"]["bold"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_bold.py"


rule format_caledna_database:
    input:
        file = lambda wildcards: next(f for f in get_caledna_files() 
                                     if os.path.splitext(os.path.basename(f))[0] == wildcards.caledna_file),
        db_ready = rules.create_database.output.report
    output:
        report = "data/processed/caledna_{caledna_file}_formatted_report.json"
    log:
        "logs/format_databases/caledna_{caledna_file}_formatting.log"
    benchmark:
        "benchmarks/format_databases/format_caledna_database_{caledna_file}.benchmark.txt"
    conda:
        "../envs/format_databases.yaml"
    resources: db_lock=GLOBAL_DB_LOCK
    params:
        database_name = lambda wildcards: f"caledna_{wildcards.caledna_file}",
        database_config = config["databases"]["caledna"],
        sequence_validation = config["sequence_validation"]
    script:
        "../scripts/format_caledna.py"


# Rule de déduplication des taxonomies
rule deduplicate_taxonomies:
    input:
        formatted = formatted_database_reports()

    output:
        report = "data/processed/taxonomy_deduplication_report.json"
    params:
        db = get_db_path(),
    log:
        "logs/format_databases/taxonomy_deduplication.log"
    benchmark:
        "benchmarks/format_databases/deduplicate_taxonomies.benchmark.txt"
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/deduplicate_taxonomies.py"

# Rule finale
rule format_all_databases:
    input:
        deduplication = "data/processed/taxonomy_deduplication_report.json"
    output:
        summary = "data/processed/formatting_summary.json"
    log:
        "logs/format_databases/all_formatting.log"
    benchmark:
        "benchmarks/format_databases/format_all_databases.benchmark.txt"
    params:
        formatted = formatted_database_reports()
    conda:
        "../envs/format_databases.yaml"
    script:
        "../scripts/format_all_databases.py"
