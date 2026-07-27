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
        return 1  # Sequential required (SQLite = single writer)
    else:  # postgresql
        return 999  # Parallelism allowed (PostgreSQL MVCC)

GLOBAL_DB_LOCK = get_db_lock()

def has_database(name):
    """Check whether a simple (path) database is configured."""
    db_config = config.get("databases", {}).get(name, {})
    return bool(db_config.get("path"))

def has_caledna_database():
    """Check whether the CALeDNA database is configured (list of files)."""
    db_config = config.get("databases", {}).get("caledna", {})
    return bool(db_config.get("files"))

def get_caledna_files():
    """Return the list of CALeDNA files from the config."""
    return config.get("databases", {}).get("caledna", {}).get("files", [])

def get_caledna_basenames():
    """Return the base names of the CALeDNA files for the wildcards."""
    files = get_caledna_files()
    return [os.path.splitext(os.path.basename(f))[0] for f in files]

def has_personal_database():
    """Check whether the PERSONAL database is configured."""
    personal_config = config.get("databases", {}).get("personal", {})
    return bool(personal_config.get("path") or personal_config.get("files"))

def get_personal_files():
    """
    Retourne la liste des fichiers PERSONAL depuis la config.
    Supporte soit un chemin unique (path) soit une liste de fichiers (files).
    """
    personal_config = config.get("databases", {}).get("personal", {})
    if "files" in personal_config:
        return personal_config["files"]
    elif "path" in personal_config and personal_config["path"]:
        return [personal_config["path"]]
    return []

def get_personal_basenames():
    """Return the base names of the PERSONAL files for the wildcards."""
    files = get_personal_files()
    return [os.path.splitext(os.path.basename(f))[0] for f in files]


def formatted_database_reports():
    """Paths of the formatted reports (configured databases + CALeDNA + PERSONAL)."""
    reports = []
    if has_database("unite"):
        reports.append("data/processed/unite_formatted_report.json")
    if has_database("planits"):
        reports.append("data/processed/planits_formatted_report.json")
    if has_database("bold"):
        reports.append("data/processed/bold_formatted_report.json")
    if has_caledna_database():
        reports.extend(
            expand(
                "data/processed/caledna_{caledna_file}_formatted_report.json",
                caledna_file=get_caledna_basenames(),
            )
        )
    if has_personal_database():
        reports.extend(
            expand(
                "data/processed/personal_{personal_file}_formatted_report.json",
                personal_file=get_personal_basenames(),
            )
        )
    return reports

if has_database("unite"):
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
            sequence_validation = config["sequence_validation"],
            database_formatting = config["database_formatting"]
        script:
            "../scripts/format_unite.py"

if has_database("planits"):
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
            sequence_validation = config["sequence_validation"],
            database_formatting = config["database_formatting"]
        script:
            "../scripts/format_planits.py"

if has_database("bold"):
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
            sequence_validation = config["sequence_validation"],
            database_formatting = config["database_formatting"]
        script:
            "../scripts/format_bold.py"


if has_caledna_database():
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
            sequence_validation = config["sequence_validation"],
            database_formatting = config["database_formatting"]
        script:
            "../scripts/format_caledna.py"


if has_personal_database():
    # Rule PERSONAL - User-provided custom database
    # Expected FASTA header format:
    # >ID;barcode=ITS;kingdom=Fungi;phylum=X;class=X;order=X;family=X;genus=X;species=X
    rule format_personal_database:
        input:
            source = lambda wildcards: next(f for f in get_personal_files()
                                            if os.path.splitext(os.path.basename(f))[0] == wildcards.personal_file),
            db_ready = rules.create_database.output.report
        output:
            report = "data/processed/personal_{personal_file}_formatted_report.json"
        log:
            "logs/format_databases/personal_{personal_file}_formatting.log"
        benchmark:
            "benchmarks/format_databases/format_personal_database_{personal_file}.benchmark.txt"
        conda:
            "../envs/format_databases.yaml"
        resources: db_lock=GLOBAL_DB_LOCK
        params:
            database_name = lambda wildcards: f"personal_{wildcards.personal_file}",
            database_config = lambda wildcards: config["databases"]["personal"],
            sequence_validation = config["sequence_validation"],
            database_formatting = config["database_formatting"]
        script:
            "../scripts/format_personal.py"


# Taxonomy deduplication rule
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
