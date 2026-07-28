# workflow/rules/database.smk
"""
Snakemake rules for managing the bioinformatics database
Handles initialisation, maintenance and statistics (SQLite or PostgreSQL).
"""

# Rule to create/initialise the database
rule create_database:
    """Initialise the database (SQLite or PostgreSQL) with the full schema."""
    output:
        db = get_db_path(),
        report = "results/processed/database_creation_report.json"
    log:
        "logs/database/creation.log"
    benchmark:
        "benchmarks/database/create_database.benchmark.txt"
    conda:
        "../envs/database.yaml"
    params:
        db_config = config["database"]
    script:
        "../scripts/create_database.py"
