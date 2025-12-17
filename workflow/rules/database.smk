# workflow/rules/database.smk
"""
Règles Snakemake pour la gestion de la base de données bioinformatique
Gère l'initialisation, la maintenance et les statistiques (SQLite ou PostgreSQL).
"""

# Règle pour créer/initialiser la base de données
rule create_database:
    """Initialise la base de données (SQLite ou PostgreSQL) avec le schéma complet."""
    output:
        db = get_db_path(),
        report = "data/processed/database_creation_report.json"
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
