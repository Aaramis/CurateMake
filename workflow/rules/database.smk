# workflow/rules/database.smk
"""
Règles Snakemake pour la gestion de la base de données bioinformatique
Gère l'initialisation, la maintenance et les statistiques de la base SQLite.
"""

# Règle pour créer/initialiser la base de données
rule create_database:
    """Initialise la base de données SQLite avec le schéma complet."""
    output:
        db = os.path.join(config["database"]["path"], config["database"]["name"]),
        report = "data/database_creation_report.json"
    log:
        "logs/database/creation.log"
    conda:
        "../envs/database.yaml"
    params:
        db_config = config["database"]
    script:
        "../scripts/create_database.py"

# Règle pour valider l'intégrité de la base de données
rule validate_database:
    """Valide l'intégrité de la base de données et génère un rapport."""
    input:
        db = "data/pipeline.db"
    output:
        report = "data/database_integrity_report.json"
    log:
        "logs/database/validation.log"
    conda:
        "../envs/database.yaml"
    script:
        "../scripts/validate_database.py"