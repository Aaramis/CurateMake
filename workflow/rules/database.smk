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

# # Règle pour optimiser la base de données
# rule optimize_database:
#     """Optimise la base de données (VACUUM, ANALYZE)."""
#     input:
#         db = "data/pipeline.db"
#     output:
#         report = "data/database_optimization_report.json"
#     log:
#         "logs/database/optimization.log"
#     conda:
#         "../envs/database.yaml"
#     script:
#         "../scripts/optimize_database.py"

# # Règle pour générer des statistiques détaillées
# rule database_stats:
#     """Génère des statistiques détaillées sur le contenu de la base."""
#     input:
#         db = "data/pipeline.db"
#     output:
#         stats = "results/database_statistics.json",
#         html_report = "results/database_report.html"
#     log:
#         "logs/database/statistics.log"
#     conda:
#         "../envs/database.yaml"
#     script:
#         "../scripts/generate_database_stats.py"

# # Règle pour sauvegarder la base de données
# rule backup_database:
#     """Crée une sauvegarde de la base de données."""
#     input:
#         db = "data/pipeline.db"
#     output:
#         backup = "backups/pipeline_backup_{timestamp}.db"
#     log:
#         "logs/database/backup_{timestamp}.log"
#     conda:
#         "../envs/database.yaml"
#     script:
#         "../scripts/backup_database.py"

# # Règle pour exporter des données spécifiques
# rule export_sequences:
#     """Exporte des séquences selon des critères spécifiques."""
#     input:
#         db = "data/pipeline.db"
#     output:
#         fasta = "exports/{export_name}.fasta",
#         metadata = "exports/{export_name}_metadata.json"
#     params:
#         filters = lambda wildcards: config.get("exports", {}).get(wildcards.export_name, {})
#     log:
#         "logs/database/export_{export_name}.log"
#     conda:
#         "../envs/database.yaml"
#     script:
#         "../scripts/export_sequences.py"