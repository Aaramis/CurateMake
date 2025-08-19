#!/usr/bin/env python3
"""
Script pour créer et initialiser la base de données bioinformatique
workflow/scripts/create_database.py

Initialise la base SQLite avec le schéma complet et les configurations.
"""

import json
import time
import logging
import sys
from pathlib import Path

# Import du gestionnaire de base
sys.path.append(str(Path(__file__).parent / "utils"))
from db_manager import DatabaseManager, initialize_database

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    """Point d'entrée principal."""
    
    try:
        # Récupération des paramètres Snakemake
        output_db = snakemake.output.db
        output_report = snakemake.output.report
        config = snakemake.config
        
        start_time = time.time()
        logger.info(f"Creating database: {output_db}")

        # Créer le dossier si nécessaire
        Path(output_db).parent.mkdir(parents=True, exist_ok=True)

        # Initialiser la base de données
        db_manager = initialize_database(config)
        
        # Vérifier l'intégrité après création
        integrity_report = db_manager.validate_integrity()
        
        # Ajouter des métadonnées initiales si spécifiées dans la config
        metadata_inserted = 0
        if 'initial_sources' in config:
            for source_config in config['initial_sources']:
                try:
                    source_id = db_manager.add_source(**source_config)
                    metadata_inserted += 1
                    logger.info(f"Added initial source: {source_config['name']} (ID: {source_id})")
                except Exception as e:
                    logger.warning(f"Failed to add source {source_config['name']}: {e}")
        
        # Optimiser la base après création
        db_manager.vacuum_database()
        
        end_time = time.time()
        creation_time = end_time - start_time
        
        # Générer le rapport
        report_data = {
            'database_path': str(output_db),
            'creation_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'creation_time_seconds': creation_time,
            'initial_sources_added': metadata_inserted,
            'integrity_check': integrity_report,
            'schema_version': '1.0',
            'status': 'success'
        }
        
        # Sauvegarder le rapport
        with open(output_report, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        logger.info(f"Database created successfully in {creation_time:.2f}s")
        logger.info(f"Database path: {output_db}")
        logger.info(f"Report saved to: {output_report}")
        logger.info(f"Integrity check: {'PASSED' if integrity_report['is_valid'] else 'FAILED'}")
        
    except Exception as e:
        logger.error(f"Database creation failed: {e}")
        
        # Créer un rapport d'erreur
        error_report = {
            'database_path': str(snakemake.output.db),
            'creation_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'status': 'failed',
            'error_message': str(e)
        }
        
        with open(snakemake.output.report, 'w') as f:
            json.dump(error_report, f, indent=2)
        
        raise

if __name__ == "__main__":
    main()