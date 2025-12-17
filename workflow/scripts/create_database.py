#!/usr/bin/env python3
"""
Script to create and initialize the bioinformatics database
workflow/scripts/create_database.py

Initializes the SQLite database with the full schema and configuration.
"""

import json
import time
import logging
import sys
import os
from pathlib import Path
from datetime import datetime

# Import database utilities
sys.path.append(str(Path(__file__).parent / "utils"))
from db_manager import initialize_database

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    """Script entry point."""
    
    try:
        # Read Snakemake parameters
        output_db = snakemake.output.db
        output_report = snakemake.output.report
        config = snakemake.config
        log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

        start_time = time.time()
        logger.info(f"Creating database: {output_db}")

        # Create the parent folder when needed
        Path(output_db).parent.mkdir(parents=True, exist_ok=True)

        # Initialize the database
        db_manager = initialize_database(config)

        # Determine the database type
        db_type = config.get("database", {}).get("type", "sqlite")

        # For PostgreSQL, create the sentinel file
        if db_type == "postgresql":
            Path(output_db).touch()
            logger.info(f"Created PostgreSQL sentinel file: {output_db}")

        # Validate integrity right after creation
        integrity_report = db_manager.validate_integrity()

        # Insert initial metadata if provided in the config
        metadata_inserted = 0
        if 'initial_sources' in config:
            for source_config in config['initial_sources']:
                try:
                    source_id = db_manager.add_source(**source_config)
                    metadata_inserted += 1
                    logger.info(f"Added initial source: {source_config['name']} (ID: {source_id})")
                except Exception as e:
                    logger.warning(f"Failed to add source {source_config['name']}: {e}")

        # Optimize the database right after creation
        db_manager.vacuum_database()

        end_time = time.time()
        creation_time = end_time - start_time

        # Compute database size
        if db_type == "sqlite":
            db_size_mb = os.path.getsize(output_db) / (1024 * 1024)
        else:
            # For PostgreSQL we could query pg_database_size
            db_size_mb = 0  # Placeholder, or query the database for the real size

        # Build the report
        report_data = {
            "timestamp": datetime.now().isoformat(),
            "database": {
                "path": output_db,
                "type": db_type,
                "size_mb": db_size_mb,
                "tables_created": db_manager.list_tables(),
                "schema_version": "1.0"  # Keep your schema versioned
            },
            "execution": {
                "duration_seconds": round(time.time() - start_time, 2),
                "sources_added": metadata_inserted
            },
            "validation": integrity_report  # Integrity report details
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
        
        # Create an error report
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
