#!/usr/bin/env python3
"""
Database Manager for Bioinformatics Pipeline
workflow/scripts/utils/db_manager.py

Professional-grade SQLite database management for sequence data pipeline.
Handles schema creation, data insertion, and integrity validation.
"""

import sqlite3
import json
import hashlib
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from contextlib import contextmanager
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DatabaseManager:
    """Professional database manager for bioinformatics sequence data."""
    
    def __init__(self, db_path: str):
        """Initialize database manager with path to SQLite database."""
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
    @contextmanager
    def get_connection(self):
        """Context manager for database connections with proper error handling."""
        conn = None
        try:
            conn = sqlite3.connect(
                self.db_path,
                timeout=30.0,
                isolation_level=None  # Autocommit mode
            )
            conn.execute("PRAGMA foreign_keys = ON")  # Enable FK constraints
            conn.execute("PRAGMA journal_mode = WAL")  # Better concurrency
            conn.row_factory = sqlite3.Row  # Dict-like access
            yield conn
        except Exception as e:
            if conn:
                conn.rollback()
            logger.error(f"Database error: {e}")
            raise
        finally:
            if conn:
                conn.close()
    
    def create_schema(self) -> None:
        """Create complete database schema with proper constraints and indices."""
        
        schema_sql = """
        -- Sources table: metadata for different data sources
        CREATE TABLE IF NOT EXISTS sources (
            source_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL UNIQUE,
            version TEXT NOT NULL,
            url TEXT,
            description TEXT,
            date_download TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            file_path TEXT,
            file_size INTEGER,
            record_count INTEGER,
            UNIQUE(name, version)
        );
        
        -- Sequences table: core sequence data with deduplication
        CREATE TABLE IF NOT EXISTS sequences (
            seq_id TEXT PRIMARY KEY,
            source_id INTEGER NOT NULL,
            original_id TEXT NOT NULL,
            sequence TEXT NOT NULL,
            length INTEGER NOT NULL,
            md5 CHAR(32) NOT NULL UNIQUE,
            gc_content REAL,
            status TEXT NOT NULL DEFAULT 'raw',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (source_id) REFERENCES sources(source_id) ON DELETE CASCADE,
            CHECK (length > 0),
            CHECK (gc_content IS NULL OR (gc_content >= 0 AND gc_content <= 1)),
            CHECK (status IN ('raw', 'itsx', 'aligned', 'augmented', 'filtered'))
        );
        
        -- Taxonomy table: normalized taxonomic hierarchy
        CREATE TABLE IF NOT EXISTS taxonomy (
            tax_id INTEGER PRIMARY KEY AUTOINCREMENT,
            kingdom TEXT,
            phylum TEXT,
            class TEXT,
            "order" TEXT,
            family TEXT,
            genus TEXT,
            species TEXT,
            scientific_name TEXT UNIQUE NOT NULL,
            ncbi_taxid INTEGER,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        
        -- Sequence-taxonomy junction table
        CREATE TABLE IF NOT EXISTS sequence_taxonomy (
            seq_id TEXT NOT NULL,
            tax_id INTEGER NOT NULL,
            method TEXT NOT NULL,
            confidence REAL,
            is_primary BOOLEAN DEFAULT FALSE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (seq_id, tax_id, method),
            FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE,
            FOREIGN KEY (tax_id) REFERENCES taxonomy(tax_id) ON DELETE CASCADE,
            CHECK (confidence IS NULL OR (confidence >= 0 AND confidence <= 1))
        );
        
        -- Transformations log: complete audit trail
        CREATE TABLE IF NOT EXISTS transformations (
            transform_id INTEGER PRIMARY KEY AUTOINCREMENT,
            seq_id TEXT NOT NULL,
            step TEXT NOT NULL,
            tool_name TEXT NOT NULL,
            tool_version TEXT,
            parameters TEXT, -- JSON string
            input_file TEXT,
            output_file TEXT,
            success BOOLEAN NOT NULL,
            error_message TEXT,
            execution_time REAL, -- seconds
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE
        );
        
        -- Datasets for ML model training
        CREATE TABLE IF NOT EXISTS datasets (
            dataset_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL UNIQUE,
            description TEXT,
            dataset_type TEXT DEFAULT 'general',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            created_by TEXT,
            CHECK (dataset_type IN ('train', 'test', 'validation', 'general'))
        );
        
        -- Dataset-sequence relationships
        CREATE TABLE IF NOT EXISTS dataset_sequences (
            dataset_id INTEGER NOT NULL,
            seq_id TEXT NOT NULL,
            split_type TEXT DEFAULT 'general',
            added_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (dataset_id, seq_id),
            FOREIGN KEY (dataset_id) REFERENCES datasets(dataset_id) ON DELETE CASCADE,
            FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE,
            CHECK (split_type IN ('train', 'test', 'validation', 'general'))
        );
        
        -- Quality metrics table
        CREATE TABLE IF NOT EXISTS sequence_quality (
            seq_id TEXT PRIMARY KEY,
            n_content REAL,
            ambiguous_bases INTEGER,
            homopolymer_max INTEGER,
            entropy REAL,
            quality_score REAL,
            is_valid BOOLEAN DEFAULT TRUE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE,
            CHECK (quality_score IS NULL OR (quality_score >= 0 AND quality_score <= 1))
        );
        """
        
        indices_sql = """
        -- Performance indices
        CREATE INDEX IF NOT EXISTS idx_sequences_md5 ON sequences(md5);
        CREATE INDEX IF NOT EXISTS idx_sequences_source ON sequences(source_id);
        CREATE INDEX IF NOT EXISTS idx_sequences_status ON sequences(status);
        CREATE INDEX IF NOT EXISTS idx_sequences_length ON sequences(length);
        
        CREATE INDEX IF NOT EXISTS idx_taxonomy_scientific_name ON taxonomy(scientific_name);
        CREATE INDEX IF NOT EXISTS idx_taxonomy_genus_species ON taxonomy(genus, species);
        
        CREATE INDEX IF NOT EXISTS idx_sequence_taxonomy_seq ON sequence_taxonomy(seq_id);
        CREATE INDEX IF NOT EXISTS idx_sequence_taxonomy_tax ON sequence_taxonomy(tax_id);
        CREATE INDEX IF NOT EXISTS idx_sequence_taxonomy_primary ON sequence_taxonomy(is_primary) WHERE is_primary = TRUE;
        
        CREATE INDEX IF NOT EXISTS idx_transformations_seq ON transformations(seq_id);
        CREATE INDEX IF NOT EXISTS idx_transformations_step ON transformations(step);
        CREATE INDEX IF NOT EXISTS idx_transformations_success ON transformations(success);
        
        CREATE INDEX IF NOT EXISTS idx_dataset_sequences_dataset ON dataset_sequences(dataset_id);
        CREATE INDEX IF NOT EXISTS idx_dataset_sequences_seq ON dataset_sequences(seq_id);
        """
        
        triggers_sql = """
        -- Triggers for data integrity and automation
        CREATE TRIGGER IF NOT EXISTS update_sequence_timestamp 
        AFTER UPDATE ON sequences
        BEGIN
            UPDATE sequences SET updated_at = CURRENT_TIMESTAMP WHERE seq_id = NEW.seq_id;
        END;
        
        CREATE TRIGGER IF NOT EXISTS validate_sequence_insert
        BEFORE INSERT ON sequences
        FOR EACH ROW
        BEGIN
            SELECT CASE
                WHEN LENGTH(NEW.sequence) != NEW.length THEN
                    RAISE(ABORT, 'Sequence length mismatch')
                WHEN NEW.sequence GLOB '*[^ACGTURYSWKMBDHVN-]*' THEN
                    RAISE(ABORT, 'Invalid nucleotide characters in sequence')
            END;
        END;
        """
        
        with self.get_connection() as conn:
            try:
                # Execute schema creation
                conn.executescript(schema_sql)
                conn.executescript(indices_sql)
                conn.executescript(triggers_sql)
                logger.info("Database schema created successfully")
            except Exception as e:
                logger.error(f"Failed to create schema: {e}")
                raise
    
    def add_source(self, name: str, version: str, url: str = None, 
                   description: str = None, file_path: str = None) -> int:
        """Add or update a data source."""
        
        with self.get_connection() as conn:
            try:
                cursor = conn.execute("""
                    INSERT OR REPLACE INTO sources 
                    (name, version, url, description, file_path, date_download)
                    VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
                """, (name, version, url, description, file_path))
                
                source_id = cursor.lastrowid
                logger.info(f"Added source: {name} v{version} (ID: {source_id})")
                return source_id
                
            except Exception as e:
                logger.error(f"Failed to add source {name}: {e}")
                raise
    
    def add_taxonomy(self, scientific_name: str, kingdom: str = None,
                     phylum: str = None, class_: str = None, order: str = None,
                     family: str = None, genus: str = None, species: str = None,
                     ncbi_taxid: int = None) -> int:
        """Add taxonomy record with deduplication."""
        
        logger.info(f"Adding taxonomy: {scientific_name}")

        with self.get_connection() as conn:
            try:
                # Check if taxonomy already exists
                cursor = conn.execute(
                    "SELECT tax_id FROM taxonomy WHERE scientific_name = ?",
                    (scientific_name,)
                )
                existing = cursor.fetchone()
                
                if existing:
                    return existing['tax_id']
                
                # Insert new taxonomy
                cursor = conn.execute("""
                    INSERT INTO taxonomy 
                    (kingdom, phylum, class, "order", family, genus, species, 
                     scientific_name, ncbi_taxid)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (kingdom, phylum, class_, order, family, genus, species,
                      scientific_name, ncbi_taxid))
                
                tax_id = cursor.lastrowid
                logger.debug(f"Added taxonomy: {scientific_name} (ID: {tax_id})")
                return tax_id
                
            except Exception as e:
                logger.error(f"Failed to add taxonomy {scientific_name}: {e}")
                raise
    
    @staticmethod
    def calculate_md5(sequence: str) -> str:
        """Calculate MD5 hash for sequence deduplication."""
        return hashlib.md5(sequence.upper().encode()).hexdigest()
    
    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculate GC content percentage."""
        sequence = sequence.upper().replace('-', '')
        if not sequence:
            return 0.0
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence)
    
    def add_sequence(self, source_id: int, original_id: str, sequence: str,
                     taxonomy_info: Dict[str, Any] = None) -> str:
        """Add sequence with automatic deduplication and quality metrics."""
        
        # Generate unique sequence ID
        seq_id = str(uuid.uuid4())
        
        # Clean and validate sequence
        clean_sequence = sequence.upper().strip()
        sequence_length = len(clean_sequence.replace('-', ''))
        md5_hash = self.calculate_md5(clean_sequence)
        gc_content = self.calculate_gc_content(clean_sequence)
        
        with self.get_connection() as conn:
            try:
                # Check for existing sequence by MD5
                cursor = conn.execute(
                    "SELECT seq_id FROM sequences WHERE md5 = ?",
                    (md5_hash,)
                )
                existing = cursor.fetchone()
                
                if existing:
                    logger.debug(f"Sequence already exists: {existing['seq_id']}")
                    return existing['seq_id']
                
                # Insert new sequence
                conn.execute("""
                    INSERT INTO sequences 
                    (seq_id, source_id, original_id, sequence, length, md5, gc_content)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (seq_id, source_id, original_id, clean_sequence, 
                      sequence_length, md5_hash, gc_content))
                
                # Add taxonomy if provided
                if taxonomy_info:
                    tax_id = self.add_taxonomy(**taxonomy_info)
                    self.link_sequence_taxonomy(seq_id, tax_id, method="original")
                
                logger.debug(f"Added sequence: {original_id} -> {seq_id}")
                return seq_id
                
            except sqlite3.IntegrityError as e:
                if "UNIQUE constraint failed" in str(e):
                    logger.warning(f"Duplicate sequence detected: {original_id}")
                    # Return existing sequence ID
                    cursor = conn.execute(
                        "SELECT seq_id FROM sequences WHERE md5 = ?",
                        (md5_hash,)
                    )
                    return cursor.fetchone()['seq_id']
                raise
            except Exception as e:
                logger.error(f"Failed to add sequence {original_id}: {e}")
                raise
    
    def link_sequence_taxonomy(self, seq_id: str, tax_id: int, method: str,
                              confidence: float = None, is_primary: bool = False) -> None:
        """Link sequence to taxonomy with method and confidence."""
        
        with self.get_connection() as conn:
            try:
                conn.execute("""
                    INSERT OR REPLACE INTO sequence_taxonomy
                    (seq_id, tax_id, method, confidence, is_primary)
                    VALUES (?, ?, ?, ?, ?)
                """, (seq_id, tax_id, method, confidence, is_primary))
                
                logger.debug(f"Linked sequence {seq_id} to taxonomy {tax_id}")
                
            except Exception as e:
                logger.error(f"Failed to link sequence to taxonomy: {e}")
                raise
    
    def log_transformation(self, seq_id: str, step: str, tool_name: str,
                          tool_version: str = None, parameters: Dict = None,
                          success: bool = True, error_message: str = None,
                          execution_time: float = None) -> None:
        """Log a transformation step for audit trail."""
        
        params_json = json.dumps(parameters) if parameters else None
        
        with self.get_connection() as conn:
            try:
                conn.execute("""
                    INSERT INTO transformations
                    (seq_id, step, tool_name, tool_version, parameters,
                     success, error_message, execution_time)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (seq_id, step, tool_name, tool_version, params_json,
                      success, error_message, execution_time))
                
                logger.debug(f"Logged transformation: {step} for {seq_id}")
                
            except Exception as e:
                logger.error(f"Failed to log transformation: {e}")
                raise
    
    def get_source_stats(self) -> List[Dict]:
        """Get statistics for all data sources."""
        
        with self.get_connection() as conn:
            cursor = conn.execute("""
                SELECT s.name, s.version, s.date_download,
                       COUNT(seq.seq_id) as sequence_count,
                       AVG(seq.length) as avg_length,
                       MIN(seq.length) as min_length,
                       MAX(seq.length) as max_length
                FROM sources s
                LEFT JOIN sequences seq ON s.source_id = seq.source_id
                GROUP BY s.source_id, s.name, s.version
                ORDER BY s.date_download DESC
            """)
            return [dict(row) for row in cursor.fetchall()]
    
    def vacuum_database(self) -> None:
        """Optimize database by reclaiming space and updating statistics."""
        
        with self.get_connection() as conn:
            conn.execute("VACUUM")
            conn.execute("ANALYZE")
            logger.info("Database optimized successfully")
    
    def validate_integrity(self) -> Dict[str, Any]:
        """Validate database integrity and return report."""
        
        checks = {}
        
        with self.get_connection() as conn:
            # Check foreign key integrity
            cursor = conn.execute("PRAGMA foreign_key_check")
            fk_violations = cursor.fetchall()
            checks['foreign_key_violations'] = len(fk_violations)
            
            # Check for orphaned sequences
            cursor = conn.execute("""
                SELECT COUNT(*) as count FROM sequences s
                LEFT JOIN sources src ON s.source_id = src.source_id
                WHERE src.source_id IS NULL
            """)
            checks['orphaned_sequences'] = cursor.fetchone()['count']
            
            # Check duplicate sequences by MD5
            cursor = conn.execute("""
                SELECT md5, COUNT(*) as count FROM sequences
                GROUP BY md5 HAVING COUNT(*) > 1
            """)
            checks['duplicate_sequences'] = len(cursor.fetchall())
            
            # Database size
            cursor = conn.execute("SELECT page_count * page_size as size FROM pragma_page_count(), pragma_page_size()")
            checks['database_size_bytes'] = cursor.fetchone()['size']
            
        checks['is_valid'] = (checks['foreign_key_violations'] == 0 and 
                             checks['orphaned_sequences'] == 0)
        
        return checks


# Convenience functions for Snakemake integration
def get_database_manager(config: Dict) -> DatabaseManager:
    """Get configured database manager from Snakemake config."""
    db_dir = Path(config.get("database", {}).get("path", "data"))
    db_name = config.get("database", {}).get("name", "pipeline.db")
    db_path = db_dir / db_name

    return DatabaseManager(db_path)

def initialize_database(config: Dict) -> DatabaseManager:
    """Initialize database with schema creation."""
    db_manager = get_database_manager(config)
    db_manager.create_schema()
    return db_manager