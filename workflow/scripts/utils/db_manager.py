#!/usr/bin/env python3
"""
Database Manager for Bioinformatics Pipeline - Simplified Version
workflow/scripts/utils/db_manager.py

Simplified SQLite database management for ITS sequence data pipeline.
Handles schema creation, data insertion, and integrity validation.
"""

import sqlite3
import uuid
from pathlib import Path
from typing import Dict, List, Any
from contextlib import contextmanager
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DatabaseManager:
    """Simplified database manager for ITS sequence data."""
    
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
        """Create simplified database schema with proper constraints and indices."""
        
        schema_sql = """
        -- Sources table: simplified metadata for different data sources
        CREATE TABLE IF NOT EXISTS sources (
            source_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL UNIQUE,
            file_path TEXT,
            date_download TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        
        -- Sequences table: core sequence data with hierarchical relationships
        CREATE TABLE IF NOT EXISTS sequences (
            seq_id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_id INTEGER NOT NULL,
            original_id TEXT NOT NULL,
            sequence TEXT NOT NULL,
            length INTEGER NOT NULL,
            percentage_n REAL DEFAULT 0.0,
            barcode TEXT DEFAULT 'ITS',
            barcode_status TEXT DEFAULT 'Pending',
            parent_seq_id TEXT,
            FOREIGN KEY (source_id) REFERENCES sources(source_id) ON DELETE CASCADE,
            FOREIGN KEY (parent_seq_id) REFERENCES sequences(seq_id) ON DELETE SET NULL,
            CHECK (length > 0),
            CHECK (percentage_n >= 0 AND percentage_n <= 100),
            CHECK (barcode_status IN ('Pending', 'Duplicate', 'ITSx')),
            CHECK (barcode IN ('ITS', 'ITS1', 'ITS2', 'LSU', 'SSU', '5.8S'))
        );
        
        -- Taxonomy table: hierarchical taxonomic classification
        CREATE TABLE IF NOT EXISTS taxonomy (
            tax_id INTEGER PRIMARY KEY AUTOINCREMENT,
            kingdom TEXT,
            phylum TEXT,
            class_ TEXT,
            order_ TEXT,
            family TEXT,
            genus TEXT,
            species TEXT,
            parent_tax_id INTEGER,
            FOREIGN KEY (parent_tax_id) REFERENCES taxonomy(tax_id) ON DELETE SET NULL
        );
        
        -- Sequence-taxonomy junction table (simplified)
        CREATE TABLE IF NOT EXISTS seq_tax (
            seq_id TEXT NOT NULL,
            tax_id INTEGER NOT NULL,
            backbone TEXT DEFAULT 'Raw',
            status TEXT DEFAULT 'Pending',
            PRIMARY KEY (seq_id, tax_id, backbone),
            FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE,
            FOREIGN KEY (tax_id) REFERENCES taxonomy(tax_id) ON DELETE CASCADE,
            CHECK (status IN ('Pending', 'Check', 'Correct', 'Failed'))
            CHECK (backbone IN ('Raw', 'CoL', 'GBIF'))
        );
        """
        
        indices_sql = """
        -- Performance indices
        CREATE INDEX IF NOT EXISTS idx_sequences_source ON sequences(source_id);
        CREATE INDEX IF NOT EXISTS idx_sequences_parent ON sequences(parent_seq_id);
        CREATE INDEX IF NOT EXISTS idx_sequences_barcode ON sequences(barcode);
        CREATE INDEX IF NOT EXISTS idx_sequences_length ON sequences(length);
        
        CREATE INDEX IF NOT EXISTS idx_taxonomy_parent ON taxonomy(parent_tax_id);
        CREATE INDEX IF NOT EXISTS idx_taxonomy_genus_species ON taxonomy(genus, species);
        
        CREATE INDEX IF NOT EXISTS idx_seq_tax_seq ON seq_tax(seq_id);
        CREATE INDEX IF NOT EXISTS idx_seq_tax_tax ON seq_tax(tax_id);
        CREATE INDEX IF NOT EXISTS idx_seq_tax_backbone ON seq_tax(backbone);
        CREATE INDEX IF NOT EXISTS idx_seq_tax_status ON seq_tax(status);
        """
        
        triggers_sql = """
        -- Triggers for data integrity and automation
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
                logger.info("Simplified database schema created successfully")
            except Exception as e:
                logger.error(f"Failed to create schema: {e}")
                raise
    
    def add_source(self, name: str, file_path: str = None) -> int:
        """Add or update a data source."""
        
        with self.get_connection() as conn:
            try:
                cursor = conn.execute("""
                    INSERT OR REPLACE INTO sources 
                    (name, file_path, date_download)
                    VALUES (?, ?, CURRENT_TIMESTAMP)
                """, (name, file_path))
                
                source_id = cursor.lastrowid
                logger.info(f"Added source: {name} (ID: {source_id})")
                return source_id
                
            except Exception as e:
                logger.error(f"Failed to add source {name}: {e}")
                raise
    
    def add_taxonomy(self, kingdom: str = None, phylum: str = None, 
                     class_: str = None, order_: str = None, family: str = None,
                     genus: str = None, species: str = None, 
                     parent_tax_id: int = None) -> int:
        """Add taxonomy record with hierarchical support."""
        
        with self.get_connection() as conn:
            try:
                # Check if taxonomy already exists (based on full lineage)
                cursor = conn.execute("""
                    SELECT tax_id FROM taxonomy 
                    WHERE kingdom = ? AND phylum = ? AND class_ = ? 
                    AND order_ = ? AND family = ? AND genus = ? AND species = ?
                """, (kingdom, phylum, class_, order_, family, genus, species))
                
                existing = cursor.fetchone()
                if existing:
                    return existing['tax_id']
                
                # Insert new taxonomy
                cursor = conn.execute("""
                    INSERT INTO taxonomy 
                    (kingdom, phylum, class_, order_, family, genus, species, parent_tax_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (kingdom, phylum, class_, order_, family, genus, species, parent_tax_id))
                
                tax_id = cursor.lastrowid
                logger.debug(f"Added taxonomy: {genus} {species} (ID: {tax_id})")
                return tax_id
                
            except Exception as e:
                logger.error(f"Failed to add taxonomy: {e}")
                raise
    
    def add_seq_tax(self, seq_id: str, tax_id: int, backbone: str,
                    status: str = 'Pending') -> None:
        """Link sequence to taxonomy with backbone information."""
        
        with self.get_connection() as conn:
            try:
                conn.execute("""
                    INSERT OR REPLACE INTO seq_tax
                    (seq_id, tax_id, backbone, status)
                    VALUES (?, ?, ?, ?)
                """, (seq_id, tax_id, backbone, status))
                
                logger.debug(f"Linked sequence {seq_id} to taxonomy {tax_id} via {backbone}")
                
            except Exception as e:
                logger.error(f"Failed to link sequence to taxonomy: {e}")
                raise
    
    @staticmethod
    def calculate_percentage_n(sequence: str) -> float:
        """Calculate percentage of N characters in sequence."""
        if not sequence:
            return 0.0
        n_count = sequence.upper().count('N')
        return (n_count / len(sequence)) * 100
    
    def add_sequence(self, source_id: int, original_id: str, sequence: str,
                     barcode: str = None, parent_seq_id: str = None, barcode_status: str = 'Pending') -> str:
        """Add sequence with automatic quality metrics calculation."""
        
        # Generate unique sequence ID
        seq_id = str(uuid.uuid4())
        
        # Clean and validate sequence
        clean_sequence = sequence.upper().strip()
        sequence_length = len(clean_sequence.replace('-', ''))
        percentage_n = self.calculate_percentage_n(clean_sequence)
        
        with self.get_connection() as conn:
            try:
                # Insert new sequence
                conn.execute("""
                    INSERT INTO sequences 
                    (seq_id, source_id, original_id, sequence, length, 
                     percentage_n, barcode, barcode_status, parent_seq_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (seq_id, source_id, original_id, clean_sequence, 
                      sequence_length, percentage_n, barcode, barcode_status, parent_seq_id))
                
                logger.debug(f"Added sequence: {original_id} -> {seq_id}")
                return seq_id
                
            except Exception as e:
                logger.error(f"Failed to add sequence {original_id}: {e}")
                raise
    
    def link_sequence_taxonomy(self, seq_id: str, tax_id: int, backbone: str,
                              status: str = 'pending', confidence: float = None) -> None:
        """Link sequence to taxonomy with backbone information."""
        
        with self.get_connection() as conn:
            try:
                conn.execute("""
                    INSERT OR REPLACE INTO seq_tax
                    (seq_id, tax_id, backbone, status, confidence)
                    VALUES (?, ?, ?, ?, ?)
                """, (seq_id, tax_id, backbone, status, confidence))
                
                logger.debug(f"Linked sequence {seq_id} to taxonomy {tax_id} via {backbone}")
                
            except Exception as e:
                logger.error(f"Failed to link sequence to taxonomy: {e}")
                raise
    
    def get_sequence_hierarchy(self, seq_id: str) -> List[Dict]:
        """Get the complete hierarchy for a sequence (parent -> children)."""
        
        with self.get_connection() as conn:
            # Get all sequences in the hierarchy
            cursor = conn.execute("""
                WITH RECURSIVE seq_hierarchy AS (
                    -- Base case: find the root sequence
                    SELECT seq_id, original_id, parent_seq_id, 0 as level
                    FROM sequences 
                    WHERE seq_id = ? OR seq_id IN (
                        SELECT parent_seq_id FROM sequences WHERE seq_id = ?
                    )
                    UNION ALL
                    -- Recursive case: find children
                    SELECT s.seq_id, s.original_id, s.parent_seq_id, h.level + 1
                    FROM sequences s
                    INNER JOIN seq_hierarchy h ON s.parent_seq_id = h.seq_id
                )
                SELECT * FROM seq_hierarchy ORDER BY level, seq_id
            """, (seq_id, seq_id))
            
            return [dict(row) for row in cursor.fetchall()]
    
    def get_taxonomic_lineage(self, tax_id: int) -> List[Dict]:
        """Get complete taxonomic lineage from root to specified taxonomy."""
        
        with self.get_connection() as conn:
            cursor = conn.execute("""
                WITH RECURSIVE tax_lineage AS (
                    -- Base case: start from the specified taxonomy
                    SELECT tax_id, kingdom, phylum, class_, order_, family, genus, species, 
                           parent_tax_id, 0 as level
                    FROM taxonomy 
                    WHERE tax_id = ?
                    UNION ALL
                    -- Recursive case: traverse up the hierarchy
                    SELECT t.tax_id, t.kingdom, t.phylum, t.class_, t.order_, t.family, 
                           t.genus, t.species, t.parent_tax_id, l.level + 1
                    FROM taxonomy t
                    INNER JOIN tax_lineage l ON t.tax_id = l.parent_tax_id
                )
                SELECT * FROM tax_lineage ORDER BY level DESC
            """, (tax_id,))
            
            return [dict(row) for row in cursor.fetchall()]
    
    def get_sequences_by_taxonomy(self, tax_id: int, backbone: str = None) -> List[Dict]:
        """Get all sequences associated with a specific taxonomy."""
        
        query = """
            SELECT s.seq_id, s.original_id, s.length, s.barcode,
                   st.backbone, st.status, st.confidence
            FROM sequences s
            INNER JOIN seq_tax st ON s.seq_id = st.seq_id
            WHERE st.tax_id = ?
        """
        params = [tax_id]
        
        if backbone:
            query += " AND st.backbone = ?"
            params.append(backbone)
            
        query += " ORDER BY s.original_id"
        
        with self.get_connection() as conn:
            cursor = conn.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]
    
    def get_source_stats(self) -> List[Dict]:
        """Get statistics for all data sources."""
        
        with self.get_connection() as conn:
            cursor = conn.execute("""
                SELECT s.name, s.date_download,
                       COUNT(seq.seq_id) as sequence_count,
                       AVG(seq.length) as avg_length,
                       MIN(seq.length) as min_length,
                       MAX(seq.length) as max_length,
                       AVG(seq.percentage_n) as avg_n_content
                FROM sources s
                LEFT JOIN sequences seq ON s.source_id = seq.source_id
                GROUP BY s.source_id, s.name
                ORDER BY s.date_download DESC
            """)
            return [dict(row) for row in cursor.fetchall()]
    
    def get_backbone_stats(self) -> List[Dict]:
        """Get statistics by backbone."""
        
        with self.get_connection() as conn:
            cursor = conn.execute("""
                SELECT st.backbone,
                       COUNT(DISTINCT st.seq_id) as sequence_count,
                       COUNT(DISTINCT st.tax_id) as taxonomy_count,
                       COUNT(CASE WHEN st.status = 'validated' THEN 1 END) as validated_count,
                       AVG(st.confidence) as avg_confidence
                FROM seq_tax st
                GROUP BY st.backbone
                ORDER BY sequence_count DESC
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
            
            # Check for circular references in taxonomy
            cursor = conn.execute("""
                WITH RECURSIVE tax_check AS (
                    SELECT tax_id, parent_tax_id, 1 as depth
                    FROM taxonomy 
                    WHERE parent_tax_id IS NOT NULL
                    UNION ALL
                    SELECT t.tax_id, t.parent_tax_id, tc.depth + 1
                    FROM taxonomy t
                    INNER JOIN tax_check tc ON t.parent_tax_id = tc.tax_id
                    WHERE tc.depth < 20  -- Prevent infinite loops
                )
                SELECT COUNT(*) as count FROM tax_check WHERE depth >= 20
            """)
            checks['circular_taxonomy_refs'] = cursor.fetchone()['count']
            
            # Database size
            cursor = conn.execute("SELECT page_count * page_size as size FROM pragma_page_count(), pragma_page_size()")
            checks['database_size_bytes'] = cursor.fetchone()['size']
            
        checks['is_valid'] = (checks['foreign_key_violations'] == 0 and 
                             checks['orphaned_sequences'] == 0 and
                             checks['circular_taxonomy_refs'] == 0)
        
        return checks

    def list_tables(self) -> List[str]:
        """List all tables in the database."""
        with self.get_connection() as conn:
            cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
            return [row['name'] for row in cursor.fetchall()]
    
    def optimize_for_bulk_insert_advanced(self):
        """Version avancée des optimisations DB"""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                logger.info("Applying advanced bulk insert optimizations...")
                
                # Configuration ultra-agressive pour vitesse maximale
                cursor.execute("PRAGMA foreign_keys = OFF")          # Pas de vérif FK
                cursor.execute("PRAGMA cache_size = -10000000")      # 10GB de cache
                cursor.execute("PRAGMA locking_mode = EXCLUSIVE")    # Verrous exclusifs
                cursor.execute("PRAGMA temp_store = MEMORY")         # Tables temp en RAM
                cursor.execute("PRAGMA mmap_size = 4294967296")      # 4GB de memory mapping
                cursor.execute("PRAGMA page_size = 65536")           # Pages plus grandes
                
                # Désactiver l'autocommit pour contrôler les transactions
                conn.isolation_level = "DEFERRED"
                
                logger.info("Advanced optimizations applied successfully")
                
        except Exception as e:
            logger.warning(f"Could not apply all advanced optimizations: {e}")

    def restore_safe_settings(self):
        """
        Restaure les paramètres de sécurité et de durabilité par défaut de SQLite.
        Cette méthode doit TOUJOURS être appelée après la fin des insertions en masse.
        """
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()

                logger.info("Restoring safe database settings...")

                # Réactiver la vérification des clés étrangères.
                cursor.execute("PRAGMA foreign_keys = ON;")
                
                # Revenir au mode de journalisation par défaut (Write-Ahead Logging),
                # qui est robuste et performant pour un usage normal.
                cursor.execute("PRAGMA journal_mode = WAL;")

                # Rétablir la synchronisation complète pour garantir la durabilité des données.
                cursor.execute("PRAGMA synchronous = FULL;")
                
                # Optionnel : lancer une vérification de l'intégrité des clés étrangères
                # integrity_check = cursor.execute("PRAGMA foreign_key_check;").fetchall()
                # if integrity_check:
                #     logger.error(f"Foreign key integrity errors found: {integrity_check}")

        except Exception as e:
            logger.warning(f"Could not restore safe database settings: {e}")


# Convenience functions for Snakemake integration
def get_database_manager(config: Dict) -> DatabaseManager:
    """Get configured database manager from Snakemake config."""
    db_dir = Path(config.get("database", {}).get("path", "data"))
    db_name = config.get("database", {}).get("name", "its_pipeline.db")
    db_path = db_dir / db_name

    return DatabaseManager(db_path)

def initialize_database(config: Dict) -> DatabaseManager:
    """Initialize database with schema creation."""
    db_manager = get_database_manager(config)
    db_manager.create_schema()
    return db_manager