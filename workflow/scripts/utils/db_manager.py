#!/usr/bin/env python3
"""
Database Manager - SQLAlchemy-based with dual backend support.

This version uses SQLAlchemy ORM for transparent SQLite/PostgreSQL support
while maintaining backward compatibility with the original db_manager.py API.

Key improvements:
- Dual backend support (SQLite for dev/small datasets, PostgreSQL for production)
- ORM-based queries for better maintainability
- Type safety and relationship handling
- Transaction management
- Connection pooling (PostgreSQL)

Usage:
    from db_manager import DatabaseManager

    # SQLite (default)
    db = DatabaseManager(db_type="sqlite", path="results/db/its.db")

    # PostgreSQL
    db = DatabaseManager(
        db_type="postgresql",
        host="localhost",
        user="its_user",
        password="secret",
        database="its_pipeline"
    )
"""

from pathlib import Path
from typing import Dict, List, Any, Optional, Callable
from contextlib import contextmanager
import logging
import uuid
import time

from sqlalchemy import create_engine, text, insert, or_, and_
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.dialects.sqlite import insert as sqlite_insert

# Import ORM models
from db_models import Base, Source, Sequence, Taxonomy, SeqTax
from db_models import get_engine, get_session, init_database

# Import profiler
from db_profiler import DBProfiler, profile_db_operation

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DatabaseManager:
    """
    Database manager with SQLAlchemy ORM and dual backend support.

    Maintains backward compatibility with original db_manager.py API while
    adding support for both SQLite and PostgreSQL backends.
    """

    def __init__(self, db_type: str = "sqlite", **kwargs):
        """
        Initialize database manager.

        Args:
            db_type: Database backend ("sqlite" or "postgresql")
            **kwargs: Backend-specific connection parameters

        SQLite kwargs:
            path (str): Database file path (default: "results/db/its.db")

        PostgreSQL kwargs:
            host (str): Database host
            port (int): Database port (default: 5432)
            user (str): Database user
            password (str): Database password
            database (str): Database name
            pool_size (int): Connection pool size (default: 10)
            max_overflow (int): Max overflow connections (default: 20)
        """
        connection_kwargs = kwargs.copy()
        retry_cfg = connection_kwargs.pop("writer_retry", {}) or {}

        # Extract batch_size from kwargs (default 1000 for backward compatibility)
        self.bulk_insert_batch_size = int(connection_kwargs.pop("batch_size", 1000))

        # Extract taxonomy_fetch_chunk_size for complex taxonomy queries (default 10)
        self.taxonomy_fetch_chunk_size = int(connection_kwargs.pop("taxonomy_fetch_chunk_size", 10))

        self.db_type = db_type
        self.connection_params = connection_kwargs.copy()

        self.writer_retry_attempts = max(0, int(retry_cfg.get("max_attempts", 3)))
        self.writer_retry_enabled = self.writer_retry_attempts > 0

        # Define schema-qualified table names for raw SQL queries
        # This avoids repeating schema logic in every script
        if db_type == "postgresql":
            schema = connection_kwargs.get('schema', 'public')
            self.sequences_table = f"{schema}.sequences"
            self.taxonomy_table = f"{schema}.taxonomy"
            self.seq_tax_table = f"{schema}.seq_tax"
            self.sources_table = f"{schema}.sources"

            # Also set schema on ORM models for ORM queries
            if schema != 'public':
                from db_models import Base
                logger.info(f"Setting schema {schema} for all ORM models")
                for table_name, table in Base.metadata.tables.items():
                    table.schema = schema
        else:  # SQLite
            self.sequences_table = "sequences"
            self.taxonomy_table = "taxonomy"
            self.seq_tax_table = "seq_tax"
            self.sources_table = "sources"

        # Create engine
        self.engine = get_engine(db_type, **connection_kwargs)

        # Create session factory with appropriate configuration
        if db_type == "postgresql":
            # For PostgreSQL, use engine with READ COMMITTED isolation
            # This provides better concurrency for parallel writes
            engine_with_options = self.engine.execution_options(
                isolation_level="READ COMMITTED"
            )
            self.SessionFactory = sessionmaker(bind=engine_with_options)
        else:
            self.SessionFactory = sessionmaker(bind=self.engine)

        # Store db_path for backward compatibility
        if db_type == "sqlite":
            self.db_path = Path(connection_kwargs.get('path', 'results/db/its.db'))
        else:
            self.db_path = None

        logger.info(f"DatabaseManager initialized with {db_type} backend")

    @contextmanager
    def get_connection(self):
        """
        Context manager for raw database connections (backward compatibility).

        Returns a raw connection object (not a Session) for scripts that use
        raw SQL. For SQLite, this is a sqlite3.Connection. For PostgreSQL,
        this is a psycopg2 connection.

        For new code, prefer get_session() which returns a SQLAlchemy Session.
        """
        # Get raw connection from engine
        conn = self.engine.raw_connection()
        try:
            yield conn
            conn.commit()
        except Exception as e:
            conn.rollback()
            logger.error(f"Database error: {e}")
            raise
        finally:
            conn.close()

    @contextmanager
    def get_session(self) -> Session:
        """
        Context manager for SQLAlchemy sessions (new API).

        Preferred over get_connection() for new code.
        """
        session = self.SessionFactory()
        try:
            yield session
            session.commit()
        except Exception as e:
            session.rollback()
            logger.error(f"Session error: {e}")
            raise
        finally:
            session.close()

    def _execute_with_retry(self, operation_name: str, func: Callable[[], Any]):
        """Execute a DB operation with retry logic for SQLAlchemy errors."""
        attempts = 0

        while True:
            try:
                return func()
            except SQLAlchemyError as exc:
                attempts += 1
                should_retry = self.writer_retry_enabled

                if not should_retry or attempts > self.writer_retry_attempts:
                    logger.error(
                        "DB operation %s failed after %s attempt(s): %s",
                        operation_name,
                        attempts,
                        exc
                    )
                    raise

                logger.info(
                    "DB operation %s: retry %s/%s...",
                    operation_name,
                    attempts,
                    self.writer_retry_attempts
                )

    @profile_db_operation("create_schema")
    def create_schema(self) -> None:
        """Create database schema using SQLAlchemy models."""
        try:
            # Pass schema for PostgreSQL if provided
            schema = self.connection_params.get('schema') if self.db_type == 'postgresql' else None
            init_database(self.engine, schema=schema)
            logger.info("Database schema created successfully")
        except Exception as e:
            logger.error(f"Failed to create schema: {e}")
            raise

    @profile_db_operation("add_source")
    def add_source(self, name: str, file_path: str = None) -> str:
        """
        Add or update a data source.

        Args:
            name: Source name (must be unique)
            file_path: Optional file path

        Returns:
            source_id: Primary key of the source (UUID string)
        """
        with self.get_session() as session:
            try:
                # Check if source already exists
                existing = session.query(Source).filter_by(name=name).first()

                if existing:
                    # Update existing source
                    existing.file_path = file_path
                    session.flush()
                    source_id = existing.source_id
                    logger.info(f"Updated source: {name} (ID: {source_id})")
                else:
                    # Create new source
                    source = Source(name=name, file_path=file_path)
                    session.add(source)
                    session.flush()
                    source_id = source.source_id
                    logger.info(f"Added source: {name} (ID: {source_id})")

                return source_id

            except Exception as e:
                logger.error(f"Failed to add source {name}: {e}")
                raise

    @profile_db_operation("add_taxonomy")
    def add_taxonomy(self, kingdom: str = None, phylum: str = None,
                     class_: str = None, order_: str = None, family: str = None,
                     genus: str = None, species: str = None,
                     parent_tax_id: str = None,
                     col_id: str = None, col_dataset_key: str = None) -> str:
        """
        Add taxonomy record with hierarchical support.

        Args:
            kingdom, phylum, class_, order_, family, genus, species: Taxonomic ranks
            parent_tax_id: Optional parent taxonomy ID (UUID string)
            col_id, col_dataset_key: Optional CoL provenance (accepted usage id + pinned
                release key). When the lineage already exists, they are stamped onto the
                existing row only if it has no col_id yet.

        Returns:
            tax_id: Primary key of the taxonomy (UUID string)
        """
        def _stamp(row):
            if col_id and not row.col_id:
                row.col_id = col_id
                row.col_dataset_key = col_dataset_key

        lineage = dict(
            kingdom=kingdom,
            phylum=phylum,
            class_=class_,
            order_=order_,
            family=family,
            genus=genus,
            species=species,
        )

        with self.get_session() as session:
            try:
                # Check if taxonomy already exists (based on full lineage)
                existing = session.query(Taxonomy).filter_by(**lineage).first()

                if existing:
                    _stamp(existing)
                    return existing.tax_id

                # Insert new taxonomy
                taxonomy = Taxonomy(parent_tax_id=parent_tax_id, col_id=col_id,
                                    col_dataset_key=col_dataset_key, **lineage)
                session.add(taxonomy)
                try:
                    session.flush()
                except IntegrityError:
                    # Another writer inserted the same lineage between our existence
                    # check and this flush (concurrent workers can request the same
                    # lineage). Roll back and return the row the other writer created.
                    # Portable to SQLite and PostgreSQL/YugabyteDB (both raise
                    # IntegrityError on the UNIQUE lineage constraint).
                    session.rollback()
                    existing = session.query(Taxonomy).filter_by(**lineage).first()
                    if existing:
                        logger.debug(
                            "add_taxonomy: lost insert race for lineage %s, reusing %s",
                            lineage, existing.tax_id
                        )
                        _stamp(existing)
                        return existing.tax_id
                    # Constraint violation was not the concurrent-duplicate case.
                    raise
                tax_id = taxonomy.tax_id

                logger.debug(f"Added taxonomy: {genus or ''} {species or ''} (ID: {tax_id})")
                return tax_id

            except Exception as e:
                logger.error(f"Failed to add taxonomy: {e}")
                raise

    def add_seq_tax(self, seq_id: str, tax_id: str, backbone: str,
                    status: str = 'Pending', confidence: float = None) -> None:
        """
        Link sequence to taxonomy with backbone information.

        Args:
            seq_id: Sequence ID (UUID string)
            tax_id: Taxonomy ID (UUID string)
            backbone: Data source ('Raw', 'CoL', 'GBIF', 'Sativa')
            status: Validation status
            confidence: Optional confidence score
        """
        with self.get_session() as session:
            try:
                # Check if link already exists
                existing = session.query(SeqTax).filter_by(
                    seq_id=seq_id,
                    tax_id=tax_id,
                    backbone=backbone
                ).first()

                if existing:
                    # Update existing link
                    existing.status = status
                    existing.confidence = confidence
                else:
                    # Create new link
                    seq_tax = SeqTax(
                        seq_id=seq_id,
                        tax_id=tax_id,
                        backbone=backbone,
                        status=status,
                        confidence=confidence
                    )
                    session.add(seq_tax)

                logger.debug(f"Linked sequence {seq_id} to taxonomy {tax_id} via {backbone}")

            except Exception as e:
                logger.error(f"Failed to link sequence to taxonomy: {e}")
                raise

    def link_sequence_taxonomy(self, seq_id: str, tax_id: str, backbone: str,
                              status: str = 'pending', confidence: float = None) -> None:
        """Alias for add_seq_tax() for backward compatibility."""
        self.add_seq_tax(seq_id, tax_id, backbone, status, confidence)

    @profile_db_operation("get_sequence_hierarchy")
    def get_sequence_hierarchy(self, seq_id: str) -> List[Dict]:
        """
        Get the complete hierarchy for a sequence (parent -> children).

        Args:
            seq_id: Sequence ID (UUID string)

        Returns:
            List of dictionaries with sequence hierarchy
        """
        with self.get_session() as session:
            # Use raw SQL for recursive CTE (SQLAlchemy doesn't have great CTE support)
            query = text("""
                WITH RECURSIVE seq_hierarchy AS (
                    -- Base case: find the root sequence
                    SELECT seq_id, original_id, parent_seq_id, 0 as level
                    FROM sequences
                    WHERE seq_id = :seq_id OR seq_id IN (
                        SELECT parent_seq_id FROM sequences WHERE seq_id = :seq_id
                    )
                    UNION ALL
                    -- Recursive case: find children
                    SELECT s.seq_id, s.original_id, s.parent_seq_id, h.level + 1
                    FROM sequences s
                    INNER JOIN seq_hierarchy h ON s.parent_seq_id = h.seq_id
                )
                SELECT * FROM seq_hierarchy ORDER BY level, seq_id
            """)

            result = session.execute(query, {"seq_id": seq_id})
            return [dict(row._mapping) for row in result]

    @profile_db_operation("get_taxonomic_lineage")
    def get_taxonomic_lineage(self, tax_id: str) -> List[Dict]:
        """
        Get complete taxonomic lineage from root to specified taxonomy.

        Args:
            tax_id: Taxonomy ID (UUID string)

        Returns:
            List of dictionaries with taxonomic lineage
        """
        with self.get_session() as session:
            query = text("""
                WITH RECURSIVE tax_lineage AS (
                    -- Base case: start from the specified taxonomy
                    SELECT tax_id, kingdom, phylum, class_, order_, family, genus, species,
                           parent_tax_id, 0 as level
                    FROM taxonomy
                    WHERE tax_id = :tax_id
                    UNION ALL
                    -- Recursive case: traverse up the hierarchy
                    SELECT t.tax_id, t.kingdom, t.phylum, t.class_, t.order_, t.family,
                           t.genus, t.species, t.parent_tax_id, l.level + 1
                    FROM taxonomy t
                    INNER JOIN tax_lineage l ON t.tax_id = l.parent_tax_id
                )
                SELECT * FROM tax_lineage ORDER BY level DESC
            """)

            result = session.execute(query, {"tax_id": tax_id})
            return [dict(row._mapping) for row in result]

    @profile_db_operation("get_sequences_by_taxonomy")
    def get_sequences_by_taxonomy(self, tax_id: str, backbone: str = None) -> List[Dict]:
        """
        Get all sequences associated with a specific taxonomy.

        Args:
            tax_id: Taxonomy ID (UUID string)
            backbone: Optional backbone filter

        Returns:
            List of dictionaries with sequence information
        """
        with self.get_session() as session:
            query = session.query(
                Sequence.seq_id,
                Sequence.original_id,
                Sequence.length,
                Sequence.barcode,
                SeqTax.backbone,
                SeqTax.status,
                SeqTax.confidence
            ).join(SeqTax).filter(SeqTax.tax_id == tax_id)

            if backbone:
                query = query.filter(SeqTax.backbone == backbone)

            query = query.order_by(Sequence.original_id)

            return [
                {
                    'seq_id': row.seq_id,
                    'original_id': row.original_id,
                    'length': row.length,
                    'barcode': row.barcode,
                    'backbone': row.backbone,
                    'status': row.status,
                    'confidence': row.confidence
                }
                for row in query.all()
            ]

    @profile_db_operation("get_source_stats")
    def get_source_stats(self) -> List[Dict]:
        """Get statistics for all data sources."""
        with self.get_session() as session:
            query = text("""
                SELECT s.name, s.date_download,
                       COUNT(seq.seq_id) as sequence_count,
                       AVG(seq.length) as avg_length,
                       MIN(seq.length) as min_length,
                       MAX(seq.length) as max_length
                FROM sources s
                LEFT JOIN sequences seq ON s.source_id = seq.source_id
                GROUP BY s.source_id, s.name
                ORDER BY s.date_download DESC
            """)

            result = session.execute(query)
            return [dict(row._mapping) for row in result]

    @profile_db_operation("get_backbone_stats")
    def get_backbone_stats(self) -> List[Dict]:
        """Get statistics by backbone."""
        with self.get_session() as session:
            query = text("""
                SELECT st.backbone,
                       COUNT(DISTINCT st.seq_id) as sequence_count,
                       COUNT(DISTINCT st.tax_id) as taxonomy_count,
                       COUNT(CASE WHEN st.status = 'validated' THEN 1 END) as validated_count,
                       AVG(st.confidence) as avg_confidence
                FROM seq_tax st
                GROUP BY st.backbone
                ORDER BY sequence_count DESC
            """)

            result = session.execute(query)
            return [dict(row._mapping) for row in result]

    # ========== SATIVA Integration Methods ==========

    def get_seq_id_from_original_id(self, original_id: str) -> str:
        """Get seq_id from original_id (returns UUID string)."""
        with self.get_session() as session:
            sequence = session.query(Sequence).filter_by(original_id=original_id).first()

            if not sequence:
                raise ValueError(f"Sequence not found for original_id: {original_id}")

            return sequence.seq_id

    @profile_db_operation("get_priority_taxonomy_id")
    def get_priority_taxonomy_id(self, seq_id: str) -> str:
        """
        Get the priority taxonomy ID for a sequence following the logic:
        1. CoL backbone if exists
        2. Otherwise, prioritize by status: Correct > Check > others
        3. If no taxonomy found, raise an error

        Args:
            seq_id: Sequence ID (UUID string)

        Returns:
            tax_id: Taxonomy ID (UUID string)
        """
        with self.get_session() as session:
            # 1. Try CoL first
            col_tax = session.query(SeqTax).filter_by(
                seq_id=seq_id,
                backbone='CoL'
            ).first()

            if col_tax:
                logger.debug(f"Found CoL taxonomy for seq_id {seq_id}: tax_id={col_tax.tax_id}")
                return col_tax.tax_id

            # 2. If no CoL, search in other backbones with status priority
            query = text("""
                SELECT tax_id, status FROM seq_tax
                WHERE seq_id = :seq_id AND backbone IN ('Raw', 'GBIF')
                ORDER BY
                    CASE status
                        WHEN 'Correct' THEN 1
                        WHEN 'Check' THEN 2
                        ELSE 3
                    END
                LIMIT 1
            """)

            result = session.execute(query, {"seq_id": seq_id}).first()

            if result:
                logger.debug(f"Found {result.status} taxonomy for seq_id {seq_id}: tax_id={result.tax_id}")
                return result.tax_id

            # 3. No taxonomy found
            raise ValueError(f"No taxonomy found for seq_id: {seq_id}")

    def check_existing_taxonomy(self, kingdom: str = None, phylum: str = None,
                               class_: str = None, order_: str = None, family: str = None,
                               genus: str = None, species: str = None) -> Optional[str]:
        """Check if a taxonomy already exists in the database (returns UUID string)."""
        with self.get_session() as session:
            taxonomy = session.query(Taxonomy).filter_by(
                kingdom=kingdom,
                phylum=phylum,
                class_=class_,
                order_=order_,
                family=family,
                genus=genus,
                species=species
            ).first()

            if taxonomy:
                logger.debug(f"Found existing taxonomy: tax_id={taxonomy.tax_id}")
                return taxonomy.tax_id

            return None

    def add_seq_tax_sativa(self, seq_id: str, tax_id: str) -> None:
        """Link sequence to SATIVA-corrected taxonomy."""
        self.add_seq_tax(seq_id, tax_id, backbone='Sativa', status=None)
        logger.debug(f"Added Sativa link: seq_id={seq_id} -> tax_id={tax_id}")

    @profile_db_operation("bulk_insert_taxonomies")
    def bulk_insert_taxonomies(self, taxonomies: List[Dict[str, Any]]) -> Dict[tuple, str]:
        """
        Insert taxonomies in bulk and return mapping of taxonomy tuples to tax_ids.

        Args:
            taxonomies: List of dicts with keys: kingdom, phylum, class_, order_, family, genus, species

        Returns:
            Dict mapping (kingdom, phylum, class_, order_, family, genus, species) -> tax_id (UUID string)
        """
        if not taxonomies:
            return {}

        # Deduplicate: keep only first occurrence of each unique taxonomy
        seen = {}
        unique_taxonomies = []
        for tax in taxonomies:
            key = (tax.get('kingdom'), tax.get('phylum'), tax.get('class_'),
                   tax.get('order_'), tax.get('family'), tax.get('genus'), tax.get('species'))
            if key not in seen:
                if 'tax_id' not in tax:
                    tax['tax_id'] = str(uuid.uuid4())
                seen[key] = tax['tax_id']
                unique_taxonomies.append(tax)

        if len(unique_taxonomies) < len(taxonomies):
            logger.debug(f"Deduplicated {len(taxonomies)} -> {len(unique_taxonomies)} taxonomies within batch")

        # Phase 1: Insert taxonomies using UPSERT (PostgreSQL) or try/except (SQLite)
        def insert_phase():
            with self.get_session() as session:
                if self.db_type == "postgresql":
                    stmt = pg_insert(Taxonomy).values(unique_taxonomies)
                    stmt = stmt.on_conflict_do_nothing(
                        index_elements=['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']
                    )
                    session.execute(stmt)
                    logger.debug(f"PostgreSQL UPSERT: inserted {len(unique_taxonomies)} taxonomies (duplicates ignored)")
                else:
                    stmt = sqlite_insert(Taxonomy).values(unique_taxonomies)
                    stmt = stmt.on_conflict_do_nothing(
                        index_elements=['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']
                    )
                    session.execute(stmt)
                    logger.debug(f"SQLite UPSERT: inserted {len(unique_taxonomies)} taxonomies (duplicates ignored)")

        self._execute_with_retry("bulk_insert_taxonomies_insert", insert_phase)

        # Phase 2: resolve the authoritative tax_id for every unique taxonomy.
        #
        # tax_ids were pre-generated BEFORE the UPSERT, so any taxonomy we actually inserted
        # keeps our pre-generated id. We therefore:
        #   (2a) fetch by PRIMARY KEY (indexed, cheap, large chunks) to confirm the ones we
        #        inserted — the vast majority on a fresh ingestion;
        #   (2b) fall back to the content lookup ONLY for the few that already existed
        #        (ON CONFLICT DO NOTHING kept another row).
        # This keeps the number of queries low, which matters most on YugabyteDB (one
        # round-trip per query).
        unique_taxonomy_tuples = {
            (t.get('kingdom'), t.get('phylum'), t.get('class_'),
             t.get('order_'), t.get('family'), t.get('genus'), t.get('species'))
            for t in taxonomies
        }

        tax_map = {}

        # (2a) Confirm inserted taxonomies via their pre-generated tax_id (PK lookup). Covers:
        #   - every genuinely new taxonomy (our insert succeeded -> pregen id present);
        #   - every taxonomy with a NULL field: the UNIQUE constraint treats NULL as distinct
        #     (NULL != NULL), so ON CONFLICT does not fire and the UPSERT inserts a row with
        #     our pregen id (a duplicate when the taxonomy already existed). These duplicates
        #     are merged later by `deduplicate_taxonomies`, which remaps seq_tax to the keeper.
        # Only fully-set keys that ALREADY existed fall through to (2b): ON CONFLICT skipped
        # our insert, so the pregen id is absent and we look up the existing id by content.
        pregen_to_key = {
            t['tax_id']: (t.get('kingdom'), t.get('phylum'), t.get('class_'),
                          t.get('order_'), t.get('family'), t.get('genus'), t.get('species'))
            for t in unique_taxonomies
        }
        pregen_ids = list(pregen_to_key.keys())
        pk_chunk = self.bulk_insert_batch_size  # same chunking as bulk_insert_sequences
        for i in range(0, len(pregen_ids), pk_chunk):
            chunk = pregen_ids[i:i + pk_chunk]

            def fetch_pk(chunk=chunk):
                with self.get_session() as session:
                    return [row[0] for row in session.query(Taxonomy.tax_id).filter(
                        Taxonomy.tax_id.in_(chunk)).all()]

            for tid in self._execute_with_retry(f"bulk_insert_taxonomies_pk_{i}", fetch_pk):
                key = pregen_to_key.get(tid)
                if key is not None:
                    tax_map[key] = tid

        # (2b) Only pre-existing taxonomies (insert conflicts) remain -> content lookup.
        missing = [k for k in unique_taxonomy_tuples if k not in tax_map]
        chunk_size = self.taxonomy_fetch_chunk_size  # small: complex OR-of-ANDs predicate
        for i in range(0, len(missing), chunk_size):
            chunk = missing[i:i + chunk_size]

            conditions = []
            for tax_tuple in chunk:
                kingdom, phylum, class_, order_, family, genus, species = tax_tuple
                conditions.append(and_(
                    (Taxonomy.kingdom == kingdom) if kingdom is not None else Taxonomy.kingdom.is_(None),
                    (Taxonomy.phylum == phylum) if phylum is not None else Taxonomy.phylum.is_(None),
                    (Taxonomy.class_ == class_) if class_ is not None else Taxonomy.class_.is_(None),
                    (Taxonomy.order_ == order_) if order_ is not None else Taxonomy.order_.is_(None),
                    (Taxonomy.family == family) if family is not None else Taxonomy.family.is_(None),
                    (Taxonomy.genus == genus) if genus is not None else Taxonomy.genus.is_(None),
                    (Taxonomy.species == species) if species is not None else Taxonomy.species.is_(None)
                ))

            def fetch_chunk(conditions=conditions):
                with self.get_session() as session:
                    return session.query(
                        Taxonomy.kingdom, Taxonomy.phylum, Taxonomy.class_,
                        Taxonomy.order_, Taxonomy.family, Taxonomy.genus,
                        Taxonomy.species, Taxonomy.tax_id
                    ).filter(or_(*conditions)).all()

            for t in self._execute_with_retry(f"bulk_insert_taxonomies_conflict_{i}", fetch_chunk):
                tax_map[(t.kingdom, t.phylum, t.class_, t.order_, t.family, t.genus, t.species)] = t.tax_id

        if len(tax_map) < len(unique_taxonomy_tuples):
            n_missing = len(unique_taxonomy_tuples) - len(tax_map)
            logger.warning(f"Could not find {n_missing}/{len(unique_taxonomy_tuples)} taxonomies in database")

        return tax_map

    @profile_db_operation("bulk_insert_sequences")
    def bulk_insert_sequences(self, sequences: List[Dict[str, Any]]) -> Dict[str, str]:
        """
        Insert sequences in bulk and return mapping of original_ids to seq_ids.

        Args:
            sequences: List of dicts with keys: source_id, original_id, sequence, length, barcode, barcode_status

        Returns:
            Dict mapping original_id -> seq_id (UUID string)
        """
        if not sequences:
            return {}

        def insert_phase():
            with self.get_session() as session:
                session.bulk_insert_mappings(Sequence, sequences)

        self._execute_with_retry("bulk_insert_sequences_insert", insert_phase)

        # Fetch phase: get seq_ids for all original_ids
        source_id = sequences[0]['source_id']
        original_ids = [s['original_id'] for s in sequences if s.get('original_id')]

        if not original_ids:
            return {}

        # De-duplicate while preserving order to reduce IN list size.
        seen = set()
        unique_ids = []
        for original_id in original_ids:
            if original_id in seen:
                continue
            seen.add(original_id)
            unique_ids.append(original_id)

        # Chunk IN queries to avoid large RPC/statement timeouts on big batches.
        # Use batch_size from config (same as bulk insert operations)
        fetch_batch_size = self.bulk_insert_batch_size
        seq_id_map: Dict[str, str] = {}

        for i in range(0, len(unique_ids), fetch_batch_size):
            chunk = unique_ids[i:i + fetch_batch_size]

            # Execute each chunk query with retry
            def fetch_chunk():
                with self.get_session() as session:
                    results = session.query(Sequence.original_id, Sequence.seq_id).filter(
                        Sequence.source_id == source_id,
                        Sequence.original_id.in_(chunk)
                    ).all()
                    return results

            results = self._execute_with_retry(
                f"bulk_insert_sequences_fetch_chunk_{i}",
                fetch_chunk
            )
            seq_id_map.update({original_id: seq_id for original_id, seq_id in results})

        return seq_id_map

    @profile_db_operation("bulk_insert_seq_tax")
    def bulk_insert_seq_tax(self, relations: List[Dict[str, Any]]) -> None:
        """
        Safe bulk insert for seq_tax with automatic batching.
        Uses PostgreSQL ON CONFLICT DO NOTHING to avoid PK violations in parallel workloads.
        Falls back to try/except for SQLite.

        Batches large inserts to avoid parameter limits and YugabyteDB RPC timeouts.
        Batch size is configurable via database_formatting.batch_size in config.
        """
        if not relations:
            return

        # Use batch_size from config (default 1000 if not specified)
        batch_size = self.bulk_insert_batch_size

        for i in range(0, len(relations), batch_size):
            batch = relations[i:i + batch_size]

            def insert_batch():
                with self.get_session() as session:
                    if self.db_type == "postgresql":
                        stmt = pg_insert(SeqTax).values(batch)
                        stmt = stmt.on_conflict_do_nothing(
                            index_elements=['seq_id', 'tax_id', 'backbone']
                        )
                        session.execute(stmt)

                    else:
                        stmt = sqlite_insert(SeqTax).values(batch)
                        stmt = stmt.on_conflict_do_nothing(
                            index_elements=['seq_id', 'tax_id', 'backbone']
                        )
                        session.execute(stmt)

            self._execute_with_retry("bulk_insert_seq_tax_insert", insert_batch)

            if (i // batch_size) % 10 == 0 and i > 0:
                logger.debug(f"bulk_insert_seq_tax: {i}/{len(relations)} rows processed")


    def get_or_create_taxonomy(self, taxonomy_dict: Dict[str, Any], original_tax_id: str = None) -> str:
        """
        Return tax_id if the taxonomy already exists, otherwise create a new one.

        Args:
            taxonomy_dict: Dictionary with taxonomic ranks
            original_tax_id: Optional original taxonomy ID (UUID string)

        Returns:
            tax_id: Taxonomy ID (UUID string)

        Example:
            taxonomy_dict = {
                "kingdom": "Fungi",
                "phylum": "Ascomycota",
                ...
            }
        """

        cleaned = {k: (v if v not in ["", None, "None", "null", "NA"] else None)
                for k, v in taxonomy_dict.items()}

        with self.get_session() as session:
            query = session.query(Taxonomy)

            # Build the WHERE clause properly
            for field, value in cleaned.items():
                col = getattr(Taxonomy, field)

                if value is None:
                    query = query.filter(col.is_(None))
                else:
                    query = query.filter(col == value)

            existing = query.first()

            if existing:
                return existing.tax_id

            # Create a brand-new taxonomy entry
            new_tax = Taxonomy(**cleaned)

            session.add(new_tax)
            session.flush()  # obtient tax_id
            tax_id = new_tax.tax_id

            # Historisation (optionnel)
            if original_tax_id:
                logger.debug(f"Created new taxonomy {tax_id} from original {original_tax_id}")

            session.commit()
            return tax_id


    @profile_db_operation("vacuum_database")
    def vacuum_database(self) -> None:
        """
        VACUUM is automatically handled by PostgreSQL autovacuum.
        Manual VACUUM during parallel jobs causes Catalog Version Mismatch
        and MUST be skipped for safety.

        For SQLite, VACUUM is safe and useful.
        """
        if self.db_type == "sqlite":
            with self.get_session() as session:
                session.execute(text("VACUUM"))
                session.execute(text("ANALYZE"))
                logger.info("SQLite database vacuumed and analyzed.")
        else:
            logger.info("Skipping manual VACUUM on PostgreSQL (autovacuum handles it).")


    @profile_db_operation("validate_integrity")
    def validate_integrity(self) -> Dict[str, Any]:
        """Validate database integrity and return report."""
        checks = {}

        with self.get_session() as session:
            # Check foreign key integrity
            if self.db_type == "sqlite":
                result = session.execute(text("PRAGMA foreign_key_check"))
                fk_violations = result.fetchall()
                checks['foreign_key_violations'] = len(fk_violations)
            else:
                # PostgreSQL has built-in FK checking
                checks['foreign_key_violations'] = 0

            # Check for orphaned sequences
            query = text("""
                SELECT COUNT(*) as count FROM sequences s
                LEFT JOIN sources src ON s.source_id = src.source_id
                WHERE src.source_id IS NULL
            """)
            result = session.execute(query).first()
            checks['orphaned_sequences'] = result.count

            # Check for circular references in taxonomy
            query = text("""
                WITH RECURSIVE tax_check AS (
                    SELECT tax_id, parent_tax_id, 1 as depth
                    FROM taxonomy
                    WHERE parent_tax_id IS NOT NULL
                    UNION ALL
                    SELECT t.tax_id, t.parent_tax_id, tc.depth + 1
                    FROM taxonomy t
                    INNER JOIN tax_check tc ON t.parent_tax_id = tc.tax_id
                    WHERE tc.depth < 20
                )
                SELECT COUNT(*) as count FROM tax_check WHERE depth >= 20
            """)
            result = session.execute(query).first()
            checks['circular_taxonomy_refs'] = result.count

            # Database size
            if self.db_type == "sqlite":
                query = text("SELECT page_count * page_size as size FROM pragma_page_count(), pragma_page_size()")
                result = session.execute(query).first()
                checks['database_size_bytes'] = result.size
            else:
                query = text("SELECT pg_database_size(current_database()) as size")
                result = session.execute(query).first()
                checks['database_size_bytes'] = result.size

        checks['is_valid'] = (checks['foreign_key_violations'] == 0 and
                             checks['orphaned_sequences'] == 0 and
                             checks['circular_taxonomy_refs'] == 0)

        return checks

    def list_tables(self) -> List[str]:
        """List all tables in the database."""
        return [table.name for table in Base.metadata.sorted_tables]

    # ========== Bulk Operations ==========

    def optimize_for_bulk_insert_advanced(self):
        """Apply optimizations for bulk insert operations."""
        logger.warning("bulk insert optimizations not implemented in V2 - use SQLAlchemy bulk_insert_mappings()")

    def restore_safe_settings(self):
        """Restore safe database settings after bulk operations."""
        logger.info("Safe settings are default in SQLAlchemy")

    # ========== Profiling Methods ==========

    @staticmethod
    def enable_profiling(enabled: bool = True):
        """Enable or disable database profiling."""
        DBProfiler.enable(enabled)

    @staticmethod
    def set_slow_query_threshold(threshold: float):
        """Set threshold in seconds for slow query detection."""
        DBProfiler.set_slow_query_threshold(threshold)

    @staticmethod
    def get_profiling_stats(operation_name: str = None) -> Dict[str, Any]:
        """Get profiling statistics."""
        return DBProfiler.get_statistics(operation_name)

    @staticmethod
    def get_slow_queries(limit: int = 10) -> List[Dict[str, Any]]:
        """Get the slowest database operations."""
        return DBProfiler.get_slow_queries(limit)

    @staticmethod
    def generate_profiling_report(output_path: str = None) -> Dict[str, Any]:
        """Generate comprehensive profiling report."""
        return DBProfiler.generate_report(output_path)

    def get_database_size_metrics(self) -> Dict[str, Any]:
        """Get detailed database size metrics."""
        if self.db_path:
            return DBProfiler.get_db_size(str(self.db_path))
        else:
            return {'backend': 'postgresql', 'note': 'Use PostgreSQL tools for size metrics'}

    @staticmethod
    def reset_profiling():
        """Clear all profiling data."""
        DBProfiler.reset()


# ========== Convenience Functions ==========

def get_database_manager(config: Dict) -> DatabaseManager:
    """
    Get configured database manager from config dictionary.

    Args:
        config: Configuration dictionary with 'database' section

    Returns:
        DatabaseManager instance

    Example:
        config = {
            'database': {
                'type': 'sqlite',
                'sqlite': {'path': 'results/db', 'name': 'its.db'}
            },
            'database_formatting': {
                'batch_size': 5000
            }
        }
        db = get_database_manager(config)
    """
    db_config = config.get("database", {})
    db_type = db_config.get("type", "sqlite")

    # Read batch_size from database_formatting section
    formatting_config = config.get("database_formatting", {})
    writer_retry_config = db_config.get("writer_retry", {})
    formatting_retry = formatting_config.get("writer_retry")
    if formatting_retry:
        writer_retry_config = {**writer_retry_config, **formatting_retry}
    batch_size = formatting_config.get("batch_size", 1000)
    taxonomy_fetch_chunk_size = formatting_config.get("taxonomy_fetch_chunk_size", 10)

    if db_type == "sqlite":
        sqlite_config = db_config.get("sqlite", {})
        db_dir = Path(sqlite_config.get("path", "results/db"))
        db_name = sqlite_config.get("name", "its.db")
        db_path = db_dir / db_name

        return DatabaseManager(
            db_type="sqlite",
            path=str(db_path),
            writer_retry=writer_retry_config,
            batch_size=batch_size,
            taxonomy_fetch_chunk_size=taxonomy_fetch_chunk_size
        )

    elif db_type == "postgresql":
        pg_config = db_config.get("postgresql", {})
        return DatabaseManager(
            db_type="postgresql",
            host=pg_config.get("host", "localhost"),
            port=pg_config.get("port", 5432),
            user=pg_config.get("user", "its_user"),
            password=pg_config["password"],  # Required
            database=pg_config.get("database", "its_pipeline"),
            schema=pg_config.get("schema", "public"),  # Schema support
            pool_size=pg_config.get("pool_size", 10),
            max_overflow=pg_config.get("max_overflow", 20),
            writer_retry=writer_retry_config,
            batch_size=batch_size,
            taxonomy_fetch_chunk_size=taxonomy_fetch_chunk_size
        )

    else:
        raise ValueError(f"Unknown database type: {db_type}")


def initialize_database(config: Dict) -> DatabaseManager:
    """Initialize database with schema creation."""
    db_manager = get_database_manager(config)
    db_manager.create_schema()
    return db_manager
