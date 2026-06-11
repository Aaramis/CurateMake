#!/usr/bin/env python3
"""
SQLAlchemy ORM models for bioinformatics pipeline database.

This module provides unified database models that work transparently with both
SQLite and PostgreSQL backends. The models match the existing database schema
while enabling ORM features and better type safety.

Compatible backends:
- SQLite 3.x (development, small datasets)
- PostgreSQL 12+ (production, large datasets)

Usage:
    from db_models import Base, Source, Sequence, Taxonomy, SeqTax
    from db_models import get_engine, get_session, init_database

    # Create engine (SQLite or PostgreSQL)
    engine = get_engine("sqlite", path="data/db/its.db")
    # or
    engine = get_engine("postgresql", host="localhost", user="its_user", ...)

    # Initialize schema
    init_database(engine)

    # Create session
    session = get_session(engine)

    # Use ORM
    source = Source(name="UNITE", file_path="/path/to/unite.fasta")
    session.add(source)
    session.commit()
"""

from sqlalchemy import create_engine, Column, Integer, String, Float, Text, DateTime
from sqlalchemy import ForeignKey, CheckConstraint, Index, UniqueConstraint, event
from sqlalchemy.orm import declarative_base, relationship, sessionmaker, Session
from sqlalchemy.sql import func
from datetime import datetime
from typing import Optional, Dict, Any
from pathlib import Path
import logging
import uuid

logger = logging.getLogger(__name__)

# Declarative base for all models
Base = declarative_base()


class Source(Base):
    """
    Data source metadata table.

    Tracks different input databases (UNITE, BOLD, GBIF, etc.) with their
    download timestamps and file paths.
    """
    __tablename__ = 'sources'

    source_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, nullable=False, unique=True)
    file_path = Column(String)
    date_download = Column(DateTime, default=datetime.utcnow, server_default=func.now())

    # Relationships
    sequences = relationship("Sequence", back_populates="source", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Source(id={self.source_id}, name='{self.name}')>"


class Sequence(Base):
    """
    Core sequence data table with hierarchical relationships.

    Stores DNA/RNA sequences with quality metrics, barcode information, and
    parent-child relationships for sequence variants.
    """
    __tablename__ = 'sequences'

    seq_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    source_id = Column(String(36), ForeignKey('sources.source_id', ondelete='CASCADE'), nullable=False)
    original_id = Column(String, nullable=False)
    sequence = Column(Text, nullable=False)
    length = Column(Integer, nullable=False)
    percentage_n = Column(Float, default=0.0)
    barcode = Column(String, default='ITS')
    old_barcode = Column(String)
    barcode_status = Column(String, default='Pending')
    msa_id = Column(String(36))
    parent_seq_id = Column(String(36), ForeignKey('sequences.seq_id', ondelete='SET NULL'))

    # Relationships
    source = relationship("Source", back_populates="sequences")
    parent = relationship("Sequence", remote_side=[seq_id], backref="children")
    seq_tax_links = relationship("SeqTax", back_populates="sequence", cascade="all, delete-orphan")

    # Constraints
    __table_args__ = (
        CheckConstraint('length > 0', name='check_length_positive'),
        CheckConstraint('percentage_n >= 0 AND percentage_n <= 100', name='check_percentage_n_range'),
        CheckConstraint(
            "barcode_status IN ('Pending', 'Duplicate', 'ITSx', 'Processed', 'no_detection')",
            name='check_barcode_status_values'
        ),
        CheckConstraint(
            "barcode IN ('ITS', 'ITS1', 'ITS2', 'LSU', 'SSU', '5.8S')",
            name='check_barcode_values'
        ),
        # Performance indices
        Index('idx_sequences_source', 'source_id'),
        Index('idx_sequences_parent', 'parent_seq_id'),
        Index('idx_sequences_barcode', 'barcode'),
        Index('idx_sequences_length', 'length'),
        Index('idx_sequences_status', 'barcode_status'),
        Index('idx_sequences_original_id', 'original_id'),  # For bulk_insert_sequences fetch queries
    )

    def __repr__(self):
        return f"<Sequence(id={self.seq_id}, original_id='{self.original_id}', length={self.length})>"


class Taxonomy(Base):
    """
    Hierarchical taxonomic classification table.

    Stores taxonomic lineages (kingdom -> species) with support for hierarchical
    relationships via parent_tax_id self-reference.
    """
    __tablename__ = 'taxonomy'

    tax_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    kingdom = Column(String)
    phylum = Column(String)
    class_ = Column('class_', String)  # 'class' is reserved in Python
    order_ = Column('order_', String)  # 'order' is reserved in Python
    family = Column(String)
    genus = Column(String)
    species = Column(String)
    parent_tax_id = Column(String(36), ForeignKey('taxonomy.tax_id', ondelete='SET NULL'))

    # Relationships
    parent = relationship("Taxonomy", remote_side=[tax_id], backref="children")
    seq_tax_links = relationship("SeqTax", back_populates="taxonomy", cascade="all, delete-orphan")

    # Indices and constraints
    __table_args__ = (
        UniqueConstraint('kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species',
                         name='uq_taxonomy_lineage'),
        Index('idx_taxonomy_parent', 'parent_tax_id'),
        Index('idx_taxonomy_genus_species', 'genus', 'species'),
    )

    def __repr__(self):
        return f"<Taxonomy(id={self.tax_id}, genus='{self.genus}', species='{self.species}')>"


class SeqTax(Base):
    """
    Sequence-taxonomy junction table.

    Links sequences to taxonomies with backbone (data source) and validation
    status information. Composite primary key (seq_id, tax_id, backbone) allows
    multiple taxonomic assignments per sequence from different sources.
    """
    __tablename__ = 'seq_tax'

    seq_id = Column(String(36), ForeignKey('sequences.seq_id', ondelete='CASCADE'), primary_key=True)
    tax_id = Column(String(36), ForeignKey('taxonomy.tax_id', ondelete='CASCADE'), primary_key=True)
    backbone = Column(String, primary_key=True, server_default='Raw')
    status = Column(String, server_default='Pending')
    confidence = Column(Float)  # Optional confidence score (e.g., from SATIVA)

    # Relationships
    sequence = relationship("Sequence", back_populates="seq_tax_links")
    taxonomy = relationship("Taxonomy", back_populates="seq_tax_links")

    # Constraints
    __table_args__ = (
        CheckConstraint(
            "status IN ('Pending', 'Check', 'Correct', 'Failed')",
            name='check_status_values'
        ),
        CheckConstraint(
            "backbone IN ('Raw', 'CoL', 'GBIF', 'Sativa')",
            name='check_backbone_values'
        ),
        # Performance indices
        Index('idx_seq_tax_seq', 'seq_id'),
        Index('idx_seq_tax_tax', 'tax_id'),
        Index('idx_seq_tax_backbone', 'backbone'),
        Index('idx_seq_tax_status', 'status'),
    )

    def __repr__(self):
        return f"<SeqTax(seq_id={self.seq_id}, tax_id={self.tax_id}, backbone='{self.backbone}')>"


class SativaWarning(Base):
    """
    SATIVA warnings table for tracking modifications that require attention.

    Stores warnings for:
    - Simple deletions (not integrated as seq_tax)
    - CoL validation discrepancies
    - Low confidence modifications
    """
    __tablename__ = 'sativa_warnings'

    warning_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    seq_id = Column(String(36), ForeignKey('sequences.seq_id', ondelete='CASCADE'), nullable=False)
    seq_original_id = Column(String, nullable=False)
    original_tax_id = Column(String(36), ForeignKey('taxonomy.tax_id', ondelete='CASCADE'), nullable=False)

    # Warning type
    warning_type = Column(String, nullable=False)  # DELETION_SIMPLE, DELETION_MULTIPLE, COL_MATCHED, COL_MODIFIED_SATIVA, COL_NO_MATCH

    # SATIVA information
    sativa_modification_type = Column(String)  # DELETION_SIMPLE, DELETION_MULTIPLE, SUBSTITUTION, COMPLETION
    sativa_confidence = Column(Float)
    sativa_original_taxonomy = Column(Text)  # Semicolon-separated path
    sativa_proposed_taxonomy = Column(Text)  # Semicolon-separated path
    deleted_levels = Column(Text)  # JSON array of deleted levels

    # CoL validation information
    col_validation_attempted = Column(String, server_default='0')  # '0' or '1' for compatibility
    col_returned_taxonomy = Column(Text)  # Semicolon-separated path
    col_match_quality = Column(String)  # exact, modified, no_match

    # Metadata
    barcode = Column(String, nullable=False)
    tax_level = Column(String, nullable=False)
    mislabeled_level = Column(String)
    per_rank_confidence = Column(Text)
    created_at = Column(DateTime, default=datetime.utcnow, server_default=func.now())

    # Relationships
    sequence = relationship("Sequence", foreign_keys=[seq_id])
    taxonomy = relationship("Taxonomy", foreign_keys=[original_tax_id])

    # Constraints
    __table_args__ = (
        CheckConstraint(
            "warning_type IN ('DELETION_SIMPLE', 'DELETION_MULTIPLE', 'COL_MATCHED', 'COL_MODIFIED_SATIVA', 'COL_NO_MATCH')",
            name='check_warning_type_values'
        ),
        CheckConstraint(
            "sativa_modification_type IS NULL OR sativa_modification_type IN ('DELETION_SIMPLE', 'DELETION_MULTIPLE', 'SUBSTITUTION', 'COMPLETION')",
            name='check_sativa_modification_type_values'
        ),
        # Performance indices
        Index('idx_sativa_warnings_seq', 'seq_id'),
        Index('idx_sativa_warnings_type', 'warning_type'),
        Index('idx_sativa_warnings_barcode', 'barcode'),
        Index('idx_sativa_warnings_original_id', 'seq_original_id'),
    )

    def __repr__(self):
        return f"<SativaWarning(id={self.warning_id}, seq_id={self.seq_id}, type='{self.warning_type}')>"


# ========== Engine and Session Factory Functions ==========

def get_engine(db_type: str = "sqlite", echo: bool = False, **kwargs):
    """
    Create SQLAlchemy engine for specified database backend.

    Args:
        db_type: Database type ("sqlite" or "postgresql")
        echo: Enable SQL query logging (useful for debugging)
        **kwargs: Backend-specific connection parameters

    SQLite kwargs:
        path (str): Database file path (default: "data/db/its.db")

    PostgreSQL kwargs:
        host (str): Database host (default: "localhost")
        port (int): Database port (default: 5432)
        user (str): Database user (default: "its_user")
        password (str): Database password (required)
        database (str): Database name (default: "its_pipeline")
        pool_size (int): Connection pool size (default: 10)
        max_overflow (int): Max overflow connections (default: 20)

    Returns:
        SQLAlchemy Engine instance

    Raises:
        ValueError: If db_type is not supported

    Examples:
        # SQLite
        engine = get_engine("sqlite", path="data/db/its.db")

        # PostgreSQL
        engine = get_engine(
            "postgresql",
            host="localhost",
            user="its_user",
            password="secret",
            database="its_pipeline"
        )
    """
    if db_type == "sqlite":
        path = kwargs.get('path', 'data/db/its.db')
        db_path = Path(path)
        db_path.parent.mkdir(parents=True, exist_ok=True)

        url = f"sqlite:///{path}"
        engine = create_engine(
            url,
            echo=echo,
            # SQLite-specific optimizations
            connect_args={
                'timeout': 30.0,
                'check_same_thread': False  # Allow multi-threading
            }
        )

        # Enable SQLite pragmas for better performance
        @event.listens_for(engine, "connect")
        def set_sqlite_pragma(dbapi_conn, connection_record):
            cursor = dbapi_conn.cursor()
            cursor.execute("PRAGMA foreign_keys=ON")
            cursor.execute("PRAGMA journal_mode=WAL")
            cursor.close()

        logger.info(f"Created SQLite engine: {path}")

    elif db_type == "postgresql":
        host = kwargs.get('host', 'localhost')
        port = kwargs.get('port', 5432)
        user = kwargs.get('user', 'its_user')
        password = kwargs.get('password')
        database = kwargs.get('database', 'its_pipeline')
        schema = kwargs.get('schema', 'public')
        pool_size = kwargs.get('pool_size', 10)
        max_overflow = kwargs.get('max_overflow', 20)

        if not password:
            raise ValueError("PostgreSQL password is required")

        url = f"postgresql://{user}:{password}@{host}:{port}/{database}"

        # PostgreSQL connection options with timeout configuration
        connect_args = {
            'connect_timeout': 30,  # Connection timeout
            'options': '-c statement_timeout=600000 -c idle_in_transaction_session_timeout=600000'  # 10min query, 10min idle
        }

        engine = create_engine(
            url,
            echo=echo,
            pool_size=pool_size,
            max_overflow=max_overflow,
            pool_pre_ping=True,  # Verify connections before using
            pool_recycle=300,  # Recycle connections after 5 minutes to avoid stale connections
            pool_timeout=30,  # Timeout for getting connection from pool
            connect_args=connect_args,
            isolation_level="READ COMMITTED"  # Best for concurrent writes (parallel rules)
        )

        # Set search_path for schema
        # Note: isolation level is already set via isolation_level parameter above
        @event.listens_for(engine, "connect")
        def set_schema_search_path(dbapi_conn, connection_record):
            cursor = dbapi_conn.cursor()
            cursor.execute(f"SET search_path TO {schema}")
            cursor.close()

        logger.info(f"Created PostgreSQL engine: {host}:{port}/{database} (schema: {schema})")

    else:
        raise ValueError(f"Unsupported database type: {db_type}. Use 'sqlite' or 'postgresql'.")

    return engine


def get_session(engine) -> Session:
    """
    Create a new SQLAlchemy session from engine.

    Args:
        engine: SQLAlchemy engine instance

    Returns:
        Session instance

    Example:
        engine = get_engine("sqlite", path="data/db/its.db")
        session = get_session(engine)

        # Use session
        sources = session.query(Source).all()

        # Always close when done
        session.close()
    """
    SessionFactory = sessionmaker(bind=engine)
    return SessionFactory()


def init_database(engine, schema=None):
    """
    Initialize database schema by creating all tables.

    This function is idempotent - safe to call multiple times.
    Uses CREATE TABLE IF NOT EXISTS internally.

    Args:
        engine: SQLAlchemy engine instance
        schema: Optional schema name for PostgreSQL

    Example:
        engine = get_engine("sqlite", path="data/db/its.db")
        init_database(engine)
    """
    # For PostgreSQL with custom schema, modify table definitions before create_all
    if schema and schema != 'public':
        logger.info(f"Setting schema {schema} for all tables")
        # Directly modify the schema attribute of each table object
        for table_name, table in Base.metadata.tables.items():
            table.schema = schema
            logger.debug(f"Table {table_name} schema set to {schema}")

    # Create all tables
    Base.metadata.create_all(engine)

    logger.info(f"Database schema initialized successfully{f' in schema {schema}' if schema else ''}")


def drop_all_tables(engine):
    """
    Drop all tables from the database.

    ⚠️  WARNING: This will delete ALL data! Use with extreme caution.

    Args:
        engine: SQLAlchemy engine instance
    """
    Base.metadata.drop_all(engine)
    logger.warning("All tables dropped from database")


# ========== Utility Functions ==========

def get_table_names():
    """
    Get list of all table names defined in models.

    Returns:
        List of table name strings
    """
    return [table.name for table in Base.metadata.sorted_tables]


def get_model_by_tablename(tablename: str):
    """
    Get model class by table name.

    Args:
        tablename: Table name string

    Returns:
        Model class or None if not found
    """
    for mapper in Base.registry.mappers:
        model = mapper.class_
        if model.__tablename__ == tablename:
            return model
    return None


# ========== Database Info Function ==========

def get_database_info(engine) -> Dict[str, Any]:
    """
    Get information about the database and schema.

    Args:
        engine: SQLAlchemy engine instance

    Returns:
        Dictionary with database information
    """
    info = {
        'dialect': engine.dialect.name,
        'driver': engine.driver,
        'tables': get_table_names(),
        'table_count': len(Base.metadata.tables)
    }

    # Add backend-specific info
    if engine.dialect.name == 'sqlite':
        info['database_path'] = str(engine.url.database)
    elif engine.dialect.name == 'postgresql':
        info['host'] = engine.url.host
        info['port'] = engine.url.port
        info['database'] = engine.url.database
        info['user'] = engine.url.username

    return info


# ========== Example Usage ==========

if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.INFO)

    print("=" * 60)
    print("SQLAlchemy Models - Example Usage")
    print("=" * 60)

    # Example 1: SQLite
    print("\n--- Example 1: SQLite Backend ---")
    sqlite_engine = get_engine("sqlite", path=":memory:", echo=True)
    init_database(sqlite_engine)

    session = get_session(sqlite_engine)

    # Create test data
    source = Source(name="TEST_DB", file_path="/test/path.fasta")
    session.add(source)
    session.commit()

    print(f"Created source: {source}")
    print(f"Database info: {get_database_info(sqlite_engine)}")

    session.close()

    print("\n✓ Examples completed successfully")
