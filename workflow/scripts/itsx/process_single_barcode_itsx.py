#!/usr/bin/env python3
"""
ITSx Results Processor - Production Version 2.0

Processes ITSx extraction results for a specific barcode (ITS, ITS1, or ITS2),
extracting subregions and reconstructing full ITS sequences.

Architecture:
- One transaction per chunk for PostgreSQL compatibility
- Streaming parent metadata fetching to minimize memory usage
- Single-pass FASTA file reading for optimal I/O
- Compatible with both SQLite and PostgreSQL backends

Usage (Snakemake):
    Automatically invoked by Snakemake workflow

Usage (Standalone CLI for testing):
    python process_single_barcode_itsx.py <barcode> <barcode_dir> <config_yaml> <output_json>

Author: Refactored for PostgreSQL compatibility
Date: 2026-01-05
"""

import sys
import json
import logging
import time
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Iterator
from datetime import datetime
from dataclasses import dataclass, asdict

from sqlalchemy import text
from sqlalchemy import insert as sa_insert
from sqlalchemy.exc import OperationalError, IntegrityError, SQLAlchemyError

try:
    from psycopg2.errors import SerializationFailure as PGSerializationFailure
except ImportError:  # SQLite-only environments
    PGSerializationFailure = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ProcessingStats:
    """Statistics for tracking processing progress."""
    barcode: str
    status: str = "Running"
    start_time: str = ""
    end_time: str = ""
    initial_sequences_processed: int = 0
    subregions: Dict[str, Dict[str, int]] = None
    reconstructed_its: Dict[str, int] = None
    no_detection: int = 0
    total_new_sequences: int = 0
    error_message: str = ""

    def __post_init__(self):
        if self.subregions is None:
            self.subregions = {}
        if self.reconstructed_its is None:
            self.reconstructed_its = {}
        if not self.start_time:
            self.start_time = datetime.now().isoformat()

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


class BarcodeITSxProcessor:
    """
    Processes ITSx results for a single barcode with PostgreSQL-compatible architecture.

    Key design principles:
    - One transaction per chunk (not one transaction for all chunks)
    - Streaming data fetching to avoid memory bloat
    - Single-pass file reading for efficiency
    - Comprehensive error handling with retries
    """

    def __init__(self, barcode: str, barcode_dir: Path, data_dir: Path, config: dict, db_manager):
        """
        Initialize the processor.

        Args:
            barcode: Target barcode (ITS, ITS1, or ITS2)
            barcode_dir: Directory containing ITSx output files
            data_dir: Directory containing original FASTA files
            config: Configuration dictionary from YAML
            db_manager: Database manager instance (from db_manager.get_database_manager)
        """
        self.barcode = barcode
        self.barcode_dir = barcode_dir
        self.data_dir = data_dir
        self.config = config
        self.db = db_manager

        # Extract configuration parameters (already flattened in main())
        # All batch sizes use database_formatting.batch_size from config (default 1000)
        self.insert_batch_size = config.get('insert_batch_size', 1000)
        self.status_update_chunk_size = config.get('status_update_chunk_size', 1000)

        # Parent metadata configuration
        self.parent_metadata_streaming = config.get('parent_metadata_streaming', True)
        self.parent_metadata_stream_fetch_size = config.get('parent_metadata_stream_fetch_size', 1000)
        self.parent_metadata_chunk_size = config.get('parent_metadata_chunk_size', 5000)
        self.parent_metadata_fallback_chunk_size = config.get('parent_metadata_fallback_chunk_size', 10000)
        self.parent_metadata_filter_barcode = config.get('parent_metadata_stream_filter_barcode', True)

        # Sequence validation configuration
        self.min_length = config.get('min_length', 50)
        self.max_n_fraction = config.get('max_n_fraction', 0.5)
        valid_chars_str = config.get('valid_chars', 'ACGTURYSWKMBDHVN-')
        self.valid_chars = set(valid_chars_str) if isinstance(valid_chars_str, str) else set(valid_chars_str)

        # Statistics
        self.stats = ProcessingStats(barcode=barcode)

        logger.info(f"Initialized BarcodeITSxProcessor for {barcode}")
        logger.info(f"  Batch configuration (from database_formatting.batch_size):")
        logger.info(f"    - Insert batch size: {self.insert_batch_size}")
        logger.info(f"    - Status update chunk size: {self.status_update_chunk_size}")
        logger.info(f"    - Parent metadata chunk size: {self.parent_metadata_chunk_size}")
        logger.info(f"  Retry configuration (from database.writer_retry):")
        logger.info(f"    - Max retry attempts: {self.db.writer_retry_attempts}")
        logger.info(f"  Other settings:")
        logger.info(f"    - Parent metadata streaming: {self.parent_metadata_streaming}")
        logger.info(f"    - Min sequence length: {self.min_length}")

    def _validate_sequence(self, sequence: str) -> bool:
        """
        Validate sequence quality.

        Args:
            sequence: DNA sequence string

        Returns:
            True if sequence passes validation, False otherwise
        """
        if not sequence:
            return False

        seq_len = len(sequence)

        # Check minimum length
        if seq_len < self.min_length:
            return False

        # Check N content
        n_count = sum(1 for c in sequence if c.upper() == 'N')
        n_fraction = n_count / seq_len
        if n_fraction > self.max_n_fraction:
            return False

        # Check valid characters (sample first 100 chars for speed)
        sample = sequence[:min(100, seq_len)]
        if not all(c.upper() in self.valid_chars for c in sample):
            return False

        return True

    def _parse_fasta_streaming(self, file_path: Path):
        """
        Parse FASTA file in streaming mode (generator).

        Yields (header, sequence) tuples without loading entire file into memory.

        Args:
            file_path: Path to FASTA file (can be gzipped)

        Yields:
            Tuple of (header, sequence)
        """
        if not file_path.exists():
            logger.warning(f"FASTA file not found: {file_path}")
            return

        open_func = gzip.open if file_path.suffix == '.gz' else open

        try:
            with open_func(file_path, 'rt', encoding='utf-8') as f:
                header = None
                sequence_parts = []

                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        # Yield previous record
                        if header and sequence_parts:
                            sequence = ''.join(sequence_parts)
                            yield header, sequence

                        # Start new record
                        header = line[1:]  # Remove '>'
                        sequence_parts = []
                    else:
                        sequence_parts.append(line)

                # Yield final record
                if header and sequence_parts:
                    sequence = ''.join(sequence_parts)
                    yield header, sequence

        except Exception as e:
            logger.error(f"Error reading FASTA file {file_path}: {e}")
            raise

    def _parse_header_for_parent_id(self, header: str) -> Optional[str]:
        """
        Extract parent seq_id from ITSx header.

        ITSx headers format: "parent_seq_id|region|coords"

        Args:
            header: FASTA header line (without '>')

        Returns:
            Parent seq_id or None
        """
        try:
            # ITSx format: seq_id|ITS1|start-end
            parts = header.split('|')
            if parts:
                return parts[0].strip()
        except Exception as e:
            logger.debug(f"Failed to parse header: {header}: {e}")

        return None

    def _bulk_update_status(self, updates: List[Tuple[str, str]]) -> int:
        """
        Bulk update sequence statuses using one transaction per chunk.

        CRITICAL: This method uses ONE TRANSACTION PER CHUNK (not one transaction
        for all chunks) to avoid PostgreSQL serialization errors.

        Args:
            updates: List of (new_status, seq_id) tuples

        Returns:
            Number of sequences updated
        """
        if not updates:
            return 0

        # Imported here, not at module scope: sys.path only gains scripts/utils in
        # main(), and Snakemake runs this file from a temp copy under .snakemake/scripts.
        from db_manager import pg_hint

        # Deduplicate updates (keep last status for each seq_id)
        deduped = {}
        for status, seq_id in updates:
            deduped[seq_id] = status

        unique_updates = [(status, seq_id) for seq_id, status in deduped.items()]

        if len(unique_updates) < len(updates):
            logger.warning(
                f"Deduped {len(updates) - len(unique_updates)} duplicate seq_id updates "
                f"for barcode {self.barcode}."
            )

        # Process in chunks
        total_updated = 0
        chunk_size = self.status_update_chunk_size
        chunks = [unique_updates[i:i + chunk_size] for i in range(0, len(unique_updates), chunk_size)]

        logger.info(f"Updating {len(unique_updates)} sequence statuses in {len(chunks)} chunks")

        # ONE TRANSACTION PER CHUNK (PostgreSQL compatibility)
        #
        # `chunks` is consumed as a QUEUE, not a fixed list: a chunk that fails is halved and
        # both halves are pushed back. _get_parent_metadata below already worked this way;
        # this method did not, so a single transient failure killed the whole rule. It is a
        # real scenario -- on 2026-08-27/28 the cluster returned "Rpc timeout, passed:
        # 15.062s, timeout: 15.000s" on a 9 946-id UPDATE while being perfectly healthy
        # minutes later (measured: the same query over 10 000 ids plans as an Index Scan and
        # runs in 0.04 s). Retrying the identical statement is NOT the answer -- a deadline
        # is deterministic for the size that hit it -- so shrink and continue instead.
        pending = list(chunks)
        chunk_idx = 0
        while pending:
            chunk = pending.pop(0)
            chunk_idx += 1
            seq_ids = [seq_id for _, seq_id in chunk]
            status_value = chunk[0][0]  # Assume same status for chunk

            # Use db_manager's retry logic
            def update_chunk():
                with self.db.get_session() as session:
                    # Use database-specific syntax for array parameters
                    if self.db.db_type == "postgresql":
                        query = text(f"""
                            {pg_hint([("s", "sequences_pkey")])}
                            UPDATE {self.db.sequences_table} AS s
                            SET barcode_status = :status
                            WHERE s.seq_id = ANY(:seq_ids)
                        """)
                        result = session.execute(query, {"status": status_value, "seq_ids": seq_ids})
                    else:  # SQLite
                        # SQLite doesn't support arrays, use IN clause
                        placeholders_named = ','.join([f':param_{i}' for i in range(len(seq_ids))])
                        query = text(f"""
                            UPDATE {self.db.sequences_table}
                            SET barcode_status = :status
                            WHERE seq_id IN ({placeholders_named})
                        """)
                        params = {"status": status_value}
                        for i, seq_id in enumerate(seq_ids):
                            params[f'param_{i}'] = seq_id
                        result = session.execute(query, params)

                    return result.rowcount

            # Execute with automatic retry handling from db_manager
            try:
                rowcount = self.db._execute_with_retry(
                    f"bulk_update_status_chunk_{chunk_idx}",
                    update_chunk,
                    retry_deadline=False,  # the except below resplits; see db_manager
                )
            except Exception as e:
                if len(chunk) == 1:
                    logger.error(
                        f"Failed to update status for single seq_id {chunk[0][1]}: {e}")
                    raise
                mid = len(chunk) // 2
                logger.warning(
                    f"Status update chunk size {len(chunk)} failed; "
                    f"splitting into {mid} and {len(chunk) - mid}")
                pending.insert(0, chunk[mid:])
                pending.insert(0, chunk[:mid])
                continue
            total_updated += rowcount

        logger.info(f"Updated {total_updated} sequences to status '{status_value}'")
        return total_updated

    def _get_parent_metadata(self, parent_ids: List[str]) -> Dict[str, Dict[str, str]]:
        """
        Fetch parent sequence metadata (source_id, original_id, tax_id).

        Uses streaming mode if enabled, otherwise falls back to chunked queries.

        Args:
            parent_ids: List of parent seq_ids

        Returns:
            Dict mapping seq_id -> {source_id, original_id, tax_id}
        """
        if not parent_ids:
            return {}

        logger.info(f"Fetching metadata for {len(parent_ids)} parent sequences")

        if self.parent_metadata_streaming:
            return self._get_parent_metadata_streaming(parent_ids)
        else:
            return self._get_parent_metadata_chunked(parent_ids)

    def _get_parent_metadata_streaming(self, parent_ids: List[str]) -> Dict[str, Dict[str, str]]:
        """
        Fetch parent metadata using server-side cursor (streaming mode).

        Memory-efficient for large datasets.

        Uses a session with READ COMMITTED isolation (requires transaction for server-side cursors).
        Wrapped with db_manager retry logic for automatic error recovery.
        """
        parent_ids_set = set(parent_ids)

        # Build query with optional barcode filter
        barcode_condition = ""
        params = {}

        if self.parent_metadata_filter_barcode:
            if self.barcode == "ITS":
                # For ITS, accept ITS barcode or NULL (original sequences)
                barcode_condition = "AND (s.barcode = :barcode OR s.barcode IS NULL)"
            else:
                # For ITS1/ITS2, exact match
                barcode_condition = "AND s.barcode = :barcode"
            params['barcode'] = self.barcode

        # Imported here, not at module scope: sys.path only gains scripts/utils in
        # main(), and Snakemake runs this file from a temp copy under .snakemake/scripts.
        from db_manager import pg_hint

        # INNER JOIN to get only sequences with CoL taxonomy
        query = text(f"""
            {pg_hint([("s", "sequences_pkey"), ("st", "seq_tax_pkey")], batched=("s", "st"), leading=("s", "st"))}
            SELECT s.seq_id, s.source_id, s.original_id, st.tax_id
            FROM {self.db.sequences_table} s
            JOIN {self.db.seq_tax_table} st ON s.seq_id = st.seq_id
            WHERE st.backbone = 'CoL'
            {barcode_condition}
        """)

        # Streaming fetch with retry logic
        def streaming_fetch():
            metadata = {}
            remaining = parent_ids_set.copy()

            # Use get_session for proper transaction management
            with self.db.get_session() as session:
                # For PostgreSQL, enable server-side cursor with yield_per
                if self.db.db_type == "postgresql":
                    # Use execution_options with stream_results for server-side cursor
                    result = session.connection().execution_options(
                        stream_results=True
                    ).execute(query, params)
                else:
                    result = session.execute(query, params)

                # Fetch in chunks using server-side cursor
                while remaining:
                    rows = result.fetchmany(self.parent_metadata_stream_fetch_size)
                    if not rows:
                        break

                    for row in rows:
                        seq_id = row[0]
                        if seq_id not in remaining:
                            continue
                        metadata[seq_id] = {
                            'source_id': row[1],
                            'original_id': row[2],
                            'tax_id': row[3]
                        }
                        remaining.remove(seq_id)

                if remaining:
                    logger.warning(
                        f"Parent metadata streaming missed {len(remaining)} ids for barcode {self.barcode}"
                    )

            return metadata

        # Execute with automatic retry handling
        try:
            metadata = self.db._execute_with_retry(
                "parent_metadata_streaming",
                streaming_fetch,
                retry_deadline=False,  # falls back to chunked mode below
            )
            logger.info(f"Fetched metadata for {len(metadata)} parent sequences (streaming mode)")
            return metadata
        except Exception as e:
            logger.warning(f"Streaming fetch failed after retries: {e}, falling back to chunked mode")
            return self._get_parent_metadata_chunked(parent_ids)

    def _get_parent_metadata_chunked(self, parent_ids: List[str]) -> Dict[str, Dict[str, str]]:
        """
        Fetch parent metadata using chunked IN queries (fallback mode).
        Uses db_manager retry logic for transient errors and chunk splitting for oversized chunks.
        """
        metadata = {}
        chunk_size = self.parent_metadata_fallback_chunk_size

        # Deduplicate parent IDs
        unique_parent_ids = list(dict.fromkeys(parent_ids))

        # Split into chunks with retry logic for failed chunks
        pending = [unique_parent_ids[i:i + chunk_size] for i in range(0, len(unique_parent_ids), chunk_size)]

        from db_manager import pg_hint  # see the note in _bulk_update_status

        logger.info(f"Fetching metadata in {len(pending)} chunks (fallback mode)")

        while pending:
            chunk = pending.pop()

            # Define fetch operation for this chunk
            def fetch_chunk():
                with self.db.get_session() as session:
                    if self.db.db_type == "postgresql":
                        # INNER JOIN with backbone='CoL' filter
                        query = text(f"""
                            {pg_hint([("s", "sequences_pkey"), ("st", "seq_tax_pkey")], batched=("s", "st"), leading=("s", "st"))}
                            SELECT s.seq_id, s.source_id, s.original_id, st.tax_id
                            FROM {self.db.sequences_table} s
                            JOIN {self.db.seq_tax_table} st ON s.seq_id = st.seq_id
                            WHERE s.seq_id = ANY(:seq_ids) AND st.backbone = 'CoL'
                        """)
                        result = session.execute(query, {"seq_ids": chunk})
                    else:  # SQLite
                        placeholders = ','.join([f':id_{i}' for i in range(len(chunk))])
                        query = text(f"""
                            SELECT s.seq_id, s.source_id, s.original_id, st.tax_id
                            FROM {self.db.sequences_table} s
                            JOIN {self.db.seq_tax_table} st ON s.seq_id = st.seq_id
                            WHERE s.seq_id IN ({placeholders}) AND st.backbone = 'CoL'
                        """)
                        params = {f'id_{i}': seq_id for i, seq_id in enumerate(chunk)}
                        result = session.execute(query, params)

                    chunk_metadata = {}
                    for row in result:
                        chunk_metadata[row[0]] = {
                            'source_id': row[1],
                            'original_id': row[2],
                            'tax_id': row[3]
                        }
                    return chunk_metadata

            try:
                # Execute with automatic retry handling
                chunk_metadata = self.db._execute_with_retry(
                    f"parent_metadata_chunk_{len(chunk)}",
                    fetch_chunk,
                    retry_deadline=False,  # the except below resplits
                )
                metadata.update(chunk_metadata)

            except Exception as e:
                # If chunk fails after retries and has more than 1 element, split it
                if len(chunk) == 1:
                    logger.error(f"Failed to fetch metadata for single ID {chunk[0]} after retries: {e}")
                    raise

                mid = len(chunk) // 2
                logger.warning(
                    f"Parent metadata chunk size {len(chunk)} failed after retries; "
                    f"splitting into {mid} and {len(chunk) - mid}"
                )
                pending.append(chunk[mid:])
                pending.append(chunk[:mid])
                continue

        logger.info(f"Fetched metadata for {len(metadata)} parent sequences (chunked mode)")
        return metadata

    def _bulk_insert_sequences(self, sequences: List[dict]) -> int:
        """
        Bulk insert sequences with their taxonomic relationships.

        Args:
            sequences: List of sequence dictionaries with keys:
                - sequence, barcode, length
                - source_id, original_id, parent_seq_id, tax_id

        Returns:
            Number of sequences inserted
        """
        if not sequences:
            return 0

        # Prepare all data first
        seq_inserts = []
        for seq_data in sequences:
            # Generate seq_id
            seq_id = seq_data.get('seq_id') or self._generate_seq_id()

            seq_inserts.append({
                'seq_id': seq_id,
                'source_id': seq_data['source_id'],
                'original_id': seq_data['original_id'],
                'sequence': seq_data['sequence'],
                'length': seq_data['length'],
                'barcode': seq_data['barcode'],
                'barcode_status': 'ITSx',  # Use 'ITSx' status for processed sequences
                'parent_seq_id': seq_data.get('parent_seq_id')
            })

            # Track for seq_tax insertion
            seq_data['seq_id'] = seq_id

        # Prepare seq_tax relationships
        seq_tax_inserts = []
        for seq_data in sequences:
            if seq_data.get('tax_id'):
                seq_tax_inserts.append({
                    'seq_id': seq_data['seq_id'],
                    'tax_id': seq_data['tax_id'],
                    'backbone': 'CoL',  # ITSx-derived sequences use CoL backbone
                    'status': 'Pending'
                })

        # Execute with retry logic
        def insert_operation():
            # Core table INSERTs, not text(): SQLAlchemy can only fold a list of parameter
            # dicts into multi-row VALUES statements for a construct it built itself. Given
            # a text() it falls back to cursor.executemany(), which psycopg2 runs as one
            # round trip PER ROW -- on a distributed database that is what this step was
            # paying. Measured on this cluster against a table carrying the same 7 indexes
            # and constraints, 10k sequences + 10k seq_tax: 26.7 s via text() vs 1.0 s here.
            # SQLAlchemy also pages the VALUES internally (1000 rows), so a 10k batch never
            # becomes the single oversized statement that closes the backend. Table-level
            # insert() stays Core-only, hence identical on SQLite; the models already carry
            # the schema db_manager set on them.
            from db_models import Sequence as SequenceModel, SeqTax as SeqTaxModel

            with self.db.get_session() as session:
                if seq_inserts:
                    session.execute(sa_insert(SequenceModel.__table__), seq_inserts)

                if seq_tax_inserts:
                    session.execute(sa_insert(SeqTaxModel.__table__), seq_tax_inserts)

                return len(seq_inserts)

        try:
            inserted_count = self.db._execute_with_retry(
                "bulk_insert_sequences",
                insert_operation
            )
            return inserted_count
        except IntegrityError as e:
            logger.warning(f"Integrity error during insert (possible duplicates): {e}")
            return 0
        except Exception as e:
            logger.error(f"Failed to insert sequences: {e}")
            raise

    def _generate_seq_id(self) -> str:
        """Generate a unique sequence ID (UUID)."""
        import uuid
        return str(uuid.uuid4())

    def _read_region_file_once(
        self,
        region: str,
        parent_metadata: Dict[str, Dict[str, str]]
    ) -> Tuple[List[dict], Dict[str, str]]:
        """
        Read a region FASTA file once, preparing both insert data and reconstruction dict.

        Args:
            region: Region name (ITS1, ITS2, 5_8S, etc.)
            parent_metadata: Parent metadata dictionary

        Returns:
            Tuple of (insert_list, sequence_dict)
            - insert_list: List of dicts ready for bulk insert
            - sequence_dict: Dict mapping parent_id -> sequence for reconstruction
        """
        file_path = self.barcode_dir / f"{self.barcode}.{region}.fasta"

        insert_list = []
        sequence_dict = {}
        found_count = 0

        for header, sequence in self._parse_fasta_streaming(file_path):
            found_count += 1

            # Extract parent ID
            parent_id = self._parse_header_for_parent_id(header)
            if not parent_id or parent_id not in parent_metadata:
                continue

            # Validate sequence
            if not self._validate_sequence(sequence):
                continue

            # Get parent metadata
            parent = parent_metadata[parent_id]

            # Prepare for insertion (only ITS1 and ITS2, not 5_8S)
            if region in ['ITS1', 'ITS2']:
                insert_list.append({
                    'sequence': sequence,
                    'barcode': region,
                    'length': len(sequence),
                    'source_id': parent['source_id'],
                    'original_id': parent['original_id'],
                    'parent_seq_id': parent_id,
                    'tax_id': parent.get('tax_id')
                })

            # Store for reconstruction
            sequence_dict[parent_id] = sequence

        logger.info(f"Region {region}: found {found_count} sequences, {len(insert_list)} valid for insert")

        return insert_list, sequence_dict

    def process(self) -> ProcessingStats:
        """
        Main processing workflow for a barcode.

        Seven-step process:
        1. Mark parent sequences as 'Processed'
        2. Collect parent IDs from ITSx output files
        3. Fetch parent metadata
        4. Read region files (ITS1, ITS2, 5_8S)
        5. Insert subregion sequences
        6. Reconstruct full ITS sequences
        7. Mark 'no_detection' sequences

        Returns:
            ProcessingStats object with results
        """
        try:
            logger.info(f"=== Processing barcode: {self.barcode} ===")

            # Step 1: Mark parent sequences as 'Processed'
            original_fasta = self.data_dir / f"{self.barcode}.fasta"
            parent_ids_from_fasta = []

            for header, _ in self._parse_fasta_streaming(original_fasta):
                parent_id = self._parse_header_for_parent_id(header)
                if parent_id:
                    parent_ids_from_fasta.append(parent_id)

            if parent_ids_from_fasta:
                updates = [("Processed", pid) for pid in parent_ids_from_fasta]
                updated = self._bulk_update_status(updates)
                self.stats.initial_sequences_processed = updated

            # Step 2: Collect parent IDs from ITSx files
            all_parent_ids = set()
            for region in ['ITS1', 'ITS2', '5_8S']:
                region_file = self.barcode_dir / f"{self.barcode}.{region}.fasta"
                for header, _ in self._parse_fasta_streaming(region_file):
                    parent_id = self._parse_header_for_parent_id(header)
                    if parent_id:
                        all_parent_ids.add(parent_id)

            logger.info(f"Found {len(all_parent_ids)} unique parent IDs in ITSx output")

            # Step 3: Fetch parent metadata
            parent_metadata = self._get_parent_metadata(list(all_parent_ids))

            # Step 4: Read region files
            its1_insert, its1_seqs = self._read_region_file_once('ITS1', parent_metadata)
            its2_insert, its2_seqs = self._read_region_file_once('ITS2', parent_metadata)
            _, its58s_seqs = self._read_region_file_once('5_8S', parent_metadata)

            # Update stats
            self.stats.subregions['ITS1'] = {
                'found': len(its1_seqs),
                'inserted': 0
            }
            self.stats.subregions['ITS2'] = {
                'found': len(its2_seqs),
                'inserted': 0
            }

            # Step 5: Insert subregions in batches
            total_inserted = 0

            for region_name, insert_list in [('ITS1', its1_insert), ('ITS2', its2_insert)]:
                for i in range(0, len(insert_list), self.insert_batch_size):
                    batch = insert_list[i:i + self.insert_batch_size]
                    inserted = self._bulk_insert_sequences(batch)
                    total_inserted += inserted
                    self.stats.subregions[region_name]['inserted'] += inserted

            # Step 6: Reconstruct full ITS sequences
            common_ids = set(its1_seqs.keys()) & set(its58s_seqs.keys()) & set(its2_seqs.keys())
            logger.info(f"Found {len(common_ids)} sequences with all three regions for ITS reconstruction")

            reconstruction_batch = []
            for parent_id in common_ids:
                full_seq = its1_seqs[parent_id] + its58s_seqs[parent_id] + its2_seqs[parent_id]

                if not self._validate_sequence(full_seq):
                    continue

                parent = parent_metadata.get(parent_id)
                if not parent:
                    continue

                reconstruction_batch.append({
                    'sequence': full_seq,
                    'barcode': 'ITS',
                    'length': len(full_seq),
                    'source_id': parent['source_id'],
                    'original_id': parent['original_id'],
                    'parent_seq_id': parent_id,
                    'tax_id': parent.get('tax_id')
                })

            # Insert reconstructed ITS sequences
            reconstructed_inserted = 0
            for i in range(0, len(reconstruction_batch), self.insert_batch_size):
                batch = reconstruction_batch[i:i + self.insert_batch_size]
                inserted = self._bulk_insert_sequences(batch)
                reconstructed_inserted += inserted

            self.stats.reconstructed_its = {
                'common_parents': len(common_ids),
                'inserted': reconstructed_inserted
            }
            total_inserted += reconstructed_inserted

            # Step 7: Mark no_detection sequences
            no_detection_file = self.barcode_dir / f"{self.barcode}_no_detections.fasta"
            no_detection_ids = []

            for header, _ in self._parse_fasta_streaming(no_detection_file):
                parent_id = self._parse_header_for_parent_id(header)
                if parent_id:
                    no_detection_ids.append(parent_id)

            if no_detection_ids:
                updates = [("no_detection", pid) for pid in no_detection_ids]
                self._bulk_update_status(updates)
                self.stats.no_detection = len(no_detection_ids)

            # Finalize stats
            self.stats.total_new_sequences = total_inserted
            self.stats.status = "Success"
            self.stats.end_time = datetime.now().isoformat()

            logger.info(f"=== Processing complete for {self.barcode} ===")
            logger.info(f"Total new sequences inserted: {total_inserted}")

        except Exception as e:
            logger.error(f"Processing failed for {self.barcode}: {e}", exc_info=True)
            self.stats.status = "Failed"
            self.stats.error_message = str(e)
            self.stats.end_time = datetime.now().isoformat()
            raise

        return self.stats


def main():
    """
    Main entry point - supports both Snakemake and standalone CLI modes.

    Snakemake mode: Uses global 'snakemake' object
    CLI mode: Uses command-line arguments
    """
    # Check if running via Snakemake
    if 'snakemake' in globals():
        # ===================================================================
        # SNAKEMAKE MODE
        # ===================================================================
        sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
        from db_manager import get_database_manager

        # Snakemake parameters
        barcode = snakemake.wildcards.barcode
        itsx_dir = Path(snakemake.params.itsx_dir)
        data_dir = Path(snakemake.params.data_dir)
        output_report = snakemake.output.report
        log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

        # Configure logging for Snakemake
        if log_file:
            handler = logging.FileHandler(log_file, mode='w')
            handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(handler)

        # Extract config
        # CENTRALIZED BATCH SIZE CONFIGURATION
        # All batch operations (inserts, updates, fetches) use the same batch_size
        # from config.yaml > database_formatting.batch_size (default: 1000)
        # This ensures consistency and makes performance tuning easier
        db_formatting = snakemake.config.get("database_formatting", {})
        seq_val_config = snakemake.config.get("sequence_validation", {})

        # Universal batch size for all database operations
        batch_size = db_formatting.get("batch_size", 1000)

        config = {
            "min_length": seq_val_config.get("min_length", 50),
            "max_n_fraction": seq_val_config.get("max_n_fraction", 0.05),
            "valid_chars": seq_val_config.get("valid_chars", "ACGTURYSWKMBDHVN-"),
            # All batch operations use the same size from database_formatting.batch_size
            "insert_batch_size": batch_size,
            "status_update_chunk_size": batch_size,
            "parent_metadata_chunk_size": batch_size,
        }

        try:
            db_manager = get_database_manager(snakemake.config)

            # Adjust barcode_dir to point to the barcode subdirectory
            barcode_subdir = itsx_dir / barcode

            processor = BarcodeITSxProcessor(barcode, barcode_subdir, data_dir, config, db_manager)
            stats = processor.process()

            # Write report
            Path(output_report).parent.mkdir(parents=True, exist_ok=True)
            with open(output_report, 'w') as f:
                json.dump(stats.to_dict(), f, indent=2)

            logger.info(f"Report written: {output_report}")

        except Exception as e:
            logger.critical(f"Critical ITSx error: {e}", exc_info=True)
            error_report = {
                "barcode": barcode,
                "status": "Critical Error",
                "error_message": str(e),
                "end_time": datetime.now().isoformat()
            }
            Path(output_report).write_text(json.dumps(error_report, indent=2))
            sys.exit(1)

    else:
        # ===================================================================
        # STANDALONE CLI MODE (for testing)
        # ===================================================================
        if len(sys.argv) != 5:
            print("Usage: python process_single_barcode_itsx.py <barcode> <barcode_dir> <config_yaml> <output_json>")
            sys.exit(1)

        import yaml
        sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
        from db_manager import get_database_manager

        barcode = sys.argv[1]
        barcode_dir = Path(sys.argv[2])
        config_yaml = sys.argv[3]
        output_json = sys.argv[4]

        # For standalone mode, data_dir is same as barcode_dir parent (../../data)
        data_dir = barcode_dir.parent.parent / "data"

        # Load configuration
        with open(config_yaml, 'r') as f:
            full_config = yaml.safe_load(f)

        # CENTRALIZED BATCH SIZE CONFIGURATION
        # All batch operations (inserts, updates, fetches) use the same batch_size
        # from config.yaml > database_formatting.batch_size (default: 1000)
        # This ensures consistency and makes performance tuning easier
        db_formatting = full_config.get("database_formatting", {})
        seq_val_config = full_config.get("sequence_validation", {})

        # Universal batch size for all database operations
        batch_size = db_formatting.get("batch_size", 1000)

        config = {
            "min_length": seq_val_config.get("min_length", 50),
            "max_n_fraction": seq_val_config.get("max_n_fraction", 0.05),
            "valid_chars": seq_val_config.get("valid_chars", "ACGTURYSWKMBDHVN-"),
            # All batch operations use the same size from database_formatting.batch_size
            "insert_batch_size": batch_size,
            "status_update_chunk_size": batch_size,
            "parent_metadata_chunk_size": batch_size,
        }

        try:
            db_manager = get_database_manager(full_config)
            processor = BarcodeITSxProcessor(barcode, barcode_dir, data_dir, config, db_manager)
            stats = processor.process()

            # Write output JSON
            output_path = Path(output_json)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, 'w') as f:
                json.dump(stats.to_dict(), f, indent=2)

            logger.info(f"Results written to {output_json}")

            # Exit with appropriate code
            sys.exit(0 if stats.status == "Success" else 1)

        except Exception as e:
            logger.critical(f"Critical error: {e}", exc_info=True)
            sys.exit(1)


if __name__ == "__main__":
    main()
