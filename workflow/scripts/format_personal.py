#!/usr/bin/env python3
"""
PERSONAL Parser - Import user-provided FASTA files into the unified SQLAlchemy DB.

Supports FASTA files with headers in the format:
>ID;barcode=ITS;kingdom=Fungi;phylum=Basidiomycota;class=Agaricomycetes;order=Russulales;family=Russulaceae;genus=Russula;species=Russula_abietina

Also handles complex headers with:
- Spaces after the ID (before the first ';')
- Additional metadata after taxonomy (e.g., strand info, length, quality scores)
- Missing taxonomic ranks
- Various separators and formatting inconsistencies

Example complex header:
>4b6ed8c0-5a11-4103-405c-b686dc74c49b ;barcode=ITS;kingdom=Fungi;phylum=Basidiomycota;class=Agaricomycetes;order=Russulales;family=Russulaceae;genus=Russula;species=Russula_abietina,+strand,0-702 length=718 error-free_length=728 read_identity=94.681%
"""

import os
import sys
import re
import time
import logging
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from queue import Queue, Empty
import multiprocessing as mp
from dataclasses import dataclass

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# Shared taxonomy normaliser (single source of truth for placeholder/species cleaning)
for _u in (Path(__file__).resolve().parent / "utils",
           Path.cwd() / "workflow" / "scripts" / "utils"):
    if (_u / "taxonomy_normalize.py").exists():
        if str(_u) not in sys.path:
            sys.path.insert(0, str(_u))
        break
from taxonomy_normalize import clean_rank, clean_species

logger = logging.getLogger(__name__)

# Constants
FASTA_CHUNK_SIZE = 10000
BATCH_SIZE = 10000
QUEUE_SIZE = 20
MIN_SEQUENCE_LENGTH = 50

# Accepted barcodes
ACCEPTED_BARCODES = {"ITS", "ITS1", "ITS2"}

# Taxonomic rank names for parsing
TAXONOMIC_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


@dataclass
class SequenceRecord:
    """Data structure representing a parsed PERSONAL sequence."""
    original_id: str
    sequence: str
    length: int
    barcode: str
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None


class PERSONALParser:
    """PERSONAL data parser and ingestion pipeline."""

    def __init__(self, db_manager, config_validation: Dict):
        self.db_manager = db_manager
        self.config_validation = config_validation
        self.batch_size = BATCH_SIZE

        # Thread-safe statistics
        self.stats = mp.Manager().dict({
            'total_sequences': 0,
            'processed_sequences': 0,
            'filtered_sequences': 0,
            'barcode_filtered': 0,
            'parsing_errors': 0,
            'processing_time': 0.0,
            'barcode_ITS': 0,
            'barcode_ITS1': 0,
            'barcode_ITS2': 0
        })

        # Precompute allowed characters
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)

        # Pre-compiled regex for cleaning species names
        # Matches trailing metadata like ",+strand,0-702" or " length=718"
        self.trailing_metadata_pattern = re.compile(r'[,\s][\+\-]?(?:strand|length|error|read)[^\;]*$', re.IGNORECASE)

    @lru_cache(maxsize=50000)
    def _parse_personal_header(self, header: str) -> Tuple[str, str, Dict[str, Optional[str]]]:
        """
        Parse a PERSONAL FASTA header.

        Expected format:
        >ID;barcode=ITS;kingdom=Fungi;phylum=Basidiomycota;...;species=Species_name

        Handles:
        - Spaces after the ID
        - Extra metadata after taxonomy fields
        - Missing taxonomic ranks

        Args:
            header: The FASTA header line (starting with '>')

        Returns:
            Tuple of (seq_id, barcode, taxonomy_dict)
        """
        try:
            # Remove '>' prefix
            header = header.lstrip('>')

            # Split by semicolon to get all fields
            parts = header.split(';')

            # First part is the ID (may have trailing spaces)
            seq_id = parts[0].strip()

            # Parse key=value pairs from remaining parts
            barcode = 'ITS'  # Default barcode
            taxonomy = {
                'kingdom': None,
                'phylum': None,
                'class_': None,
                'order_': None,
                'family': None,
                'genus': None,
                'species': None
            }

            for part in parts[1:]:
                part = part.strip()
                if not part or '=' not in part:
                    continue

                # Split on first '=' only to handle values with '=' in them
                key_value = part.split('=', 1)
                if len(key_value) != 2:
                    continue

                key = key_value[0].strip().lower()
                value = key_value[1].strip()

                # Clean up the value - remove trailing metadata
                # E.g., "Russula_abietina,+strand,0-702" -> "Russula_abietina"
                if ',' in value:
                    value = value.split(',')[0].strip()
                if ' ' in value and key in TAXONOMIC_RANKS:
                    # For taxonomy fields, take only the first part if there's a space
                    # Exception: species names can have spaces (genus + epithet)
                    if key != 'species':
                        value = value.split()[0].strip()

                # Handle empty or placeholder values
                if not value or value.lower() in ['na', 'n/a', 'none', 'unknown', 'unidentified', '']:
                    value = None

                # Map to internal field names
                if key == 'barcode':
                    barcode = value.upper() if value else 'ITS'
                elif key == 'kingdom':
                    taxonomy['kingdom'] = self._clean_taxonomic_value(value)
                elif key == 'phylum':
                    taxonomy['phylum'] = self._clean_taxonomic_value(value)
                elif key == 'class':
                    taxonomy['class_'] = self._clean_taxonomic_value(value)
                elif key == 'order':
                    taxonomy['order_'] = self._clean_taxonomic_value(value)
                elif key == 'family':
                    taxonomy['family'] = self._clean_taxonomic_value(value)
                elif key == 'genus':
                    taxonomy['genus'] = self._clean_taxonomic_value(value)
                elif key == 'species':
                    taxonomy['species'] = self._clean_species_value(value)

            return seq_id, barcode, taxonomy

        except Exception as e:
            logger.debug(f"Error parsing header {header[:100]}...: {e}")
            self.stats['parsing_errors'] += 1
            return header.lstrip('>').split(';')[0].strip(), 'ITS', {
                'kingdom': None, 'phylum': None, 'class_': None,
                'order_': None, 'family': None, 'genus': None, 'species': None
            }

    def _clean_taxonomic_value(self, value: Optional[str]) -> Optional[str]:
        """Clean a non-species taxonomic field value (shared normaliser)."""
        return clean_rank(value)

    def _clean_species_value(self, value: Optional[str]) -> Optional[str]:
        """Clean a species field value (shared normaliser: binomial or None)."""
        return clean_species(value)

    def validate_sequence(self, sequence: str) -> Tuple[bool, int, float]:
        """Early-exit sequence validation."""
        if not sequence:
            return False, 0, 1.0

        seq_len = len(sequence)
        if seq_len < self.min_length:
            return False, seq_len, 0.0

        # Count Ns
        n_count = sequence.count('N') + sequence.count('n')
        n_fraction = n_count / seq_len

        if n_fraction > self.max_n_fraction:
            return False, seq_len, n_fraction

        # Character validation (only first 100 chars for speed)
        if not all(c.upper() in self.valid_chars for c in sequence[:100]):
            return False, seq_len, n_fraction

        return True, seq_len, n_fraction

    def fasta_reader(self, file_path: str, queue: Queue):
        """FASTA reader for PERSONAL files (uncompressed FASTA)."""
        logger.info(f"Starting FASTA reader for {file_path}")
        chunk = []

        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                header, sequence_parts = None, []

                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        if header and sequence_parts:
                            sequence = ''.join(sequence_parts).replace('-', '')
                            self._add_to_chunk(header, sequence, chunk, queue)

                        header = line
                        sequence_parts = []
                    else:
                        sequence_parts.append(line)

                # Handle final record
                if header and sequence_parts:
                    sequence = ''.join(sequence_parts).replace('-', '')
                    self._add_to_chunk(header, sequence, chunk, queue)

        except Exception as e:
            logger.error(f"Error in FASTA reader: {e}", exc_info=True)
        finally:
            if chunk:
                queue.put(chunk)
            queue.put(None)  # End signal
            logger.info("FASTA reader finished")

    def _add_to_chunk(self, header: str, sequence: str, chunk: List, queue: Queue):
        """Add a validated sequence to the chunk."""
        # Parse header to extract barcode
        _, barcode, _ = self._parse_personal_header(header)

        # Filter by barcode if not in accepted list
        if barcode not in ACCEPTED_BARCODES:
            self.stats['barcode_filtered'] += 1
            return

        # Validate sequence
        is_valid, seq_len, n_frac = self.validate_sequence(sequence)

        if is_valid:
            chunk.append((header, sequence, seq_len, barcode))
            self.stats['total_sequences'] += 1
            self.stats[f'barcode_{barcode}'] += 1

            if len(chunk) >= FASTA_CHUNK_SIZE:
                queue.put(chunk[:])  # Copy the chunk
                chunk.clear()
        else:
            self.stats['filtered_sequences'] += 1

    def sequence_processor_worker(self, chunk: List[Tuple[str, str, int, str]]) -> List[SequenceRecord]:
        """Worker that processes a chunk of sequences."""
        results = []

        for header, sequence, seq_len, barcode in chunk:
            seq_id, _, taxonomy = self._parse_personal_header(header)

            # Build the SequenceRecord
            record = SequenceRecord(
                original_id=seq_id,
                sequence=sequence,
                length=seq_len,
                barcode=barcode,
                kingdom=taxonomy.get('kingdom'),
                phylum=taxonomy.get('phylum'),
                class_=taxonomy.get('class_'),
                order_=taxonomy.get('order_'),
                family=taxonomy.get('family'),
                genus=taxonomy.get('genus'),
                species=taxonomy.get('species')
            )

            results.append(record)

        return results

    def db_writer_consumer(self, source_id: str, results_queue: Queue):
        """Writer with SQLAlchemy bulk inserts (source_id is UUID string)."""
        logger.info("Starting DB writer")
        batch = []

        try:
            while True:
                try:
                    result_chunk = results_queue.get(timeout=6000)
                    if result_chunk is None:
                        break

                    batch.extend(result_chunk)

                    if len(batch) >= self.batch_size:
                        self._flush_batch_sqlalchemy(batch, source_id)
                        self.stats['processed_sequences'] += len(batch)
                        logger.info(f"Processed {self.stats['processed_sequences']} sequences")
                        batch.clear()

                except Empty:
                    logger.warning("DB writer timeout")
                    break

            # Process final batch
            if batch:
                self._flush_batch_sqlalchemy(batch, source_id)
                self.stats['processed_sequences'] += len(batch)

            logger.info("DB writer finished successfully")

        except Exception as e:
            logger.error(f"DB writer error: {e}", exc_info=True)
            raise

    def _flush_batch_sqlalchemy(self, batch: List[SequenceRecord], source_id: str):
        """Flush batch using SQLAlchemy bulk operations (source_id is UUID string)."""
        # 1. Prepare taxonomy payload
        taxonomy_data = []
        taxonomy_map = {}

        for record in batch:
            tax_dict = {
                'kingdom': record.kingdom,
                'phylum': record.phylum,
                'class_': record.class_,
                'order_': record.order_,
                'family': record.family,
                'genus': record.genus,
                'species': record.species
            }
            taxonomy_data.append(tax_dict)
            tax_tuple = (record.kingdom, record.phylum, record.class_,
                         record.order_, record.family, record.genus, record.species)
            taxonomy_map[record.original_id] = tax_tuple

        # 2. Bulk insert taxonomies
        tax_id_map = self.db_manager.bulk_insert_taxonomies(taxonomy_data)

        # 3. Prepare sequence rows
        sequence_data = []
        for record in batch:
            seq_dict = {
                'source_id': source_id,
                'original_id': record.original_id,
                'sequence': record.sequence,
                'length': record.length,
                'barcode': record.barcode,
                'barcode_status': 'Pending'
            }
            sequence_data.append(seq_dict)

        # 4. Bulk insert sequences
        seq_id_map = self.db_manager.bulk_insert_sequences(sequence_data)

        # 5. Prepare seq_tax relations
        relation_data = []
        for record in batch:
            seq_id = seq_id_map.get(record.original_id)
            tax_tuple = taxonomy_map.get(record.original_id)
            tax_id = tax_id_map.get(tax_tuple)

            if seq_id and tax_id:
                relation_dict = {
                    'seq_id': seq_id,
                    'tax_id': tax_id,
                    'backbone': 'Raw',
                    'status': 'Pending'
                }
                relation_data.append(relation_dict)

        # 6. Bulk insert relations
        self.db_manager.bulk_insert_seq_tax(relation_data)

    def process_personal_database(self, input_file: str, source_name: str = "PERSONAL") -> Dict:
        """PERSONAL processing pipeline."""
        start_time = time.time()
        logger.info(f"Starting PERSONAL processing: {input_file}")

        # Configure DB for peak performance
        self.db_manager.optimize_for_bulk_insert_advanced()
        source_id = self.db_manager.add_source(name=source_name, file_path=input_file)

        # Queues
        fasta_queue = Queue(maxsize=QUEUE_SIZE)
        results_queue = Queue(maxsize=QUEUE_SIZE)

        # Compute optimal worker count
        cpu_count = os.cpu_count()
        num_workers = max(2, cpu_count - 2)
        logger.info(f"Using {num_workers} worker threads on {cpu_count} CPUs")

        try:
            with ThreadPoolExecutor(max_workers=num_workers + 2) as executor:
                # Launch reader and writer
                writer_future = executor.submit(
                    self.db_writer_consumer, source_id, results_queue
                )
                reader_future = executor.submit(
                    self.fasta_reader, input_file, fasta_queue
                )

                # Process chunks in parallel
                worker_futures = []
                processed_chunks = 0

                while True:
                    try:
                        chunk = fasta_queue.get(timeout=120)
                        if chunk is None:
                            break

                        # Submit for processing
                        future = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(future)
                        processed_chunks += 1

                        # Stream results to avoid memory pressure
                        if len(worker_futures) >= num_workers * 2:
                            completed_future = next(as_completed(worker_futures))
                            worker_futures.remove(completed_future)
                            results_queue.put(completed_future.result())

                    except Empty:
                        logger.warning("Reader timeout - assuming completion")
                        break

                # Handle remaining futures
                for future in as_completed(worker_futures):
                    results_queue.put(future.result())

                # Signal completion and wait
                results_queue.put(None)
                writer_future.result()
                reader_future.result()

                logger.info(f"Processed {processed_chunks} chunks total")

        except Exception as e:
            logger.error(f"Critical error in processing pipeline: {e}", exc_info=True)
            raise
        finally:
            # Restore default settings
            self.db_manager.restore_safe_settings()

        # Finalization
        end_time = time.time()
        self.stats['processing_time'] = end_time - start_time

        # Final statistics
        total_seqs = self.stats['total_sequences']
        processed_seqs = self.stats['processed_sequences']
        filtered_seqs = self.stats['filtered_sequences']
        barcode_filtered = self.stats['barcode_filtered']
        processing_time = self.stats['processing_time']

        logger.info(f"Processing completed in {processing_time:.2f}s")
        logger.info(f"Total: {total_seqs}, Processed: {processed_seqs}, "
                    f"Filtered (sequence): {filtered_seqs}, Filtered (barcode): {barcode_filtered}")
        logger.info(f"Rate: {processed_seqs / processing_time:.0f} seqs/sec")
        logger.info(f"Barcodes - ITS: {self.stats['barcode_ITS']}, "
                    f"ITS1: {self.stats['barcode_ITS1']}, ITS2: {self.stats['barcode_ITS2']}")

        return dict(self.stats)


# Entry point
def main():
    """Script entry point."""
    try:
        # Snakemake configuration
        if 'snakemake' in globals():
            input_file = (
                snakemake.input.source
                if hasattr(snakemake.input, "source")
                else snakemake.input[0]
            )
            output_report = snakemake.output[0]
            config_validation = snakemake.params.sequence_validation
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            # Test mode
            logger.warning("Running in standalone test mode")
            input_file = "data/raw/PERSONAL/personal.fasta"
            output_report = "test_personal_report.json"
            config = {
                'database': {
                    'type': 'sqlite',
                    'sqlite': {'path': 'data/db', 'name': 'its.db'}
                }
            }
            config_validation = {
                'min_length': 50,
                'valid_chars': 'ATCGN',
                'max_n_fraction': 0.1
            }

        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        # Parser
        parser = PERSONALParser(db_manager, config_validation)

        # Processing
        stats = parser.process_personal_database(
            input_file=input_file,
            source_name="PERSONAL"
        )

        # Final optimization
        db_manager.vacuum_database()

        # Report generation
        processing_time = stats.get('processing_time', 1)
        processed_seqs = stats.get('processed_sequences', 0)
        rate = processed_seqs / processing_time if processing_time > 0 else 0

        report_data = {
            'source': 'PERSONAL',
            'input_file': input_file,
            'accepted_barcodes': sorted(list(ACCEPTED_BARCODES)),
            'stats': dict(stats),
            'processing_rate_per_second': rate,
            'barcode_distribution': {
                'ITS': stats.get('barcode_ITS', 0),
                'ITS1': stats.get('barcode_ITS1', 0),
                'ITS2': stats.get('barcode_ITS2', 0)
            },
            'performance_metrics': {
                'cpu_cores_used': os.cpu_count(),
                'memory_optimizations': True,
                'batch_size': BATCH_SIZE,
                'chunk_size': FASTA_CHUNK_SIZE
            },
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }

        if 'snakemake' in globals():
            Path(output_report).parent.mkdir(parents=True, exist_ok=True)
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"PERSONAL processing report saved to: {output_report}")
        else:
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Standalone PERSONAL report saved to: {output_report}")

        logger.info(f"Final processing rate: {rate:.0f} sequences/second")

    except Exception as e:
        logger.critical(f"Fatal error in PERSONAL processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
