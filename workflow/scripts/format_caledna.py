#!/usr/bin/env python3
"""
CALeDNA Parser
Parses CALeDNA archives (.zip containing .fasta and taxonomy .txt files).
"""

import os
import sys
import re
import time
import mmap
import logging
import json
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from queue import Queue, Empty
import multiprocessing as mp
from dataclasses import dataclass
import tempfile
import shutil
from sqlalchemy.exc import OperationalError

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# Shared taxonomy normaliser (single source of truth for placeholder/species cleaning)
for _u in (Path(__file__).resolve().parent / "utils",
           Path.cwd() / "workflow" / "scripts" / "utils"):
    if (_u / "taxonomy_normalize.py").exists():
        if str(_u) not in sys.path:
            sys.path.insert(0, str(_u))
        break
from taxonomy_normalize import clean_rank, clean_species, is_domain

logger = logging.getLogger(__name__)

# Constants
FASTA_CHUNK_SIZE = 10000
DEFAULT_BATCH_SIZE = 2000
DEFAULT_RETRY_ATTEMPTS = 3
DEFAULT_RETRY_DELAY = 5
DEFAULT_RETRY_BACKOFF = 2.0
QUEUE_SIZE = 20
MIN_SEQUENCE_LENGTH = 50

@dataclass
class SequenceRecord:
    """Data structure representing a normalized sequence entry."""
    original_id: str
    sequence: str
    length: int
    n_frac: float
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

class CALeDNAParser:
    """CALeDNA parser and ingestion pipeline."""

    def __init__(self, db_manager, config_validation: Dict, database_formatting: Optional[Dict] = None):
        self.db_manager = db_manager
        self.config_validation = config_validation

        # Read batch_size from database_formatting config (consistent with other parsers)
        formatting_config = database_formatting or {}
        batch_size_cfg = formatting_config.get('batch_size', DEFAULT_BATCH_SIZE)
        try:
            self.batch_size = max(1, int(batch_size_cfg))
        except (TypeError, ValueError):
            logger.warning(f"Invalid batch_size '{batch_size_cfg}', falling back to {DEFAULT_BATCH_SIZE}")
            self.batch_size = DEFAULT_BATCH_SIZE

        # Read writer_retry config
        writer_retry_cfg = formatting_config.get('writer_retry', {})
        try:
            self.flush_retry_attempts = max(0, int(writer_retry_cfg.get('max_attempts', DEFAULT_RETRY_ATTEMPTS)))
        except (TypeError, ValueError):
            logger.warning("Invalid writer_retry.max_attempts, defaulting to %s", DEFAULT_RETRY_ATTEMPTS)
            self.flush_retry_attempts = DEFAULT_RETRY_ATTEMPTS
        
        # Pre-compiled regex patterns
        self.version_pattern = re.compile(r'(\d{1,2}\.\d{1,2}|\d{4})')
        
        # Pre-compute allowed characters
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)

    def _reset_for_new_file(self):
        """Reset parser state before processing a new file."""
        # Thread-safe statistics for this file
        self.stats = mp.Manager().dict({
            'total_sequences': 0,
            'processed_sequences': 0,
            'filtered_sequences': 0,
            'parsing_errors': 0,
            'taxonomy_mapping_errors': 0,
            'processing_time': 0.0
        })
        
        # Fresh taxonomy mapping dictionary (thread-safe)
        self.taxonomy_map = {}
        
        # Clear the LRU cache
        self._parse_caledna_taxonomy.cache_clear()

    def extract_zip_files(self, zip_path: str) -> Tuple[str, str]:
        """
        Extract FASTA and taxonomy .txt files from the zip archive
        and return their paths.
        """
        logger.info(f"Extracting files from {zip_path}")
        
        # Create a temporary folder
        temp_dir = tempfile.mkdtemp(prefix="caledna_extract_")
        
        fasta_file = None
        taxonomy_file = None
        
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                # List every file in the archive
                file_list = zip_ref.namelist()
                logger.info(f"Files in archive: {file_list}")
                
                # Extract all files
                zip_ref.extractall(temp_dir)
                
                # Locate .fasta and .txt files
                for root, dirs, files in os.walk(temp_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        if file.lower().endswith(('.fasta', '.fas', '.fa')):
                            fasta_file = file_path
                            logger.info(f"Found FASTA file: {fasta_file}")
                        elif file.lower().endswith('.txt') and 'tax' in file.lower():
                            taxonomy_file = file_path
                            logger.info(f"Found taxonomy file: {taxonomy_file}")
                        elif file.lower().endswith('.txt') and not taxonomy_file:
                            # If no file contains 'tax', grab the first .txt file
                            taxonomy_file = file_path
                            logger.info(f"Found potential taxonomy file: {taxonomy_file}")
                
                if not fasta_file:
                    raise FileNotFoundError("No FASTA file found in the archive")
                if not taxonomy_file:
                    raise FileNotFoundError("No taxonomy file found in the archive")
                
        except Exception as e:
            # Clean up on error
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
            raise e
        
        return fasta_file, taxonomy_file, temp_dir

    def load_taxonomy_mapping(self, taxonomy_file: str):
        """
        Load the taxonomy mapping from the .txt file.
        Expected format: ID\tTaxonomy;separated;by;semicolons
        """
        logger.info(f"Loading taxonomy mapping from {taxonomy_file}")
        
        taxonomy_count = 0
        
        try:
            with open(taxonomy_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        # Split ID and taxonomy (tab or space separated)
                        if '\t' in line:
                            seq_id, taxonomy_str = line.split('\t', 1)
                        elif ' ' in line:
                            parts = line.split(' ', 1)
                            if len(parts) == 2:
                                seq_id, taxonomy_str = parts
                            else:
                                continue
                        else:
                            continue
                        
                        seq_id = seq_id.strip()
                        taxonomy_str = taxonomy_str.strip()
                        
                        # Parse taxonomy string
                        taxonomy = self._parse_caledna_taxonomy(taxonomy_str)
                        self.taxonomy_map[seq_id] = taxonomy
                        taxonomy_count += 1
                        
                        if taxonomy_count % 10000 == 0:
                            logger.info(f"Loaded {taxonomy_count} taxonomy mappings")
                            
                    except Exception as e:
                        logger.warning(f"Error parsing taxonomy line {line_num}: {line[:50]}... - {e}")
                        self.stats['taxonomy_mapping_errors'] += 1
                        continue
        
        except Exception as e:
            logger.error(f"Error loading taxonomy file: {e}")
            raise
        
        logger.info(f"Successfully loaded {taxonomy_count} taxonomy mappings")

    @lru_cache(maxsize=50000)
    def _parse_caledna_taxonomy(self, taxonomy_str: str) -> Dict[str, Optional[str]]:
        """
        Parse a CALeDNA taxonomy string.
        Example format: Eukaryota;Ascomycota;Sordariomycetes;Microascales;Microascaceae;Microascus;Microascus cirrosus
        """
        try:
            ranks = ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']
            taxonomy = {}
            
            if not taxonomy_str or taxonomy_str.strip() == '':
                return {rank: None for rank in ranks}
            
            # Split on semicolons. CALeDNA lineages are fixed 7-field and positional,
            # with missing intermediate ranks kept as NA in place (never omitted), so
            # positional mapping is safe. Guard against an unexpected layout.
            tax_parts = [part.strip() for part in taxonomy_str.split(';')]
            if len(tax_parts) != len(ranks):
                logger.warning(
                    f"CALeDNA lineage has {len(tax_parts)} fields (expected {len(ranks)}), "
                    f"positional mapping may be off: {taxonomy_str[:120]!r}"
                )

            # Map onto taxonomy ranks. Field 0 is a DOMAIN (always "Eukaryota"), not a
            # kingdom; storing a domain there is inconsistent with the other sources and
            # is a no-op for CoL's kingdom-based disambiguation, so we leave kingdom empty
            # and let the CoL harmonisation fill the real kingdom from genus/species.
            for i, rank in enumerate(ranks):
                raw = tax_parts[i] if i < len(tax_parts) else None
                if rank == 'kingdom':
                    taxonomy[rank] = None if is_domain(raw) else clean_rank(raw)
                elif rank == 'species':
                    taxonomy[rank] = clean_species(raw)
                else:
                    taxonomy[rank] = clean_rank(raw)

            return taxonomy
            
        except Exception as e:
            logger.warning(f"Error parsing taxonomy '{taxonomy_str}': {e}")
            return {rank: None for rank in ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']}

    def validate_sequence(self, sequence: str) -> Tuple[bool, int, float]:
        """Early-exit sequence validation routine."""
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
        
        # Validate characters (set membership is faster)
        if not all(c.upper() in self.valid_chars for c in sequence[:100]):
            return False, seq_len, n_fraction
        
        return True, seq_len, n_fraction

    def fasta_reader_mmap(self, file_path: str, queue: Queue):
        """
        FASTA reader.
        Iterate sequence by sequence and use _add_to_chunk to create batches of
        tuples (header, sequence, len, n_frac) with the expected format.
        """
        logger.info(f"Starting FASTA reader for {file_path}")
        
        batch = []
        header, sequence_parts = None, []

        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        if header and sequence_parts:
                            sequence = ''.join(sequence_parts)
                            self._add_to_chunk(header, sequence, batch, queue)
                        
                        header = line
                        sequence_parts = []
                    elif header:
                        sequence_parts.append(line)

                if header and sequence_parts:
                    sequence = ''.join(sequence_parts)
                    self._add_to_chunk(header, sequence, batch, queue)

        except Exception as e:
            logger.error(f"Error in FASTA reader: {e}", exc_info=True)
        finally:
            if batch:
                queue.put(batch)
            queue.put(None)
            logger.info("FASTA reader finished and sent sentinel.")

    def _process_fasta_stream(self, stream, chunk, queue):
        """Process a FASTA stream sequentially."""
        header, sequence_parts = None, []
        
        for line in stream:
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
        
        if header and sequence_parts:
            sequence = ''.join(sequence_parts).replace('-', '')
            self._add_to_chunk(header, sequence, chunk, queue)

    def _add_to_chunk(self, header: str, sequence: str, chunk: List, queue: Queue):
        """Add a sequence to the chunk after early validation."""
        is_valid, seq_len, n_frac = self.validate_sequence(sequence)
        
        if is_valid:
            chunk.append((header, sequence, seq_len, n_frac))
            self.stats['total_sequences'] += 1
            
            if len(chunk) >= FASTA_CHUNK_SIZE:
                queue.put(chunk[:])
                chunk.clear()
        else:
            self.stats['filtered_sequences'] += 1

    def sequence_processor_worker(self, chunk: List[Tuple[str, str, int, float]]) -> List[SequenceRecord]:
        """Worker that processes a chunk of sequences."""
        results = []
        
        for header, sequence, seq_len, n_frac in chunk:
            seq_id = header.lstrip('>').split()[0]  # First token after >
            
            taxonomy = self.taxonomy_map.get(seq_id, {})
            
            if not taxonomy:
                logger.debug(f"No taxonomy found for sequence ID: {seq_id}")
                taxonomy = {rank: None for rank in ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']}
            
            # Build the SequenceRecord
            record = SequenceRecord(
                original_id=seq_id,
                sequence=sequence,
                length=seq_len,
                n_frac=n_frac,
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
        """Writer using SQLAlchemy bulk inserts (source_id is UUID string)."""
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
                        self._flush_with_retry(batch, source_id)
                        self.stats['processed_sequences'] += len(batch)
                        logger.info(f"Processed {self.stats['processed_sequences']} sequences")
                        batch.clear()

                except Empty:
                    logger.warning("DB writer timeout")
                    break

            if batch:
                self._flush_with_retry(batch, source_id)
                self.stats['processed_sequences'] += len(batch)

            logger.info("DB writer finished successfully")

        except Exception as e:
            logger.error(f"Critical error in DB writer: {e}", exc_info=True)
            raise

    def _flush_with_retry(self, batch: List[SequenceRecord], source_id: str):
        """Flush batch with immediate retry for transient OperationalError failures."""
        attempts = 0

        while True:
            try:
                self._flush_batch_sqlalchemy(batch, source_id)
                return
            except OperationalError as exc:
                if attempts >= self.flush_retry_attempts:
                    logger.error("Failed to flush CALeDNA batch after %s attempts: %s", attempts, exc)
                    raise

                attempts += 1
                logger.warning(
                    "OperationalError flushing CALeDNA batch (attempt %s/%s). Retrying immediately... Error: %s",
                    attempts,
                    self.flush_retry_attempts,
                    exc
                )

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
                'barcode': 'ITS2',
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

    def process_caledna_database(self, zip_file: str, source_name: str = "CALeDNA") -> Dict:
        """Main CALeDNA processing pipeline."""
        
        # Reset for this file
        self._reset_for_new_file()

        start_time = time.time()
        logger.info(f"Starting CALeDNA processing: {zip_file}")
        
        temp_dir = None
        
        try:
            # Extract files
            fasta_file, taxonomy_file, temp_dir = self.extract_zip_files(zip_file)
            
            # Load taxonomy mapping
            self.load_taxonomy_mapping(taxonomy_file)
            
            # Configure DB for peak performance
            self.db_manager.optimize_for_bulk_insert_advanced()
            source_id = self.db_manager.add_source(name=source_name, file_path=zip_file)
            
            # Queues
            fasta_queue = Queue(maxsize=QUEUE_SIZE)
            results_queue = Queue(maxsize=QUEUE_SIZE)
            
            # Compute optimal worker count
            cpu_count = os.cpu_count()
            num_workers = max(2, cpu_count - 2)
            logger.info(f"Using {num_workers} worker threads on {cpu_count} CPUs")
            
            with ThreadPoolExecutor(max_workers=num_workers + 2) as executor:
                # Launch reader and writer
                writer_future = executor.submit(
                    self.db_writer_consumer, source_id, results_queue
                )
                reader_future = executor.submit(
                    self.fasta_reader_mmap, fasta_file, fasta_queue
                )
                
                # Process chunks in parallel
                worker_futures = []
                processed_chunks = 0
                
                while True:
                    try:
                        chunk = fasta_queue.get(timeout=120)
                        if chunk is None:
                            break
                        
                        future = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(future)
                        processed_chunks += 1
                        
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
            # Cleanup
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temporary directory: {temp_dir}")
        
        # Finalization
        end_time = time.time()
        self.stats['processing_time'] = end_time - start_time
        
        # Final statistics
        total_seqs = self.stats['total_sequences']
        processed_seqs = self.stats['processed_sequences']
        filtered_seqs = self.stats['filtered_sequences']
        processing_time = self.stats['processing_time']
        
        logger.info(f"Processing completed in {processing_time:.2f}s")
        logger.info(f"Total: {total_seqs}, Processed: {processed_seqs}, "
                   f"Filtered: {filtered_seqs}, Rate: {processed_seqs/processing_time:.0f} seqs/sec")
        
        return dict(self.stats)

# Entry point
def main():
    """Main entry point."""
    try:
        # Snakemake configuration
        if 'snakemake' in globals():
            input_file = (
                snakemake.input.file
                if hasattr(snakemake.input, "file")
                else snakemake.input[0]
            )
            output_report = snakemake.output[0]
            config_validation = snakemake.params.sequence_validation
            database_formatting = snakemake.params.database_formatting
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            # Test mode
            logger.warning("Running in standalone test mode")
            input_files = ["test_FITS.zip", "test_PITS.zip"] # Mettez vos fichiers de test ici
            input_db = "test_output.sqlite"
            log_file = "test_format_caledna.log"
            config = {'database': {'path': input_db}}
            config_validation = {
                'min_length': 50,
                'valid_chars': 'ATCGN',
                'max_n_fraction': 0.1
            }
            database_formatting = {}

        # Import and initialization
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)
        parser = CALeDNAParser(db_manager, config_validation, database_formatting)
        
        # Restore safe settings before processing each file
        db_manager.restore_safe_settings()
        
        
        stats = parser.process_caledna_database(
            zip_file=str(input_file), 
            source_name=f"CALeDNA_{Path(input_file).stem}"
        )
        
        # Final optimization
        db_manager.restore_safe_settings()
        db_manager.vacuum_database()
        
        # Report generation
        total_time = stats.get('processing_time', 1)
        total_processed = stats.get('processed_sequences', 0)
        rate = total_processed / total_time if total_time > 0 else 0
        
        report_data = {
            'source': 'CALeDNA_Multi',
            'input_file': input_file,
            'stats': stats,
            'processing_rate_per_second': rate,
            'performance_metrics': {
                'cpu_cores_used': os.cpu_count(),
                'memory_optimizations': True,
                'batch_size': parser.batch_size,
                'chunk_size': FASTA_CHUNK_SIZE,
            },
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        # Save the report
        if 'snakemake' in globals():
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Optimization report saved to: {output_report}")
        
        logger.info(f"Final processing rate: {rate:.0f} sequences/second")
        
    except Exception as e:
        logger.critical(f"Fatal error in CALeDNA processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
