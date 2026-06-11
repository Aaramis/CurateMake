#!/usr/bin/env python3
"""
UNITE Parser
"""

import os
import sys
import gzip
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
logger = logging.getLogger(__name__)

# Constants
FASTA_CHUNK_SIZE = 10000  # Increased to reduce overhead
QUEUE_SIZE = 20         # More buffer to avoid stalls
MIN_SEQUENCE_LENGTH = 50  # Early filtering
DEFAULT_BATCH_SIZE = 1000  # Default batch size if not provided in config

@dataclass
class SequenceRecord:
    """Data structure representing a UNITE sequence."""
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

class UNITEParser:
    """UNITE parser and ingestion pipeline."""

    def __init__(self, db_manager, config_validation: Dict, database_formatting: Optional[Dict] = None):
        self.db_manager = db_manager
        self.config_validation = config_validation

        # Read batch_size from config
        formatting_config = database_formatting or {}
        self.batch_size = formatting_config.get('batch_size', DEFAULT_BATCH_SIZE)
        
        # Thread-safe statistics
        self.stats = mp.Manager().dict({
            'total_sequences': 0,
            'processed_sequences': 0,
            'filtered_sequences': 0,
            'parsing_errors': 0,
            'processing_time': 0.0
        })
        
        # Pre-compiled regex
        self.version_pattern = re.compile(r'(\d{1,2}\.\d{1,2}\.\d{4})')
        
        # Precompute allowed characters
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)

    @lru_cache(maxsize=50000)  # Larger cache
    def _parse_taxonomy(self, header: str) -> Tuple[str, Tuple]:
        """Parse taxonomy information via string splitting rather than regex."""
        try:
            header = header.lstrip('>')
            parts = header.split('|', 2)
            seq_id = parts[0].strip()

            if len(parts) < 2:
                return seq_id, ()

            taxo_string = parts[1]
            ranks = {}
            for entry in taxo_string.split(';'):
                if '__' in entry:
                    rank, value = entry.split('__', 1)
                    rank = {
                        "k": "kingdom",
                        "p": "phylum",
                        "c": "class_",
                        "o": "order_",
                        "f": "family",
                        "g": "genus",
                        "s": "species",
                    }.get(rank[0], rank)  # map short → long

                    clean_val = value.replace('_', ' ').strip()
                    ranks[rank] = clean_val if clean_val and clean_val != "unidentified" else None

            return seq_id, tuple(sorted(ranks.items()))

        except Exception:
            self.stats['parsing_errors'] += 1
            return header.lstrip('>'), ()

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

    def fasta_reader_mmap(self, file_path: str, queue: Queue):
        """FASTA reader optimized for large files."""
        logger.info(f"Starting reader for {file_path}")
        chunk = []
        
        try:
            with gzip.open(file_path, 'rt', encoding='utf-8') as f:
                self._process_fasta_stream(f, chunk, queue)
        except Exception as e:
            logger.error(f"Error in FASTA reader: {e}", exc_info=True)
        finally:
            if chunk:
                queue.put(chunk)
            queue.put(None)  # End signal
            logger.info("Reader finished")

    def _process_fasta_stream(self, stream, chunk, queue):
        """Process a FASTA stream sequentially."""
        header, sequence_parts = None, []
        
        for line in stream:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if header and sequence_parts:
                    sequence = ''.join(sequence_parts).replace('-', '')  # Remove alignment marks
                    self._add_to_chunk(header, sequence, chunk, queue)
                
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)
        
        # Handle final record
        if header and sequence_parts:
            sequence = ''.join(sequence_parts).replace('-', '')  # Remove alignment marks
            self._add_to_chunk(header, sequence, chunk, queue)

    def _add_to_chunk(self, header: str, sequence: str, chunk: List, queue: Queue):
        """Add a validated sequence to the chunk."""
        # Validate early to skip useless work
        is_valid, seq_len, n_frac = self.validate_sequence(sequence)
        
        if is_valid:
            chunk.append((header, sequence, seq_len, n_frac))
            self.stats['total_sequences'] += 1
            
            if len(chunk) >= FASTA_CHUNK_SIZE:
                queue.put(chunk[:])  # Copy the chunk
                chunk.clear()
        else:
            self.stats['filtered_sequences'] += 1

    def sequence_processor_worker(self, chunk: List[Tuple[str, str, int, float]]) -> List[SequenceRecord]:
        """Worker that leverages efficient data structures."""
        results = []
        
        for header, sequence, seq_len, n_frac in chunk:
            seq_id, taxonomy_tuple = self._parse_taxonomy(header)
            tax_dict = dict(taxonomy_tuple) if taxonomy_tuple else {}
            
            # Build the SequenceRecord directly
            record = SequenceRecord(
                original_id=seq_id,
                sequence=sequence,
                length=seq_len,
                n_frac=n_frac,
                **{k: tax_dict.get(k) for k in ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']}
            )
            
            results.append(record)
        
        return results

    def db_writer_consumer(self, source_id: str, results_queue: Queue):
        """Writer avec SQLAlchemy bulk inserts (source_id is UUID string)."""
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

            # Traiter le dernier batch
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
        taxonomy_map = {}  # record -> tax_tuple

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

        # 2. Bulk insert taxonomies and fetch tax_ids
        tax_id_map = self.db_manager.bulk_insert_taxonomies(taxonomy_data)

        # 3. Prepare sequence rows
        sequence_data = []
        for record in batch:
            seq_dict = {
                'source_id': source_id,
                'original_id': record.original_id,
                'sequence': record.sequence,
                'length': record.length,
                'percentage_n': record.n_frac,
                'barcode': 'ITS',
                'barcode_status': 'Pending'
            }
            sequence_data.append(seq_dict)

        # 4. Bulk insert sequences and fetch seq_ids
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

    def process_unite_database(self, input_file: str, source_name: str = "UNITE", 
                                       source_version: str = "unknown") -> Dict:
        """UNITE processing pipeline."""
        start_time = time.time()
        logger.info(f"Starting UNITE processing: {input_file}")
        
        # Configure DB for peak performance
        self.db_manager.optimize_for_bulk_insert_advanced()
        source_id = self.db_manager.add_source(name=source_name, file_path=input_file)
        
        # Queues
        fasta_queue = Queue(maxsize=QUEUE_SIZE)
        results_queue = Queue(maxsize=QUEUE_SIZE)
        
        # Compute the optimal number of workers
        cpu_count = os.cpu_count()
        num_workers = max(2, cpu_count - 2)  # Keep 2 CPUs for reader/writer
        logger.info(f"Using {num_workers} worker threads on {cpu_count} CPUs")
        
        try:
            with ThreadPoolExecutor(max_workers=num_workers + 2) as executor:
                # Launch reader and writer
                writer_future = executor.submit(
                    self.db_writer_consumer, source_id, results_queue
                )
                reader_future = executor.submit(
                    self.fasta_reader_mmap, input_file, fasta_queue
                )
                
                # Process chunks in parallel
                worker_futures = []
                processed_chunks = 0
                
                while True:
                    try:
                        chunk = fasta_queue.get(timeout=120)  # Timeout plus long
                        if chunk is None:
                            break
                        
                        # Submit for processing
                        future = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(future)
                        processed_chunks += 1
                        
                        # Stream processing to avoid memory spikes
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
        processing_time = self.stats['processing_time']
        
        logger.info(f"Processing completed in {processing_time:.2f}s")
        logger.info(f"Total: {total_seqs}, Processed: {processed_seqs}, "
                   f"Filtered: {filtered_seqs}, Rate: {processed_seqs/processing_time:.0f} seqs/sec")
        
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
            database_formatting = snakemake.params.database_formatting
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            # Test mode
            logger.warning("Running in standalone test mode")
            input_file = "test_unite.fasta.gz"
            input_db = "test_output.sqlite"
            log_file = "test_format_unite.log"
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

        # Parser
        parser = UNITEParser(db_manager, config_validation, database_formatting)
        
        # Detect version
        version_match = parser.version_pattern.search(Path(input_file).name)
        version = version_match.group(1) if version_match else "unknown"
        
        # Processing
        stats = parser.process_unite_database(
            input_file=input_file, 
            source_version=version
        )
        
        # Final optimization
        db_manager.vacuum_database()
        
        # Report generation
        processing_time = stats.get('processing_time', 1)
        processed_seqs = stats.get('processed_sequences', 0)
        rate = processed_seqs / processing_time if processing_time > 0 else 0
        
        report_data = {
            'source': 'UNITE',
            'version': version,
            'input_file': input_file,
            'stats': dict(stats),
            'processing_rate_per_second': rate,
            'performance_metrics': {
                'cpu_cores_used': os.cpu_count(),
                'memory_optimizations': True,
                'batch_size': parser.batch_size,
                'chunk_size': FASTA_CHUNK_SIZE
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
        logger.critical(f"Fatal error in UNITE processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
