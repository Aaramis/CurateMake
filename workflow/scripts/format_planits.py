#!/usr/bin/env python3
"""
PLANITS Parser
"""

import os
import sys
import zipfile
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
from io import TextIOWrapper

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants
FASTA_CHUNK_SIZE = 10000    # Increased to reduce overhead
QUEUE_SIZE = 20             # More buffer to avoid stalls
MIN_SEQUENCE_LENGTH = 50    # Early filtering
DEFAULT_BATCH_SIZE = 1000   # Default batch size if not provided in config

@dataclass
class SequenceRecord:
    """Data container for PLANITS sequences."""
    original_id: str
    sequence: str
    length: int
    barcode: str  # ITS, ITS1, ou ITS2
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

class PLANITSParser:
    """PLANITS parser and ingestion orchestrator."""

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
            'processing_time': 0.0,
            'its_ITS': 0, 
            'its_ITS1': 0,
            'its_ITS2': 0
        })
        
        # Pre-compiled regex
        self.version_pattern = re.compile(r'(\d{1,2}\.\d{1,2}\.\d{4})')
        
        # Precompute allowed characters
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)

    @lru_cache(maxsize=50000)  # Larger cache
    def _parse_planits_taxonomy(self, header: str, barcode_type: str) -> Tuple[str, Tuple]:
        """Parse PLANITS headers in SINTAX format.
        
        Format: >MK696133.1;tax=p:Chlorophyta,c:Chlorophyceae,o:Chlamydomonadales,f:NA,g:NA,s:unknown Chlamydomonadales ;

        Note that the header does not contain kingdom information.
        """
        try:
            header = header.lstrip('>')
            
            # Split the ID from taxonomy
            if ';tax=' in header:
                seq_id_part, tax_part = header.split(';tax=', 1)
                seq_id = seq_id_part.strip()
                seq_id = seq_id + f"_{barcode_type}"
                
                # Extract the taxonomy (before the trailing ';')
                tax_string = tax_part.split(';')[0]
                
                # Parse the taxonomy
                ranks = {}
                for entry in tax_string.split(','):
                    if ':' in entry:
                        rank_code, value = entry.split(':', 1)
                        rank_code = rank_code.strip()
                        value = value.strip()
                        
                        # Map rank codes to full names
                        rank_mapping = {
                            'k': 'kingdom',
                            'p': 'phylum', 
                            'c': 'class_',
                            'o': 'order_',
                            'f': 'family',
                            'g': 'genus',
                            's': 'species'
                        }
                        
                        rank_name = rank_mapping.get(rank_code, rank_code)
                        
                        # Clean the value
                        if value and value.lower() not in ['na', 'unknown', '']:
                            # Handle the special case for species
                            if rank_name == 'species' and ' ' in value:
                                # Keep only the species name (after the genus)
                                parts = value.split()
                                if len(parts) >= 2:
                                    value = ' '.join(parts[:2])  # Genus + species
                            
                            ranks[rank_name] = value
                        else:
                            ranks[rank_name] = None
                
                return seq_id, tuple(sorted(ranks.items()))
            
            else:
                # No taxonomy found
                return header.strip(), ()
                
        except Exception as e:
            logger.debug(f"Error parsing header {header[:100]}...: {e}")
            self.stats['parsing_errors'] += 1
            return header.lstrip('>').strip(), ()

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
        
        # Character validation (limited to first 100 chars for speed)
        if not all(c.upper() in self.valid_chars for c in sequence[:100]):
            return False, seq_len, n_fraction
        
        return True, seq_len, n_fraction

    def extract_barcode_type(self, filename: str) -> str:
        """Extract the barcode type from the filename."""
        filename_lower = filename.lower()
        if 'its1.sintax' in filename_lower:
            return 'ITS1'
        elif 'its2.sintax' in filename_lower:
            return 'ITS2'
        elif 'its.sintax' in filename_lower:
            return 'ITS'
        else:
            logger.warning(f"Unknown barcode type for file: {filename}")
            return 'ITS'  # Default

    def fasta_reader_zip(self, zip_path: str, queue: Queue):
        """FASTA reader dedicated to PLANITS ZIP archives."""
        logger.info(f"Starting ZIP reader for {zip_path}")
        
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_file:
                # Look for *.SINTAX_format.fas files
                sintax_files = [f for f in zip_file.namelist() 
                              if f.endswith('.SINTAX_format.fas')]
                
                logger.info(f"Found SINTAX files: {sintax_files}")
                
                for sintax_file in sintax_files:
                    logger.info(f"Processing file: {sintax_file}")
                    barcode_type = self.extract_barcode_type(sintax_file)
                    
                    with zip_file.open(sintax_file) as file:
                        # Wrap as a text stream
                        text_file = TextIOWrapper(file, encoding='utf-8')
                        self._process_fasta_stream(text_file, barcode_type, queue)
                        
        except Exception as e:
            logger.error(f"Error in ZIP reader: {e}", exc_info=True)
        finally:
            queue.put(None)  # End signal
            logger.info("ZIP reader finished")

    def _process_fasta_stream(self, stream, barcode_type: str, queue: Queue):
        """Process a FASTA stream."""
        header, sequence_parts = None, []
        chunk = []
        
        for line in stream:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if header and sequence_parts:
                    sequence = ''.join(sequence_parts).replace('-', '')  # Remove alignment characters
                    self._add_to_chunk(header, sequence, barcode_type, chunk, queue)
                
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)
        
        # Handle final record
        if header and sequence_parts:
            sequence = ''.join(sequence_parts).replace('-', '')
            self._add_to_chunk(header, sequence, barcode_type, chunk, queue)
        
        # Send any remaining chunk
        if chunk:
            queue.put(chunk[:])

    def _add_to_chunk(self, header: str, sequence: str, barcode_type: str, chunk: List, queue: Queue):
        """Add a validated sequence to the chunk."""
        # Validate early to skip useless work
        is_valid, seq_len, n_frac = self.validate_sequence(sequence)
        
        if is_valid:
            chunk.append((header, sequence, seq_len, barcode_type))
            self.stats['total_sequences'] += 1
            self.stats[f'its_{barcode_type}'] += 1
            
            if len(chunk) >= FASTA_CHUNK_SIZE:
                queue.put(chunk[:])  # Copy the chunk
                chunk.clear()
        else:
            self.stats['filtered_sequences'] += 1

    def sequence_processor_worker(self, chunk: List[Tuple[str, str, int, str]]) -> List[SequenceRecord]:
        """Worker leveraging efficient data structures."""
        results = []
        
        for header, sequence, seq_len, barcode_type in chunk:
            seq_id, taxonomy_tuple = self._parse_planits_taxonomy(header, barcode_type)
            tax_dict = dict(taxonomy_tuple) if taxonomy_tuple else {}
            
            # Build the SequenceRecord directly
            record = SequenceRecord(
                original_id=seq_id,
                sequence=sequence,
                length=seq_len,
                barcode=barcode_type,
                kingdom="Plantae",  # Forced value
                phylum=tax_dict.get('phylum'),
                class_=tax_dict.get('class_'),
                order_=tax_dict.get('order_'),
                genus=tax_dict.get('genus'),
                species=tax_dict.get('species')
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

    def process_planits_database(self, input_file: str, source_name: str = "PLANITS", 
                                source_version: str = "unknown") -> Dict:
        """PLANITS processing pipeline."""
        start_time = time.time()
        logger.info(f"Starting PLANITS processing: {input_file}")
        
        # Ensure we are dealing with a ZIP file
        if not input_file.endswith('.zip'):
            raise ValueError(f"Expected ZIP file, got: {input_file}")
        
        # Configure DB for peak performance
        self.db_manager.optimize_for_bulk_insert_advanced()
        source_id = self.db_manager.add_source(name=source_name, file_path=input_file)
        
        # Queues
        fasta_queue = Queue(maxsize=QUEUE_SIZE)
        results_queue = Queue(maxsize=QUEUE_SIZE)
        
        # Compute optimal worker count
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
                    self.fasta_reader_zip, input_file, fasta_queue
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
        processing_time = self.stats['processing_time']
        its_ITS = self.stats['its_ITS']
        its_ITS1 = self.stats['its_ITS1']
        its_ITS2 = self.stats['its_ITS2']
        
        logger.info(f"Processing completed in {processing_time:.2f}s")
        logger.info(f"Total: {total_seqs}, Processed: {processed_seqs}, "
                   f"Filtered: {filtered_seqs}, Rate: {processed_seqs/processing_time:.0f} seqs/sec")
        logger.info(f"ITS1: {its_ITS1}, ITS2: {its_ITS2}, ITS: {its_ITS}")
        
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
            input_file = "PLANITS_small.zip"
            input_db = "test_output.sqlite"
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
        parser = PLANITSParser(db_manager, config_validation, database_formatting)
        
        # Detect version
        version_match = parser.version_pattern.search(Path(input_file).name)
        version = version_match.group(1) if version_match else "unknown"
        
        # Processing
        stats = parser.process_planits_database(
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
            'source': 'PLANITS',
            'version': version,
            'input_file': input_file,
            'stats': dict(stats),
            'processing_rate_per_second': rate,
            'its_ITS': stats.get('its_ITS', 0),
            'its_ITS1': stats.get('its_ITS1', 0),
            'its_ITS2': stats.get('its_ITS2', 0),
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
        logger.critical(f"Fatal error in PLANITS processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
