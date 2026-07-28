#!/usr/bin/env python3
"""
BOLD Parser - Import BOLD TSV-in-TAR archives into the unified SQLAlchemy DB.

- Read a tar.gz archive containing a BOLD TSV
- Stream-parse it with pandas chunks
- Filter on ITS / ITS1 / ITS2 markers
- Validate and clean the sequences
- Perform optimized bulk inserts through SQLAlchemy:
    - Taxonomy
    - Sequences
    - SeqTax (with silent deduplication, Option A)
- Produce processing statistics + JSON report (used by Snakemake)
"""

import os
import sys
import time
import logging
import json
import tarfile
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue, Empty
import multiprocessing as mp

import pandas as pd

# Global logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
# Shared taxonomy normaliser (single source of truth for placeholder/species cleaning)
for _u in (Path(__file__).resolve().parent / "utils",
           Path.cwd() / "workflow" / "scripts" / "utils"):
    if (_u / "taxonomy_normalize.py").exists():
        if str(_u) not in sys.path:
            sys.path.insert(0, str(_u))
        break
from taxonomy_normalize import clean_rank, clean_species

logger = logging.getLogger(__name__)

# Performance constants
QUEUE_SIZE = 20           # Queue size between reader / workers / writer
MIN_SEQUENCE_LENGTH = 50  # Default minimum sequence length
TSV_CHUNK_SIZE = 10_000   # pandas chunk size
DEFAULT_BATCH_SIZE = 1000  # Default batch size if not provided in config

# Accepted ITS markers
ACCEPTED_MARKERS = {"ITS", "ITS1", "ITS2"}


@dataclass
class SequenceRecord:
    """Data container for a normalized BOLD sequence."""
    original_id: str
    sequence: str
    length: int
    marker_code: str
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None


class BOLDParser:
    """Parsing and ingestion pipeline from BOLD into the SQLAlchemy database."""

    def __init__(self, db_manager, config_validation: Dict, database_formatting: Optional[Dict] = None):
        self.db_manager = db_manager
        self.config_validation = config_validation

        # Read batch_size from config
        formatting_config = database_formatting or {}
        self.batch_size = formatting_config.get('batch_size', DEFAULT_BATCH_SIZE)

        # Shared statistics
        self.stats = mp.Manager().dict({
            "total_sequences": 0,
            "processed_sequences": 0,
            "filtered_sequences": 0,
            "marker_filtered": 0,
            "parsing_errors": 0,
            "processing_time": 0.0,
        })

        # Sequence validation settings
        self.valid_chars = set(self.config_validation.get("valid_chars", "ATCGN"))
        self.min_length = self.config_validation.get("min_length", MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get("max_n_fraction", 0.1)

        # Mapping of BOLD columns kept for potential reuse later
        self.column_mapping = {
            "sampleid": "sampleid",
            "record_id": "record_id",
            "kingdom": "kingdom",
            "phylum": "phylum",
            "class": "class_",
            "order": "order_",
            "family": "family",
            "genus": "genus",
            "species": "species",
            "marker_code": "marker_code",
            "nuc": "sequence",
        }

    # --------------------------
    # 1. Sequence validation
    # --------------------------
    def validate_sequence(self, sequence: str) -> Tuple[bool, int, float]:
        """Fast sequence validation (length, N content, allowed characters)."""
        if not sequence:
            return False, 0, 1.0

        sequence = (
            sequence.upper()
            .strip()
            .replace("-", "")
            .replace(" ", "")
            .replace(".", "")
        )

        seq_len = len(sequence)
        if seq_len < self.min_length:
            return False, seq_len, 0.0

        n_count = sequence.count("N")
        n_fraction = n_count / seq_len if seq_len > 0 else 1.0

        if n_fraction > self.max_n_fraction:
            return False, seq_len, n_fraction

        if not all(c in self.valid_chars for c in sequence):
            return False, seq_len, n_fraction

        return True, seq_len, n_fraction

    # --------------------------
    # 2. TSV extraction
    # --------------------------
    def extract_tsv_from_tar(self, tar_path: str) -> str:
        """Extract the TSV from a BOLD .tar.gz archive and return its path."""
        logger.info(f"Extracting TSV from {tar_path}")

        with tarfile.open(tar_path, "r:gz") as tar:
            tsv_members = [
                m for m in tar.getmembers()
                if m.isfile() and m.name.lower().endswith(".tsv")
            ]

            if not tsv_members:
                raise ValueError(f"No TSV file found in {tar_path}")

            if len(tsv_members) > 1:
                logger.warning(
                    f"Multiple TSV files in archive, using first: {tsv_members[0].name}"
                )

            tsv_member = tsv_members[0]
            temp_dir = Path(tar_path).parent / "temp_bold_extract"
            temp_dir.mkdir(exist_ok=True)

            tar.extract(tsv_member, temp_dir)
            extracted_path = temp_dir / tsv_member.name

        logger.info(f"Extracted TSV to: {extracted_path}")
        return str(extracted_path)

    # --------------------------
    # 3. TSV reader → queue
    # --------------------------
    def tsv_reader(self, file_path: str, queue: Queue) -> None:
        """Robust chunked TSV reader using pandas that pushes into the queue."""
        logger.info(f"Starting pandas TSV reader for {file_path}")
        line_count = 0

        try:
            reader_iterator = pd.read_csv(
                file_path,
                sep="\t",
                engine="python",
                quoting=csv.QUOTE_NONE,
                on_bad_lines="warn",
                chunksize=TSV_CHUNK_SIZE,
                dtype=str,
            )

            for chunk_df in reader_iterator:
                line_count += len(chunk_df)

                # Casting en str
                chunk_df = chunk_df.astype(str)

                # Filtrage markers
                if "marker_code" not in chunk_df.columns:
                    logger.error("Column 'marker_code' missing in TSV.")
                    continue

                chunk_df["marker_code"] = (
                    chunk_df["marker_code"].str.strip().str.upper()
                )
                initial_count = len(chunk_df)
                filtered_chunk = chunk_df[chunk_df["marker_code"].isin(ACCEPTED_MARKERS)]

                self.stats["marker_filtered"] += initial_count - len(filtered_chunk)

                if filtered_chunk.empty:
                    continue

                valid_rows_for_queue: List[Tuple[Dict, int]] = []

                for row in filtered_chunk.itertuples(index=False):
                    row_dict = row._asdict()

                    sequence_raw = row_dict.get("nuc", "")
                    if pd.isna(sequence_raw):
                        self.stats["filtered_sequences"] += 1
                        continue

                    sequence_raw = str(sequence_raw).strip().replace("-", "")

                    is_valid, seq_len, _ = self.validate_sequence(sequence_raw)
                    if not is_valid:
                        self.stats["filtered_sequences"] += 1
                        continue

                    valid_rows_for_queue.append((row_dict, seq_len))
                    self.stats["total_sequences"] += 1

                if valid_rows_for_queue:
                    queue.put(valid_rows_for_queue)

                if line_count % (TSV_CHUNK_SIZE * 5) == 0:
                    logger.info(
                        "Processed %d lines, accepted %d sequences so far",
                        line_count,
                        self.stats["total_sequences"],
                    )

        except Exception as e:
            logger.error(f"Error in pandas TSV reader: {e}", exc_info=True)
            raise
        finally:
            queue.put(None)
            logger.info(
                "Pandas TSV reader finished. Total lines processed: %d", line_count
            )

    # --------------------------
    # 4. Field cleaning
    # --------------------------
    def _clean_sequence(self, sequence: str) -> str:
        """Deep-clean a sequence: uppercase + valid character filter."""
        if not sequence:
            return ""

        sequence = sequence.upper().strip()
        sequence = sequence.replace("-", "").replace(" ", "").replace(".", "")

        return "".join(c for c in sequence if c in self.valid_chars)


    # --------------------------
    # 5. Transformation worker
    # --------------------------
    def sequence_processor_worker(
        self, chunk: List[Tuple[Dict, int]]
    ) -> List[SequenceRecord]:
        """
        Transform a chunk of TSV rows (dict) into normalized SequenceRecords.
        """
        results: List[SequenceRecord] = []

        for row_data, seq_len in chunk:
            try:
                kingdom = clean_rank(row_data.get("kingdom"))
                phylum = clean_rank(row_data.get("phylum"))
                class_ = clean_rank(row_data.get("class"))
                order_ = clean_rank(row_data.get("order"))
                family = clean_rank(row_data.get("family"))
                genus = clean_rank(row_data.get("genus"))
                species = clean_species(row_data.get("species"))

                raw_sequence = (
                    row_data.get("nuc", "") or ""
                )
                raw_sequence = raw_sequence.strip().upper().replace("-", "")
                sequence = self._clean_sequence(raw_sequence)

                original_id = (
                    row_data.get("record_id")
                    or row_data.get("processid")
                    or row_data.get("sampleid")
                    or ""
                )
                record = SequenceRecord(
                    original_id=str(original_id).strip(),
                    sequence=sequence,
                    length=seq_len,
                    marker_code=(row_data.get("marker_code", "") or "").strip(),
                    kingdom=kingdom,
                    phylum=phylum,
                    class_=class_,
                    order_=order_,
                    family=family,
                    genus=genus,
                    species=species,
                )
                results.append(record)

            except Exception as e:
                logger.warning(f"Error processing record: {e}")
                self.stats["parsing_errors"] += 1

        return results

    # --------------------------
    # 6. Database writer (bulk SQLAlchemy)
    # --------------------------
    def db_writer_consumer(self, source_id: str, results_queue: Queue) -> None:
        """
        Consumer that pulls SequenceRecords from results_queue and flushes them
        to the database in batches through SQLAlchemy.

        Args:
            source_id: Source ID (UUID string)
        """
        logger.info("Starting BOLD DB writer")
        batch: List[SequenceRecord] = []

        try:
            while True:
                try:
                    result_chunk = results_queue.get(timeout=6000)
                    if result_chunk is None:
                        break

                    batch.extend(result_chunk)

                    if len(batch) >= self.batch_size:
                        self._flush_bold_batch_sqlalchemy(batch, source_id)
                        self.stats["processed_sequences"] += len(batch)
                        logger.info(
                            "Processed %d BOLD sequences",
                            self.stats["processed_sequences"],
                        )
                        batch.clear()

                except Empty:
                    logger.warning("DB writer timeout")
                    break

            if batch:
                self._flush_bold_batch_sqlalchemy(batch, source_id)
                self.stats["processed_sequences"] += len(batch)

            logger.info("BOLD DB writer finished successfully")

        except Exception as e:
            logger.error(f"BOLD DB writer error: {e}", exc_info=True)
            raise

    def _flush_bold_batch_sqlalchemy(
        self, batch: List[SequenceRecord], source_id: str
    ) -> None:
        """
        Flush a batch of SequenceRecords through db_manager (source_id is UUID):

        1. Prepare taxonomy records
        2. Bulk insert taxonomies → tax_id_map
        3. Prepare sequences
        4. Bulk insert sequences → seq_id_map
        5. Prepare seq_tax relations (Raw/Pending) with silent deduplication
        6. Bulk insert seq_tax entries
        """
        # 1. Taxonomies
        taxonomy_data: List[Dict] = []
        taxonomy_map: Dict[str, Tuple] = {}

        for record in batch:
            tax_tuple = (
                record.kingdom,
                record.phylum,
                record.class_,
                record.order_,
                record.family,
                record.genus,
                record.species,
            )
            taxonomy_map[record.original_id] = tax_tuple

            taxonomy_data.append(
                {
                    "kingdom": record.kingdom,
                    "phylum": record.phylum,
                    "class_": record.class_,
                    "order_": record.order_,
                    "family": record.family,
                    "genus": record.genus,
                    "species": record.species,
                }
            )

        tax_id_map = self.db_manager.bulk_insert_taxonomies(taxonomy_data)

        # 2. Sequences
        sequence_data: List[Dict] = []
        for record in batch:
            seq_dict = {
                "source_id": source_id,
                "original_id": record.original_id,
                "sequence": record.sequence,
                "length": record.length,
                "barcode": record.marker_code,
                "barcode_status": "Pending",
            }
            sequence_data.append(seq_dict)

        seq_id_map = self.db_manager.bulk_insert_sequences(sequence_data)

        # 3. Seq_tax relations (with deduplication)
        relation_set = set()  # (seq_id, tax_id, 'Raw')
        relation_data: List[Dict] = []

        for record in batch:
            seq_id = seq_id_map.get(record.original_id)
            tax_tuple = taxonomy_map.get(record.original_id)
            tax_id = tax_id_map.get(tax_tuple)

            if not seq_id or not tax_id:
                continue

            key = (seq_id, tax_id, "Raw")
            if key in relation_set:
                # Option A: silent deduplication
                continue

            relation_set.add(key)
            relation_data.append(
                {
                    "seq_id": seq_id,
                    "tax_id": tax_id,
                    "backbone": "Raw",
                    "status": "Pending",
                    "confidence": None,
                }
            )

        if relation_data:
            self.db_manager.bulk_insert_seq_tax(relation_data)

    # --------------------------
    # 7. Full orchestration
    # --------------------------
    def process_bold_database(self, input_file: str, source_name: str = "BOLD") -> Dict:
        """
        Complete BOLD processing pipeline:
        - configure DB for bulk operations
        - extract the TSV file
        - launch reader + workers + writer
        - clean up and restore DB settings
        - return processing statistics
        """
        start_time = time.time()
        logger.info(f"Starting BOLD processing: {input_file}")

        # DB optimizations (no-op in V2 but kept for compatibility)
        self.db_manager.optimize_for_bulk_insert_advanced()
        source_id = self.db_manager.add_source(name=source_name, file_path=input_file)

        temp_tsv_path = None
        try:
            temp_tsv_path = self.extract_tsv_from_tar(input_file)

            tsv_queue: Queue = Queue(maxsize=QUEUE_SIZE)
            results_queue: Queue = Queue(maxsize=QUEUE_SIZE)

            cpu_count = os.cpu_count() or 4
            num_workers = max(2, cpu_count - 2)
            logger.info(
                "Using %d worker threads on %d CPUs", num_workers, cpu_count
            )

            with ThreadPoolExecutor(max_workers=num_workers + 2) as executor:
                writer_future = executor.submit(
                    self.db_writer_consumer, source_id, results_queue
                )
                reader_future = executor.submit(
                    self.tsv_reader, temp_tsv_path, tsv_queue
                )

                worker_futures = []
                processed_chunks = 0

                while True:
                    try:
                        chunk = tsv_queue.get(timeout=6000)
                        if chunk is None:
                            break

                        fut = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(fut)
                        processed_chunks += 1

                        if len(worker_futures) >= num_workers * 2:
                            completed_future = next(as_completed(worker_futures))
                            worker_futures.remove(completed_future)
                            results_queue.put(completed_future.result())

                    except Empty:
                        logger.warning("Reader timeout - assuming completion")
                        break

                # Remaining workers
                for fut in as_completed(worker_futures):
                    results_queue.put(fut.result())

                results_queue.put(None)
                writer_future.result()
                reader_future.result()

                logger.info("Processed %d chunks total", processed_chunks)

        except Exception as e:
            logger.error(
                f"Critical error in BOLD processing pipeline: {e}", exc_info=True
            )
            raise
        finally:
            # Clean up temporary files
            if temp_tsv_path and os.path.exists(temp_tsv_path):
                try:
                    os.remove(temp_tsv_path)
                    temp_dir = Path(temp_tsv_path).parent
                    if temp_dir.name == "temp_bold_extract":
                        temp_dir.rmdir()
                except Exception as e:
                    logger.warning(f"Could not cleanup temp files: {e}")

            self.db_manager.restore_safe_settings()

        end_time = time.time()
        self.stats["processing_time"] = end_time - start_time

        return dict(self.stats)


# --------------------------
# 8. Script entry point
# --------------------------
def main():
    try:
        # --- Contexte Snakemake ---
        if "snakemake" in globals():
            input_file = snakemake.input["source"]
            output_report = snakemake.output["report"]
            config_validation = snakemake.params.sequence_validation
            database_formatting = snakemake.params.database_formatting
            config = snakemake.config

            log_file = (
                snakemake.log[0]
                if isinstance(snakemake.log, list)
                else snakemake.log
            )
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(
                logging.Formatter(
                    "%(asctime)s - %(levelname)s - %(message)s"
                )
            )
            logger.addHandler(file_handler)
        else:
            # --- Mode standalone ---
            logger.warning("Running in standalone test mode")
            input_file = "test_bold.tar.gz"
            output_report = "test_bold_report.json"
            config_validation = {
                "min_length": 50,
                "valid_chars": "ATCGN",
                "max_n_fraction": 0.1,
            }
            config = {
                "database": {
                    "type": "sqlite",
                    "sqlite": {"path": "data/db", "name": "its.db"},
                }
            }
            database_formatting = {}

        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)
        parser = BOLDParser(db_manager, config_validation, database_formatting)

        stats = parser.process_bold_database(input_file=input_file)

        # VACUUM / ANALYZE
        db_manager.vacuum_database()

        processing_time = stats.get("processing_time", 1.0) or 1.0
        processed_seqs = stats.get("processed_sequences", 0)
        rate = processed_seqs / processing_time if processing_time > 0 else 0.0

        report_data = {
            "source": "BOLD",
            "version": "latest",
            "input_file": input_file,
            "accepted_markers": sorted(list(ACCEPTED_MARKERS)),
            "stats": stats,
            "processing_rate_per_second": rate,
            "performance_metrics": {
                "cpu_cores_used": os.cpu_count(),
                "memory_optimizations": True,
                "batch_size": parser.batch_size,
                "chunk_size": TSV_CHUNK_SIZE,
            },
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }

        if "snakemake" in globals():
            Path(output_report).parent.mkdir(parents=True, exist_ok=True)
            with open(output_report, "w") as f:
                json.dump(report_data, f, indent=2)
            logger.info("BOLD processing report saved to: %s", output_report)
        else:
            with open(output_report, "w") as f:
                json.dump(report_data, f, indent=2)
            logger.info("Standalone BOLD report saved to: %s", output_report)

        logger.info("Final BOLD processing rate: %.0f sequences/second", rate)

    except Exception as e:
        logger.critical(f"Fatal error in BOLD processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
