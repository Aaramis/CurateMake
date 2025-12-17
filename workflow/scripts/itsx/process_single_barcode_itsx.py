#!/usr/bin/env python3
"""
Process a single ITSx barcode (parallelized via Snakemake).
"""
import logging
import sys
import json
import uuid
import time
from pathlib import Path
from typing import Dict, List, Tuple, Iterator, Optional
from datetime import datetime
from sqlalchemy import text
from sqlalchemy.exc import OperationalError

try:
    from psycopg2.errors import SerializationFailure as PGSerializationFailure
except ImportError:  # pragma: no cover - psycopg2 not always installed (e.g., SQLite mode)
    PGSerializationFailure = None

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# FONCTIONS UTILITAIRES
# ============================================================================

def parse_fasta_streaming(file_path: Path) -> Iterator[Tuple[str, str]]:
    """Stream a FASTA file and yield (header, sequence) tuples."""
    if not file_path.exists():
        return

    with file_path.open('r') as f:
        header, sequence = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)


def parse_header(header: str) -> Optional[str]:
    """Extract the seq_id from a FASTA header."""
    try:
        return header.split('|', 1)[0]
    except IndexError:
        return None


# ============================================================================
# CLASSE PROCESSEUR BARCODE
# ============================================================================

class BarcodeITSxProcessor:
    """Process a single ITSx barcode."""

    def __init__(self, db_manager, barcode: str, itsx_dir: Path, data_dir: Path, config: Dict):
        self.db_manager = db_manager
        self.barcode = barcode
        self.itsx_dir = Path(itsx_dir)
        self.data_dir = Path(data_dir)
        self.config = config

        self.insert_batch_size = config.get('insert_batch_size', 10000)
        self.min_length = config.get('min_length', 50)
        self.max_n_fraction = config.get('max_n_fraction', 0.05)
        self.sql_variable_limit = 900 if db_manager.db_type == "sqlite" else 10000

    def _get_parent_metadata(self, parent_ids: List[str]) -> Dict[str, Dict]:
        """Retrieve metadata for the parent sequences."""
        if not parent_ids:
            return {}

        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table
        all_metadata = {}

        with self.db_manager.get_session() as session:
            if self.db_manager.db_type == "postgresql":
                query = text(f"""
                    SELECT s.seq_id, s.source_id, s.original_id, st.tax_id
                    FROM {sequences_table} s
                    JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                    WHERE s.seq_id = ANY(:parent_ids) AND st.backbone = 'CoL'
                """)
                rows = session.execute(query, {"parent_ids": parent_ids}).fetchall()

                for row in rows:
                    all_metadata[row.seq_id] = {
                        "source_id": row.source_id,
                        "original_id": row.original_id,
                        "tax_id": row.tax_id,
                    }
            else:
                for i in range(0, len(parent_ids), self.sql_variable_limit):
                    chunk = parent_ids[i:i + self.sql_variable_limit]
                    placeholders = ",".join([f":id{j}" for j in range(len(chunk))])
                    params = {f"id{j}": chunk[j] for j in range(len(chunk))}

                    query = text(f"""
                        SELECT s.seq_id, s.source_id, s.original_id, st.tax_id
                        FROM {sequences_table} s
                        JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                        WHERE s.seq_id IN ({placeholders}) AND st.backbone = 'CoL'
                    """)
                    rows = session.execute(query, params).fetchall()

                    for row in rows:
                        all_metadata[row.seq_id] = {
                            "source_id": row.source_id,
                            "original_id": row.original_id,
                            "tax_id": row.tax_id,
                        }

        return all_metadata

    def _bulk_insert_sequences(self, sequences_batch: List[Dict]) -> int:
        """Insert sequences in bulk using executemany."""
        if not sequences_batch:
            return 0

        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table

        seq_params = []
        tax_params = []

        for seq_data in sequences_batch:
            seq_id = str(uuid.uuid4())

            seq_params.append({
                'seq_id': seq_id,
                'source_id': seq_data['source_id'],
                'original_id': seq_data['original_id'],
                'sequence': seq_data['sequence'],
                'length': len(seq_data['sequence']),
                'barcode': seq_data['barcode'],
                'barcode_status': 'ITSx',
                'parent_seq_id': seq_data['parent_seq_id']
            })

            tax_params.append({
                'seq_id': seq_id,
                'tax_id': seq_data['tax_id']
            })

        try:
            with self.db_manager.get_session() as session:
                session.execute(
                    text(f"""
                        INSERT INTO {sequences_table}
                        (seq_id, source_id, original_id, sequence, length,
                         barcode, barcode_status, parent_seq_id)
                        VALUES (:seq_id, :source_id, :original_id, :sequence, :length,
                                :barcode, :barcode_status, :parent_seq_id)
                    """),
                    seq_params
                )

                session.execute(
                    text(f"""
                        INSERT INTO {seq_tax_table} (seq_id, tax_id, backbone)
                        VALUES (:seq_id, :tax_id, 'CoL')
                    """),
                    tax_params
                )

                session.commit()
            return len(sequences_batch)

        except Exception as e:
            logger.error(f"Insert error: {e}")
            raise

    def _is_serialization_failure(self, error: OperationalError) -> bool:
        """Return True if the exception comes from a PostgreSQL serialization failure."""
        orig = getattr(error, "orig", None)
        if orig is None:
            return False
        if getattr(orig, "pgcode", None) == "40001":  # PostgreSQL serialization_failure
            return True
        if PGSerializationFailure and isinstance(orig, PGSerializationFailure):
            return True
        return False

    def _bulk_update_status(self, updates: List[Tuple[str, str]]) -> int:
        """Bulk update barcode statuses."""
        if not updates:
            return 0

        sequences_table = self.db_manager.sequences_table

        if self.db_manager.db_type == "postgresql":
            payload = [{"seq_id": seq_id, "status": status} for status, seq_id in updates]
            max_retries = self.config.get("status_update_retries", 3)

            for attempt in range(1, max_retries + 1):
                try:
                    with self.db_manager.get_session() as session:
                        session.execute(text("SET LOCAL TRANSACTION ISOLATION LEVEL READ COMMITTED"))
                        session.execute(text("DROP TABLE IF EXISTS temp_updates"))
                        session.execute(text("""
                            CREATE TEMP TABLE temp_updates (seq_id TEXT PRIMARY KEY, status TEXT)
                        """))
                        session.execute(
                            text("INSERT INTO temp_updates (seq_id, status) VALUES (:seq_id, :status)"),
                            payload
                        )
                        session.execute(text(f"""
                            UPDATE {sequences_table} AS s
                            SET barcode_status = t.status
                            FROM temp_updates AS t
                            WHERE s.seq_id = t.seq_id
                        """))
                    break
                except OperationalError as exc:
                    if self._is_serialization_failure(exc) and attempt < max_retries:
                        backoff = 0.5 * attempt
                        logger.warning(
                            "Serialization failure while updating barcode statuses for %s "
                            "(attempt %d/%d). Retrying in %.1fs...",
                            self.barcode,
                            attempt,
                            max_retries,
                            backoff
                        )
                        time.sleep(backoff)
                        continue
                    raise
        else:
            with self.db_manager.get_session() as session:
                session.execute(text("""
                    CREATE TEMPORARY TABLE IF NOT EXISTS temp_updates (
                        seq_id TEXT PRIMARY KEY, status TEXT
                    )
                """))
                session.execute(text("DELETE FROM temp_updates"))

                for i in range(0, len(updates), 10000):
                    batch = updates[i:i + 10000]
                    session.execute(
                        text("INSERT INTO temp_updates (seq_id, status) VALUES (:seq_id, :status)"),
                        [{"seq_id": seq_id, "status": status} for status, seq_id in batch]
                    )

                session.execute(text(f"""
                    UPDATE {sequences_table}
                    SET barcode_status = (
                        SELECT status FROM temp_updates
                        WHERE temp_updates.seq_id = {sequences_table}.seq_id
                    )
                    WHERE seq_id IN (SELECT seq_id FROM temp_updates)
                """))

        return len(updates)

    def _validate_sequence(self, sequence: str) -> bool:
        """Validate a sequence against length and ambiguity thresholds."""
        if len(sequence) < self.min_length:
            return False
        n_fraction = sequence.upper().count('N') / len(sequence) if sequence else 1.0
        return n_fraction <= self.max_n_fraction

    def _read_region_file_once(self, barcode_dir: Path, region: str,
                               parent_metadata: Dict) -> Tuple[List[Dict], Dict[str, str]]:
        """
        Read a region-specific ITSx FASTA file exactly once.

        Returns:
            (insert_data, region_sequences_dict)
        """
        region_fasta = barcode_dir / f"{self.barcode}.{region}.fasta"

        insert_data = []
        region_sequences = {}

        if not region_fasta.exists():
            return insert_data, region_sequences

        for header, sequence in parse_fasta_streaming(region_fasta):
            parent_id = parse_header(header)
            if not parent_id or parent_id not in parent_metadata:
                continue

            if not self._validate_sequence(sequence):
                continue

            region_sequences[parent_id] = sequence

            meta = parent_metadata[parent_id]
            insert_data.append({
                "sequence": sequence,
                "barcode": region,
                "parent_seq_id": parent_id,
                "source_id": meta["source_id"],
                "original_id": meta["original_id"],
                "tax_id": meta["tax_id"]
            })

        return insert_data, region_sequences

    def process(self) -> Dict:
        """Run the complete ITSx workflow for this barcode."""
        logger.info(f"=== Processing barcode: {self.barcode} ===")

        stats = {
            "barcode": self.barcode,
            "status": "Success",
            "start_time": datetime.now().isoformat(),
            "initial_sequences_processed": 0,
            "subregions": {},
            "reconstructed_its": {"common_parents": 0, "inserted": 0},
            "no_detection": 0,
            "total_new_sequences": 0
        }

        try:
            barcode_dir = self.itsx_dir / self.barcode

            # ================================================================
            # STEP 1: mark parent sequences as "Processed"
            # ================================================================
            original_fasta = self.data_dir / f"{self.barcode}.fasta"

            if original_fasta.exists():
                parent_ids = [
                    parse_header(h) for h, _ in parse_fasta_streaming(original_fasta)
                    if parse_header(h)
                ]

                if parent_ids:
                    updates = [("Processed", pid) for pid in parent_ids]
                    self._bulk_update_status(updates)
                    stats["initial_sequences_processed"] = len(parent_ids)
                    logger.info(f"Marked {len(parent_ids)} sequences as 'Processed'")

            # ================================================================
            # STEP 2: collect parent_ids
            # ================================================================
            all_parent_ids = set()
            for region in ["ITS1", "ITS2", "5_8S"]:
                region_fasta = barcode_dir / f"{self.barcode}.{region}.fasta"
                if region_fasta.exists():
                    for header, _ in parse_fasta_streaming(region_fasta):
                        if pid := parse_header(header):
                            all_parent_ids.add(pid)

            parent_metadata = self._get_parent_metadata(list(all_parent_ids))

            # ================================================================
            # STEP 3: single-pass parsing for each region file
            # ================================================================
            logger.info("Reading each ITSx region file once...")

            its1_insert, its1_seqs = self._read_region_file_once(barcode_dir, "ITS1", parent_metadata)
            its2_insert, its2_seqs = self._read_region_file_once(barcode_dir, "ITS2", parent_metadata)
            _, its58s_seqs = self._read_region_file_once(barcode_dir, "5_8S", parent_metadata)

            stats["subregions"]["ITS1"] = {"found": len(its1_insert), "inserted": 0}
            stats["subregions"]["ITS2"] = {"found": len(its2_insert), "inserted": 0}

            # ================================================================
            # STEP 4: insert ITS1 and ITS2 sequences
            # ================================================================
            total_inserted = 0

            if its1_insert:
                for i in range(0, len(its1_insert), self.insert_batch_size):
                    batch = its1_insert[i:i + self.insert_batch_size]
                    inserted = self._bulk_insert_sequences(batch)
                    stats["subregions"]["ITS1"]["inserted"] += inserted
                    total_inserted += inserted
                logger.info(f"ITS1: inserted {stats['subregions']['ITS1']['inserted']} sequences")

            if its2_insert:
                for i in range(0, len(its2_insert), self.insert_batch_size):
                    batch = its2_insert[i:i + self.insert_batch_size]
                    inserted = self._bulk_insert_sequences(batch)
                    stats["subregions"]["ITS2"]["inserted"] += inserted
                    total_inserted += inserted
                logger.info(f"ITS2: inserted {stats['subregions']['ITS2']['inserted']} sequences")

            # ================================================================
            # STEP 5: reconstruct full ITS sequences
            # ================================================================
            common_ids = set(its1_seqs.keys()) & set(its58s_seqs.keys()) & set(its2_seqs.keys())
            stats["reconstructed_its"]["common_parents"] = len(common_ids)

            if common_ids:
                logger.info(f"Reconstructing {len(common_ids)} ITS sequences")

                reconstruction_batch = []
                for parent_id in common_ids:
                    full_seq = its1_seqs[parent_id] + its58s_seqs[parent_id] + its2_seqs[parent_id]

                    if not self._validate_sequence(full_seq):
                        continue

                    meta = parent_metadata[parent_id]
                    reconstruction_batch.append({
                        "sequence": full_seq,
                        "barcode": "ITS",
                        "parent_seq_id": parent_id,
                        "source_id": meta["source_id"],
                        "original_id": meta["original_id"],
                        "tax_id": meta["tax_id"]
                    })

                if reconstruction_batch:
                    for i in range(0, len(reconstruction_batch), self.insert_batch_size):
                        batch = reconstruction_batch[i:i + self.insert_batch_size]
                        inserted = self._bulk_insert_sequences(batch)
                        stats["reconstructed_its"]["inserted"] += inserted
                        total_inserted += inserted

                    logger.info(f"ITS: reconstructed {stats['reconstructed_its']['inserted']} sequences")

            # ================================================================
            # STEP 6: mark "no detection" sequences
            # ================================================================
            no_detection_fasta = barcode_dir / f"{self.barcode}_no_detections.fasta"

            if no_detection_fasta.exists():
                no_detection_ids = [
                    parse_header(h) for h, _ in parse_fasta_streaming(no_detection_fasta)
                    if parse_header(h)
                ]

                if no_detection_ids:
                    updates = [("no_detection", nid) for nid in no_detection_ids]
                    self._bulk_update_status(updates)
                    stats["no_detection"] = len(no_detection_ids)
                    logger.info(f"Marked {len(no_detection_ids)} sequences as 'no_detection'")

            stats["total_new_sequences"] = total_inserted
            stats["end_time"] = datetime.now().isoformat()
            logger.info(f"✓ {self.barcode}: inserted {total_inserted} sequences")

        except Exception as e:
            stats["status"] = "Error"
            stats["error_message"] = str(e)
            stats["end_time"] = datetime.now().isoformat()
            logger.error(f"Error while processing {self.barcode}: {e}", exc_info=True)

        return stats


# ============================================================================
# MAIN (SNAKEMAKE)
# ============================================================================

def main():
    """Snakemake entry point."""
    if 'snakemake' not in globals():
        logger.error("This script must be executed via Snakemake")
        sys.exit(1)

    sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
    from db_manager import get_database_manager

    # Snakemake parameters
    barcode = snakemake.wildcards.barcode
    itsx_dir = snakemake.params.itsx_dir
    data_dir = snakemake.params.data_dir
    output_report = snakemake.output.report
    log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

    # Logging
    if log_file:
        handler = logging.FileHandler(log_file, mode='w')
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)

    config = {
        "min_length": snakemake.config.get('sequence_validation', {}).get('min_length', 50),
        "max_n_fraction": snakemake.config.get('sequence_validation', {}).get('max_n_fraction', 0.05),
        "insert_batch_size": 10000,  # Taille de batch pour insertions DB
    }

    try:
        db_manager = get_database_manager(snakemake.config)
        processor = BarcodeITSxProcessor(db_manager, barcode, itsx_dir, data_dir, config)

        report = processor.process()

        # Write report
        Path(output_report).parent.mkdir(parents=True, exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(report, f, indent=2)

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


if __name__ == "__main__":
    main()
