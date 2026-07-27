#!/usr/bin/env python3
"""
Extract sequences for downstream MSA using SQLAlchemy (SQLite + PostgreSQL).
Produces a CSV with sequences plus full taxonomy grouped by barcode.

Relies on DatabaseManager (db_manager.py).
"""

import json
import sys
import os
import time
import csv
from pathlib import Path
from collections import defaultdict
import logging
from typing import Dict, List, Tuple, Optional

from sqlalchemy import text

# Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MSASequenceExtractor:
    """SQLAlchemy-powered extractor for the MSA step."""

    def __init__(self, db_manager, db_label: str, batch_size: int = 5000,
                 create_indexes: bool = False):
        """
        Args:
            db_manager: DatabaseManager instance
            db_label: identifier shown in logs/reports only
            batch_size: number of sequences fetched per query
            create_indexes: whether to attempt index creation (can fail on YugabyteDB)
        """
        self.db_manager = db_manager
        self.db_label = str(db_label)
        self.db_type = db_manager.db_type
        self.batch_size = max(int(batch_size), 1)
        self.sql_variable_limit = 900 if self.db_type == "sqlite" else 10000
        self.create_indexes = create_indexes

    def _get_table(self, table_name: str) -> str:
        """Get schema-qualified table name for PostgreSQL, unqualified for SQLite."""
        if self.db_manager.db_type == "postgresql":
            schema = self.db_manager.connection_params.get('schema', 'public')
            return f"{schema}.{table_name}"
        return table_name

    # =====================================================================
    # INDEX CREATION (PG + SQLite)
    # =====================================================================
    def ensure_indexes(self):
        """Create indexes required to speed up downstream joins.

        Each index is created in its own transaction with retry logic.
        Failures are logged but do not stop the pipeline.
        """
        if not self.create_indexes:
            logger.info("Index creation disabled by config, skipping...")
            return

        logger.info("Checking/creating MSA indexes...")

        seq_tax_table = self._get_table('seq_tax')
        sequences_table = self._get_table('sequences')
        taxonomy_table = self._get_table('taxonomy')

        index_statements = [
            (f"idx_seq_tax_backbone_status", f"""
                CREATE INDEX IF NOT EXISTS idx_seq_tax_backbone_status
                ON {seq_tax_table}(backbone, status)
            """),
            (f"idx_sequences_barcode_status", f"""
                CREATE INDEX IF NOT EXISTS idx_sequences_barcode_status
                ON {sequences_table}(barcode_status)
            """),
            (f"idx_sequences_barcode", f"""
                CREATE INDEX IF NOT EXISTS idx_sequences_barcode
                ON {sequences_table}(barcode)
            """),
            (f"idx_sequences_barcode_seq_id", f"""
                CREATE INDEX IF NOT EXISTS idx_sequences_barcode_seq_id
                ON {sequences_table}(barcode, seq_id)
            """),
            (f"idx_seq_tax_backbone_seq_id", f"""
                CREATE INDEX IF NOT EXISTS idx_seq_tax_backbone_seq_id
                ON {seq_tax_table}(backbone, seq_id)
            """),
            (f"idx_taxonomy_tax_id", f"""
                CREATE INDEX IF NOT EXISTS idx_taxonomy_tax_id
                ON {taxonomy_table}(tax_id)
            """),
        ]

        # Create each index in its own transaction with retry
        for idx_name, sql in index_statements:
            def create_single_index(sql_stmt=sql):
                with self.db_manager.get_session() as session:
                    session.execute(text(sql_stmt))

            try:
                self.db_manager._execute_with_retry(
                    f"msa_create_index_{idx_name}",
                    create_single_index
                )
                logger.debug(f"Index {idx_name} created or already exists")
            except Exception as e:
                # Index creation failure should not stop the pipeline
                logger.warning(f"Could not create index {idx_name}: {e}")
                # Continue with next index

    # =====================================================================
    # ESTIMATE DATASET SIZE
    # =====================================================================
    def estimate_sequence_count(self) -> int:
        """Estimate non-duplicate sequences for logging purposes.

        Uses a simple COUNT without JOIN to avoid timeout on large tables.
        This is just an approximation for progress logging.
        """
        sequences_table = self._get_table('sequences')

        # Simple count without expensive JOIN - just an estimate
        query = text(f"""
            SELECT COUNT(*)
            FROM {sequences_table}
            WHERE barcode_status != 'Duplicate'
        """)

        def fetch_count():
            with self.db_manager.get_session() as session:
                return session.execute(query).scalar()

        try:
            count = self.db_manager._execute_with_retry(
                "msa_estimate_sequence_count",
                fetch_count
            )
            logger.info(f"Estimated sequences to process: ~{count:,}")
            return count
        except Exception as e:
            logger.warning(f"Unable to compute estimate (non-blocking): {e}")
            return 0

    # =====================================================================
    # MSA EXTRACTION
    # =====================================================================
    def _fetch_barcodes(self) -> List[str]:
        """Fetch distinct barcodes to partition extraction.

        Uses simple query without JOIN to avoid timeout on large tables.
        """
        sequences_table = self._get_table('sequences')

        # Simple distinct without expensive JOIN
        query = text(f"""
            SELECT DISTINCT COALESCE(barcode, 'ITS') AS barcode
            FROM {sequences_table}
            WHERE barcode_status != 'Duplicate'
            ORDER BY barcode
        """)

        def fetch():
            with self.db_manager.get_session() as session:
                return [row[0] for row in session.execute(query).fetchall()]

        barcodes = self.db_manager._execute_with_retry("msa_fetch_barcodes", fetch)
        return barcodes

    def extract_msa_sequences(self) -> List[Dict]:
        """Extract sequences in batches using SQLAlchemy.

        Uses server-side cursor streaming to avoid loading all sequences
        into memory and removes ORDER BY to avoid temp_file_limit errors.
        """
        self.ensure_indexes()
        total_estimate = self.estimate_sequence_count()

        seq_tax_table = self._get_table('seq_tax')
        taxonomy_table = self._get_table('taxonomy')

        sequences_with_keys: List[Tuple[Tuple[str, str, str, str], Dict]] = []

        logger.info("Starting MSA sequence extraction (batch_size=%d)...", self.batch_size)

        start = time.time()
        last_log = start

        barcodes = self._fetch_barcodes()
        if not barcodes:
            return []

        logger.info("Found %d barcodes for MSA extraction.", len(barcodes))

        for barcode_idx, barcode in enumerate(barcodes, start=1):
            logger.info("Barcode %d/%d: %s", barcode_idx, len(barcodes), barcode)

            # Stream sequences for this barcode (no ORDER BY to avoid temp_file_limit)
            def process_barcode():
                batches = []
                with self.db_manager.get_session() as session:
                    if self.db_type == "postgresql":
                        barcode_clause = "s.barcode = :barcode"
                        if barcode == "ITS":
                            barcode_clause = "(s.barcode = :barcode OR s.barcode IS NULL)"
                        seq_query = text(f"""
                            SELECT
                                s.seq_id,
                                s.sequence,
                                COALESCE(s.barcode, 'ITS') AS barcode,
                                s.source_id,
                                s.original_id
                            FROM {self._get_table('sequences')} s
                            WHERE s.barcode_status != 'Duplicate'
                              AND {barcode_clause}
                        """)
                        result = session.execute(
                            seq_query,
                            {"barcode": barcode}
                        ).yield_per(self.batch_size)
                    else:
                        barcode_clause = "s.barcode = :barcode"
                        if barcode == "ITS":
                            barcode_clause = "(s.barcode = :barcode OR s.barcode IS NULL)"
                        seq_query = text(f"""
                            SELECT
                                s.seq_id,
                                s.sequence,
                                COALESCE(s.barcode, 'ITS') AS barcode,
                                s.source_id,
                                s.original_id
                            FROM {self._get_table('sequences')} s
                            WHERE s.barcode_status != 'Duplicate'
                              AND {barcode_clause}
                        """)
                        result = session.execute(seq_query, {"barcode": barcode})

                    batch = []
                    for row in result:
                        batch.append(row)
                        if len(batch) >= self.batch_size:
                            batches.append(batch)
                            batch = []
                    if batch:
                        batches.append(batch)
                return batches

            seq_batches = self.db_manager._execute_with_retry(
                f"msa_extract_barcode_{barcode}",
                process_barcode
            )

            for seq_rows in seq_batches:
                seq_map = {row.seq_id: row for row in seq_rows}
                seq_ids = list(seq_map.keys())

                tax_rows = []
                if self.db_type == "postgresql":
                    tax_query = text(f"""
                        SELECT
                            st.seq_id,
                            COALESCE(t.kingdom, '') AS kingdom,
                            COALESCE(t.phylum, '') AS phylum,
                            COALESCE(t.class_, '') AS class_,
                            COALESCE(t.order_, '') AS order_,
                            COALESCE(t.family, '') AS family,
                            COALESCE(t.genus, '') AS genus,
                            COALESCE(t.species, '') AS species
                        FROM {seq_tax_table} st
                        LEFT JOIN {taxonomy_table} t ON st.tax_id = t.tax_id
                        WHERE st.backbone = 'CoL'
                          AND st.seq_id = ANY(:seq_ids)
                    """)

                    # Capture seq_ids in closure properly
                    def fetch_taxonomy(query=tax_query, ids=seq_ids):
                        with self.db_manager.get_session() as session:
                            return session.execute(query, {"seq_ids": ids}).fetchall()

                    tax_rows = self.db_manager._execute_with_retry(
                        "msa_extract_taxonomy",
                        fetch_taxonomy
                    )
                else:
                    # SQLite: chunk queries with retry
                    for i in range(0, len(seq_ids), self.sql_variable_limit):
                        chunk = seq_ids[i:i + self.sql_variable_limit]
                        placeholders = ",".join([f":id{j}" for j in range(len(chunk))])
                        params = {f"id{j}": chunk[j] for j in range(len(chunk))}

                        chunk_tax_query = text(f"""
                            SELECT
                                st.seq_id,
                                COALESCE(t.kingdom, '') AS kingdom,
                                COALESCE(t.phylum, '') AS phylum,
                                COALESCE(t.class_, '') AS class_,
                                COALESCE(t.order_, '') AS order_,
                                COALESCE(t.family, '') AS family,
                                COALESCE(t.genus, '') AS genus,
                                COALESCE(t.species, '') AS species
                            FROM {seq_tax_table} st
                            LEFT JOIN {taxonomy_table} t ON st.tax_id = t.tax_id
                            WHERE st.backbone = 'CoL'
                              AND st.seq_id IN ({placeholders})
                        """)

                        # Capture variables properly for retry
                        def fetch_tax_chunk(query=chunk_tax_query, p=params):
                            with self.db_manager.get_session() as session:
                                return session.execute(query, p).fetchall()

                        chunk_results = self.db_manager._execute_with_retry(
                            f"msa_extract_taxonomy_chunk_{i}",
                            fetch_tax_chunk
                        )
                        tax_rows.extend(chunk_results)

                for row in tax_rows:
                    seq_row = seq_map.get(row.seq_id)
                    if not seq_row:
                        continue

                    item = {
                        'seq_id': seq_row.seq_id,
                        'sequence': seq_row.sequence,
                        'barcode': seq_row.barcode,
                        'kingdom': row.kingdom,
                        'phylum': row.phylum,
                        'class': row.class_,
                        'order': row.order_,
                        'family': row.family,
                        'genus': row.genus,
                        'species': row.species,
                    }

                    # Rule 1: Streptophyta -> Tracheophyta
                    if item['phylum'].lower() == "streptophyta":
                        item['phylum'] = "Tracheophyta"

                    # Rule 2: replace "incertae sedis" with empty string
                    for k, v in item.items():
                        if isinstance(v, str) and "incertae sedis" in v.lower():
                            item[k] = ""

                    sort_key = (
                        str(seq_row.barcode),
                        str(seq_row.source_id or ""),
                        str(seq_row.original_id or ""),
                        str(seq_row.seq_id)
                    )
                    sequences_with_keys.append((sort_key, item))

                # Periodic logging
                now = time.time()
                if now - last_log >= 30:
                    elapsed = now - start
                    processed = len(sequences_with_keys)
                    rate = processed / elapsed if elapsed > 0 else 0

                    if total_estimate > 0:
                        pct = processed / total_estimate * 100
                        eta = (total_estimate - processed) / rate if rate > 0 else 0
                        logger.info(
                            f"Progress: {processed:,}/{total_estimate:,} "
                            f"({pct:.1f}%) — {rate:.0f} seq/s — ETA {eta/60:.1f} min"
                        )
                    else:
                        logger.info(f"{processed:,} sequences — {rate:.0f} seq/s")

                    last_log = now

        elapsed = time.time() - start
        sequences_with_keys.sort(key=lambda item: item[0])
        sequences = [item for _, item in sequences_with_keys]

        logger.info(
            f"MSA extraction finished: {len(sequences):,} sequences "
            f"in {elapsed:.1f}s ({len(sequences)/elapsed:.0f} seq/s)"
        )

        return sequences

    # =====================================================================
    # CSV WRITING
    # =====================================================================
    @staticmethod
    def write_csv(sequences: List[Dict], out_csv: str) -> str:
        out = Path(out_csv)
        out.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Writing {len(sequences):,} sequences to CSV → {out_csv}")

        fields = [
            'seq_id', 'sequence', 'kingdom', 'phylum', 'class', 'order',
            'family', 'genus', 'species', 'barcode'
        ]

        with out.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()

            batch = 1000
            for i in range(0, len(sequences), batch):
                writer.writerows(sequences[i:i+batch])
                if len(sequences) > 20000 and i % 10000 == 0:
                    logger.info(f"  CSV progress: {i:,}/{len(sequences):,}")

        return str(out)

    # =====================================================================
# STATISTICS
    # =====================================================================
    @staticmethod
    def generate_stats(sequences: List[Dict]) -> Dict:
        if not sequences:
            return {
                "total_sequences": 0,
                "unique_taxonomies": 0,
                "barcode_distribution": {},
                "taxonomy_coverage": {},
                "sequence_length_stats": {}
            }

        barcode_counts = defaultdict(int)
        lengths = []
        taxonomy_sets = defaultdict(set)

        for seq in sequences:
            barcode_counts[seq["barcode"]] += 1
            lengths.append(len(seq["sequence"]))

            for lvl in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
                if seq[lvl]:
                    taxonomy_sets[lvl].add(seq[lvl])

        length_stats = {
            "min": min(lengths),
            "max": max(lengths),
            "mean": sum(lengths) / len(lengths),
            "median": sorted(lengths)[len(lengths)//2]
        }

        taxonomy_coverage = {lvl: len(vals) for lvl, vals in taxonomy_sets.items()}

        return {
            "total_sequences": len(sequences),
            "unique_taxonomies": len(set(tuple(seq[k] for k in taxonomy_sets) for seq in sequences)),
            "barcode_distribution": dict(barcode_counts),
            "taxonomy_coverage": taxonomy_coverage,
            "sequence_length_stats": length_stats
        }


# =====================================================================
# MAIN
# =====================================================================
def main():
    start = time.time()

    try:
        # ------------------------
        # Snakemake / standalone
        # ------------------------
        if "snakemake" in globals():
            db_label = snakemake.params.db
            out_report = snakemake.output.report
            out_csv = getattr(snakemake.output, "sequences", None)
            if not out_csv:
                out_csv = snakemake.params.get("csv_path") or config["msa"]["csv_path"]
            config = snakemake.config

            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            fh = logging.FileHandler(log_file, mode="w")
            fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            logger.addHandler(fh)

        else:
            logger.warning("Mode standalone")
            db_label = "data/db/its.db"
            out_csv = "tmp/msa_sequences.csv"
            out_report = "tmp/msa_prepare_report.json"
            config = {
                "database": {
                    "type": "sqlite",
                    "sqlite": {"path": "data/db", "name": "its.db"},
                }
            }
            Path(out_csv).parent.mkdir(parents=True, exist_ok=True)

        # ------------------------
        # Initialisation DB
        # ------------------------
        SCRIPT_DIR = Path(__file__).resolve().parent
        PROJECT_SCRIPTS_DIR = SCRIPT_DIR
        while PROJECT_SCRIPTS_DIR.name != "scripts" and PROJECT_SCRIPTS_DIR.parent != PROJECT_SCRIPTS_DIR:
            PROJECT_SCRIPTS_DIR = PROJECT_SCRIPTS_DIR.parent
        UTILS_DIR = PROJECT_SCRIPTS_DIR / "utils"
        sys.path.insert(0, str(UTILS_DIR))
        from db_manager import get_database_manager
        db_manager = get_database_manager(config)

        # Get batch_size from database_formatting (same as db_manager uses)
        formatting_config = config.get("database_formatting", {})
        msa_config = config.get("msa", {})
        batch_size = msa_config.get("batch_size", formatting_config.get("batch_size", 5000))

        # Index creation is optional - can fail on YugabyteDB due to timeouts
        create_indexes = msa_config.get("create_indexes", False)

        extractor = MSASequenceExtractor(
            db_manager,
            db_label,
            batch_size=batch_size,
            create_indexes=create_indexes
        )

        # ------------------------
        # Extraction
        # ------------------------
        sequences = extractor.extract_msa_sequences()

        if not sequences:
            logger.warning("No sequences found for MSA")
            stats = {
                "total_sequences": 0,
                "unique_taxonomies": 0,
                "barcode_distribution": {},
                "taxonomy_coverage": {},
                "sequence_length_stats": {}
            }
            csv_file = None
        else:
            csv_file = extractor.write_csv(sequences, out_csv)
            stats = extractor.generate_stats(sequences)

        elapsed = time.time() - start

        # ------------------------
        # Rapport JSON
        # ------------------------
        report = {
            "source": "MSA_preparation",
            "database": db_label,
            "output_csv": out_csv,
            "stats": stats,
            "csv_file_created": csv_file,
            "processing_time_sec": round(elapsed, 2),
            "sequence_rate": len(sequences) / elapsed if elapsed > 0 else 0.0,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "status": "success"
        }

        Path(out_report).parent.mkdir(parents=True, exist_ok=True)
        with open(out_report, "w") as f:
            json.dump(report, f, indent=2)

        logger.info(f"Report written to: {out_report}")

    except Exception as e:
        logger.critical(f"Fatal MSA error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
