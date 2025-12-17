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

    def __init__(self, db_manager, db_label: str):
        """
        Args:
            db_manager: DatabaseManager instance
            db_label: identifier shown in logs/reports only
        """
        self.db_manager = db_manager
        self.db_label = str(db_label)
        self.db_type = db_manager.db_type

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
        """Create indexes required to speed up downstream joins."""
        logger.info("Checking/creating MSA indexes...")

        seq_tax_table = self._get_table('seq_tax')
        sequences_table = self._get_table('sequences')
        taxonomy_table = self._get_table('taxonomy')

        index_statements = []

        # Index used by MSA queries (compatible with SQLite + PostgreSQL)
        index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS idx_seq_tax_backbone_status
            ON {seq_tax_table}(backbone, status)
        """)

        index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS idx_sequences_barcode_status
            ON {sequences_table}(barcode_status)
        """)

        index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS idx_sequences_barcode
            ON {sequences_table}(barcode)
        """)

        index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS idx_seq_tax_seq_id
            ON {seq_tax_table}(seq_id)
        """)

        index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS idx_taxonomy_tax_id
            ON {taxonomy_table}(tax_id)
        """)

        with self.db_manager.get_session() as session:
            for sql in index_statements:
                try:
                    session.execute(text(sql))
                except Exception as e:
                    # Index creation failure should not stop the pipeline (maybe already exists / permission issue)
                    logger.warning(f"Could not create index: {e}")
            session.commit()

    # =====================================================================
    # ESTIMATE DATASET SIZE
    # =====================================================================
    def estimate_sequence_count(self) -> int:
        """Count non-duplicate CoL sequences for estimation purposes."""
        sequences_table = self._get_table('sequences')
        seq_tax_table = self._get_table('seq_tax')

        query = text(f"""
            SELECT COUNT(DISTINCT s.seq_id)
            FROM {sequences_table} s
            INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
            WHERE st.backbone = 'CoL'
              AND s.barcode_status != 'Duplicate'
        """)

        try:
            with self.db_manager.get_session() as session:
                count = session.execute(query).scalar()
                logger.info(f"Estimated sequences to process: ~{count:,}")
                return count
        except Exception as e:
            logger.warning(f"Unable to compute estimate: {e}")
            return 0

    # =====================================================================
    # MSA EXTRACTION
    # =====================================================================
    def extract_msa_sequences(self) -> List[Dict]:
        """Extract sequences in batches using SQLAlchemy."""
        self.ensure_indexes()
        total_estimate = self.estimate_sequence_count()

        sequences_table = self._get_table('sequences')
        seq_tax_table = self._get_table('seq_tax')
        taxonomy_table = self._get_table('taxonomy')

        # Main query template
        base_query = text(f"""
            SELECT
                s.seq_id,
                s.sequence,
                COALESCE(s.barcode, 'ITS') AS barcode,
                COALESCE(t.kingdom, '') AS kingdom,
                COALESCE(t.phylum, '') AS phylum,
                COALESCE(t.class_, '') AS class_,
                COALESCE(t.order_, '') AS order_,
                COALESCE(t.family, '') AS family,
                COALESCE(t.genus, '') AS genus,
                COALESCE(t.species, '') AS species
            FROM {sequences_table} s
            INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
            LEFT JOIN {taxonomy_table} t ON st.tax_id = t.tax_id
            WHERE st.backbone = 'CoL'
              AND s.barcode_status != 'Duplicate'
            ORDER BY s.barcode, s.seq_id
            LIMIT :limit OFFSET :offset
        """)

        sequences = []
        batch_size = 10000
        offset = 0

        logger.info("Starting MSA sequence extraction...")

        start = time.time()
        last_log = start

        with self.db_manager.get_session() as session:

            while True:
                rows = session.execute(
                    base_query,
                    {"limit": batch_size, "offset": offset}
                ).fetchall()

                if not rows:
                    break

                for row in rows:
                    item = {
                        'seq_id': row.seq_id,
                        'sequence': row.sequence,
                        'barcode': row.barcode,
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

                    sequences.append(item)

                offset += batch_size

                # Periodic logging
                now = time.time()
                if now - last_log >= 30:
                    elapsed = now - start
                    processed = len(sequences)
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

        extractor = MSASequenceExtractor(db_manager, db_label)

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
