#!/usr/bin/env python3
"""
Sequence deduplication script.
Detect and mark duplicate sequences based on sequence + barcode + tax_id.
"""

import os
import sys
import time
import logging
import json
from datetime import datetime
from typing import List, Tuple, Set, Dict, Any
from pathlib import Path
from sqlalchemy import func, text

# Import ORM models
sys.path.append(str(Path(__file__).parent / "utils"))
from db_models import Sequence

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SequenceDeduplicator:
    """Sequence deduplication manager."""

    def __init__(self, db_manager):
        self.db_manager = db_manager
        self.stats = {
            'initial_sequence_count': 0,
            'final_sequence_count': 0,
            'duplicate_groups': 0,
            'duplicates_marked': 0,
            'unique_sequences_kept': 0,
            'processing_time': 0.0
        }

    def setup_database(self):
        """Create indexes required for performant deduplication."""
        logger.info("Verifying/creating indexes for deduplication...")

        with self.db_manager.get_session() as session:
            # Essential indexes for deduplication
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sequences_seq_barcode ON {self.db_manager.sequences_table}(sequence, barcode)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sequences_barcode_status ON {self.db_manager.sequences_table}(barcode_status)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_seq_tax_seq_id ON {self.db_manager.seq_tax_table}(seq_id)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_seq_tax_tax_id ON {self.db_manager.seq_tax_table}(tax_id)"))

            # Composite index to speed grouping
            session.execute(text(f"""
                CREATE INDEX IF NOT EXISTS idx_sequences_dedup_composite
                ON {self.db_manager.sequences_table}(sequence, barcode, barcode_status)
            """))

            session.commit()
        logger.info("Indexes ready.")

    def find_duplicate_sequences(self) -> List[Tuple[str, List[str]]]:
        """
        Find duplicate sequences via a single SQL query.
        A duplicate means identical (sequence, barcode, tax_id).
        Returns a list of tuples: (seq_id_to_keep, [seq_ids_to_mark_duplicate]).
        """
        logger.info("Searching for duplicate sequences via SQL...")

        duplicate_groups = []

        sequences_table = self.db_manager.sequences_table
        seq_tax_table  = self.db_manager.seq_tax_table

        # Dialect-specific GROUP_CONCAT / STRING_AGG
        if self.db_manager.db_type == "postgresql":
            group_concat_expr = "STRING_AGG(CAST(seq_id AS TEXT), ',' ORDER BY seq_id)"
        else:
            group_concat_expr = "GROUP_CONCAT(seq_id ORDER BY seq_id)"

        query = f"""
            WITH sequence_groups AS (
                SELECT
                    s.seq_id,
                    s.sequence,
                    s.barcode,
                    st.tax_id,
                    s.barcode_status,
                    ROW_NUMBER() OVER (
                        PARTITION BY s.sequence, s.barcode, st.tax_id
                        ORDER BY s.seq_id
                    ) as row_num
                FROM {sequences_table} s
                INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                WHERE s.barcode_status != 'Duplicate'
            ),
            duplicate_candidates AS (
                SELECT
                    sequence,
                    barcode,
                    tax_id,
                    COUNT(*) as duplicate_count,
                    {group_concat_expr} as all_seq_ids
                FROM sequence_groups
                GROUP BY sequence, barcode, tax_id
                HAVING COUNT(*) > 1
            )
            SELECT
                sequence,
                barcode,
                tax_id,
                duplicate_count,
                all_seq_ids
            FROM duplicate_candidates
            ORDER BY duplicate_count DESC
        """

        with self.db_manager.get_session() as session:
            result = session.execute(text(query))

            for row in result:
                sequence, barcode, tax_id, duplicate_count, all_seq_ids_str = row

                all_seq_ids = [sid.strip() for sid in all_seq_ids_str.split(',')]

                keeper_seq_id = all_seq_ids[0]
                seq_ids_to_mark = all_seq_ids[1:]

                if seq_ids_to_mark:
                    duplicate_groups.append((keeper_seq_id, seq_ids_to_mark))

        self.stats['duplicate_groups'] = len(duplicate_groups)
        self.stats['duplicates_marked'] = sum(len(seq_ids) for _, seq_ids in duplicate_groups)
        self.stats['unique_sequences_kept'] = len(duplicate_groups)

        logger.info(
            f"Found {self.stats['duplicate_groups']} duplicate groups "
            f"({self.stats['duplicates_marked']} sequences to flag)"
        )

        return duplicate_groups

    def mark_duplicate_sequences(self, duplicate_groups: List[Tuple[str, List[str]]]):
        """
        Mark duplicate sequences by setting barcode_status = 'Duplicate'.
        """
        logger.info("Marking duplicate sequences...")

        if not duplicate_groups:
            return

        total_marked = 0
        sequences_table = self.db_manager.sequences_table

        with self.db_manager.get_session() as session:
            try:
                batch_size = 1000
                all_seq_ids_to_mark = []

                # Collect all seq_ids
                for _, seq_ids_to_mark in duplicate_groups:
                    all_seq_ids_to_mark.extend(seq_ids_to_mark)

                logger.info(f"Flagging {len(all_seq_ids_to_mark)} sequences as duplicates...")

                for i in range(0, len(all_seq_ids_to_mark), batch_size):
                    batch = all_seq_ids_to_mark[i:i + batch_size]

                    if self.db_manager.db_type == "postgresql":
                        update_query = text(f"""
                            UPDATE {sequences_table}
                            SET barcode_status = 'Duplicate'
                            WHERE seq_id = ANY(:seq_ids)
                            AND barcode_status != 'Duplicate'
                        """)
                        result = session.execute(update_query, {"seq_ids": batch})

                    else:  # SQLite
                        placeholders = ", ".join([f":id{j}" for j in range(len(batch))])
                        params = {f'id{j}': batch[j] for j in range(len(batch))}
                        update_query = text(f"""
                            UPDATE {sequences_table}
                            SET barcode_status = 'Duplicate'
                            WHERE seq_id IN ({placeholders})
                            AND barcode_status != 'Duplicate'
                        """)
                        result = session.execute(update_query, params)

                    total_marked += result.rowcount
                    logger.debug(f"Batch {i//batch_size + 1}: {result.rowcount} sequences marked")

                session.commit()
                logger.info(f"Marking complete. {total_marked} sequences marked as duplicates.")

            except Exception as e:
                session.rollback()
                logger.error(f"Error while marking duplicates: {e}", exc_info=True)
                raise

            self.stats['duplicates_marked'] = total_marked


    def verify_integrity(self) -> bool:
        """Verify integrity after deduplication."""
        logger.info("Performing final integrity checks...")

        sequences_table = self.db_manager.sequences_table
        seq_tax_table  = self.db_manager.seq_tax_table

        with self.db_manager.get_session() as session:
            # Ensure no unmarked duplicates remain
            result = session.execute(text(f"""
                SELECT COUNT(*) FROM (
                    SELECT s.sequence, s.barcode, st.tax_id, COUNT(*) as cnt
                    FROM {sequences_table} s
                    INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                    WHERE s.barcode_status != 'Duplicate'
                    GROUP BY s.sequence, s.barcode, st.tax_id
                    HAVING COUNT(*) > 1
                ) AS subquery
            """))
            remaining_duplicates = result.scalar()

            # Count sequences marked as duplicates
            result = session.execute(text(f"""
                SELECT COUNT(*) FROM {sequences_table}
                WHERE barcode_status = 'Duplicate'
            """))
            marked_duplicates = result.scalar()

            # Check for orphaned seq_tax rows
            result = session.execute(text(f"""
                SELECT COUNT(*) FROM {seq_tax_table} st
                LEFT JOIN {sequences_table} s ON st.seq_id = s.seq_id
                WHERE s.seq_id IS NULL
            """))
            orphaned_relations = result.scalar()

            logger.info(f"Integrity check: {remaining_duplicates} unmarked duplicates remaining, "
                       f"{marked_duplicates} sequences flagged as duplicates, "
                       f"{orphaned_relations} orphaned relations")

            if remaining_duplicates == 0 and orphaned_relations == 0:
                logger.info("✓ Integrity verification passed.")
                return True
            else:
                logger.warning(f"⚠ Integrity issues detected: "
                               f"{remaining_duplicates} unmarked duplicates, "
                               f"{orphaned_relations} orphaned relations.")
                return False

    def get_deduplication_summary(self) -> dict:
        """Build a detailed deduplication summary."""
        summary = {}

        sequences_table = self.db_manager.sequences_table
        seq_tax_table  = self.db_manager.seq_tax_table

        with self.db_manager.get_session() as session:
            # Stats per barcode_status
            result = session.execute(text(f"""
                SELECT barcode_status, COUNT(*) as count
                FROM {sequences_table}
                GROUP BY barcode_status
                ORDER BY count DESC
            """))
            summary['sequences_by_status'] = {row[0]: row[1] for row in result}

            # Stats per barcode
            result = session.execute(text(f"""
                SELECT barcode, COUNT(*) as count
                FROM {sequences_table}
                WHERE barcode_status = 'Duplicate'
                GROUP BY barcode
                ORDER BY count DESC
            """))
            summary['duplicates_by_barcode'] = {row[0]: row[1] for row in result}

            # Top 10 tax IDs with most duplicates
            result = session.execute(text(f"""
                SELECT st.tax_id, COUNT(*) as duplicate_count
                FROM {sequences_table} s
                INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                WHERE s.barcode_status = 'Duplicate'
                GROUP BY st.tax_id
                ORDER BY duplicate_count DESC
                LIMIT 10
            """))
            summary['top_duplicated_tax_ids'] = [
                {'tax_id': row[0], 'duplicate_count': row[1]}
                for row in result
            ]

        return summary

    def optimize_database(self):
        """Optimize the database after deduplication."""
        logger.info("Running database optimization (ANALYZE)...")

        if self.db_manager.db_type == "sqlite":
            with self.db_manager.get_session() as session:
                session.execute(text("ANALYZE"))
                session.commit()
        elif self.db_manager.db_type == "postgresql":
            # PostgreSQL requires AUTOCOMMIT for ANALYZE
            with self.db_manager.engine.connect().execution_options(isolation_level="AUTOCOMMIT") as conn:
                conn.execute(text("ANALYZE"))

        logger.info("Database optimization completed.")
    
    def deduplicate_sequences(self):
        """Full end-to-end deduplication workflow."""
        start_time = time.time()
        logger.info("Starting sequence deduplication")

        try:
            with self.db_manager.get_session() as session:
                # Use ORM query instead of raw SQL
                self.stats['initial_sequence_count'] = session.query(func.count(Sequence.seq_id)).scalar()
            logger.info(f"Initial sequence count: {self.stats['initial_sequence_count']}")

            # 0. Prepare indexes
            self.setup_database()

            # 1. Find duplicates
            duplicate_groups = self.find_duplicate_sequences()

            if not duplicate_groups:
                logger.info("No duplicates detected.")
                self.stats['processing_time'] = time.time() - start_time
                return

            # 2. Mark duplicates
            self.mark_duplicate_sequences(duplicate_groups)

            # 3. Verify integrity
            integrity_ok = self.verify_integrity()

            # 4. Update final stats
            with self.db_manager.get_session() as session:
                # Use ORM queries instead of raw SQL
                self.stats['final_sequence_count'] = session.query(func.count(Sequence.seq_id)).scalar()

                # Count unique (non-duplicate) sequences
                self.stats['unique_sequences_kept'] = session.query(
                    func.count(Sequence.seq_id)
                ).filter(
                    Sequence.barcode_status != 'Duplicate'
                ).scalar()

            # 5. Optimize database
            self.optimize_database()

        except Exception as e:
            logger.error(f"Error during deduplication: {e}", exc_info=True)
            raise
        finally:
            self.stats['processing_time'] = time.time() - start_time
            logger.info(f"Deduplication finished in {self.stats['processing_time']:.2f}s")
            logger.info(f"Sequences: {self.stats['initial_sequence_count']} total, "
                       f"{self.stats['unique_sequences_kept']} unique, "
                       f"{self.stats['duplicates_marked']} marked as duplicates")

def write_report(stats: dict, summary: dict, path: str, db_path: str):
    """
    Write a JSON report summarizing the sequence deduplication.

    Parameters
    ----------
    stats : dict
        Statistics returned by SequenceDeduplicator.
    summary : dict
        Detailed deduplication summary.
    path : str
        Output path for the JSON report.
    db_path : str
        Database path (for context).
    """
    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "database": os.path.abspath(db_path),
        "summary": {
            "initial_sequence_count": stats.get("initial_sequence_count"),
            "final_sequence_count": stats.get("final_sequence_count"),
            "unique_sequences_kept": stats.get("unique_sequences_kept"),
            "duplicate_groups": stats.get("duplicate_groups"),
            "duplicates_marked": stats.get("duplicates_marked"),
            "processing_time_sec": round(stats.get("processing_time", 0), 3),
        },
        "detailed_summary": summary
    }

    # Ensure parent directory exists
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    logger.info(f"Deduplication report written to: {path}")

def main():
    """Entry point (Snakemake or standalone)."""
    try:
        if 'snakemake' in globals():
            db_path = snakemake.params.db
            output_report = snakemake.output.report
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            logger.warning("Test mode - use Snakemake for production runs")
            db_path = "test_database.sqlite"
            output_report = "sequence_deduplication_report.json"
            config = {'database': {'path': db_path}}

        # Initialize database manager
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        deduplicator = SequenceDeduplicator(db_manager)
        deduplicator.deduplicate_sequences()

        # Generate detailed summary
        summary = deduplicator.get_deduplication_summary()

        print(f"\n{'='*60}")
        print("SEQUENCE DEDUPLICATION SUMMARY")
        print(f"{'='*60}")
        print(f"Initial sequences:       {deduplicator.stats['initial_sequence_count']}")
        print(f"Final sequences:         {deduplicator.stats['final_sequence_count']}")
        print(f"Unique sequences kept:   {deduplicator.stats['unique_sequences_kept']}")
        print(f"Duplicate groups:        {deduplicator.stats['duplicate_groups']}")
        print(f"Duplicates flagged:      {deduplicator.stats['duplicates_marked']}")
        print(f"Processing time:         {deduplicator.stats['processing_time']:.2f}s")
        print(f"{'='*60}")

        # Display additional stats
        if summary.get('sequences_by_status'):
            print("\nDistribution by barcode_status:")
            for status, count in summary['sequences_by_status'].items():
                print(f"  {status}: {count}")

        if summary.get('duplicates_by_barcode'):
            print("\nDuplicates by barcode:")
            for barcode, count in summary['duplicates_by_barcode'].items():
                print(f"  {barcode}: {count}")

        write_report(deduplicator.stats, summary, output_report, db_path)

    except Exception as e:
        logger.critical(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
