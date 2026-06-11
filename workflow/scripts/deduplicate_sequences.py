#!/usr/bin/env python3
"""
Sequence deduplication script.
Marks duplicate sequences based on (sequence, barcode, tax_id).
"""

import os
import sys
import time
import logging
import json
from datetime import datetime
from typing import List, Tuple, Set, Dict, Any, Optional
from pathlib import Path
from sqlalchemy import func, text

# Import ORM models
sys.path.append(str(Path(__file__).parent / "utils"))
from db_models import Sequence

# Logging configuration
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


class SequenceDeduplicator:
    """Sequence deduplication manager."""

    def __init__(self, db_manager, config: Optional[Dict[str, Any]] = None):
        self.db_manager = db_manager
        config = config or {}
        dedup_config = config.get("deduplication", {})
        sequence_config = config.get("sequence_validation", {})

        self.create_indexes = dedup_config.get("create_indexes", True)
        self.skip_verify_integrity = dedup_config.get("skip_verify_integrity", False)
        self.barcode_log_every = max(1, int(dedup_config.get("barcode_log_every", 20)))
        self.update_batch_size = max(1, int(dedup_config.get("update_batch_size", 500)))
        self.tax_id_batch_size = max(1, int(dedup_config.get("tax_id_batch_size", 100)))

        self.stats = {
            "initial_sequence_count": 0,
            "final_sequence_count": 0,
            "duplicate_groups": 0,
            "duplicates_marked": 0,
            "duplicates_to_mark": 0,
            "unique_sequences_kept": 0,
            "processing_time": 0.0,
        }

    def setup_database(self):
        """Create indexes required for performant deduplication."""
        logger.info("Verifying/creating indexes for deduplication...")

        index_statements = [
            ("seq_barcode", f"CREATE INDEX IF NOT EXISTS idx_sequences_seq_barcode ON {self.db_manager.sequences_table}(sequence, barcode)"),
            ("seq_tax_seq", f"CREATE INDEX IF NOT EXISTS idx_seq_tax_seq_id ON {self.db_manager.seq_tax_table}(seq_id)"),
            ("seq_tax_tax", f"CREATE INDEX IF NOT EXISTS idx_seq_tax_tax_id ON {self.db_manager.seq_tax_table}(tax_id)"),
            ("dedup_composite", f"""
                CREATE INDEX IF NOT EXISTS idx_sequences_dedup_composite
                ON {self.db_manager.sequences_table}(sequence, barcode, barcode_status)
            """),
        ]

        for name, sql in index_statements:
            def create_index(stmt=sql):
                with self.db_manager.get_session() as session:
                    session.execute(text(stmt))
                    session.commit()

            self.db_manager._execute_with_retry(
                f"dedup_sequences_create_index_{name}",
                create_index,
            )

        logger.info("Indexes ready.")

    def _fetch_barcodes(self) -> List[Optional[str]]:
        sequences_table = self.db_manager.sequences_table

        def fetch():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT DISTINCT barcode
                    FROM {sequences_table}
                    WHERE barcode_status != 'Duplicate'
                    ORDER BY barcode
                """))
                return [row[0] for row in result]

        return self.db_manager._execute_with_retry(
            "dedup_sequences_fetch_barcodes",
            fetch,
        )

    def _update_sequences_by_ids(self, seq_ids: List[str]) -> int:
        if not seq_ids:
            return 0

        sequences_table = self.db_manager.sequences_table
        if self.db_manager.db_type == "postgresql":
            update_sql = text(f"""
                UPDATE {sequences_table}
                SET barcode_status = 'Duplicate'
                WHERE seq_id = ANY(:seq_ids)
                  AND barcode_status != 'Duplicate'
            """)
            params = {"seq_ids": list(seq_ids)}

            def apply_update():
                with self.db_manager.get_session() as session:
                    result = session.execute(update_sql, params)
                    session.commit()
                    return result.rowcount

            return self.db_manager._execute_with_retry(
                "dedup_sequences_mark_batch",
                apply_update,
            )

        placeholders = ", ".join([f":id{i}" for i in range(len(seq_ids))])
        params = {f"id{i}": seq_id for i, seq_id in enumerate(seq_ids)}
        update_sql = text(f"""
            UPDATE {sequences_table}
            SET barcode_status = 'Duplicate'
            WHERE seq_id IN ({placeholders})
              AND barcode_status != 'Duplicate'
        """)

        def apply_update():
            with self.db_manager.get_session() as session:
                result = session.execute(update_sql, params)
                session.commit()
                return result.rowcount

        return self.db_manager._execute_with_retry(
            "dedup_sequences_mark_batch",
            apply_update,
        )

    def _mark_duplicates_for_barcode_streaming(self, barcode: Optional[str]) -> Tuple[int, int, int]:
        """
        Mark duplicate sequences for a barcode using seq_id batching:
        1. Fetch seq_ids for the barcode in batches (no JOIN, avoids temp_file_limit)
        2. For each seq_id batch, fetch sequences with their tax_ids
        3. Accumulate all data, sort in Python, mark duplicates

        This completely avoids large JOINs that exceed YugabyteDB's temp_file_limit.
        """
        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table

        if barcode is None:
            barcode_filter = "barcode IS NULL"
            params = {}
        else:
            barcode_filter = "barcode = :barcode"
            params = {"barcode": barcode}

        # Step 1: Fetch seq_ids for this barcode in batches (no JOIN)
        # We use OFFSET/LIMIT to stream through seq_ids without loading all at once
        seq_id_batch_size = self.tax_id_batch_size * 10  # Larger batches for seq_ids
        offset = 0
        all_rows = []

        while True:
            seq_id_query = f"""
                SELECT seq_id
                FROM {sequences_table}
                WHERE barcode_status != 'Duplicate'
                  AND {barcode_filter}
                LIMIT :limit OFFSET :offset
            """

            def fetch_seq_ids():
                with self.db_manager.get_session() as session:
                    result = session.execute(
                        text(seq_id_query),
                        {**params, "limit": seq_id_batch_size, "offset": offset}
                    )
                    return [row[0] for row in result]

            try:
                seq_ids = self.db_manager._execute_with_retry(
                    f"dedup_fetch_seq_ids_{barcode}_offset_{offset}",
                    fetch_seq_ids
                )
            except Exception as e:
                logger.error(f"Failed to fetch seq_ids at offset {offset}: {e}")
                break

            if not seq_ids:
                break

            logger.info(f"Fetching data for {len(seq_ids)} sequences (offset {offset}) for barcode {barcode}...")

            # Step 2: Fetch sequences and tax_ids for this batch of seq_ids
            if self.db_manager.db_type == "postgresql":
                seq_filter = "s.seq_id = ANY(:seq_ids)"
                seq_params = {"seq_ids": seq_ids}
            else:
                placeholders = ", ".join([f":sid{i}" for i in range(len(seq_ids))])
                seq_filter = f"s.seq_id IN ({placeholders})"
                seq_params = {f"sid{i}": sid for i, sid in enumerate(seq_ids)}

            data_query = f"""
                SELECT s.seq_id, s.sequence, st.tax_id
                FROM {sequences_table} s
                INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                WHERE {seq_filter}
            """

            def fetch_sequences():
                with self.db_manager.get_session() as session:
                    result = session.execute(text(data_query), seq_params)
                    return list(result)

            try:
                rows = self.db_manager._execute_with_retry(
                    f"dedup_fetch_data_{barcode}_offset_{offset}",
                    fetch_sequences
                )
                all_rows.extend(rows)
            except Exception as e:
                logger.error(f"Failed to fetch data for seq_id batch at offset {offset}: {e}")

            offset += seq_id_batch_size

        if not all_rows:
            logger.info(f"No sequences found for barcode {barcode}")
            return 0, 0, 0

        logger.info(f"Collected {len(all_rows)} sequence records for barcode {barcode}, sorting...")

        # Step 3: Sort all data in Python by (tax_id, sequence, seq_id)
        all_rows.sort(key=lambda r: (r[2], r[1], r[0]))

        logger.info(f"Sorted, now detecting duplicates...")

        # Step 4: Detect duplicates and mark them
        batch_ids: List[str] = []
        batch_seen: Set[str] = set()
        total_marked = 0
        total_groups = 0
        total_to_mark = 0

        def flush_batch():
            nonlocal total_marked, batch_ids, batch_seen
            if not batch_ids:
                return
            total_marked += self._update_sequences_by_ids(batch_ids)
            batch_ids = []
            batch_seen = set()

        current_key = None
        keep_seq_id = None
        group_size = 0

        def close_group():
            nonlocal total_groups, total_to_mark, group_size
            if group_size > 1:
                total_groups += 1
                total_to_mark += group_size - 1
            group_size = 0

        for seq_id, sequence, tax_id in all_rows:
            key = (tax_id, sequence)
            if key != current_key:
                if current_key is not None:
                    close_group()
                current_key = key
                keep_seq_id = seq_id
                group_size = 1
                continue

            group_size += 1
            if seq_id == keep_seq_id:
                continue
            if seq_id in batch_seen:
                continue
            batch_seen.add(seq_id)
            batch_ids.append(seq_id)
            batch_limit = self.update_batch_size
            if self.db_manager.db_type == "sqlite":
                batch_limit = min(batch_limit, 900)
            if len(batch_ids) >= batch_limit:
                flush_batch()

        if current_key is not None:
            close_group()

        flush_batch()
        return total_marked, total_groups, total_to_mark

    def mark_duplicate_sequences(self) -> None:
        """Mark duplicate sequences by streaming through sorted groups."""
        logger.info("Marking duplicate sequences...")

        barcodes = self._fetch_barcodes()
        total_barcodes = len(barcodes)
        total_marked = 0
        total_groups = 0
        total_to_mark = 0

        for idx, barcode in enumerate(barcodes, start=1):
            if idx == 1 or idx % self.barcode_log_every == 0 or idx == total_barcodes:
                logger.info(f"Barcode stream {idx}/{total_barcodes}...")
            marked, groups, to_mark = self._mark_duplicates_for_barcode_streaming(barcode)
            total_marked += marked
            total_groups += groups
            total_to_mark += to_mark

        self.stats["duplicates_marked"] = total_marked
        self.stats["duplicate_groups"] = total_groups
        self.stats["duplicates_to_mark"] = total_to_mark
        logger.info(f"Marking complete. {total_marked} sequences marked as duplicates.")

    def _count_remaining_duplicates_for_barcode_streaming(self, barcode: Optional[str]) -> int:
        """
        Count remaining duplicate groups using seq_id batching (same as marking).
        """
        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table

        if barcode is None:
            barcode_filter = "barcode IS NULL"
            params = {}
        else:
            barcode_filter = "barcode = :barcode"
            params = {"barcode": barcode}

        # Step 1: Fetch seq_ids in batches and accumulate all data
        seq_id_batch_size = self.tax_id_batch_size * 10
        offset = 0
        all_rows = []

        while True:
            seq_id_query = f"""
                SELECT seq_id
                FROM {sequences_table}
                WHERE barcode_status != 'Duplicate'
                  AND {barcode_filter}
                LIMIT :limit OFFSET :offset
            """

            def fetch_seq_ids():
                with self.db_manager.get_session() as session:
                    result = session.execute(
                        text(seq_id_query),
                        {**params, "limit": seq_id_batch_size, "offset": offset}
                    )
                    return [row[0] for row in result]

            try:
                seq_ids = self.db_manager._execute_with_retry(
                    f"dedup_verify_seq_ids_{barcode}_offset_{offset}",
                    fetch_seq_ids
                )
            except Exception as e:
                logger.error(f"Failed to fetch seq_ids for verification at offset {offset}: {e}")
                break

            if not seq_ids:
                break

            # Fetch sequences and tax_ids for this batch
            if self.db_manager.db_type == "postgresql":
                seq_filter = "s.seq_id = ANY(:seq_ids)"
                seq_params = {"seq_ids": seq_ids}
            else:
                placeholders = ", ".join([f":sid{i}" for i in range(len(seq_ids))])
                seq_filter = f"s.seq_id IN ({placeholders})"
                seq_params = {f"sid{i}": sid for i, sid in enumerate(seq_ids)}

            data_query = f"""
                SELECT s.sequence, st.tax_id
                FROM {sequences_table} s
                INNER JOIN {seq_tax_table} st ON s.seq_id = st.seq_id
                WHERE {seq_filter}
            """

            def fetch_sequences():
                with self.db_manager.get_session() as session:
                    result = session.execute(text(data_query), seq_params)
                    return list(result)

            try:
                rows = self.db_manager._execute_with_retry(
                    f"dedup_verify_data_{barcode}_offset_{offset}",
                    fetch_sequences
                )
                all_rows.extend(rows)
            except Exception as e:
                logger.error(f"Failed to fetch data for verification at offset {offset}: {e}")

            offset += seq_id_batch_size

        if not all_rows:
            return 0

        # Sort in Python
        all_rows.sort(key=lambda r: (r[1], r[0]))  # (tax_id, sequence)

        # Count duplicate groups
        total_groups = 0
        current_key = None
        group_size = 0

        def close_group():
            nonlocal total_groups, group_size
            if group_size > 1:
                total_groups += 1
            group_size = 0

        for sequence, tax_id in all_rows:
            key = (tax_id, sequence)
            if key != current_key:
                if current_key is not None:
                    close_group()
                current_key = key
                group_size = 1
                continue
            group_size += 1

        if current_key is not None:
            close_group()

        return total_groups

    def verify_integrity(self) -> bool:
        """Verify integrity after deduplication using streaming checks."""
        logger.info("Performing final integrity checks...")

        barcodes = self._fetch_barcodes()
        total_barcodes = len(barcodes)
        remaining_duplicates = 0
        for idx, barcode in enumerate(barcodes, start=1):
            if idx == 1 or idx % self.barcode_log_every == 0 or idx == total_barcodes:
                logger.info(f"Barcode verify {idx}/{total_barcodes}...")
            remaining_duplicates += self._count_remaining_duplicates_for_barcode_streaming(barcode)

        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table

        def count_marked_duplicates():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT COUNT(*) FROM {sequences_table}
                    WHERE barcode_status = 'Duplicate'
                """))
                return result.scalar()

        def count_orphaned_relations():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT COUNT(*) FROM {seq_tax_table} st
                    LEFT JOIN {sequences_table} s ON st.seq_id = s.seq_id
                    WHERE s.seq_id IS NULL
                """))
                return result.scalar()

        marked_duplicates = self.db_manager._execute_with_retry(
            "dedup_sequences_verify_marked",
            count_marked_duplicates,
        )
        orphaned_relations = self.db_manager._execute_with_retry(
            "dedup_sequences_verify_orphans",
            count_orphaned_relations,
        )

        logger.info(
            "Integrity check: %s unmarked duplicates remaining, %s sequences flagged as duplicates, %s orphaned relations",
            remaining_duplicates,
            marked_duplicates,
            orphaned_relations,
        )

        if remaining_duplicates == 0 and orphaned_relations == 0:
            logger.info("✓ Integrity verification passed.")
            return True

        logger.warning(
            "⚠ Integrity issues detected: %s unmarked duplicates, %s orphaned relations.",
            remaining_duplicates,
            orphaned_relations,
        )
        return False

    def _get_top_duplicated_tax_ids_streaming(self, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Get top duplicated tax_ids using batched streaming to avoid temp_file_limit.
        Fetches duplicate seq_ids in batches, then aggregates in Python.
        """
        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table

        # Count duplicates per tax_id using batched approach
        tax_id_counts: Dict[str, int] = {}
        batch_size = self.tax_id_batch_size * 10
        offset = 0

        while True:
            # Fetch duplicate seq_ids in batches
            seq_id_query = f"""
                SELECT seq_id
                FROM {sequences_table}
                WHERE barcode_status = 'Duplicate'
                LIMIT :limit OFFSET :offset
            """

            def fetch_seq_ids():
                with self.db_manager.get_session() as session:
                    result = session.execute(
                        text(seq_id_query),
                        {"limit": batch_size, "offset": offset}
                    )
                    return [row[0] for row in result]

            try:
                seq_ids = self.db_manager._execute_with_retry(
                    f"dedup_summary_seq_ids_offset_{offset}",
                    fetch_seq_ids
                )
            except Exception as e:
                logger.warning(f"Failed to fetch seq_ids for summary at offset {offset}: {e}")
                break

            if not seq_ids:
                break

            # Fetch tax_ids for these seq_ids
            if self.db_manager.db_type == "postgresql":
                tax_query = f"""
                    SELECT tax_id FROM {seq_tax_table}
                    WHERE seq_id = ANY(:seq_ids)
                """
                tax_params = {"seq_ids": seq_ids}
            else:
                placeholders = ", ".join([f":sid{i}" for i in range(len(seq_ids))])
                tax_query = f"""
                    SELECT tax_id FROM {seq_tax_table}
                    WHERE seq_id IN ({placeholders})
                """
                tax_params = {f"sid{i}": sid for i, sid in enumerate(seq_ids)}

            def fetch_tax_ids():
                with self.db_manager.get_session() as session:
                    result = session.execute(text(tax_query), tax_params)
                    return [row[0] for row in result]

            try:
                tax_ids = self.db_manager._execute_with_retry(
                    f"dedup_summary_tax_ids_offset_{offset}",
                    fetch_tax_ids
                )
                for tax_id in tax_ids:
                    tax_id_counts[tax_id] = tax_id_counts.get(tax_id, 0) + 1
            except Exception as e:
                logger.warning(f"Failed to fetch tax_ids for summary at offset {offset}: {e}")

            offset += batch_size

        # Sort and return top N
        sorted_counts = sorted(tax_id_counts.items(), key=lambda x: x[1], reverse=True)[:limit]
        return [{"tax_id": tax_id, "duplicate_count": count} for tax_id, count in sorted_counts]

    def get_deduplication_summary(self) -> dict:
        """Build a detailed deduplication summary."""
        summary = {}

        sequences_table = self.db_manager.sequences_table

        def build_basic_summary():
            summary_local = {}
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT barcode_status, COUNT(*) as count
                    FROM {sequences_table}
                    GROUP BY barcode_status
                    ORDER BY count DESC
                """))
                summary_local["sequences_by_status"] = {row[0]: row[1] for row in result}

                result = session.execute(text(f"""
                    SELECT barcode, COUNT(*) as count
                    FROM {sequences_table}
                    WHERE barcode_status = 'Duplicate'
                    GROUP BY barcode
                    ORDER BY count DESC
                """))
                summary_local["duplicates_by_barcode"] = {row[0]: row[1] for row in result}

            return summary_local

        summary.update(self.db_manager._execute_with_retry(
            "dedup_sequences_summary",
            build_basic_summary,
        ))

        # Get top duplicated tax_ids using streaming approach to avoid temp_file_limit
        logger.info("Computing top duplicated tax_ids (streaming)...")
        summary["top_duplicated_tax_ids"] = self._get_top_duplicated_tax_ids_streaming(limit=10)

        return summary

    def optimize_database(self):
        """Optimize the database after deduplication."""
        logger.info("Running database optimization (ANALYZE)...")

        if self.db_manager.db_type == "sqlite":
            def analyze_op():
                with self.db_manager.get_session() as session:
                    session.execute(text("ANALYZE"))
                    session.commit()

            self.db_manager._execute_with_retry("dedup_sequences_analyze", analyze_op)
        elif self.db_manager.db_type == "postgresql":
            with self.db_manager.engine.connect().execution_options(isolation_level="AUTOCOMMIT") as conn:
                conn.execute(text("ANALYZE"))

        logger.info("Database optimization completed.")

    def deduplicate_sequences(self):
        """Full end-to-end deduplication workflow."""
        start_time = time.time()
        logger.info("Starting sequence deduplication")

        try:
            def fetch_initial_count():
                with self.db_manager.get_session() as session:
                    return session.query(func.count(Sequence.seq_id)).scalar()

            self.stats["initial_sequence_count"] = self.db_manager._execute_with_retry(
                "dedup_sequences_initial_count",
                fetch_initial_count,
            )
            logger.info(f"Initial sequence count: {self.stats['initial_sequence_count']}")

            if self.create_indexes:
                self.setup_database()
            else:
                logger.info("Skipping index creation for deduplication (deduplication.create_indexes=false)")

            self.mark_duplicate_sequences()
            if self.skip_verify_integrity:
                logger.info("Skipping integrity verification (deduplication.skip_verify_integrity=true)")
            else:
                self.verify_integrity()

            def fetch_final_stats():
                with self.db_manager.get_session() as session:
                    final_count = session.query(func.count(Sequence.seq_id)).scalar()
                    unique_kept = session.query(
                        func.count(Sequence.seq_id)
                    ).filter(
                        Sequence.barcode_status != "Duplicate"
                    ).scalar()
                    return final_count, unique_kept

            final_count, unique_kept = self.db_manager._execute_with_retry(
                "dedup_sequences_final_stats",
                fetch_final_stats,
            )
            self.stats["final_sequence_count"] = final_count
            self.stats["unique_sequences_kept"] = unique_kept

            self.optimize_database()

        except Exception as e:
            logger.error(f"Error during deduplication: {e}", exc_info=True)
            raise
        finally:
            self.stats["processing_time"] = time.time() - start_time
            logger.info(f"Deduplication finished in {self.stats['processing_time']:.2f}s")
            logger.info(
                "Sequences: %s total, %s unique, %s marked as duplicates",
                self.stats["initial_sequence_count"],
                self.stats["unique_sequences_kept"],
                self.stats["duplicates_marked"],
            )


def write_report(stats: dict, summary: dict, path: str, db_path: str):
    """Write a JSON report summarizing the sequence deduplication."""
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
        "detailed_summary": summary,
    }

    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    logger.info(f"Deduplication report written to: {path}")


def main():
    """Entry point (Snakemake or standalone)."""
    try:
        if "snakemake" in globals():
            db_path = snakemake.params.db
            output_report = snakemake.output.report
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            logger.addHandler(file_handler)
        else:
            logger.warning("Test mode - use Snakemake for production runs")
            db_path = "test_database.sqlite"
            output_report = "sequence_deduplication_report.json"
            config = {"database": {"path": db_path}}

        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)
        deduplicator = SequenceDeduplicator(db_manager, config=config)
        deduplicator.deduplicate_sequences()

        summary = deduplicator.get_deduplication_summary()
        write_report(deduplicator.stats, summary, output_report, db_path)

    except Exception as e:
        logger.critical(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
