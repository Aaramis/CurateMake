#!/usr/bin/env python3
"""
Sequence deduplication script.
Marks duplicate sequences based on (sequence, barcode, tax_id).
"""

import os
import sys
import threading
import time
import logging
import json
from datetime import datetime
from typing import List, Tuple, Set, Dict, Any, Optional
from pathlib import Path
from sqlalchemy import func, text

sys.path.append(str(Path(__file__).parent / "utils"))
from db_models import Sequence

# Logging configuration
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


class SequenceDeduplicator:
    """Sequence deduplication manager."""

    def __init__(self, db_manager, config: Optional[Dict[str, Any]] = None,
                 allow_resume_shortcut: bool = True):
        self.db_manager = db_manager
        # Whether "duplicates already exist -> skip re-marking" may be taken. Safe for the
        # pre-ITSx pass (nothing else has marked anything yet), WRONG for the post-ITSx pass:
        # the pre pass always leaves Duplicate rows behind, so the post pass would conclude it
        # had already run and silently skip the several million ITSx children it alone can see.
        self.allow_resume_shortcut = allow_resume_shortcut
        config = config or {}
        dedup_config = config.get("deduplication", {})
        sequence_config = config.get("sequence_validation", {})

        self.create_indexes = dedup_config.get("create_indexes", True)
        self.skip_verify_integrity = dedup_config.get("skip_verify_integrity", False)
        self.barcode_log_every = max(1, int(dedup_config.get("barcode_log_every", 20)))
        self.update_batch_size = max(1, int(dedup_config.get("update_batch_size", 500)))
        self.tax_id_batch_size = max(1, int(dedup_config.get("tax_id_batch_size", 100)))

        # Duplicates per tax_id, accumulated while marking. Recomputing this from the
        # database afterwards means walking every duplicate row a second time.
        self._dup_by_tax: Dict[str, int] = {}

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
        """Check the indexes deduplication relies on, and rebuild them only if missing.

        Both are declared on the Sequence model, so create_database builds them while the
        tables are still empty. This method is therefore a no-op on any database the
        current schema created; the build path below exists only for databases predating
        that change, where the index has to be backfilled over rows that already exist.

        seq_tax lookups by seq_id/tax_id are already served by the composite PK (seq_id
        leading) and the ORM index idx_seq_tax_tax, so no seq_tax index is needed here.
        """
        logger.info("Verifying indexes for deduplication...")

        index_statements = [
            ("seq_barcode", f"CREATE INDEX IF NOT EXISTS idx_sequences_seq_barcode ON {self.db_manager.sequences_table}(sequence, barcode)"),
            ("dedup_composite", f"""
                CREATE INDEX IF NOT EXISTS idx_sequences_dedup_composite
                ON {self.db_manager.sequences_table}(sequence, barcode, barcode_status)
            """),
        ]

        for name, sql in index_statements:
            index_name = f"idx_sequences_{name}"
            if self._index_exists(index_name):
                logger.info("Index %s already present, skipping.", index_name)
                continue

            logger.warning("Index %s missing — this database predates its declaration on "
                           "the schema; backfilling it over existing rows may be slow.",
                           index_name)

            def create_index(stmt=sql):
                with self.db_manager.get_session() as session:
                    self._execute_ddl_bounded(session, stmt)
                    session.commit()

            self.db_manager._execute_with_retry(
                f"dedup_sequences_create_index_{name}",
                create_index,
            )

        logger.info("Indexes ready.")

    def _index_exists(self, index_name: str) -> bool:
        """Whether the index is already in the catalog, so a rebuild is not attempted.

        Cheap, and it is what makes an index created out-of-band (see below) enough to keep
        this step from ever blocking again.
        """
        if self.db_manager.db_type == "postgresql":
            sql = "SELECT 1 FROM pg_class WHERE relname = :n"
        else:
            sql = "SELECT 1 FROM sqlite_master WHERE type='index' AND name = :n"

        def check():
            with self.db_manager.get_session() as session:
                return session.execute(text(sql), {"n": index_name}).first() is not None

        try:
            return bool(self.db_manager._execute_with_retry(
                f"dedup_sequences_index_exists_{index_name}", check))
        except Exception:
            return False

    def _execute_ddl_bounded(self, session, stmt: str, budget_seconds: int = 1800):
        """Run a DDL statement with a wall-clock budget, cancelling it if it overruns.

        On YugabyteDB a CREATE INDEX is carried out by the master and is NOT bound by
        statement_timeout, and during a cluster stall the socket stays healthy so TCP
        keepalives never fire either. A blocked statement therefore hangs the client
        forever (observed: 62 minutes at 0 s of CPU). Cancelling on the raw connection is
        the only lever that works; raising afterwards lets _execute_with_retry try again.
        """
        if self.db_manager.db_type != "postgresql":
            session.execute(text(stmt))
            return

        raw = session.connection().connection
        box: Dict[str, Any] = {}

        def work():
            try:
                session.execute(text(stmt))
                box["ok"] = True
            except Exception as exc:
                box["err"] = exc

        worker = threading.Thread(target=work, daemon=True)
        worker.start()
        worker.join(budget_seconds)

        if worker.is_alive():
            logger.warning("DDL exceeded %ss, cancelling: %s", budget_seconds, stmt.strip()[:80])
            try:
                raw.cancel()
            except Exception:
                pass
            worker.join(120)
            raise TimeoutError(f"DDL cancelled after {budget_seconds}s")
        if "err" in box:
            raise box["err"]

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

    # yb_hash_code() spans this range; it is fixed by YugabyteDB.
    YB_HASH_SPACE = 65536
    # Starting slice count. Few wide slices beat many narrow ones: each slice pays a fixed
    # index-scan cost, so over-slicing multiplies it (measured: 64 slices -> 0.7 ms/id,
    # 512 slices -> 5.6 ms/id). Failures subdivide from here.
    HASH_SLICES = 64

    def _iter_seq_ids(self, extra_where: str, params: Dict[str, Any], label: str):
        """Yield lists of seq_id matching `extra_where`, in slices.

        Paginating with `ORDER BY seq_id LIMIT n` has no index to satisfy the ordering, so
        every page re-sorts the whole subset -- 25 ms per id, 16 h for a single barcode.
        Walking the table's own hash space instead needs no ordering, and each range is an
        independent, retryable query. A range the backend refuses is halved rather than
        abandoned. SQLite has neither yb_hash_code nor the problem, so it reads in one go.
        """
        sequences_table = self.db_manager.sequences_table

        if self.db_manager.db_type != "postgresql":
            query = f"SELECT seq_id FROM {sequences_table} WHERE {extra_where}"

            def fetch_all():
                with self.db_manager.get_session() as session:
                    return [row[0] for row in session.execute(text(query), params)]

            ids = self.db_manager._execute_with_retry(f"{label}_all", fetch_all)
            if ids:
                yield ids
            return

        def read_range(lo: int, hi: int, depth: int = 0):
            query = f"""
                SELECT seq_id FROM {sequences_table}
                WHERE {extra_where}
                  AND yb_hash_code(seq_id) >= :hash_lo
                  AND yb_hash_code(seq_id) < :hash_hi
            """

            def fetch():
                with self.db_manager.get_session() as session:
                    result = session.execute(text(query),
                                             {**params, "hash_lo": lo, "hash_hi": hi})
                    return [row[0] for row in result]

            try:
                # Two attempts, then split: a range that keeps failing is more likely too
                # large than transiently unlucky.
                last_exc = None
                for _ in range(2):
                    try:
                        return fetch()
                    except Exception as exc:  # noqa: BLE001
                        last_exc = exc
                        time.sleep(5)
                raise last_exc
            except Exception:
                if hi - lo <= 1:
                    raise
                mid = (lo + hi) // 2
                logger.warning("%s: hash range [%s,%s) refused, splitting", label, lo, hi)
                return read_range(lo, mid, depth + 1) + read_range(mid, hi, depth + 1)

        step = max(1, self.YB_HASH_SPACE // self.HASH_SLICES)
        for start in range(0, self.YB_HASH_SPACE, step):
            ids = read_range(start, min(start + step, self.YB_HASH_SPACE))
            if ids:
                yield ids

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

        # Step 1: enumerate the seq_ids of this barcode, sliced over the table's own
        # hash space (see _iter_seq_ids: an ORDER BY page would re-sort the whole subset).
        all_rows = []
        n_slices = 0

        for seq_ids in self._iter_seq_ids(
            f"barcode_status != 'Duplicate' AND {barcode_filter}", params,
            f"dedup_seq_ids_{barcode}"
        ):
            n_slices += 1
            logger.info(f"Slice {n_slices}: fetching data for {len(seq_ids)} sequences of barcode {barcode}...")

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
                    f"dedup_fetch_data_{barcode}_slice_{n_slices}",
                    fetch_sequences
                )
                all_rows.extend(rows)
            except Exception as e:
                logger.error(f"Failed to fetch data for slice {n_slices}: {e}")

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
            self._dup_by_tax[tax_id] = self._dup_by_tax.get(tax_id, 0) + 1
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
        """Most duplicated tax_ids, taken from the counts gathered while marking.

        This used to re-read every duplicate row from the database with an ordered keyset
        walk -- the same pattern that makes enumeration quadratic here, and it cost more
        than a third of the whole step for a ten-row summary. The marking pass already
        sees each duplicate exactly once, so the tally is free.
        """
        top = sorted(self._dup_by_tax.items(), key=lambda kv: kv[1], reverse=True)[:limit]
        return [{"tax_id": tax_id, "duplicate_count": count} for tax_id, count in top]

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

    def _duplicates_already_marked(self) -> bool:
        """True if a prior run already marked duplicates (a 'Duplicate' barcode_status exists)."""
        def check():
            with self.db_manager.get_session() as session:
                row = session.execute(text(
                    f"SELECT 1 FROM {self.db_manager.sequences_table} "
                    f"WHERE barcode_status = 'Duplicate' LIMIT 1"
                )).first()
                return row is not None
        return self.db_manager._execute_with_retry("dedup_check_already_marked", check)

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

            # Resumable restart: if a previous run already marked duplicates (persisted in the DB),
            # skip the expensive re-marking and go straight to summary/report. Lets a rule that died
            # in post-processing (e.g. the summary) be re-run cheaply without redoing the whole mark.
            if self.allow_resume_shortcut and self._duplicates_already_marked():
                logger.info("Duplicates already marked in the DB — skipping re-marking (resumable restart).")
            else:
                if not self.allow_resume_shortcut:
                    logger.info("Resume shortcut disabled for this pass — marking duplicates over "
                                "every unmarked row, including the ones ITSx has just added.")
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
            if not self.stats.get("duplicates_marked"):
                # re-marking was skipped (resumable restart) -> derive the count for the report
                self.stats["duplicates_marked"] = (self.stats.get("initial_sequence_count", 0) or 0) - unique_kept

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
            allow_resume_shortcut = bool(getattr(snakemake.params, "allow_resume_shortcut", True))
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            logger.addHandler(file_handler)
        else:
            logger.warning("Test mode - use Snakemake for production runs")
            db_path = "test_database.sqlite"
            output_report = "sequence_deduplication_report.json"
            config = {"database": {"path": db_path}}
            allow_resume_shortcut = True

        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)
        deduplicator = SequenceDeduplicator(db_manager, config=config,
                                            allow_resume_shortcut=allow_resume_shortcut)
        deduplicator.deduplicate_sequences()

        summary = deduplicator.get_deduplication_summary()
        write_report(deduplicator.stats, summary, output_report, db_path)

    except Exception as e:
        logger.critical(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
