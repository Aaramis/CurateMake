#!/usr/bin/env python3
"""
Taxonomy deduplication script.
Detect duplicate taxonomy rows and update seq_tax references accordingly.
"""

import os
import sys
import time
import logging
import json
from datetime import datetime
from typing import List, Tuple, Dict, Any
from pathlib import Path
from sqlalchemy import func, text, Index

# Import ORM models
sys.path.append(str(Path(__file__).parent / "utils"))
from db_models import Taxonomy, SeqTax

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TaxonomyDeduplicator:
    """Taxonomy deduplication manager."""

    def __init__(self, db_manager):
        self.db_manager = db_manager
        self.stats = {
            'initial_taxonomy_count': 0,
            'final_taxonomy_count': 0,
            'duplicate_groups': 0,
            'duplicates_removed': 0,
            'sequences_updated': 0,
            'processing_time': 0.0
        }

    def setup_database(self):
        """Ensure necessary indexes exist (mostly handled by the ORM)."""
        logger.info("Checking indexes for deduplication...")
        # Note: Indexes are already defined in ORM models via __table_args__
        logger.info("Indexes already managed by ORM.")

    def find_duplicate_taxonomies(self) -> List[Tuple[int, List[int]]]:
        """
        Find duplicate taxonomies via a single SQL query.
        Returns tuples: (tax_id_to_keep, [ids_to_remove]).
        """
        logger.info("Searching taxonomy duplicates via SQL...")

        duplicate_groups = []

        # PostgreSQL uses STRING_AGG instead of GROUP_CONCAT
        if self.db_manager.db_type == "postgresql":
            group_concat_expr = "STRING_AGG(CAST(tax_id AS TEXT), ',')"
        else:
            group_concat_expr = "GROUP_CONCAT(tax_id)"

        # Group on COALESCE(col, '') so that a missing rank stored as NULL and the same
        # rank stored as '' collapse into ONE group. The two write paths disagree on the
        # convention (bulk_insert_taxonomies -> NULL, add_taxonomy -> ''), and a plain
        # GROUP BY treats NULL and '' as distinct, leaving such same-lineage rows unmerged.
        # COALESCE closes that gap and is portable (SQLite + PostgreSQL/YugabyteDB).
        group_by_expr = (
            "COALESCE(kingdom, ''), COALESCE(phylum, ''), COALESCE(class_, ''), "
            "COALESCE(order_, ''), COALESCE(family, ''), COALESCE(genus, ''), "
            "COALESCE(species, '')"
        )

        query = f"""
            SELECT
                MIN(tax_id) as keeper_id,
                {group_concat_expr} as all_ids,
                COUNT(*) as count
            FROM {self.db_manager.taxonomy_table}
            GROUP BY {group_by_expr}
            HAVING COUNT(*) > 1
        """

        def fetch_duplicates():
            groups = []
            with self.db_manager.get_session() as session:
                result = session.execute(text(query))

                for row in result:
                    keeper_id = row[0]
                    # Convert "uuid1,uuid2,uuid3" string into a list
                    all_ids = [id_str.strip() for id_str in row[1].split(',')]
                    # Remove every ID except the keeper
                    ids_to_remove = [id_val for id_val in all_ids if id_val != keeper_id]

                    if ids_to_remove:
                        groups.append((keeper_id, ids_to_remove))

            return groups

        duplicate_groups = self.db_manager._execute_with_retry(
            "dedup_taxonomies_find_duplicates",
            fetch_duplicates
        )

        self.stats['duplicate_groups'] = len(duplicate_groups)
        self.stats['duplicates_removed'] = sum(len(ids) for _, ids in duplicate_groups)

        logger.info(f"Found {self.stats['duplicate_groups']} duplicate groups "
                   f"({self.stats['duplicates_removed']} rows slated for removal)")

        return duplicate_groups


    def merge_duplicate_taxonomies(self, duplicate_groups: List[Tuple[str, List[str]]]):
        """
        Merge duplicate taxonomies by updating seq_tax references while
        respecting the unique constraint (seq_id, tax_id, backbone).
        Each UPDATE is retried individually to handle YugabyteDB RPC timeouts.
        """

        logger.info("Merging duplicate taxonomies...")

        if not duplicate_groups:
            return

        seq_tax_table = self.db_manager.seq_tax_table

        # 1. Build the global remapping dict old -> new
        id_remap = {
            old_id: keep_id
            for keep_id, ids_to_remove in duplicate_groups
            for old_id in ids_to_remove
        }

        all_ids_to_remove = list(id_remap.keys())
        total_sequences_updated = 0
        total_remaps = len(id_remap)

        logger.info(f"Updating seq_tax references for {total_remaps} taxonomies...")

        # 2. UPDATE each old tax_id -> new tax_id safely with individual retry
        for idx, (old_tax_id, new_tax_id) in enumerate(id_remap.items(), start=1):

            # Log progress every 100 updates or at the end
            if idx % 100 == 0 or idx == total_remaps:
                logger.info(f"  Processing taxonomy {idx}/{total_remaps}...")

            def merge_single_taxonomy():
                """Merge a single taxonomy with its own transaction."""
                with self.db_manager.get_session() as session:
                    # ⚠ Count affected rows before updating
                    count = session.query(func.count(SeqTax.seq_id)).filter(
                        SeqTax.tax_id == old_tax_id
                    ).scalar()

                    # ⚠ SAFE UPDATE to avoid PK collisions (NOT EXISTS)
                    safe_update_sql = text(f"""
                        UPDATE {seq_tax_table} AS st
                        SET tax_id = :new_tax_id
                        WHERE st.tax_id = :old_tax_id
                        AND NOT EXISTS (
                            SELECT 1
                            FROM {seq_tax_table} AS st2
                            WHERE st2.seq_id   = st.seq_id
                            AND st2.backbone = st.backbone
                            AND st2.tax_id   = :new_tax_id
                        )
                    """)

                    session.execute(
                        safe_update_sql,
                        {"new_tax_id": new_tax_id, "old_tax_id": old_tax_id}
                    )

                    # ⚠ Delete any leftover rows blocked by constraints
                    delete_remaining_sql = text(f"""
                        DELETE FROM {seq_tax_table}
                        WHERE tax_id = :old_tax_id
                    """)

                    session.execute(delete_remaining_sql, {"old_tax_id": old_tax_id})

                    return count

            # Execute with retry for this specific taxonomy
            count = self.db_manager._execute_with_retry(
                f"dedup_taxonomies_merge_{idx}",
                merge_single_taxonomy
            )
            total_sequences_updated += count

        logger.info("All seq_tax references updated safely.")
        logger.info("Seq_tax reference updates committed. Proceeding with batched taxonomy deletions...")
        self._delete_taxonomies_in_batches(all_ids_to_remove)
        logger.info("Taxonomy merge completed successfully.")

        self.stats['sequences_updated'] = total_sequences_updated

    def _delete_taxonomies_in_batches(self, tax_ids: List[str], batch_size: int = 2000):
        """
        Delete taxonomy rows in manageable batches to avoid long-running transactions.

        Parameters
        ----------
        tax_ids : List[str]
            Taxonomy identifiers slated for removal.
        batch_size : int
            Number of rows to delete per transaction.
        """
        if not tax_ids:
            logger.info("No duplicate taxonomy rows to delete.")
            return

        total = len(tax_ids)
        logger.info(f"Deleting {total} taxonomy rows in batches of {batch_size}...")

        for batch_index, start in enumerate(range(0, total, batch_size), start=1):
            batch = tax_ids[start:start + batch_size]

            def delete_op(batch_ids=batch):
                with self.db_manager.get_session() as session:
                    session.query(Taxonomy).filter(
                        Taxonomy.tax_id.in_(batch_ids)
                    ).delete(synchronize_session=False)

            if hasattr(self.db_manager, "_execute_with_retry"):
                operation_name = f"taxonomy_delete_batch_{batch_index}"
                self.db_manager._execute_with_retry(operation_name, delete_op)
            else:
                delete_op()

            processed = min(start + batch_size, total)
            if batch_index % 10 == 0 or processed == total:
                logger.info("  - Deleted %d/%d taxonomy rows", processed, total)

        logger.info("Duplicate taxonomy rows deleted in batches.")

    
    def verify_integrity(self) -> bool:
        """Verify integrity after deduplication."""
        logger.info("Performing final integrity checks...")

        def check_op():
            with self.db_manager.get_session() as session:
                # Ensure no duplicate groups remain (COALESCE so NULL and '' collapse,
                # matching the grouping used by find_duplicate_taxonomies).
                result = session.execute(text(f"""
                    SELECT COUNT(*) FROM (
                        SELECT 1 FROM {self.db_manager.taxonomy_table}
                        GROUP BY COALESCE(kingdom, ''), COALESCE(phylum, ''),
                                 COALESCE(class_, ''), COALESCE(order_, ''),
                                 COALESCE(family, ''), COALESCE(genus, ''),
                                 COALESCE(species, '')
                        HAVING COUNT(*) > 1
                    ) AS subquery
                """))
                remaining_duplicates = result.scalar()

                # Detect orphan seq_tax relations
                result = session.execute(text(f"""
                    SELECT COUNT(*) FROM {self.db_manager.seq_tax_table} st
                    LEFT JOIN {self.db_manager.taxonomy_table} t ON st.tax_id = t.tax_id
                    WHERE t.tax_id IS NULL
                """))
                orphaned_relations = result.scalar()

            return remaining_duplicates, orphaned_relations

        remaining_duplicates, orphaned_relations = self.db_manager._execute_with_retry(
            "dedup_taxonomies_verify",
            check_op
        )

        if remaining_duplicates == 0 and orphaned_relations == 0:
            logger.info("✓ Integrity verification passed.")
            return True
        else:
            logger.warning(f"⚠ Integrity issues detected: "
                           f"{remaining_duplicates} duplicate groups remaining, "
                           f"{orphaned_relations} orphaned relations.")
            return False

    def optimize_database(self):
        """Optimize the database after deduplication."""
        logger.info("Running database optimization (VACUUM)...")
        # Use the db_manager helper which handles SQLite and PostgreSQL
        self.db_manager._execute_with_retry(
            "dedup_taxonomies_vacuum",
            self.db_manager.vacuum_database
        )
        logger.info("Database optimization completed.")
    
    def deduplicate_taxonomies(self):
        """Execute the full taxonomy deduplication workflow."""
        start_time = time.time()
        logger.info("Starting taxonomy deduplication")
        self.stats['final_unique_tax_in_seq'] = 0

        try:
            def fetch_initial_count():
                with self.db_manager.get_session() as session:
                    # Use ORM query instead of raw SQL
                    return session.query(func.count(Taxonomy.tax_id)).scalar()

            self.stats['initial_taxonomy_count'] = self.db_manager._execute_with_retry(
                "dedup_taxonomies_initial_count",
                fetch_initial_count
            )
            logger.info(f"Initial taxonomy count: {self.stats['initial_taxonomy_count']}")

            # 0. Prepare indexes
            self.setup_database()

            # 1. Find duplicates
            duplicate_groups = self.find_duplicate_taxonomies()

            if not duplicate_groups:
                logger.info("No duplicate taxonomies detected.")
                self.stats['processing_time'] = time.time() - start_time
                return

            # 2. Merge duplicates
            self.merge_duplicate_taxonomies(duplicate_groups)

            # 3. Verify integrity
            integrity_ok = self.verify_integrity()

            # 4. Update final stats
            def fetch_final_stats():
                with self.db_manager.get_session() as session:
                    final_tax_count = session.query(func.count(Taxonomy.tax_id)).scalar()
                    final_unique_seq = session.query(
                        func.count(func.distinct(SeqTax.tax_id))
                    ).scalar()
                    return final_tax_count, final_unique_seq

            final_tax_count, final_unique_seq = self.db_manager._execute_with_retry(
                "dedup_taxonomies_final_stats",
                fetch_final_stats
            )
            self.stats['final_taxonomy_count'] = final_tax_count
            self.stats['final_unique_tax_in_seq'] = final_unique_seq
            logger.info(f"Final unique tax IDs referenced in seq_tax: {self.stats['final_unique_tax_in_seq']}")

            # 5. Optimize storage
            self.optimize_database()

        except Exception as e:
            logger.error(f"Error during taxonomy deduplication: {e}", exc_info=True)
            raise
        finally:
            self.stats['processing_time'] = time.time() - start_time
            logger.info(f"Deduplication finished in {self.stats['processing_time']:.2f}s")
            logger.info(f"Taxonomies: {self.stats['initial_taxonomy_count']} -> {self.stats['final_taxonomy_count']}")
            logger.info(f"Removed {self.stats['duplicates_removed']} duplicates, "
                       f"updated {self.stats['sequences_updated']} seq_tax relations.")

def write_report(stats: dict, path: str, db_path: str):
    """
    Write a JSON report summarizing the taxonomy deduplication.

    Parameters
    ----------
    stats : dict
        Statistics returned by TaxonomyDeduplicator.
    path : str
        Output file path.
    db_path : str
        Database path for context.
    """
    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "database": os.path.abspath(db_path),
        "summary": {
            "initial_taxonomy_count": stats.get("initial_taxonomy_count"),
            "final_taxonomy_count":   stats.get("final_taxonomy_count"),
            "final_unique_tax_id_in_seq_tax":   stats.get("final_unique_tax_in_seq"),
            "duplicate_groups":       stats.get("duplicate_groups"),
            "duplicates_removed":     stats.get("duplicates_removed"),
            "sequences_updated":      stats.get("sequences_updated"),
            "processing_time_sec":    round(stats.get("processing_time", 0), 3),
        }
    }

    # Ensure parent directory exists
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    logger.info(f"Taxonomy deduplication report written to: {path}")

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
            logger.warning("Test mode - Use Snakemake for production runs")
            db_path = "test_database.sqlite"
            output_report = "taxonomy_deduplication_report.json"
            config = {'database': {'path': db_path}}

        # Import et initialisation du database manager
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        deduplicator = TaxonomyDeduplicator(db_manager)
        deduplicator.deduplicate_taxonomies()

        print(f"\n{'='*60}")
        print("TAXONOMY DEDUPLICATION SUMMARY")
        print(f"{'='*60}")
        print(f"Initial taxonomies: {deduplicator.stats['initial_taxonomy_count']}")
        print(f"Final taxonomies:   {deduplicator.stats['final_taxonomy_count']}")
        print(f"Unique tax IDs in seq_tax: {deduplicator.stats['final_unique_tax_in_seq']}")
        print(f"Duplicate groups:   {deduplicator.stats['duplicate_groups']}")
        print(f"Duplicates removed: {deduplicator.stats['duplicates_removed']}")
        print(f"Relations updated:  {deduplicator.stats['sequences_updated']}")
        print(f"Processing time:    {deduplicator.stats['processing_time']:.2f}s")
        print(f"{'='*60}")

        write_report(deduplicator.stats, output_report, db_path)

    except Exception as e:
        logger.critical(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
