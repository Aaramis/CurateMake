#!/usr/bin/env python3
"""
Extract CoL sequences from the ITS database for Snakemake runs
and emit separate FASTA files per barcode.

SQLAlchemy / DatabaseManager flavor:
- Relies on DatabaseManager (db_manager.py) with SQLAlchemy sessions.
- Works with both PostgreSQL and SQLite depending on Snakemake config.
"""

import json
import sys
import os
import time
from pathlib import Path
from collections import defaultdict
import logging
from typing import Dict, List, Tuple

from sqlalchemy import text

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SequenceExtractor:
    """Extract CoL sequences from the ITS database via DatabaseManager."""

    def __init__(self, db_manager, db_label: str):
        """
        Args:
            db_manager: DatabaseManager instance (PostgreSQL or SQLite).
            db_label: database label/path used purely for logs and reports.
        """
        self.db_manager = db_manager
        # For context in the report / logs only
        self.db_label = str(db_label)

    def extract_col_sequences(self) -> Dict[str, List[Tuple[int, str, str, int]]]:
        """
        Extract every sequence with backbone 'CoL' and barcode_status 'Pending'
        using backend-agnostic SQL (PostgreSQL/SQLite).
        """

        sequences_table = self.db_manager.sequences_table
        seq_tax_table = self.db_manager.seq_tax_table

        query = text(f"""
            SELECT s.seq_id,
                s.sequence,
                s.original_id,
                s.barcode,
                st.tax_id
            FROM {sequences_table} AS s
            INNER JOIN {seq_tax_table} AS st ON s.seq_id = st.seq_id
            WHERE st.backbone = 'CoL'
            AND s.barcode_status = 'Pending'
        """)

        sequences_by_barcode: Dict[str, List[Tuple[int, str, str, int]]] = defaultdict(list)

        def fetch_rows():
            with self.db_manager.get_session() as session:
                return list(session.execute(query))

        try:
            rows = self.db_manager._execute_with_retry(
                "itsx_extract_col_sequences",
                fetch_rows
            )

            for row in rows:
                seq_id = row.seq_id
                sequence = row.sequence
                original_id = row.original_id
                barcode = row.barcode
                tax_id = row.tax_id

                sequences_by_barcode[barcode].append(
                    (seq_id, sequence, original_id, tax_id)
                )

            logger.info("Sequences extracted per barcode:")
            for barcode, sequences in sequences_by_barcode.items():
                logger.info(f"  {barcode}: {len(sequences)} sequences")

            return dict(sequences_by_barcode)

        except Exception as e:
            logger.error(f"Error while extracting CoL sequences: {e}", exc_info=True)
            raise

    @staticmethod
    def write_fasta_files(
        sequences_by_barcode: Dict[str, List[Tuple[int, str, str, int]]],
        output_dir: str = "output"
    ) -> List[str]:
        """
        Write one FASTA file per barcode.

        Args:
            sequences_by_barcode: Mapping of barcode -> list of sequences.
            output_dir: Destination directory for FASTA files.

        Returns:
            List of created FASTA paths as strings.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        created_files: List[str] = []

        for barcode in sorted(sequences_by_barcode):
            sequences = sequences_by_barcode[barcode]
            # Output: [barcode].fasta
            filename = f"{barcode}.fasta"
            filepath = output_path / filename

            try:
                with filepath.open('w') as f:
                    for seq_id, sequence, original_id, tax_id in sorted(
                        sequences,
                        key=lambda row: (str(row[2]), str(row[1]))
                    ):
                        # Header uses seq_id and tax_id as requested
                        header = f">{seq_id}|tax_id:{tax_id}|original_id:{original_id}"
                        f.write(f"{header}\n")
                        f.write(f"{sequence}\n")

                logger.info(f"Created FASTA {filepath} with {len(sequences)} sequences")
                created_files.append(str(filepath))

            except Exception as e:
                logger.error(f"Failed to write {filename}: {e}", exc_info=True)
                raise

        return created_files

    def generate_stats(
        self,
        sequences_by_barcode: Dict[str, List[Tuple[int, str, str, int]]]
    ) -> Dict:
        """
        Build extraction statistics.

        Args:
            sequences_by_barcode: Mapping of barcode -> sequences.

        Returns:
            Dict containing global metrics and per-barcode breakdowns.
        """
        total_sequences = 0
        all_tax_ids: List[int] = []
        barcode_stats: Dict[str, Dict] = {}

        for barcode, sequences in sequences_by_barcode.items():
            total_sequences += len(sequences)

            lengths = [len(seq) for _, seq, _, _ in sequences]
            barcode_tax_ids = {tax_id for _, _, _, tax_id in sequences}
            all_tax_ids.extend(tax_id for _, _, _, tax_id in sequences)

            barcode_stats[barcode] = {
                'sequence_count': len(sequences),
                'unique_taxonomies': len(barcode_tax_ids),
                'length_stats': {
                    'min': min(lengths) if lengths else 0,
                    'max': max(lengths) if lengths else 0,
                    'mean': (sum(lengths) / len(lengths)) if lengths else 0.0,
                }
            }

        unique_tax_ids = len(set(all_tax_ids))

        return {
            'total_sequences': total_sequences,
            'unique_taxonomies': unique_tax_ids,
            'barcode_types': len(sequences_by_barcode),
            'barcode_details': barcode_stats
        }


def main():
    """Main entry point (Snakemake or standalone)."""
    start_time = time.time()

    try:
        # ================== Configuration via Snakemake or standalone ==================
        if 'snakemake' in globals():
            input_db = snakemake.params.db          # label ou chemin (pour le rapport)
            output_dir = str(Path(snakemake.params.output_fodler))
            output_report = snakemake.output.report
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            # File logging for Snakemake
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)

        else:
            # Mode test standalone
            logger.warning("Running in standalone test mode")
            input_db = "data/db/its.db"
            output_dir = "tmp/output"
            output_report = "tmp/output/itsx_report.json"

            # SQLite config compatible with DatabaseManager V2
            config = {
                'database': {
                    'type': 'sqlite',
                    'sqlite': {
                        'path': 'data/db',
                        'name': 'its.db'
                    }
                }
            }

            Path(output_dir).mkdir(parents=True, exist_ok=True)

        logger.info("Starting CoL sequence extraction...")

        # ================== Initialisation du DatabaseManager ==================
        sys.path.append(str(Path(__file__).parent.parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        # ================== Extraction ==================
        extractor = SequenceExtractor(db_manager, input_db)

        sequences_by_barcode = extractor.extract_col_sequences()

        if not sequences_by_barcode:
            logger.warning("No CoL sequences were found in the database")

            processing_time = time.time() - start_time

            report_data = {
                'source': 'CoL_extraction',
                'input_database': str(input_db),
                'output_directory': output_dir,
                'stats': {
                    'total_sequences': 0,
                    'unique_taxonomies': 0,
                    'barcode_types': 0,
                    'barcode_details': {}
                },
                'files_created': [],
                'processing_time': processing_time,
                'processing_rate_per_second': 0.0,
                'performance_metrics': {
                    'cpu_cores_used': os.cpu_count(),
                    'database_label': str(input_db)
                },
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'status': 'no_sequences_found'
            }

        else:
            # FASTA generation
            logger.info("Generating FASTA files...")
            created_files = extractor.write_fasta_files(sequences_by_barcode, output_dir)

            # Ensure every expected barcode has at least an empty FASTA placeholder, so
            # downstream ITSx barcode rules remain optional (no crash when a barcode
            # yields no sequences).
            expected_barcodes = config.get("msa", {}).get("barcodes", [])
            if expected_barcodes:
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)
                for barcode in expected_barcodes:
                    filepath = output_path / f"{barcode}.fasta"
                    if not filepath.exists():
                        filepath.touch()

            # Statistics generation
            logger.info("Computing extraction statistics...")
            stats = extractor.generate_stats(sequences_by_barcode)

            processing_time = time.time() - start_time
            total_sequences = stats['total_sequences']
            processing_rate = total_sequences / processing_time if processing_time > 0 else 0.0

            report_data = {
                'source': 'CoL_extraction',
                'input_database': str(input_db),
                'output_directory': output_dir,
                'stats': stats,
                'files_created': created_files,
                'processing_time': processing_time,
                'processing_rate_per_second': processing_rate,
                'performance_metrics': {
                    'cpu_cores_used': os.cpu_count(),
                    'database_label': str(input_db)
                },
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'status': 'success'
            }

            logger.info("Extraction completed successfully!")
            logger.info(f"Total sequences extracted: {total_sequences}")
            logger.info(f"FASTA files generated: {len(created_files)}")
            logger.info(f"Processing rate: {processing_rate:.0f} sequences/second")

        # ================== Sauvegarde du rapport ==================
        if 'snakemake' in globals() or output_report:
            Path(output_report).parent.mkdir(parents=True, exist_ok=True)
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Report saved to: {output_report}")

        # Standalone summary
        if 'snakemake' not in globals():
            print("\n=== SUMMARY ===")
            print(f"Total sequences extracted: {report_data['stats']['total_sequences']}")
            print(f"FASTA files generated: {len(report_data['files_created'])}")
            print(f"Output directory: {output_dir}")
            print(f"Processing time: {processing_time:.2f}s")

    except Exception as e:
        logger.critical(f"Fatal error during CoL extraction: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
