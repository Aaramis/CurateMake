#!/usr/bin/env python3
"""
Script d'extraction des séquences CoL depuis la base de données ITS pour Snakemake
Génère des fichiers FASTA séparés par type de barcode.

Version SQLAlchemy / DatabaseManager :
- Utilise DatabaseManager (db_manager.py) et SQLAlchemy.
- Compatible PostgreSQL et SQLite via la config Snakemake.
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

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SequenceExtractor:
    """Extracteur de séquences depuis la base de données ITS via DatabaseManager."""

    def __init__(self, db_manager, db_label: str):
        """
        Args:
            db_manager: instance de DatabaseManager (PostgreSQL ou SQLite).
            db_label: chemin ou label de la base (uniquement pour le rapport).
        """
        self.db_manager = db_manager
        # Juste pour le contexte dans le rapport / logs
        self.db_label = str(db_label)

    def extract_col_sequences(self) -> Dict[str, List[Tuple[int, str, str, int]]]:
        """
        Extrait toutes les séquences avec backbone 'CoL' et barcode_status 'Pending',
        en utilisant du SQL brut compatible PostgreSQL/SQLite et les bons noms de tables.
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

        try:
            with self.db_manager.get_session() as session:
                result = session.execute(query)

                for row in result:
                    seq_id = row.seq_id
                    sequence = row.sequence
                    original_id = row.original_id
                    barcode = row.barcode
                    tax_id = row.tax_id

                    sequences_by_barcode[barcode].append(
                        (seq_id, sequence, original_id, tax_id)
                    )

            logger.info("Séquences extraites par barcode :")
            for barcode, sequences in sequences_by_barcode.items():
                logger.info(f"  {barcode}: {len(sequences)} séquences")

            return dict(sequences_by_barcode)

        except Exception as e:
            logger.error(f"Erreur lors de l'extraction : {e}", exc_info=True)
            raise

    @staticmethod
    def write_fasta_files(
        sequences_by_barcode: Dict[str, List[Tuple[int, str, str, int]]],
        output_dir: str = "output"
    ) -> List[str]:
        """
        Écrit les fichiers FASTA pour chaque barcode.

        Args
        ----
        sequences_by_barcode :
            Dict des séquences groupées par barcode.
        output_dir :
            Répertoire de sortie.

        Returns
        -------
        List[str]
            Liste des fichiers créés.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        created_files: List[str] = []

        for barcode, sequences in sequences_by_barcode.items():
            # Nom du fichier : [barcode].fasta
            filename = f"{barcode}.fasta"
            filepath = output_path / filename

            try:
                with filepath.open('w') as f:
                    for seq_id, sequence, original_id, tax_id in sequences:
                        # Header avec seq_id et tax_id comme demandé
                        header = f">{seq_id}|tax_id:{tax_id}|original_id:{original_id}"
                        f.write(f"{header}\n")
                        f.write(f"{sequence}\n")

                logger.info(f"Fichier créé: {filepath} ({len(sequences)} séquences)")
                created_files.append(str(filepath))

            except Exception as e:
                logger.error(f"Erreur lors de l'écriture de {filename} : {e}", exc_info=True)
                raise

        return created_files

    def generate_stats(
        self,
        sequences_by_barcode: Dict[str, List[Tuple[int, str, str, int]]]
    ) -> Dict:
        """
        Génère les statistiques d'extraction.

        Args
        ----
        sequences_by_barcode :
            Dict des séquences groupées par barcode.

        Returns
        -------
        Dict
            Statistiques globales + détaillées par barcode.
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
    """Point d'entrée principal."""
    start_time = time.time()

    try:
        # ================== Configuration via Snakemake ou standalone ==================
        if 'snakemake' in globals():
            input_db = snakemake.params.db          # label ou chemin (pour le rapport)
            output_dir = str(Path(snakemake.params.output_fodler))
            output_report = snakemake.output.report
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            # Logging fichier pour Snakemake
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)

        else:
            # Mode test standalone
            logger.warning("Running in standalone test mode")
            input_db = "data/db/its.db"
            output_dir = "tmp/output"
            output_report = "tmp/output/itsx_report.json"

            # Config SQLite compatible avec DatabaseManager V2
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

        logger.info("Début de l'extraction des séquences CoL...")

        # ================== Initialisation du DatabaseManager ==================
        sys.path.append(str(Path(__file__).parent.parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        # ================== Extraction ==================
        extractor = SequenceExtractor(db_manager, input_db)

        sequences_by_barcode = extractor.extract_col_sequences()

        if not sequences_by_barcode:
            logger.warning("Aucune séquence CoL trouvée dans la base de données")

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
            # Génération des fichiers FASTA
            logger.info("Génération des fichiers FASTA...")
            created_files = extractor.write_fasta_files(sequences_by_barcode, output_dir)

            # Génération des statistiques
            logger.info("Calcul des statistiques...")
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

            logger.info("Extraction terminée avec succès !")
            logger.info(f"Total des séquences extraites : {total_sequences}")
            logger.info(f"Fichiers FASTA générés : {len(created_files)}")
            logger.info(f"Taux de traitement : {processing_rate:.0f} séquences/seconde")

        # ================== Sauvegarde du rapport ==================
        if 'snakemake' in globals() or output_report:
            Path(output_report).parent.mkdir(parents=True, exist_ok=True)
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Rapport sauvegardé : {output_report}")

        # Résumé pour mode standalone
        if 'snakemake' not in globals():
            print("\n=== RÉSUMÉ ===")
            print(f"Total des séquences extraites : {report_data['stats']['total_sequences']}")
            print(f"Fichiers FASTA générés : {len(report_data['files_created'])}")
            print(f"Répertoire de sortie : {output_dir}")
            print(f"Temps de traitement : {processing_time:.2f}s")

    except Exception as e:
        logger.critical(f"Erreur fatale lors de l'extraction CoL: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
