#!/usr/bin/env python3
"""
Script d'extraction des séquences CoL depuis la base de données ITS pour Snakemake
Génère des fichiers FASTA séparés par type de barcode
"""

import sqlite3
import json
import sys
import os
import time
from pathlib import Path
from collections import defaultdict
import logging
from typing import Dict, List, Tuple

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SequenceExtractor:
    """Extracteur de séquences depuis la base de données ITS"""
    
    def __init__(self, db_path: str):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Base de données introuvable: {db_path}")
    
    def extract_col_sequences(self) -> Dict[str, List[Tuple[str, str, str, int]]]:
        """
        Extrait toutes les séquences avec backbone 'CoL'
        
        Returns:
            Dict avec barcode comme clé et liste de tuples (seq_id, sequence, original_id, tax_id)
        """
        query = """
        SELECT DISTINCT 
            s.seq_id,
            s.sequence,
            s.original_id,
            s.barcode,
            st.tax_id
        FROM sequences s
        INNER JOIN seq_tax st ON s.seq_id = st.seq_id
        WHERE st.backbone = 'CoL'
        AND s.sequence IS NOT NULL
        AND s.sequence != ''
        ORDER BY s.barcode, s.seq_id
        """
        
        sequences_by_barcode = defaultdict(list)
        
        try:
            with sqlite3.connect(self.db_path) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.execute(query)
                
                for row in cursor:
                    barcode = row['barcode'] or 'ITS'  # Valeur par défaut si NULL
                    sequences_by_barcode[barcode].append((
                        row['seq_id'],
                        row['sequence'],
                        row['original_id'],
                        row['tax_id']
                    ))
                
                logger.info(f"Séquences extraites par barcode:")
                for barcode, sequences in sequences_by_barcode.items():
                    logger.info(f"  {barcode}: {len(sequences)} séquences")
                
                return dict(sequences_by_barcode)
                
        except Exception as e:
            logger.error(f"Erreur lors de l'extraction: {e}")
            raise
    
    def write_fasta_files(self, sequences_by_barcode: Dict[str, List[Tuple]], 
                         output_dir: str = "output") -> List[str]:
        """
        Écrit les fichiers FASTA pour chaque barcode
        
        Args:
            sequences_by_barcode: Dict des séquences groupées par barcode
            output_dir: Répertoire de sortie
            
        Returns:
            Liste des fichiers créés
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        created_files = []
        
        for barcode, sequences in sequences_by_barcode.items():
            # Modification du nom de fichier : juste [barcode].fasta
            filename = f"{barcode}.fasta"
            filepath = output_path / filename
            
            try:
                with open(filepath, 'w') as f:
                    for seq_id, sequence, original_id, tax_id in sequences:
                        # Header avec seq_id et tax_id comme demandé
                        header = f">{seq_id}|tax_id:{tax_id}|original_id:{original_id}"
                        f.write(f"{header}\n")
                        
                        # Séquence formatée
                        f.write(f"{sequence}\n")
                
                logger.info(f"Fichier créé: {filepath} ({len(sequences)} séquences)")
                created_files.append(str(filepath))
                
            except Exception as e:
                logger.error(f"Erreur lors de l'écriture de {filename}: {e}")
                raise
        
        return created_files
    
    def get_taxonomy_info(self, tax_ids: List[int]) -> Dict[int, Dict]:
        """
        Récupère les informations taxonomiques pour les tax_ids donnés
        
        Args:
            tax_ids: Liste des identifiants taxonomiques
            
        Returns:
            Dict avec tax_id comme clé et info taxonomique comme valeur
        """
        if not tax_ids:
            return {}
        
        placeholders = ','.join(['?' for _ in tax_ids])
        query = f"""
        SELECT tax_id, kingdom, phylum, class_, order_, family, genus, species
        FROM taxonomy
        WHERE tax_id IN ({placeholders})
        """
        
        taxonomy_info = {}
        
        try:
            with sqlite3.connect(self.db_path) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.execute(query, tax_ids)
                
                for row in cursor:
                    taxonomy_info[row['tax_id']] = {
                        'kingdom': row['kingdom'],
                        'phylum': row['phylum'],
                        'class': row['class_'],
                        'order': row['order_'],
                        'family': row['family'],
                        'genus': row['genus'],
                        'species': row['species']
                    }
                    
            return taxonomy_info
            
        except Exception as e:
            logger.error(f"Erreur lors de la récupération des taxonomies: {e}")
            return {}
    
    def generate_stats(self, sequences_by_barcode: Dict[str, List[Tuple]]) -> Dict:
        """
        Génère les statistiques d'extraction
        
        Args:
            sequences_by_barcode: Dict des séquences groupées par barcode
            
        Returns:
            Dict contenant les statistiques
        """
        # Collecte des statistiques
        total_sequences = sum(len(sequences) for sequences in sequences_by_barcode.values())
        all_tax_ids = []
        barcode_stats = {}
        
        for barcode, sequences in sequences_by_barcode.items():
            # Statistiques sur les longueurs de séquences
            lengths = [len(seq) for _, seq, _, _ in sequences]
            barcode_tax_ids = set(tax_id for _, _, _, tax_id in sequences)
            all_tax_ids.extend([tax_id for _, _, _, tax_id in sequences])
            
            barcode_stats[barcode] = {
                'sequence_count': len(sequences),
                'unique_taxonomies': len(barcode_tax_ids),
                'length_stats': {
                    'min': min(lengths) if lengths else 0,
                    'max': max(lengths) if lengths else 0,
                    'mean': sum(lengths) / len(lengths) if lengths else 0
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
        # Configuration Snakemake
        if 'snakemake' in globals():
            input_db = snakemake.params.db
            output_dir = str(Path(snakemake.params.output_fodler))
            output_report = snakemake.output.report
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            
            # Configuration du logging pour Snakemake
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
            
        else:
            # Mode test standalone
            logger.warning("Running in standalone test mode")
            input_db = "data/db/its.db"
            output_dir = "tmp/output"
            output_report = "tmp/output/itsx_report.json"
            config = {}
            
            # Créer le répertoire de test
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        logger.info("Début de l'extraction des séquences CoL...")
        
        # Initialisation de l'extracteur
        extractor = SequenceExtractor(input_db)
        
        # Extraction des séquences
        sequences_by_barcode = extractor.extract_col_sequences()
        
        if not sequences_by_barcode:
            logger.warning("Aucune séquence CoL trouvée dans la base de données")
            
            # Créer un rapport vide pour Snakemake
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
                'processing_time': time.time() - start_time,
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
            processing_rate = total_sequences / processing_time if processing_time > 0 else 0
            
            # Génération du rapport final
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
                    'database_path': str(input_db)
                },
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'status': 'success'
            }
            
            logger.info("Extraction terminée avec succès!")
            logger.info(f"Total des séquences extraites: {total_sequences}")
            logger.info(f"Fichiers FASTA générés: {len(created_files)}")
            logger.info(f"Taux de traitement: {processing_rate:.0f} séquences/seconde")
        
        # Sauvegarde du rapport
        if 'snakemake' in globals() or output_report:
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Rapport sauvegardé: {output_report}")
        
        # Affichage du résumé pour le mode standalone
        if 'snakemake' not in globals():
            print(f"\n=== RÉSUMÉ ===")
            print(f"Total des séquences extraites: {report_data['stats']['total_sequences']}")
            print(f"Fichiers FASTA générés: {len(report_data['files_created'])}")
            print(f"Répertoire de sortie: {output_dir}")
            print(f"Temps de traitement: {processing_time:.2f}s")
        
    except Exception as e:
        logger.critical(f"Erreur fatale lors de l'extraction CoL: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()