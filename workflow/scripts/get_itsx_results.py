#!/usr/bin/env python3
"""
Script pour traiter les sorties de ITSx, valider les séquences, 
mettre à jour la base de données et générer un rapport JSON.
"""

import sqlite3
import logging
import argparse
import time
import sys
import json
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Iterator, Optional
from datetime import datetime

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)
logger = logging.getLogger(__name__)

# --- Fonctions utilitaires ---

def parse_fasta(file_path: Path) -> Iterator[Tuple[str, str]]:
    """Parse un fichier FASTA de manière efficace."""
    if not file_path.exists():
        logger.warning(f"Fichier FASTA non trouvé, ignoré : {file_path}")
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

def parse_header(header: str) -> Optional[int]:
    """Extrait le seq_id du header."""
    try:
        return int(header.split('|', 1)[0])
    except (ValueError, IndexError):
        logger.error(f"Impossible de parser le seq_id du header : {header}")
        return None

# --- Classe principale de traitement ---

class ITSxProcessor:
    """Gère le traitement des résultats ITSx et la mise à jour de la BDD."""
    
    def __init__(self, db_path: str, itsx_main_dir: str, config: Dict):
        self.db_path = Path(db_path)
        self.itsx_main_dir = Path(itsx_main_dir)
        self.config = config
        self.conn = None
        self.report = {
            "run_summary": {},
            "parameters": {
                "database": str(self.db_path),
                "itsx_directory": str(self.itsx_main_dir),
                **config
            },
            "barcodes": {}
        }

    def __enter__(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            self.conn.close()

    def _get_parent_metadata(self, parent_ids: List[int]) -> Dict[int, Dict]:
        """Récupère les métadonnées des séquences parentes."""
        if not parent_ids:
            return {}
        
        placeholders = ','.join('?' for _ in parent_ids)
        query = f"""
            SELECT s.seq_id, s.source_id, s.original_id, st.tax_id
            FROM sequences s
            JOIN seq_tax st ON s.seq_id = st.seq_id
            WHERE s.seq_id IN ({placeholders}) AND st.backbone = 'CoL'
        """
        cursor = self.conn.execute(query, list(parent_ids))
        return {row['seq_id']: dict(row) for row in cursor.fetchall()}

    def _validate_sequence(self, sequence: str) -> bool:
        """Valide une séquence selon la configuration."""
        if len(sequence) < self.config.get('min_length', 1):
            return False
        n_fraction = sequence.upper().count('N') / len(sequence)
        if n_fraction > self.config.get('max_n_fraction', 1.0):
            return False
        return True

    def run(self):
        """Lance le traitement pour tous les barcodes trouvés."""
        start_time = time.time()
        self.report["run_summary"]["start_time"] = datetime.now().isoformat()
        
        barcode_dirs = [d for d in self.itsx_main_dir.iterdir() if d.is_dir()]
        logger.info(f"Début du traitement pour {len(barcode_dirs)} barcodes...")
        
        for barcode_dir in barcode_dirs:
            barcode = barcode_dir.name
            self.process_barcode(barcode, barcode_dir)
        
        end_time = time.time()
        self.report["run_summary"]["end_time"] = datetime.now().isoformat()
        self.report["run_summary"]["total_duration_seconds"] = round(end_time - start_time, 2)
        self.report["run_summary"]["barcodes_processed"] = len(barcode_dirs)
        
        # Calcul des totaux
        total_inserted = sum(
            stats.get("total_new_sequences", 0)
            for stats in self.report["barcodes"].values()
        )
        self.report["run_summary"]["total_new_sequences_inserted"] = total_inserted
        self.report["run_summary"]["status"] = "Success"
        logger.info("Traitement de tous les barcodes terminé.")

    def process_barcode(self, barcode: str, barcode_dir: Path):
        """Exécute les 4 étapes de traitement pour un seul barcode."""
        logger.info(f"--- Traitement du barcode : {barcode} ---")
        stats = {
            "status": "Success",
            "initial_sequences_processed": 0,
            "subregions": {},
            "reconstructed_its": {"common_parents": 0, "inserted": 0},
            "no_detection": 0,
            "total_new_sequences": 0
        }
        cursor = self.conn.cursor()

        try:
            # Étape 1: Initialisation
            original_fasta_path = self.itsx_main_dir.parent / f"{barcode}.fasta"
            if original_fasta_path.exists():
                parent_ids = [pid for h, s in parse_fasta(original_fasta_path) if (pid := parse_header(h))]
                stats["initial_sequences_processed"] = len(parent_ids)
                if parent_ids:
                    cursor.executemany(
                        "UPDATE sequences SET barcode_status = ? WHERE seq_id = ?",
                        [("Processed", pid) for pid in parent_ids]
                    )

            # Étape 2: Sous-régions
            new_sequences_to_insert = []
            for region in ["ITS1", "ITS2"]:
                region_fasta = barcode_dir / f"{barcode}.{region}.fasta"
                sequences = list(parse_fasta(region_fasta))
                inserted_count = 0
                
                if sequences:
                    parent_ids = [pid for h, s in sequences if (pid := parse_header(h))]
                    parent_metadata = self._get_parent_metadata(parent_ids)
                    
                    for header, sequence in sequences:
                        if self._validate_sequence(sequence):
                            parent_id = parse_header(header)
                            if parent_id in parent_metadata:
                                meta = parent_metadata[parent_id]
                                new_sequences_to_insert.append({
                                    "sequence": sequence, "barcode": region, "parent_seq_id": parent_id,
                                    "original_id": meta["original_id"], "tax_id": meta["tax_id"],
                                    "barcode_status": "ITSx", "source_id": meta['source_id']
                                })
                                inserted_count += 1
                stats["subregions"][region] = {"found": len(sequences), "inserted": inserted_count}

            # Étape 3: Reconstruction
            its1_seqs = {parse_header(h): s for h, s in parse_fasta(barcode_dir / f"{barcode}.ITS1.fasta")}
            its58s_seqs = {parse_header(h): s for h, s in parse_fasta(barcode_dir / f"{barcode}.5_8S.fasta")}
            its2_seqs = {parse_header(h): s for h, s in parse_fasta(barcode_dir / f"{barcode}.ITS2.fasta")}
            common_ids = set(its1_seqs.keys()) & set(its58s_seqs.keys()) & set(its2_seqs.keys())
            
            stats["reconstructed_its"]["common_parents"] = len(common_ids)
            inserted_reconstructed = 0
            if common_ids:
                parent_metadata = self._get_parent_metadata(list(common_ids))
                for parent_id in common_ids:
                    full_seq = its1_seqs[parent_id] + its58s_seqs[parent_id] + its2_seqs[parent_id]
                    if self._validate_sequence(full_seq) and parent_id in parent_metadata:
                        meta = parent_metadata[parent_id]
                        new_sequences_to_insert.append({
                            "sequence": full_seq, "barcode": "ITS", "parent_seq_id": parent_id,
                            "original_id": meta["original_id"], "tax_id": meta["tax_id"],
                            "barcode_status": "ITSx", "source_id": meta['source_id']
                        })
                        inserted_reconstructed += 1
            stats["reconstructed_its"]["inserted"] = inserted_reconstructed

            # Insertion BDD
            if new_sequences_to_insert:
                for new_seq in new_sequences_to_insert:
                    cursor.execute(
                        """INSERT INTO sequences (source_id, original_id, sequence, length, barcode, barcode_status, parent_seq_id)
                           VALUES (?, ?, ?, ?, ?, ?, ?)""",
                        (new_seq['source_id'], new_seq["original_id"], new_seq["sequence"],
                         len(new_seq["sequence"]), new_seq["barcode"], new_seq["barcode_status"],
                         new_seq["parent_seq_id"])
                    )
                    new_seq_id = cursor.lastrowid
                    cursor.execute(
                        "INSERT INTO seq_tax (seq_id, tax_id, backbone) VALUES (?, ?, 'CoL')",
                        (new_seq_id, new_seq["tax_id"])
                    )
            stats["total_new_sequences"] = len(new_sequences_to_insert)

            # Étape 4: Non-détections
            no_detection_fasta = barcode_dir / f"{barcode}_no_detections.fasta"
            ids_no_detection = [pid for h, s in parse_fasta(no_detection_fasta) if (pid := parse_header(h))]
            stats["no_detection"] = len(ids_no_detection)
            if ids_no_detection:
                cursor.executemany(
                    "UPDATE sequences SET barcode_status = ? WHERE seq_id = ?",
                    [("no_detection", pid) for pid in ids_no_detection]
                )
            
            self.conn.commit()
            logger.info(f"Traitement du barcode '{barcode}' terminé et validé.")

        except Exception as e:
            self.conn.rollback()
            stats["status"] = "Error"
            stats["error_message"] = str(e)
            logger.error(f"Transaction annulée pour '{barcode}' : {e}", exc_info=True)
        
        self.report["barcodes"][barcode] = stats

    def write_report(self, output_path: str):
        """Écrit le dictionnaire du rapport dans un fichier JSON."""
        logger.info(f"Rédaction du rapport JSON : {output_path}")
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with output_file.open('w') as f:
            json.dump(self.report, f, indent=2)

def main():
    parser = argparse.ArgumentParser(description="Traite les résultats de ITSx et génère un rapport.")
    parser.add_argument("--database", required=True, help="Chemin vers la base de données SQLite.")
    parser.add_argument("--itsx-dir", required=True, help="Répertoire principal des sorties ITSx.")
    parser.add_argument("--output-report", required=True, help="Chemin du fichier de rapport JSON en sortie.")
    parser.add_argument("--min-length", type=int, default=50)
    parser.add_argument("--max-n-fraction", type=float, default=0.05)
    parser.add_argument("--log-file", help="Fichier de log où enregistrer la sortie.")
    
    args = parser.parse_args()

    if args.log_file:
        file_handler = logging.FileHandler(args.log_file, mode='w')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

    config = {
        "min_length": args.min_length,
        "max_n_fraction": args.max_n_fraction
    }
    
    try:
        with ITSxProcessor(args.database, args.itsx_dir, config) as processor:
            processor.run()
            processor.write_report(args.output_report)
    except Exception as e:
        logger.critical(f"Une erreur critique est survenue : {e}", exc_info=True)
        report = {"run_summary": {"status": "Critical Error", "error_message": str(e)}}
        Path(args.output_report).write_text(json.dumps(report, indent=2))
        sys.exit(1)

if __name__ == "__main__":
    main()