#!/usr/bin/env python3
"""
BOLD Parser - Adaptation du parser UNITE pour BOLD
"""

import os
import sys
import gzip
import re
import time
import logging
import json
import sqlite3
import tarfile
import csv
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from queue import Queue, Empty
import multiprocessing as mp
from dataclasses import dataclass

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constantes
BATCH_SIZE = 10000       # Plus gros batches pour moins d'I/O DB
QUEUE_SIZE = 20         # Plus de buffer pour éviter les blocages
MIN_SEQUENCE_LENGTH = 50  # Filtrage précoce
TSV_CHUNK_SIZE = 10000   # Taille des chunks pour le TSV

# Markers ITS acceptés
ACCEPTED_MARKERS = {'ITS', 'ITS1', 'ITS2'}

@dataclass
class SequenceRecord:
    """Structure de données pour les séquences BOLD."""
    original_id: str
    sequence: str
    length: int
    marker_code: str
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

class BOLDParser:
    """Parser BOLD pour données TSV."""
    
    def __init__(self, db_manager, config_validation: Dict):
        self.db_manager = db_manager
        self.config_validation = config_validation
        self.batch_size = BATCH_SIZE
        
        # Statistiques thread-safe
        self.stats = mp.Manager().dict({
            'total_sequences': 0,
            'processed_sequences': 0,
            'filtered_sequences': 0,
            'marker_filtered': 0,
            'parsing_errors': 0,
            'processing_time': 0.0
        })
        
        # Précompilation des caractères valides
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)
        
        # Mapping des colonnes BOLD (basé sur votre exemple)
        self.column_mapping = {
            'sampleid': 'sampleid', 
            'record_id': 'record_id',
            'kingdom': 'kingdom',
            'phylum': 'phylum',
            'class': 'class_',
            'order': 'order_',
            'family': 'family',
            'genus': 'genus',
            'species': 'species',
            'marker_code': 'marker_code',
            'nuc': 'sequence'
        }

    def validate_sequence(self, sequence: str) -> Tuple[bool, int, float]:
        """Validation de séquence avec early exit."""
        if not sequence:
            return False, 0, 1.0
        
        sequence = sequence.upper().strip().replace('-', '').replace(' ', '').replace('.', '')

        seq_len = len(sequence)
        if seq_len < self.min_length:
            return False, seq_len, 0.0
        
        # Comptage des N
        n_count = sequence.count('N') + sequence.count('n')
        n_fraction = n_count / seq_len
        
        if n_fraction > self.max_n_fraction:
            return False, seq_len, n_fraction
        
        # Validation des caractères (plus rapide avec set)
        if not all(c.upper() in self.valid_chars for c in sequence):
            return False, seq_len, n_fraction
        
        return True, seq_len, n_fraction


    def extract_tsv_from_tar(self, tar_path: str) -> str:
        """Extrait le fichier TSV du tar.gz et retourne son chemin."""
        logger.info(f"Extracting TSV from {tar_path}")
        
        with tarfile.open(tar_path, 'r:gz') as tar:
            # Chercher les fichiers TSV
            tsv_files = [member for member in tar.getmembers() 
                        if member.isfile() and member.name.lower().endswith('.tsv')]
            
            if not tsv_files:
                raise ValueError(f"No TSV file found in {tar_path}")
            
            if len(tsv_files) > 1:
                logger.warning(f"Multiple TSV files found, using first one: {tsv_files[0].name}")
            
            # Extraire dans un répertoire temporaire
            tsv_member = tsv_files[0]
            temp_dir = Path(tar_path).parent / "temp_bold_extract"
            temp_dir.mkdir(exist_ok=True)
            
            tar.extract(tsv_member, temp_dir)
            extracted_path = temp_dir / tsv_member.name
            
            logger.info(f"Extracted TSV to: {extracted_path}")
            return str(extracted_path)

    def tsv_reader(self, file_path: str, queue: Queue):
            """Lecteur TSV robuste utilisant pandas avec chunking."""
            logger.info(f"Starting pandas TSV reader for {file_path}")
            line_count = 0

            try:
                # Création d'un itérateur de DataFrames (chunks)
                reader_iterator = pd.read_csv(
                    file_path,
                    sep='\t',
                    engine='python',        # Moteur plus lent mais plus tolérant
                    quoting=csv.QUOTE_NONE, # Gère correctement les guillemets non standard
                    on_bad_lines='warn',    # Avertit des lignes mal formées au lieu de planter
                    chunksize=TSV_CHUNK_SIZE, # La clé pour ne pas tout charger en mémoire
                    dtype=str
                )

                for chunk_df in reader_iterator:
                    line_count += len(chunk_df)
                    
                    # --- 1. Filtrage précoce par marker (optimisé avec pandas) ---
                    initial_count = len(chunk_df)
                    
                    chunk_df = chunk_df.astype(str)

                    # Assurer que la colonne est de type string et normaliser
                    chunk_df['marker_code'] = chunk_df['marker_code'].str.strip().str.upper()
                    
                    # Masque booléen pour le filtrage
                    mask_marker = chunk_df['marker_code'].isin(ACCEPTED_MARKERS)
                    
                    # Appliquer le filtre
                    filtered_chunk = chunk_df[mask_marker]
                    
                    # Mettre à jour les statistiques
                    self.stats['marker_filtered'] += (initial_count - len(filtered_chunk))

                    if filtered_chunk.empty:
                        continue

                    # --- 2. Validation de séquence (plus simple de le faire en itérant sur le chunk filtré) ---
                    valid_rows_for_queue = []
                    
                    # itertuples est beaucoup plus rapide que iterrows
                    for row in filtered_chunk.itertuples(index=False):
                        sequence = getattr(row, 'nuc', '')
                        
                        # Nettoyage initial de la séquence
                        if pd.isna(sequence):
                            self.stats['filtered_sequences'] += 1
                            continue
                        
                        sequence = str(sequence).strip().replace('-', '')

                        # Validation via votre fonction existante
                        is_valid, seq_len, _ = self.validate_sequence(sequence)
                        
                        if is_valid:
                            # On convertit la NamedTuple de itertuples en dictionnaire
                            # pour garder la compatibilité avec le reste du code
                            valid_rows_for_queue.append((row._asdict(), seq_len))
                            self.stats['total_sequences'] += 1
                        else:
                            self.stats['filtered_sequences'] += 1

                    # --- 3. Envoyer le chunk traité dans la file d'attente ---
                    if valid_rows_for_queue:
                        queue.put(valid_rows_for_queue)

                    if line_count % (TSV_CHUNK_SIZE * 5) == 0:
                        logger.info(f"Processed {line_count} lines, "
                                f"accepted {self.stats['total_sequences']} sequences so far")

            except Exception as e:
                logger.error(f"Error in pandas TSV reader: {e}", exc_info=True)
                raise
            finally:
                queue.put(None)  # Signal de fin pour les consommateurs
                logger.info(f"Pandas TSV reader finished. Total lines processed: {line_count}")

    def _clean_sequence(self, sequence: str) -> str:
        """Nettoyage complet de la séquence - supprime tous les caractères non valides."""
        if not sequence:
            return ""
        
        # Conversion en majuscules et suppression des caractères spéciaux
        sequence = sequence.upper().strip()
        
        # Suppression des tirets, espaces, etc.
        sequence = sequence.replace('-', '').replace(' ', '').replace('.', '')
        
        # Filtrage des caractères valides seulement
        valid_chars = self.valid_chars
        cleaned_sequence = ''.join(c for c in sequence if c in valid_chars)
        
        return cleaned_sequence

    def sequence_processor_worker(self, chunk: List[Tuple[Dict, int]]) -> List[SequenceRecord]:
        """Worker pour traitement des enregistrements BOLD."""
        results = []
        
        for row_data, seq_len in chunk:

            # logger.info(f"Processing record: {row_data}")
            try:
                # Extraction des champs taxonomiques
                kingdom = self._clean_taxonomic_field(row_data.get('kingdom'))
                phylum = self._clean_taxonomic_field(row_data.get('phylum'))
                class_ = self._clean_taxonomic_field(row_data.get('class'))
                order_ = self._clean_taxonomic_field(row_data.get('order'))
                family = self._clean_taxonomic_field(row_data.get('family'))
                genus = self._clean_taxonomic_field(row_data.get('genus'))
                species = self._clean_taxonomic_field(row_data.get('species'))
                
                # Nettoyage de la séquence
                raw_sequence = row_data.get('nuc', '').strip().upper().replace('-', '')
                sequence = self._clean_sequence(raw_sequence)
                
                # Création de l'enregistrement
                record = SequenceRecord(
                    original_id=row_data.get('sampleid', '').strip(),
                    sequence=sequence,
                    length=seq_len,
                    marker_code=row_data.get('marker_code', '').strip(),
                    kingdom=kingdom,
                    phylum=phylum,
                    class_=class_,
                    order_=order_,
                    family=family,
                    genus=genus,
                    species=species,
                )
                
                results.append(record)
                
            except Exception as e:
                logger.warning(f"Error processing record: {e}")
                self.stats['parsing_errors'] += 1
                continue
        
        return results

    def _clean_taxonomic_field(self, value: str) -> Optional[str]:
        """Nettoyage des champs taxonomiques."""
        if not value or value.strip().lower() in ['none', '', 'null', 'na']:
            return None
        
        cleaned = value.strip().replace('_', ' ')
        return cleaned if cleaned else None

    def _clean_field(self, value: str) -> Optional[str]:
        """Nettoyage général des champs."""
        val = str(value).strip()
        if not val or val.lower() in ['none', '', 'null', 'na']:
            return None
        return val

    def db_writer_consumer(self, source_id: int, results_queue: Queue):
        """Writer avec transactions optimisées pour BOLD."""
        logger.info("Starting BOLD DB writer")
        batch = []
        
        with self.db_manager.get_connection() as conn:
            conn.execute("BEGIN IMMEDIATE TRANSACTION")
            
            try:
                while True:
                    try:
                        result_chunk = results_queue.get(timeout=6000)
                        if result_chunk is None:
                            break
                        
                        batch.extend(result_chunk)
                        
                        if len(batch) >= self.batch_size:
                            self._flush_bold_batch(conn, batch, source_id)
                            self.stats['processed_sequences'] += len(batch)
                            logger.info(f"Processed {self.stats['processed_sequences']} BOLD sequences")
                            batch.clear()
                            
                    except Empty:
                        logger.warning("DB writer timeout")
                        break
                
                # Traiter le dernier batch
                if batch:
                    self._flush_bold_batch(conn, batch, source_id)
                    self.stats['processed_sequences'] += len(batch)
                
                conn.execute("COMMIT")
                logger.info("BOLD DB writer finished successfully")
                
            except Exception as e:
                conn.execute("ROLLBACK")
                logger.error(f"BOLD DB writer error: {e}", exc_info=True)
                raise

    def _flush_bold_batch(self, conn: sqlite3.Connection, batch: List[SequenceRecord], source_id: int):
        """Flush optimisé pour les données BOLD."""
        cursor = conn.cursor()
        
        # 1. Préparer les données taxonomiques
        taxonomy_data = []
        sequence_data = []
        
        for record in batch:
            # Taxonomie étendue pour BOLD
            tax_tuple = (
                record.kingdom, record.phylum, record.class_, record.order_,
                record.family, record.genus, record.species
            )
            taxonomy_data.append(tax_tuple)
            
            # Données de séquence avec informations BOLD spécifiques
            sequence_data.append((
                source_id, record.original_id, record.sequence, 
                record.length, record.marker_code
            ))
        
        # 2. Insertion des taxonomies avec structure étendue
        cursor.executemany("""
            INSERT OR IGNORE INTO taxonomy 
            (kingdom, phylum, class_, order_, family, genus, species) 
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, taxonomy_data)
        
        # 3. Insertion des séquences avec métadonnées BOLD
        cursor.executemany("""
            INSERT INTO sequences 
            (source_id, original_id, sequence, length, barcode) 
            VALUES (?, ?, ?, ?, ?)
        """, sequence_data)
        
        # 4. Création des relations taxonomiques
        unique_taxonomies = list(set(taxonomy_data))
        tax_id_map = {}
        
        for tax_tuple in unique_taxonomies:
            columns = ["kingdom", "phylum", "class_", "order_", "family", "genus", "species"]
            
            query_parts = []
            values = []
            for col, val in zip(columns, tax_tuple):
                if val is None:
                    query_parts.append(f"{col} IS NULL")
                else:
                    query_parts.append(f"{col}=?")
                    values.append(val)
            
            sql = f"SELECT tax_id FROM taxonomy WHERE {' AND '.join(query_parts)}"
            cursor.execute(sql, values)
            result = cursor.fetchone()
            if result:
                tax_id_map[tax_tuple] = result[0]
        
        # 5. Relations seq_tax
        sequences_for_mapping = [record.sequence for record in batch]

        cursor.execute(
            f"""SELECT sequence, seq_id FROM sequences 
                WHERE source_id = ? AND sequence IN ({','.join(['?'] * len(batch))})""",
            (source_id, *sequences_for_mapping)
        )
        seq_id_map = {row[0]: row[1] for row in cursor.fetchall()}
        
        relation_data = []
        for record in batch:
            seq_id = seq_id_map.get(record.sequence)
            tax_tuple = (
                record.kingdom, record.phylum, record.class_, record.order_,
                record.family, record.genus, record.species
            )
            tax_id = tax_id_map.get(tax_tuple)
            
            if seq_id and tax_id:
                relation_data.append((seq_id, tax_id, 'Raw', 'Pending'))
        
        cursor.executemany(
            "INSERT OR IGNORE INTO seq_tax (seq_id, tax_id, backbone, status) VALUES (?, ?, ?, ?)",
            relation_data
        )

    def process_bold_database(self, input_file: str, source_name: str = "BOLD") -> Dict:
        """Pipeline de traitement BOLD."""
        start_time = time.time()
        logger.info(f"Starting BOLD processing: {input_file}")
        
        # Configuration DB pour performance maximale
        self.db_manager.optimize_for_bulk_insert_advanced()
        source_id = self.db_manager.add_source(name=source_name, file_path=input_file)
        
        # Extraction du TSV
        temp_tsv_path = None
        try:
            temp_tsv_path = self.extract_tsv_from_tar(input_file)
            
            # Files d'attente
            tsv_queue = Queue(maxsize=QUEUE_SIZE)
            results_queue = Queue(maxsize=QUEUE_SIZE)
            
            # Calcul optimal du nombre de workers
            cpu_count = os.cpu_count()
            num_workers = max(2, cpu_count - 2)
            logger.info(f"Using {num_workers} worker threads on {cpu_count} CPUs")
            
            with ThreadPoolExecutor(max_workers=num_workers + 2) as executor:
                # Lancer reader et writer
                writer_future = executor.submit(
                    self.db_writer_consumer, source_id, results_queue
                )
                reader_future = executor.submit(
                    self.tsv_reader, temp_tsv_path, tsv_queue
                )
                
                # Traitement parallèle des chunks
                worker_futures = []
                processed_chunks = 0
                
                while True:
                    try:
                        chunk = tsv_queue.get(timeout=120)
                        if chunk is None:
                            break
                        
                        future = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(future)
                        processed_chunks += 1
                        
                        # Traitement en streaming
                        if len(worker_futures) >= num_workers * 2:
                            completed_future = next(as_completed(worker_futures))
                            worker_futures.remove(completed_future)
                            results_queue.put(completed_future.result())
                            
                    except Empty:
                        logger.warning("Reader timeout - assuming completion")
                        break
                
                # Traiter les futures restants
                for future in as_completed(worker_futures):
                    results_queue.put(future.result())
                
                # Signaler la fin et attendre
                results_queue.put(None)
                writer_future.result()
                reader_future.result()
                
                logger.info(f"Processed {processed_chunks} chunks total")
                
        except Exception as e:
            logger.error(f"Critical error in BOLD processing pipeline: {e}", exc_info=True)
            raise
        finally:
            # Nettoyage des fichiers temporaires
            if temp_tsv_path and os.path.exists(temp_tsv_path):
                try:
                    os.remove(temp_tsv_path)
                    temp_dir = Path(temp_tsv_path).parent
                    if temp_dir.name == "temp_bold_extract":
                        temp_dir.rmdir()
                except Exception as e:
                    logger.warning(f"Could not cleanup temp files: {e}")
            
            # Restaurer les paramètres normaux
            self.db_manager.restore_safe_settings()
        
        # Finalisation
        end_time = time.time()
        self.stats['processing_time'] = end_time - start_time
        
        # Statistiques finales
        total_seqs = self.stats['total_sequences']
        processed_seqs = self.stats['processed_sequences']
        filtered_seqs = self.stats['filtered_sequences']
        marker_filtered = self.stats['marker_filtered']
        processing_time = self.stats['processing_time']
        
        logger.info(f"BOLD processing completed in {processing_time:.2f}s")
        logger.info(f"Total accepted: {total_seqs}, Processed: {processed_seqs}")
        logger.info(f"Filtered by sequence: {filtered_seqs}, Filtered by marker: {marker_filtered}")
        logger.info(f"Processing rate: {processed_seqs/processing_time:.0f} seqs/sec")
        
        return dict(self.stats)

# Point d'entrée
def main():
    """Point d'entrée principal."""
    try:
        # Configuration Snakemake
        if 'snakemake' in globals():
            input_db = snakemake.input[0]
            input_file = snakemake.input[1]
            output_report = snakemake.output[0]
            config_validation = snakemake.params.sequence_validation
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            # Mode test
            logger.warning("Running in standalone test mode")
            input_file = "test_bold.tar.gz"
            input_db = "test_output.sqlite"
            log_file = "test_format_bold.log"
            config = {'database': {'path': input_db}}
            config_validation = {
                'min_length': 50, 
                'valid_chars': 'ATCGN', 
                'max_n_fraction': 0.1
            }

        # Import et initialisation
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager
        
        db_manager = get_database_manager(config)
        
        # Parser
        parser = BOLDParser(db_manager, config_validation)
        
        # Traitement
        stats = parser.process_bold_database(input_file=input_file)
        
        # Optimisation finale
        db_manager.vacuum_database()
        
        # Génération du rapport
        processing_time = stats.get('processing_time', 1)
        processed_seqs = stats.get('processed_sequences', 0)
        rate = processed_seqs / processing_time if processing_time > 0 else 0
        
        report_data = {
            'source': 'BOLD',
            'version': 'latest',
            'input_file': input_file,
            'accepted_markers': list(ACCEPTED_MARKERS),
            'stats': dict(stats),
            'processing_rate_per_second': rate,
            'performance_metrics': {
                'cpu_cores_used': os.cpu_count(),
                'memory_optimizations': True,
                'batch_size': BATCH_SIZE,
                'chunk_size': TSV_CHUNK_SIZE
            },
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        # Sauvegarde du rapport
        if 'snakemake' in globals():
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"BOLD processing report saved to: {output_report}")
        
        logger.info(f"Final BOLD processing rate: {rate:.0f} sequences/second")
        
    except Exception as e:
        logger.critical(f"Fatal error in BOLD processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()