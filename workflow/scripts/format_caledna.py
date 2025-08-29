#!/usr/bin/env python3
"""
CALeDNA Parser
Parser pour les fichiers CALeDNA (.zip contenant .fasta et .txt de taxonomie)
"""

import os
import sys
import re
import time
import mmap
import logging
import json
import sqlite3
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from queue import Queue, Empty
import multiprocessing as mp
from dataclasses import dataclass
import tempfile
import shutil

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constantes
FASTA_CHUNK_SIZE = 10000
BATCH_SIZE = 10000
QUEUE_SIZE = 20
MIN_SEQUENCE_LENGTH = 50

@dataclass
class SequenceRecord:
    """Structure de données pour les séquences."""
    original_id: str
    sequence: str
    length: int
    n_frac: float
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

class CALeDNAParser:
    """Parser CALeDNA."""
    
    def __init__(self, db_manager, config_validation: Dict):
        self.db_manager = db_manager
        self.config_validation = config_validation
        self.batch_size = BATCH_SIZE
        
        # Patterns regex compilés
        self.version_pattern = re.compile(r'(\d{1,2}\.\d{1,2}|\d{4})')
        
        # Précompilation des caractères valides
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)

    def _reset_for_new_file(self):
        """Réinitialise le parser pour un nouveau fichier."""
        # Statistiques thread-safe pour ce fichier
        self.stats = mp.Manager().dict({
            'total_sequences': 0,
            'processed_sequences': 0,
            'filtered_sequences': 0,
            'parsing_errors': 0,
            'taxonomy_mapping_errors': 0,
            'processing_time': 0.0
        })
        
        # Nouveau dictionnaire pour mapping taxonomie (thread-safe)
        self.taxonomy_map = {}
        
        # Vider le cache LRU
        self._parse_caledna_taxonomy.cache_clear()

    def extract_zip_files(self, zip_path: str) -> Tuple[str, str]:
        """
        Extrait les fichiers .fasta et .txt du .zip
        Retourne les chemins vers le fichier FASTA et le fichier taxonomie
        """
        logger.info(f"Extracting files from {zip_path}")
        
        # Créer un dossier temporaire
        temp_dir = tempfile.mkdtemp(prefix="caledna_extract_")
        
        fasta_file = None
        taxonomy_file = None
        
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                # Lister tous les fichiers dans l'archive
                file_list = zip_ref.namelist()
                logger.info(f"Files in archive: {file_list}")
                
                # Extraire tous les fichiers
                zip_ref.extractall(temp_dir)
                
                # Trouver les fichiers .fasta et .txt
                for root, dirs, files in os.walk(temp_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        if file.lower().endswith(('.fasta', '.fas', '.fa')):
                            fasta_file = file_path
                            logger.info(f"Found FASTA file: {fasta_file}")
                        elif file.lower().endswith('.txt') and 'tax' in file.lower():
                            taxonomy_file = file_path
                            logger.info(f"Found taxonomy file: {taxonomy_file}")
                        elif file.lower().endswith('.txt') and not taxonomy_file:
                            # Si aucun fichier avec 'tax' dans le nom, prendre le premier .txt
                            taxonomy_file = file_path
                            logger.info(f"Found potential taxonomy file: {taxonomy_file}")
                
                if not fasta_file:
                    raise FileNotFoundError("No FASTA file found in the archive")
                if not taxonomy_file:
                    raise FileNotFoundError("No taxonomy file found in the archive")
                
        except Exception as e:
            # Nettoyer en cas d'erreur
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
            raise e
        
        return fasta_file, taxonomy_file, temp_dir

    def load_taxonomy_mapping(self, taxonomy_file: str):
        """
        Charge le mapping taxonomique depuis le fichier .txt
        Format attendu: ID\tTaxonomie;séparée;par;points-virgules
        """
        logger.info(f"Loading taxonomy mapping from {taxonomy_file}")
        
        taxonomy_count = 0
        
        try:
            with open(taxonomy_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        # Séparer ID et taxonomie (par tabulation ou espace)
                        if '\t' in line:
                            seq_id, taxonomy_str = line.split('\t', 1)
                        elif ' ' in line:
                            parts = line.split(' ', 1)
                            if len(parts) == 2:
                                seq_id, taxonomy_str = parts
                            else:
                                continue
                        else:
                            continue
                        
                        seq_id = seq_id.strip()
                        taxonomy_str = taxonomy_str.strip()
                        
                        # Parser la taxonomie
                        taxonomy = self._parse_caledna_taxonomy(taxonomy_str)
                        self.taxonomy_map[seq_id] = taxonomy
                        taxonomy_count += 1
                        
                        if taxonomy_count % 10000 == 0:
                            logger.info(f"Loaded {taxonomy_count} taxonomy mappings")
                            
                    except Exception as e:
                        logger.warning(f"Error parsing taxonomy line {line_num}: {line[:50]}... - {e}")
                        self.stats['taxonomy_mapping_errors'] += 1
                        continue
        
        except Exception as e:
            logger.error(f"Error loading taxonomy file: {e}")
            raise
        
        logger.info(f"Successfully loaded {taxonomy_count} taxonomy mappings")

    @lru_cache(maxsize=50000)
    def _parse_caledna_taxonomy(self, taxonomy_str: str) -> Dict[str, Optional[str]]:
        """
        Parse la chaîne taxonomique CALeDNA
        Format: Eukaryota;Ascomycota;Sordariomycetes;Microascales;Microascaceae;Microascus;Microascus cirrosus
        """
        try:
            ranks = ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']
            taxonomy = {}
            
            if not taxonomy_str or taxonomy_str.strip() == '':
                return {rank: None for rank in ranks}
            
            # Diviser par point-virgule
            tax_parts = [part.strip() for part in taxonomy_str.split(';')]
            
            # Mapper aux rangs taxonomiques
            for i, rank in enumerate(ranks):
                if i < len(tax_parts):
                    value = tax_parts[i].strip()
                    # Nettoyer les valeurs vides ou "NA"
                    if value and value.upper() not in ['', 'NA', 'UNKNOWN', 'UNIDENTIFIED']:
                        taxonomy[rank] = value.replace('_', ' ')
                    else:
                        taxonomy[rank] = None
                else:
                    taxonomy[rank] = None
            
            return taxonomy
            
        except Exception as e:
            logger.warning(f"Error parsing taxonomy '{taxonomy_str}': {e}")
            return {rank: None for rank in ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']}

    def validate_sequence(self, sequence: str) -> Tuple[bool, int, float]:
        """Validation de séquence avec early exit."""
        if not sequence:
            return False, 0, 1.0
        
        seq_len = len(sequence)
        if seq_len < self.min_length:
            return False, seq_len, 0.0
        
        # Comptage des N
        n_count = sequence.count('N') + sequence.count('n')
        n_fraction = n_count / seq_len
        
        if n_fraction > self.max_n_fraction:
            return False, seq_len, n_fraction
        
        # Validation des caractères (plus rapide avec set)
        if not all(c.upper() in self.valid_chars for c in sequence[:100]):
            return False, seq_len, n_fraction
        
        return True, seq_len, n_fraction

    def fasta_reader_mmap(self, file_path: str, queue: Queue):
        """
        Lecteur FASTA.
        Lit le fichier séquence par séquence et utilise _add_to_chunk pour créer
        des "batchs" de tuples (header, sequence, len, n_frac) au format attendu.
        """
        logger.info(f"Starting FASTA reader for {file_path}")
        
        batch = []
        header, sequence_parts = None, []

        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        if header and sequence_parts:
                            sequence = ''.join(sequence_parts)
                            self._add_to_chunk(header, sequence, batch, queue)
                        
                        header = line
                        sequence_parts = []
                    elif header:
                        sequence_parts.append(line)

                if header and sequence_parts:
                    sequence = ''.join(sequence_parts)
                    self._add_to_chunk(header, sequence, batch, queue)

        except Exception as e:
            logger.error(f"Error in FASTA reader: {e}", exc_info=True)
        finally:
            if batch:
                queue.put(batch)
            queue.put(None)
            logger.info("FASTA reader finished and sent sentinel.")

    def _process_fasta_stream(self, stream, chunk, queue):
        """Traitement d'un stream FASTA."""
        header, sequence_parts = None, []
        
        for line in stream:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if header and sequence_parts:
                    sequence = ''.join(sequence_parts).replace('-', '')
                    self._add_to_chunk(header, sequence, chunk, queue)
                
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)
        
        if header and sequence_parts:
            sequence = ''.join(sequence_parts).replace('-', '')
            self._add_to_chunk(header, sequence, chunk, queue)

    def _add_to_chunk(self, header: str, sequence: str, chunk: List, queue: Queue):
        """Ajoute une séquence au chunk avec validation précoce."""
        is_valid, seq_len, n_frac = self.validate_sequence(sequence)
        
        if is_valid:
            chunk.append((header, sequence, seq_len, n_frac))
            self.stats['total_sequences'] += 1
            
            if len(chunk) >= FASTA_CHUNK_SIZE:
                queue.put(chunk[:])
                chunk.clear()
        else:
            self.stats['filtered_sequences'] += 1

    def sequence_processor_worker(self, chunk: List[Tuple[str, str, int, float]]) -> List[SequenceRecord]:
        """Worker qui traite un chunk de séquences."""
        results = []
        
        for header, sequence, seq_len, n_frac in chunk:
            seq_id = header.lstrip('>').split()[0]  # Premier token après >
            
            taxonomy = self.taxonomy_map.get(seq_id, {})
            
            if not taxonomy:
                logger.debug(f"No taxonomy found for sequence ID: {seq_id}")
                taxonomy = {rank: None for rank in ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']}
            
            # Création du SequenceRecord
            record = SequenceRecord(
                original_id=seq_id,
                sequence=sequence,
                length=seq_len,
                n_frac=n_frac,
                kingdom=taxonomy.get('kingdom'),
                phylum=taxonomy.get('phylum'),
                class_=taxonomy.get('class_'),
                order_=taxonomy.get('order_'),
                family=taxonomy.get('family'),
                genus=taxonomy.get('genus'),
                species=taxonomy.get('species')
            )
            
            results.append(record)
        
        return results

    def db_writer_consumer(self, source_id: int, results_queue: Queue):
        """Writer avec transactions plus efficaces."""
        logger.info("Starting DB writer")
        batch = []
        
        try:
            # Créer une nouvelle connexion pour ce thread
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
                                self._flush_batch(conn, batch, source_id)
                                self.stats['processed_sequences'] += len(batch)
                                logger.info(f"Processed {self.stats['processed_sequences']} sequences")
                                batch.clear()
                                
                        except Empty:
                            logger.warning("DB writer timeout")
                            break
                    
                    if batch:
                        self._flush_batch(conn, batch, source_id)
                        self.stats['processed_sequences'] += len(batch)
                    
                    conn.execute("COMMIT")
                    logger.info("DB writer finished successfully")
                    
                except Exception as e:
                    conn.execute("ROLLBACK")
                    logger.error(f"DB writer error: {e}", exc_info=True)
                    raise
                    
        except Exception as e:
            logger.error(f"Critical error in DB writer: {e}", exc_info=True)
            raise

    def _flush_batch(self, conn: sqlite3.Connection, batch: List[SequenceRecord], source_id: int):
        """Flush avec prepared statements."""
        cursor = conn.cursor()
        
        taxonomy_data = []
        sequence_data = []
        relation_data = []
        
        for record in batch:
            tax_tuple = (record.kingdom, record.phylum, record.class_, 
                        record.order_, record.family, record.genus, record.species)
            taxonomy_data.append(tax_tuple)
            sequence_data.append((source_id, record.original_id, record.sequence, record.length, record.n_frac))
        
        # Insertion des taxonomies
        cursor.executemany(
            "INSERT OR IGNORE INTO taxonomy (kingdom, phylum, class_, order_, family, genus, species) VALUES (?, ?, ?, ?, ?, ?, ?)",
            taxonomy_data
        )
        
        # Récupération des tax_ids
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

        # Insertion des séquences
        cursor.executemany(
            "INSERT INTO sequences (source_id, original_id, sequence, length, percentage_n) VALUES (?, ?, ?, ?, ?)",
            sequence_data
        )
        
        # Relations seq_tax
        cursor.execute(
            f"SELECT original_id, seq_id FROM sequences WHERE source_id = ? AND original_id IN ({','.join(['?'] * len(batch))})",
            (source_id, *[record.original_id for record in batch])
        )
        seq_id_map = {row[0]: row[1] for row in cursor.fetchall()}
        
        for record in batch:
            seq_id = seq_id_map.get(record.original_id)
            tax_tuple = (record.kingdom, record.phylum, record.class_, 
                        record.order_, record.family, record.genus, record.species)
            tax_id = tax_id_map.get(tax_tuple)
            
            if seq_id and tax_id:
                relation_data.append((seq_id, tax_id, 'Raw', 'Pending'))
        
        cursor.executemany(
            "INSERT INTO seq_tax (seq_id, tax_id, backbone, status) VALUES (?, ?, ?, ?)",
            relation_data
        )

    def process_caledna_database(self, zip_file: str, source_name: str = "CALeDNA") -> Dict:
        """Pipeline de traitement principal."""
        
        # Réinitialiser pour ce nouveau fichier
        self._reset_for_new_file()

        start_time = time.time()
        logger.info(f"Starting CALeDNA processing: {zip_file}")
        
        temp_dir = None
        
        try:
            # Extraction des fichiers
            fasta_file, taxonomy_file, temp_dir = self.extract_zip_files(zip_file)
            
            # Chargement du mapping taxonomique
            self.load_taxonomy_mapping(taxonomy_file)
            
            # Configuration DB pour performance maximale
            self.db_manager.optimize_for_bulk_insert_advanced()
            source_id = self.db_manager.add_source(name=source_name, file_path=zip_file)
            
            # Files d'attente
            fasta_queue = Queue(maxsize=QUEUE_SIZE)
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
                    self.fasta_reader_mmap, fasta_file, fasta_queue
                )
                
                # Traitement parallèle des chunks
                worker_futures = []
                processed_chunks = 0
                
                while True:
                    try:
                        chunk = fasta_queue.get(timeout=120)
                        if chunk is None:
                            break
                        
                        future = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(future)
                        processed_chunks += 1
                        
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
            logger.error(f"Critical error in processing pipeline: {e}", exc_info=True)
            raise
        finally:
            # Nettoyage
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temporary directory: {temp_dir}")
        
        # Finalisation
        end_time = time.time()
        self.stats['processing_time'] = end_time - start_time
        
        # Statistiques finales
        total_seqs = self.stats['total_sequences']
        processed_seqs = self.stats['processed_sequences']
        filtered_seqs = self.stats['filtered_sequences']
        processing_time = self.stats['processing_time']
        
        logger.info(f"Processing completed in {processing_time:.2f}s")
        logger.info(f"Total: {total_seqs}, Processed: {processed_seqs}, "
                   f"Filtered: {filtered_seqs}, Rate: {processed_seqs/processing_time:.0f} seqs/sec")
        
        return dict(self.stats)

# Point d'entrée
def main():
    """Point d'entrée principal."""
    try:
        # Configuration Snakemake
        if 'snakemake' in globals():
            input_db = snakemake.input[0]
            input_file = snakemake.input.file
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
            input_files = ["test_FITS.zip", "test_PITS.zip"] # Mettez vos fichiers de test ici
            input_db = "test_output.sqlite"
            log_file = "test_format_caledna.log"
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
        parser = CALeDNAParser(db_manager, config_validation)
        
        # Restaurer les paramètres sûrs avant chaque fichier
        db_manager.restore_safe_settings()
        
        
        stats = parser.process_caledna_database(
            zip_file=str(input_file), 
            source_name=f"CALeDNA_{Path(input_file).stem}"
        )
        
        # Optimisation finale
        db_manager.restore_safe_settings()
        db_manager.vacuum_database()
        
        # Génération du rapport
        total_time = stats.get('processing_time', 1)
        total_processed = stats.get('processed_sequences', 0)
        rate = total_processed / total_time if total_time > 0 else 0
        
        report_data = {
            'source': 'CALeDNA_Multi',
            'input_file': input_file,
            'stats': stats,
            'processing_rate_per_second': rate,
            'performance_metrics': {
                'cpu_cores_used': os.cpu_count(),
                'memory_optimizations': True,
                'batch_size': BATCH_SIZE,
                'chunk_size': FASTA_CHUNK_SIZE,
            },
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        # Sauvegarde du rapport
        if 'snakemake' in globals():
            with open(output_report, 'w') as f:
                json.dump(report_data, f, indent=2)
            logger.info(f"Optimization report saved to: {output_report}")
        
        logger.info(f"Final processing rate: {rate:.0f} sequences/second")
        
    except Exception as e:
        logger.critical(f"Fatal error in CALeDNA processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()