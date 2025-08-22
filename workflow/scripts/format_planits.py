#!/usr/bin/env python3
"""
PLANITS Parser
"""

import os
import sys
import zipfile
import re
import time
import logging
import json
import sqlite3
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from queue import Queue, Empty
import multiprocessing as mp
from dataclasses import dataclass
from io import TextIOWrapper

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constantes
FASTA_CHUNK_SIZE = 10000  # Augmenté pour réduire l'overhead
BATCH_SIZE = 10000       # Plus gros batches pour moins d'I/O DB
QUEUE_SIZE = 20         # Plus de buffer pour éviter les blocages
MIN_SEQUENCE_LENGTH = 50  # Filtrage précoce

@dataclass
class SequenceRecord:
    """Structure de données pour les séquences."""
    original_id: str
    sequence: str
    length: int
    barcode: str  # ITS, ITS1, ou ITS2
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

class PLANITSParser:
    """Parser PLANITS."""
    
    def __init__(self, db_manager, config_validation: Dict):
        self.db_manager = db_manager
        self.config_validation = config_validation
        self.batch_size = BATCH_SIZE
        
        # Statistiques thread-safe
        self.stats = mp.Manager().dict({
            'total_sequences': 0,
            'processed_sequences': 0,
            'filtered_sequences': 0,
            'parsing_errors': 0,
            'processing_time': 0.0,
            'its_ITS': 0, 
            'its_ITS1': 0,
            'its_ITS2': 0
        })
        
        # Patterns regex compilés
        self.version_pattern = re.compile(r'(\d{1,2}\.\d{1,2}\.\d{4})')
        
        # Précompilation des caractères valides
        self.valid_chars = set(self.config_validation.get('valid_chars', 'ATCGN'))
        self.min_length = self.config_validation.get('min_length', MIN_SEQUENCE_LENGTH)
        self.max_n_fraction = self.config_validation.get('max_n_fraction', 0.1)

    @lru_cache(maxsize=50000)  # Cache plus grand
    def _parse_planits_taxonomy(self, header: str, barcode_type: str) -> Tuple[str, Tuple]:
        """Parse des headers PLANITS format SINTAX.
        
        Format: >MK696133.1;tax=p:Chlorophyta,c:Chlorophyceae,o:Chlamydomonadales,f:NA,g:NA,s:unknown Chlamydomonadales ;

        Attention il n'y a pas de kingdom dans le header
        """
        try:
            header = header.lstrip('>')
            
            # Séparer l'ID de la taxonomie
            if ';tax=' in header:
                seq_id_part, tax_part = header.split(';tax=', 1)
                seq_id = seq_id_part.strip()
                seq_id = seq_id + f"_{barcode_type}"
                
                # Extraire la partie taxonomique (avant le dernier ';')
                tax_string = tax_part.split(';')[0]  # Prendre seulement la partie avant le dernier ';'
                
                # Parser la taxonomie
                ranks = {}
                for entry in tax_string.split(','):
                    if ':' in entry:
                        rank_code, value = entry.split(':', 1)
                        rank_code = rank_code.strip()
                        value = value.strip()
                        
                        # Mapping des codes de rang vers les noms complets
                        rank_mapping = {
                            'k': 'kingdom',
                            'p': 'phylum', 
                            'c': 'class_',
                            'o': 'order_',
                            'f': 'family',
                            'g': 'genus',
                            's': 'species'
                        }
                        
                        rank_name = rank_mapping.get(rank_code, rank_code)
                        
                        # Nettoyer la valeur
                        if value and value.lower() not in ['na', 'unknown', '']:
                            # Traiter les cas spéciaux pour species
                            if rank_name == 'species' and ' ' in value:
                                # Garder seulement le nom d'espèce (après le genre)
                                parts = value.split()
                                if len(parts) >= 2:
                                    value = ' '.join(parts[:2])  # Genre + espèce
                            
                            ranks[rank_name] = value
                        else:
                            ranks[rank_name] = None
                
                return seq_id, tuple(sorted(ranks.items()))
            
            else:
                # Pas de taxonomie trouvée
                return header.strip(), ()
                
        except Exception as e:
            logger.debug(f"Error parsing header {header[:100]}...: {e}")
            self.stats['parsing_errors'] += 1
            return header.lstrip('>').strip(), ()

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
        if not all(c.upper() in self.valid_chars for c in sequence[:100]):  # Check sur 100 premiers chars seulement
            return False, seq_len, n_fraction
        
        return True, seq_len, n_fraction

    def extract_barcode_type(self, filename: str) -> str:
        """Extrait le type de barcode depuis le nom de fichier."""
        filename_lower = filename.lower()
        if 'its1.sintax' in filename_lower:
            return 'ITS1'
        elif 'its2.sintax' in filename_lower:
            return 'ITS2'
        elif 'its.sintax' in filename_lower:
            return 'ITS'
        else:
            logger.warning(f"Unknown barcode type for file: {filename}")
            return 'ITS'  # Default

    def fasta_reader_zip(self, zip_path: str, queue: Queue):
        """Lecteur FASTA pour fichiers ZIP PLANITS."""
        logger.info(f"Starting ZIP reader for {zip_path}")
        
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_file:
                # Chercher les fichiers SINTAX_format.fas
                sintax_files = [f for f in zip_file.namelist() 
                              if f.endswith('.SINTAX_format.fas')]
                
                logger.info(f"Found SINTAX files: {sintax_files}")
                
                for sintax_file in sintax_files:
                    logger.info(f"Processing file: {sintax_file}")
                    barcode_type = self.extract_barcode_type(sintax_file)
                    
                    with zip_file.open(sintax_file) as file:
                        # Wrapper pour traiter comme fichier texte
                        text_file = TextIOWrapper(file, encoding='utf-8')
                        self._process_fasta_stream(text_file, barcode_type, queue)
                        
        except Exception as e:
            logger.error(f"Error in ZIP reader: {e}", exc_info=True)
        finally:
            queue.put(None)  # Signal de fin
            logger.info("ZIP reader finished")

    def _process_fasta_stream(self, stream, barcode_type: str, queue: Queue):
        """Traitement d'un stream FASTA."""
        header, sequence_parts = None, []
        chunk = []
        
        for line in stream:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if header and sequence_parts:
                    sequence = ''.join(sequence_parts).replace('-', '')  # Enlever les alignements
                    self._add_to_chunk(header, sequence, barcode_type, chunk, queue)
                
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)
        
        # Dernier enregistrement
        if header and sequence_parts:
            sequence = ''.join(sequence_parts).replace('-', '')
            self._add_to_chunk(header, sequence, barcode_type, chunk, queue)
        
        # Envoyer le chunk restant
        if chunk:
            queue.put(chunk[:])

    def _add_to_chunk(self, header: str, sequence: str, barcode_type: str, chunk: List, queue: Queue):
        """Ajoute une séquence au chunk avec validation précoce."""
        # Validation précoce pour éviter le traitement inutile
        is_valid, seq_len, n_frac = self.validate_sequence(sequence)
        
        if is_valid:
            chunk.append((header, sequence, seq_len, barcode_type))
            self.stats['total_sequences'] += 1
            self.stats[f'its_{barcode_type}'] += 1
            
            if len(chunk) >= FASTA_CHUNK_SIZE:
                queue.put(chunk[:])  # Copie du chunk
                chunk.clear()
        else:
            self.stats['filtered_sequences'] += 1

    def sequence_processor_worker(self, chunk: List[Tuple[str, str, int, str]]) -> List[SequenceRecord]:
        """Worker avec structures de données plus efficaces."""
        results = []
        
        for header, sequence, seq_len, barcode_type in chunk:
            seq_id, taxonomy_tuple = self._parse_planits_taxonomy(header, barcode_type)
            tax_dict = dict(taxonomy_tuple) if taxonomy_tuple else {}
            
            # Création directe du SequenceRecord
            record = SequenceRecord(
                original_id=seq_id,
                sequence=sequence,
                length=seq_len,
                barcode=barcode_type,
                kingdom="Plantae",  # Valeur forcée
                phylum=tax_dict.get('phylum'),
                class_=tax_dict.get('class_'),
                order_=tax_dict.get('order_'),
                genus=tax_dict.get('genus'),
                species=tax_dict.get('species')
            )
            
            results.append(record)
        
        return results

    def db_writer_consumer(self, source_id: int, results_queue: Queue):
        """Writer avec transactions plus efficaces."""
        logger.info("Starting DB writer")
        batch = []
        
        # Préparer les statements SQL
        with self.db_manager.get_connection() as conn:
            conn.execute("BEGIN IMMEDIATE TRANSACTION")
            
            try:
                while True:
                    try:
                        result_chunk = results_queue.get(timeout=60)
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
                
                # Traiter le dernier batch
                if batch:
                    self._flush_batch(conn, batch, source_id)
                    self.stats['processed_sequences'] += len(batch)
                
                conn.execute("COMMIT")
                logger.info("DB writer finished successfully")
                
            except Exception as e:
                conn.execute("ROLLBACK")
                logger.error(f"DB writer error: {e}", exc_info=True)
                raise

    def _flush_batch(self, conn: sqlite3.Connection, batch: List[SequenceRecord], source_id: int):
        """Flush avec prepared statements."""
        cursor = conn.cursor()
        
        # 1. Préparer toutes les données d'un coup
        taxonomy_data = []
        sequence_data = []
        relation_data = []

        for record in batch:
            # Taxonomie
            tax_tuple = (record.kingdom, record.phylum, record.class_, 
                        record.order_, record.family, record.genus, record.species)
            taxonomy_data.append(tax_tuple)
            
            # Séquence
            sequence_data.append((source_id, record.original_id, record.sequence, record.length))
        
        # 2. Insertion des taxonomies (plus efficace)
        cursor.executemany(
            "INSERT OR IGNORE INTO taxonomy (kingdom, phylum, class_, order_, family, genus, species) VALUES (?, ?, ?, ?, ?, ?, ?)",
            taxonomy_data
        )
        
        # 3. Récupération des tax_ids en batch
        unique_taxonomies = list(set(taxonomy_data))
        tax_id_map = {}
        
        for tax_tuple in unique_taxonomies:
            # On enlève kingdom (1er élément)
            query_tuple = tax_tuple[1:]  

            # Colonnes de la table correspondantes
            columns = ["phylum", "class_", "order_", "family", "genus", "species"]

            # Construction dynamique de la requête
            query_parts = []
            values = []
            for col, val in zip(columns, query_tuple):
                if val is None:
                    query_parts.append(f"{col} IS NULL")
                else:
                    query_parts.append(f"{col}=?")
                    values.append(val)

            sql = f"SELECT tax_id FROM taxonomy WHERE {' AND '.join(query_parts)}"
            cursor.execute(sql, values)

            result = cursor.fetchone()

            if result:
                tax_id_map[query_tuple] = result[0]
        
        
        # 4. Insertion des séquences
        cursor.executemany(
            "INSERT INTO sequences (source_id, original_id, sequence, length) VALUES (?, ?, ?, ?)",
            sequence_data
        )
        
        # 5. Récupération des seq_ids et création des relations
        cursor.execute(
            f"SELECT original_id, seq_id FROM sequences WHERE source_id = ? AND original_id IN ({','.join(['?'] * len(batch))})",
            (source_id, *[record.original_id for record in batch])
        )
        seq_id_map = {row[0]: row[1] for row in cursor.fetchall()}
        
        # 6. Relations seq_tax
        for record in batch:
            seq_id = seq_id_map.get(record.original_id)
            tax_tuple = (record.phylum, record.class_, record.order_, record.family, record.genus, record.species)
            tax_id = tax_id_map.get(tax_tuple)
            
            if seq_id and tax_id:
                relation_data.append((seq_id, tax_id, 'Raw', 'Pending'))
        

        cursor.executemany(
            "INSERT INTO seq_tax (seq_id, tax_id, backbone, status) VALUES (?, ?, ?, ?)",
            relation_data   
        )

        logger.info(f"Creating {len(relation_data)} seq_tax relations")

    def process_planits_database(self, input_file: str, source_name: str = "PLANITS", 
                                source_version: str = "unknown") -> Dict:
        """Pipeline de traitement PLANITS."""
        start_time = time.time()
        logger.info(f"Starting PLANITS processing: {input_file}")
        
        # Vérification que c'est un fichier ZIP
        if not input_file.endswith('.zip'):
            raise ValueError(f"Expected ZIP file, got: {input_file}")
        
        # Configuration DB pour performance maximale
        self.db_manager.optimize_for_bulk_insert_advanced()
        source_id = self.db_manager.add_source(name=source_name, file_path=input_file)
        
        # Files d'attente
        fasta_queue = Queue(maxsize=QUEUE_SIZE)
        results_queue = Queue(maxsize=QUEUE_SIZE)
        
        # Calcul optimal du nombre de workers
        cpu_count = os.cpu_count()
        num_workers = max(2, cpu_count - 2)  # Garder 2 CPUs pour reader/writer
        logger.info(f"Using {num_workers} worker threads on {cpu_count} CPUs")
        
        try:
            with ThreadPoolExecutor(max_workers=num_workers + 2) as executor:
                # Lancer reader et writer
                writer_future = executor.submit(
                    self.db_writer_consumer, source_id, results_queue
                )
                reader_future = executor.submit(
                    self.fasta_reader_zip, input_file, fasta_queue
                )
                
                # Traitement parallèle des chunks
                worker_futures = []
                processed_chunks = 0
                
                while True:
                    try:
                        chunk = fasta_queue.get(timeout=120)  # Timeout plus long
                        if chunk is None:
                            break
                        
                        # Soumettre le traitement
                        future = executor.submit(self.sequence_processor_worker, chunk)
                        worker_futures.append(future)
                        processed_chunks += 1
                        
                        # Traitement en streaming pour éviter la surcharge mémoire
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
            # Restaurer les paramètres normaux
            self.db_manager.restore_safe_settings()
        
        # Finalisation
        end_time = time.time()
        self.stats['processing_time'] = end_time - start_time
        
        # Statistiques finales
        total_seqs = self.stats['total_sequences']
        processed_seqs = self.stats['processed_sequences']
        filtered_seqs = self.stats['filtered_sequences']
        processing_time = self.stats['processing_time']
        its_ITS = self.stats['its_ITS']
        its_ITS1 = self.stats['its_ITS1']
        its_ITS2 = self.stats['its_ITS2']
        
        logger.info(f"Processing completed in {processing_time:.2f}s")
        logger.info(f"Total: {total_seqs}, Processed: {processed_seqs}, "
                   f"Filtered: {filtered_seqs}, Rate: {processed_seqs/processing_time:.0f} seqs/sec")
        logger.info(f"ITS1: {its_ITS1}, ITS2: {its_ITS2}, ITS: {its_ITS}")
        
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
        else:
            # Mode test
            logger.warning("Running in standalone test mode")
            input_file = "PLANITS_small.zip"
            input_db = "test_output.sqlite"
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
        db_manager.create_schema()
        
        # Parser
        parser = PLANITSParser(db_manager, config_validation)
        
        # Détection de version
        version_match = parser.version_pattern.search(Path(input_file).name)
        version = version_match.group(1) if version_match else "unknown"
        
        # Traitement
        stats = parser.process_planits_database(
            input_file=input_file, 
            source_version=version
        )
        
        # Optimisation finale
        db_manager.vacuum_database()
        
        # Génération du rapport
        processing_time = stats.get('processing_time', 1)
        processed_seqs = stats.get('processed_sequences', 0)
        rate = processed_seqs / processing_time if processing_time > 0 else 0
        
        report_data = {
            'source': 'PLANITS',
            'version': version,
            'input_file': input_file,
            'stats': dict(stats),
            'processing_rate_per_second': rate,
            'its_ITS': stats.get('its_ITS', 0),
            'its_ITS1': stats.get('its_ITS1', 0),
            'its_ITS2': stats.get('its_ITS2', 0),
            'optimizations_applied': [
                'zip_file_processing',
                'sintax_format_parsing', 
                'barcode_type_detection',
                'early_sequence_validation',
                'optimized_regex_parsing',
                'lru_cache_taxonomy_parsing',
                'bulk_database_operations',
                'streaming_chunk_processing',
                'thread_pool_parallelism',
                'prepared_sql_statements',
                'transaction_batching',
                'dataclass_optimization'
            ],
            'performance_metrics': {
                'cpu_cores_used': os.cpu_count(),
                'memory_optimizations': True,
                'batch_size': BATCH_SIZE,
                'chunk_size': FASTA_CHUNK_SIZE
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
        logger.critical(f"Fatal error in PLANITS processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()