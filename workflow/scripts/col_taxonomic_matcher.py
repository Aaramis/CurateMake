#!/usr/bin/env python3
"""
CoL Taxonomic Matcher
workflow/scripts/CoL_taxonomic_matcher.py

Script pour matcher les taxonomies avec l'API Catalogue of Life (CoL).
Utilise une approche hiérarchique, concurrente et par lots pour améliorer les performances.
"""

import sys
import logging
import json
import time
import requests
from requests.adapters import HTTPAdapter
import os # NOUVEAU: Pour déterminer le nombre de workers
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import re

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class TaxonomyRecord:
    """Structure pour les données taxonomiques."""
    tax_id: int
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

class CoLMatcher:
    """Classe principale pour le matching taxonomique avec CoL."""
    
    # MODIFIÉ: Ajout de max_workers pour le ThreadPoolExecutor
    def __init__(self, db_manager, max_requests_per_second: int = 100, max_workers: int = 10):
        self.db_manager = db_manager
        self.base_url = "https://api.checklistbank.org/dataset/3LR/match/nameusage"
        self.session = requests.Session()
    
        # Créer un adaptateur avec une taille de pool plus grande
        adapter = HTTPAdapter(
            pool_connections=max_workers,
            pool_maxsize=max_workers,
            max_retries=3,
            pool_block=False
        )
       
        # Monter l'adaptateur pour les requêtes HTTP et HTTPS
        self.session.mount('http://', adapter)
        self.session.mount('https://', adapter)
        
        self.session.headers.update({'User-Agent': 'BioinformaticsPipeline/1.0'})
        
        # Rate limiting
        self.max_requests_per_second = max_requests_per_second
        self.request_interval = 1.0 / max_requests_per_second
        self.last_request_time = 0
        
        # NOUVEAU: Configuration pour la concurrence
        self.max_workers = max_workers
        
        # Statistiques
        self.stats = {
            'total_taxonomies': 0,
            'matched_taxonomies': 0,
            'failed_taxonomies': 0,
            'api_calls': 0,
            'processing_time': 0.0
        }
        
        # Ordre des niveaux taxonomiques (du plus spécifique au plus général)
        self.taxonomy_levels = ["species", "genus", "family", "order_", "class_", "phylum", "kingdom"]

    def rate_limit(self):
        """Implemente le rate limiting pour les appels API."""
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        if time_since_last_request < self.request_interval:
            sleep_time = self.request_interval - time_since_last_request
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()

    def get_taxonomies(self, backbone: str = 'Raw', status:str = 'Pending') -> List[TaxonomyRecord]:
        """
        Récupère les taxonomies de la base de données qui nécessitent un matching CoL.
        """
        query = """
            SELECT DISTINCT t.tax_id, t.kingdom, t.phylum, t.class_, t.order_, 
                   t.family, t.genus, t.species
            FROM taxonomy t
            INNER JOIN seq_tax st ON t.tax_id = st.tax_id
            WHERE st.backbone = ? 
            AND st.status = ?
            ORDER BY t.tax_id
        """
        
        with self.db_manager.get_connection() as conn:
            cursor = conn.execute(query, (backbone, status,))
            taxonomies = []
            
            for row in cursor.fetchall():
                taxonomy = TaxonomyRecord(
                    tax_id=row['tax_id'],
                    kingdom=row['kingdom'],
                    phylum=row['phylum'],
                    class_=row['class_'],
                    order_=row['order_'],
                    family=row['family'],
                    genus=row['genus'],
                    species=row['species']
                )
                taxonomies.append(taxonomy)
            
            logger.info(f"Retrieved {len(taxonomies)} taxonomies for CoL matching")
            return taxonomies

    def clean_taxonomies(self, taxonomy: TaxonomyRecord) -> Dict[str, Optional[str]]:
        """
        Nettoie et prétraite les données taxonomiques.
        """
        cleaned = {}
        
        for level in self.taxonomy_levels:
            value = getattr(taxonomy, level, None)
            
            if value:
                cleaned_value = value.replace('_', ' ').strip()
                if 'unknown' in cleaned_value.lower() or cleaned_value.lower() in ['n/a', 'null', 'none', 'na', '']:
                    cleaned_value = None
                
                if cleaned_value and cleaned_value.endswith(' sp'):
                    cleaned_value = cleaned_value[:-3].strip()
                
                cleaned[level] = cleaned_value if cleaned_value else None
            else:
                cleaned[level] = None
        
        return cleaned

    def match_taxon(self, cleaned_tax: Dict[str, Optional[str]]) -> Optional[Dict]:
        """Logique interne de matching taxonomique."""
        for i, level in enumerate(self.taxonomy_levels):
            if not cleaned_tax.get(level):
                continue
            
            q_name = cleaned_tax[level]
            params = {"q": q_name, "fuzzy": "true"}
            
            for context_level in self.taxonomy_levels[i+1:]:
                if cleaned_tax.get(context_level):
                    api_param = context_level.replace('_', '')
                    params[api_param] = cleaned_tax[context_level]
            
            try:
                self.rate_limit()
                response = self.session.get(self.base_url, params=params, timeout=2)
                self.stats['api_calls'] += 1
                
                if response.status_code != 200:
                    logger.warning(f"API error {response.status_code} for {q_name}")
                    continue
                
                data = response.json()
                
                if data.get("match") and "usage" in data:
                    usage = data["usage"]
                    classification = {}
                    
                    if "classification" in usage:
                        for taxon in usage["classification"]:
                            if taxon.get('rank') and taxon.get('name'):
                                rank = taxon['rank'].lower().replace('class', 'class_').replace('order', 'order_')
                                classification[rank] = taxon['name']

                    if usage.get("rank") and usage.get("name"):
                         rank = usage["rank"].lower().replace('class', 'class_').replace('order', 'order_')
                         classification[rank] = usage["name"]

                    found_level_index = next((idx for idx, tax_level in enumerate(self.taxonomy_levels) if tax_level in classification), None)

                    if found_level_index is not None:
                        for idx in range(found_level_index):
                            specific_level = self.taxonomy_levels[idx]
                            if cleaned_tax.get(specific_level) and specific_level not in classification:
                                classification[specific_level] = cleaned_tax[specific_level]

                    logger.debug(f"Found match for {q_name}: {classification}")
                    return classification
                
            except requests.RequestException as e:
                logger.warning(f"Request failed for {q_name}: {e}")
                continue
        
        return None
    
    # NOUVEAU: Tâche exécutée par chaque thread pour un seul item taxonomique
    def _match_single_taxonomy_task(self, taxonomy: TaxonomyRecord) -> Tuple[int, Optional[Dict]]:
        """
        Tâche isolée qui nettoie et matche une seule taxonomie.
        Retourne l'ID original et le résultat du match.
        """
        try:
            cleaned_tax = self.clean_taxonomies(taxonomy)
            CoL_taxonomy = self.match_taxon(cleaned_tax)
            return (taxonomy.tax_id, CoL_taxonomy)
        except Exception as e:
            logger.error(f"Error in thread for tax_id {taxonomy.tax_id}: {e}")
            return (taxonomy.tax_id, None)

    # NOUVEAU: Gère les écritures en base de données par lot
    def _process_batch_results(self, successful_matches: List[Tuple[int, Dict]]):
        """
        Traite les résultats d'un lot, en groupant les écritures DB.
        C'est comme apporter une liste de courses au supermarché et tout payer en une fois.
        """
        if not successful_matches:
            return

        with self.db_manager.get_connection() as conn:
            cursor = conn.cursor()
            try:
                # Étape 1: Insérer les nouvelles taxonomies et mapper les IDs
                tax_id_map = {} # {original_tax_id: new_CoL_tax_id}
                for original_tax_id, CoL_taxonomy in successful_matches:
                    tax_values = [CoL_taxonomy.get(lvl) for lvl in self.taxonomy_levels]
                    
                    # Vérifier si la taxonomie existe déjà
                    cursor.execute("""
                        SELECT tax_id FROM taxonomy WHERE kingdom IS ? AND phylum IS ? AND class_ IS ? AND
                        order_ IS ? AND family IS ? AND genus IS ? AND species IS ?
                    """, tax_values)
                    existing = cursor.fetchone()

                    if existing:
                        new_tax_id = existing[0] # SQLite row est un tuple
                    else:
                        cursor.execute("""
                            INSERT INTO taxonomy (kingdom, phylum, class_, order_, family, genus, species, parent_tax_id)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        """, (*tax_values, original_tax_id))
                        new_tax_id = cursor.lastrowid
                    tax_id_map[original_tax_id] = new_tax_id

                # Étape 2: Préparer les données pour les mises à jour par lot (executemany)
                raw_updates_params = [] # Pour passer Raw -> 'Correct'
                col_inserts_params = [] # Pour insérer les nouvelles relations CoL

                # Récupérer tous les seq_ids nécessaires en une seule grosse requête
                original_ids_to_query = list(tax_id_map.keys())
                placeholders = ','.join('?' for _ in original_ids_to_query)
                cursor.execute(f"""
                    SELECT seq_id, tax_id FROM seq_tax
                    WHERE tax_id IN ({placeholders}) AND backbone = 'Raw'
                """, original_ids_to_query)

                # Construire les listes de paramètres
                for seq_id, original_tax_id in cursor.fetchall():
                    new_tax_id = tax_id_map.get(original_tax_id)
                    if new_tax_id:
                         # Paramètres pour INSERT OR REPLACE: (seq_id, new_tax_id)
                        col_inserts_params.append((seq_id, new_tax_id))

                # La liste des tax_id à mettre à jour est simplement les clés de tax_id_map
                raw_updates_params = [(original_tax_id,) for original_tax_id in original_ids_to_query]

                # Étape 3: Exécuter les mises à jour en lots
                if raw_updates_params:
                    cursor.executemany("""
                        UPDATE seq_tax SET status = 'Correct'
                        WHERE tax_id = ? AND backbone = 'Raw'
                    """, raw_updates_params)
                    logger.debug(f"Batched update of {len(raw_updates_params)} 'Raw' statuses.")

                if col_inserts_params:
                    cursor.executemany("""
                        INSERT OR REPLACE INTO seq_tax (seq_id, tax_id, backbone, status)
                        VALUES (?, ?, 'CoL', 'Correct')
                    """, col_inserts_params)
                    logger.debug(f"Batched insert of {len(col_inserts_params)} 'CoL' relations.")

                conn.commit() # Valider la transaction pour ce lot

            except Exception as e:
                logger.error(f"Failed to process batch results in DB: {e}")
                conn.rollback()
                raise

    def remove_duplicate_taxonomies(self):
        """
        Supprime les doublons dans la table taxonomy et met à jour les relations.
        """
        with self.db_manager.get_connection() as conn:
            try:
                # Identifier les doublons (même classification, parent_tax_id différent)
                cursor = conn.execute("""
                    SELECT kingdom, phylum, class_, order_, family, genus, species,
                           MIN(tax_id) as keep_id, 
                           GROUP_CONCAT(tax_id) as all_ids
                    FROM taxonomy 
                    WHERE parent_tax_id IS NOT NULL
                    GROUP BY kingdom, phylum, class_, order_, family, genus, species
                    HAVING COUNT(*) > 1
                """)
                
                duplicates = cursor.fetchall()
                if not duplicates:
                    logger.info("No duplicate taxonomies to remove.")
                    return

                removed_count = 0
                for dup in duplicates:
                    keep_id = dup['keep_id']
                    all_ids = [int(x) for x in dup['all_ids'].split(',')]
                    remove_ids = [x for x in all_ids if x != keep_id]
                    
                    placeholders = ','.join('?' for _ in remove_ids)
                    conn.execute(f"UPDATE seq_tax SET tax_id = ? WHERE tax_id IN ({placeholders})", (keep_id, *remove_ids))
                    conn.execute(f"DELETE FROM taxonomy WHERE tax_id IN ({placeholders})", remove_ids)
                    removed_count += len(remove_ids)
                
                logger.info(f"Removed {removed_count} duplicate taxonomies")
                
            except Exception as e:
                logger.error(f"Failed to remove duplicates: {e}")
                raise

    def process_taxonomy_matching(self, batch_size: int = 100) -> Dict[str, Any]:
        """
        Pipeline principal de matching taxonomique.
        """
        start_time = time.time()
        logger.info("Starting CoL taxonomic matching...")
        
        taxonomies = self.get_taxonomies()
        self.stats['total_taxonomies'] = len(taxonomies)
        
        if not taxonomies:
            logger.info("No taxonomies found for CoL matching")
            return self.stats
        
        # MODIFIÉ: Utilise ThreadPoolExecutor pour traiter les lots en parallèle
        matched_count = 0
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for i in range(0, len(taxonomies), batch_size):
                batch = taxonomies[i:i + batch_size]
                
                # Soumettre toutes les tâches du lot au pool de threads
                future_to_tax = {executor.submit(self._match_single_taxonomy_task, tax): tax for tax in batch}
                
                successful_matches_in_batch = []
                
                for future in as_completed(future_to_tax):
                    try:
                        original_tax_id, CoL_result = future.result()
                        if CoL_result:
                            successful_matches_in_batch.append((original_tax_id, CoL_result))
                    except Exception as e:
                        tax = future_to_tax[future]
                        logger.error(f"Error processing future for tax_id {tax.tax_id}: {e}")

                # Une fois tous les résultats du lot reçus, les traiter en base de données
                if successful_matches_in_batch:
                    self._process_batch_results(successful_matches_in_batch)
                    matched_count += len(successful_matches_in_batch)

                progress = (i + len(batch)) / len(taxonomies) * 100
                logger.info(f"Progress: {progress:.1f}% - Batch matched: {len(successful_matches_in_batch)}/{len(batch)}")

        self.remove_duplicate_taxonomies()
        
        failed_count = self.stats['total_taxonomies'] - matched_count
        end_time = time.time()
        self.stats.update({
            'matched_taxonomies': matched_count,
            'failed_taxonomies': failed_count,
            'processing_time': end_time - start_time
        })
        
        self._mark_failed_taxonomies_as_failed()
        
        logger.info(f"CoL matching completed in {self.stats['processing_time']:.2f}s")
        logger.info(f"Matched: {matched_count}, Failed: {failed_count}, API calls: {self.stats['api_calls']}")
        
        return self.stats

    # MODIFIÉ: Renommage pour plus de clarté
    def _mark_failed_taxonomies_as_failed(self):
        """Marque les taxonomies Raw qui n'ont pas pu être matchées comme 'Failed'."""
        with self.db_manager.get_connection() as conn:
            try:
                # Une taxonomie a échoué si elle est 'Raw', toujours 'Pending', et qu'aucune
                # relation 'CoL' n'a été créée pour les séquences associées.
                cursor = conn.execute("""
                    UPDATE seq_tax
                    SET status = 'Failed'
                    WHERE backbone = 'Raw' AND status = 'Pending'
                """)
                failed_count = cursor.rowcount
                logger.info(f"Marked {failed_count} remaining 'Pending' taxonomies as 'Failed'")
            except Exception as e:
                logger.error(f"Failed to mark failed taxonomies: {e}")

    # Le reste du code (find_official_taxonomy, generate_summary_report, main) reste inchangé
    # car il n'est pas directement affecté par les optimisations de performance.
    # ... (Collez ici les fonctions find_official_taxonomy, generate_summary_report, et main de votre code original)
    def find_official_taxonomy(self, tax_hierarchy: Dict[str, str]) -> Optional[Dict]:
        """
        Fonction adaptée de votre exemple pour trouver la taxonomie officielle via CoL.
        """
        # Ordre des niveaux à tester pour q
        levels = ["species", "genus", "family", "order_", "class_", "phylum", "kingdom"]
        
        # Créer une copie modifiable de la hiérarchie
        hierarchy = {k: v for k, v in tax_hierarchy.items() if v}
        
        for i, lvl in enumerate(levels):
            if lvl not in hierarchy:
                continue
            
            q_name = hierarchy[lvl]
            
            # Ne garder que les niveaux supérieurs ou égaux à celui choisi
            params = {"q": q_name, "fuzzy": "true"}
            for key, val in hierarchy.items():
                if levels.index(key) >= i and key != lvl:
                    # Mapper order_ et class_ pour l'API
                    api_key = key.replace('_', '') if key.endswith('_') else key
                    params[api_key] = val
            
            try:
                # Rate limiting et appel API
                self.rate_limit()
                response = self.session.get(self.base_url, params=params, timeout=30)
                self.stats['api_calls'] += 1
                
                if response.status_code != 200:
                    continue
                
                data = response.json()
                
                if data.get("match") and "usage" in data:
                    usage = data["usage"]
                    classification = {}
                    
                    # Extraire la classification du résultat
                    for t in usage.get("classification", []):
                        if t['rank'] and t['name']:
                            rank = t['rank'].lower()
                            if rank == "class":
                                rank = "class_"
                            elif rank == "order":
                                rank = "order_"
                            if rank in levels:
                                classification[rank] = t['name']
                    
                    # Ajouter le taxon principal
                    if usage.get("rank") and usage.get("name"):
                        rank = usage["rank"].lower()
                        if rank == "class":
                            rank = "class_"
                        elif rank == "order":
                            rank = "order_"
                        classification[rank] = usage["name"]
                    
                    return classification
                
            except requests.RequestException as e:
                logger.warning(f"Request failed for {q_name}: {e}")
                continue
            
            # Supprimer le niveau actuel si pas trouvé pour remonter
            hierarchy.pop(lvl)
        
        return None

    def generate_summary_report(self) -> Dict[str, Any]:
        """Génère un rapport de résumé final."""
        
        with self.db_manager.get_connection() as conn:
            # Statistiques CoL
            cursor = conn.execute("""
                SELECT 
                    COUNT(DISTINCT CASE WHEN st.backbone = 'CoL' THEN st.seq_id END) as CoL_matches,
                    COUNT(DISTINCT CASE WHEN st.backbone = 'Raw' AND st.status = 'Failed' THEN st.seq_id END) as failed_matches,
                    COUNT(DISTINCT CASE WHEN st.backbone = 'Raw' THEN st.seq_id END) as total_raw
                FROM seq_tax st
            """)
            db_stats = cursor.fetchone()
            
            # Statistiques par niveau taxonomique pour les entrées CoL uniques
            cursor = conn.execute("""
                SELECT 
                    SUM(CASE WHEN t.species IS NOT NULL THEN 1 ELSE 0 END) as species_level,
                    SUM(CASE WHEN t.genus IS NOT NULL THEN 1 ELSE 0 END) as genus_level,
                    SUM(CASE WHEN t.family IS NOT NULL THEN 1 ELSE 0 END) as family_level
                FROM taxonomy t
                WHERE t.tax_id IN (SELECT DISTINCT tax_id FROM seq_tax WHERE backbone = 'CoL')
            """)
            level_stats = cursor.fetchone()
        
        report = {
            'processing_stats': self.stats,
            'database_stats': {
                'CoL_matches': db_stats['CoL_matches'],
                'failed_matches': db_stats['failed_matches'], 
                'total_raw_taxonomies': db_stats['total_raw'],
                'success_rate': (db_stats['CoL_matches'] / db_stats['total_raw'] * 100) if db_stats['total_raw'] > 0 else 0
            },
            'taxonomy_resolution': {
                'species_level': level_stats['species_level'],
                'genus_level': level_stats['genus_level'],
                'family_level': level_stats['family_level']
            },
            'api_performance': {
                'total_api_calls': self.stats['api_calls'],
                'avg_time_per_call': self.stats['processing_time'] / max(self.stats['api_calls'], 1),
            }
        }
        
        return report

def main():
    """Point d'entrée principal."""
    try:
        if 'snakemake' in globals():
            input_db = snakemake.params.db
            output_report = snakemake.output[0]
            config = snakemake.config
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            logger.warning("Running in standalone test mode")
            input_db = "test_output.sqlite"
            output_report = "CoL_matching_report.json"
            config = {'database': {'path': input_db}}

        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager
        
        db_manager = get_database_manager(config)
        
        rate_limit = config.get('CoL_api', {}).get('max_requests_per_second', 10)
        # NOUVEAU: Définition dynamique du nombre de workers
        # max_workers = config.get('CoL_api', {}).get('max_workers', min(32, (os.cpu_count() or 1) + 4))
        max_workers = config.get('CoL_api', {}).get('max_workers', min(10, (os.cpu_count() or 1) + 2))  # Réduit à 10 maximum

        logger.info(f"Using up to {max_workers} workers for API calls.")
        matcher = CoLMatcher(db_manager, max_requests_per_second=rate_limit, max_workers=max_workers)
        
        matcher.process_taxonomy_matching(batch_size=rate_limit)
        final_report = matcher.generate_summary_report()
        
        if 'snakemake' in globals():
            with open(output_report, 'w', encoding='utf-8') as f:
                json.dump(final_report, f, indent=2, ensure_ascii=False)
            logger.info(f"CoL matching report saved to: {output_report}")
        
        db_stats = final_report['database_stats']
        logger.info("\n" + "="*50)
        logger.info("CoL MATCHING SUMMARY")
        logger.info("="*50)
        logger.info(f"Total taxonomies processed: {db_stats['total_raw_taxonomies']}")
        logger.info(f"Successfully matched with CoL: {db_stats['CoL_matches']}")
        logger.info(f"Failed to match: {db_stats['failed_matches']}")
        logger.info(f"Success rate: {db_stats['success_rate']:.1f}%")
        logger.info(f"API calls made: {final_report['api_performance']['total_api_calls']}")
        logger.info("="*50)
        
    except Exception as e:
        logger.critical(f"Fatal error in CoL matching: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()