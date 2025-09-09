#!/usr/bin/env python3
"""
CoL Taxonomic Matcher - Version améliorée
workflow/scripts/CoL_taxonomic_matcher.py

Script pour matcher les taxonomies avec l'API Catalogue of Life (CoL).
Utilise une approche modulaire avec séparation des responsabilités.
"""

import sys
import logging
import json
import time
import requests
from requests.adapters import HTTPAdapter
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Iterator
from dataclasses import dataclass
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)
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

@dataclass
class MatchResult:
    """Résultat d'un matching taxonomique."""
    original_tax_id: int
    col_taxonomy: Optional[Dict[str, str]] = None
    success: bool = False
    error_message: Optional[str] = None

class CoLAPIClient:
    """Client pour l'API Catalogue of Life avec gestion du rate limiting."""
    
    def __init__(self, max_requests_per_second: int = 100, timeout: int = 30):
        self.base_url = "https://api.checklistbank.org/dataset/3LR/match/nameusage"
        self.session = requests.Session()
        
        # Configuration de la session
        adapter = HTTPAdapter(
            pool_connections=20,
            pool_maxsize=20,
            max_retries=3,
            pool_block=False
        )
        self.session.mount('http://', adapter)
        self.session.mount('https://', adapter)
        self.session.headers.update({'User-Agent': 'BioinformaticsPipeline/1.0'})
        
        # Rate limiting
        self.max_requests_per_second = max_requests_per_second
        self.request_interval = 1.0 / max_requests_per_second if max_requests_per_second > 0 else 0
        self.last_request_time = 0
        self.timeout = timeout
        
        # Statistiques
        self.api_calls = 0
        self.failed_calls = 0

    def rate_limit(self):
        """Implémente le rate limiting pour les appels API."""
        if self.request_interval <= 0:
            return
            
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        if time_since_last_request < self.request_interval:
            sleep_time = self.request_interval - time_since_last_request
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()

    def query_taxon(self, params: Dict[str, str]) -> Optional[Dict]:
        """Effectue une requête vers l'API CoL."""
        try:
            self.rate_limit()
            response = self.session.get(self.base_url, params=params, timeout=self.timeout)
            self.api_calls += 1
            
            if response.status_code != 200:
                logger.warning(f"API error {response.status_code} for params: {params}")
                self.failed_calls += 1
                return None
            
            data = response.json()
            
            if data.get("match") and "usage" in data:
                return data["usage"]
            
            return None
            
        except requests.RequestException as e:
            logger.warning(f"Request failed for params {params}: {e}")
            self.failed_calls += 1
            return None

    def close(self):
        """Ferme la session."""
        self.session.close()

class TaxonomyMatcher:
    """Logique de matching taxonomique."""
    
    def __init__(self, api_client: CoLAPIClient):
        self.api_client = api_client
        self.taxonomy_levels = ["species", "genus", "family", "order_", "class_", "phylum", "kingdom"]

    def clean_taxonomy(self, taxonomy: TaxonomyRecord) -> Dict[str, Optional[str]]:
        """Nettoie et prétraite les données taxonomiques."""
        cleaned = {}
        
        for level in self.taxonomy_levels:
            value = getattr(taxonomy, level, None)
            
            if value:
                cleaned_value = value.replace('_', ' ').strip()
                if any(term in cleaned_value.lower() for term in ['unknown', 'n/a', 'null', 'none', 'na']):
                    cleaned_value = None
                elif cleaned_value.endswith(' sp'):
                    cleaned_value = cleaned_value[:-3].strip()
                
                cleaned[level] = cleaned_value if cleaned_value else None
            else:
                cleaned[level] = None
        
        return cleaned

    def match_taxonomy(self, taxonomy: TaxonomyRecord) -> MatchResult:
        """Effectue le matching pour une taxonomie donnée."""
        try:
            cleaned_tax = self.clean_taxonomy(taxonomy)
            col_taxonomy = self._hierarchical_match(cleaned_tax)
            
            return MatchResult(
                original_tax_id=taxonomy.tax_id,
                col_taxonomy=col_taxonomy,
                success=col_taxonomy is not None
            )
            
        except Exception as e:
            logger.error(f"Error matching taxonomy {taxonomy.tax_id}: {e}")
            return MatchResult(
                original_tax_id=taxonomy.tax_id,
                success=False,
                error_message=str(e)
            )

    def _hierarchical_match(self, cleaned_tax: Dict[str, Optional[str]]) -> Optional[Dict[str, str]]:
        """Effectue le matching hiérarchique."""
        for i, level in enumerate(self.taxonomy_levels):
            if not cleaned_tax.get(level):
                continue
            
            q_name = cleaned_tax[level]
            params = {"q": q_name, "fuzzy": "true"}
            
            # Ajouter le contexte taxonomique
            for context_level in self.taxonomy_levels[i+1:]:
                if cleaned_tax.get(context_level):
                    api_param = context_level.replace('_', '')
                    params[api_param] = cleaned_tax[context_level]
            
            usage = self.api_client.query_taxon(params)
            if usage:
                classification = self._extract_classification(usage)
                if classification:
                    # Compléter avec les niveaux plus spécifiques
                    self._fill_missing_levels(classification, cleaned_tax, i)
                    return classification
        
        return None

    def _extract_classification(self, usage: Dict) -> Dict[str, str]:
        """Extrait la classification du résultat API."""
        classification = {}
        
        # Extraire de la classification
        if "classification" in usage:
            for taxon in usage["classification"]:
                if taxon.get('rank') and taxon.get('name'):
                    rank = taxon['rank'].lower().replace('class', 'class_').replace('order', 'order_')
                    classification[rank] = taxon['name']

        # Ajouter le taxon principal
        if usage.get("rank") and usage.get("name"):
            rank = usage["rank"].lower().replace('class', 'class_').replace('order', 'order_')
            classification[rank] = usage["name"]

        return classification

    def _fill_missing_levels(self, classification: Dict[str, str], original_tax: Dict[str, Optional[str]], found_index: int):
        """Complète les niveaux manquants avec les données originales."""
        for idx in range(found_index):
            specific_level = self.taxonomy_levels[idx]
            if original_tax.get(specific_level) and specific_level not in classification:
                classification[specific_level] = original_tax[specific_level]

class CoLProcessor:
    """Processeur principal pour le matching CoL."""
    
    def __init__(self, db_path, config: Dict):
        self.db_path = Path(db_path)
        # self.db_manager = db_manager
        self.config = config
        
        # Configuration API
        api_config = config.get('CoL_api', {})
        max_requests = api_config.get('max_requests_per_second', 1000)
        max_workers=min(20, os.cpu_count() * 5)

        self.api_client = CoLAPIClient(max_requests_per_second=max_requests)
        self.matcher = TaxonomyMatcher(self.api_client)
        self.max_workers = max_workers
        
        # Rapport
        self.report = {
            "run_summary": {},
            "parameters": {
                "max_requests_per_second": max_requests,
                "max_workers": max_workers,
                **config
            },
            "statistics": {}
        }

    def __enter__(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.api_client.close()
        if self.conn:
            self.conn.close()

    def get_pending_taxonomies(self, backbone: str = 'Raw', status: str = 'Pending') -> List[TaxonomyRecord]:
        """Récupère les taxonomies à traiter."""
        query = """
            SELECT DISTINCT t.tax_id, t.kingdom, t.phylum, t.class_, t.order_, 
                   t.family, t.genus, t.species
            FROM taxonomy t
            INNER JOIN seq_tax st ON t.tax_id = st.tax_id
            WHERE st.backbone = ? AND st.status = ?
            ORDER BY t.tax_id
        """
        cursor = self.conn.execute(query, (backbone, status))
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

    def process_batch(self, taxonomies: List[TaxonomyRecord]) -> List[MatchResult]:
        """Traite un lot de taxonomies avec threading."""
        if not taxonomies:
            return []
        
        results = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_taxonomy = {
                executor.submit(self.matcher.match_taxonomy, taxonomy): taxonomy
                for taxonomy in taxonomies
            }
            
            for future in as_completed(future_to_taxonomy):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    taxonomy = future_to_taxonomy[future]
                    logger.error(f"Future failed for taxonomy {taxonomy.tax_id}: {e}")
                    results.append(MatchResult(
                        original_tax_id=taxonomy.tax_id,
                        success=False,
                        error_message=str(e)
                    ))
        
        return results

    def update_database(self, results: List[MatchResult]):
        """Met à jour la base de données avec les résultats."""
        successful_matches = [r for r in results if r.success and r.col_taxonomy]
        
        if not successful_matches:
            return 0
        
        cursor = self.conn.cursor()
        
        try:
            # Mapper les IDs
            tax_id_map = {}
            for result in successful_matches:
                col_taxonomy = result.col_taxonomy
                tax_values = [col_taxonomy.get(lvl) for lvl in 
                                ["kingdom", "phylum", "class_", "order_", "family", "genus", "species"]]
                
                # Vérifier si existe
                cursor.execute("""
                    SELECT tax_id FROM taxonomy 
                    WHERE kingdom IS ? AND phylum IS ? AND class_ IS ? AND
                            order_ IS ? AND family IS ? AND genus IS ? AND species IS ?
                """, tax_values)
                existing = cursor.fetchone()

                if existing:
                    new_tax_id = existing[0]
                else:
                    cursor.execute("""
                        INSERT INTO taxonomy (kingdom, phylum, class_, order_, family, genus, species, parent_tax_id)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, (*tax_values, result.original_tax_id))
                    new_tax_id = cursor.lastrowid
                
                tax_id_map[result.original_tax_id] = new_tax_id

            # Mise à jour des relations
            original_ids = list(tax_id_map.keys())
            placeholders = ','.join('?' for _ in original_ids)
            cursor.execute(f"""
                SELECT seq_id, tax_id FROM seq_tax
                WHERE tax_id IN ({placeholders}) AND backbone = 'Raw'
            """, original_ids)

            # Préparer les mises à jour
            col_inserts = []
            for seq_id, original_tax_id in cursor.fetchall():
                new_tax_id = tax_id_map.get(original_tax_id)
                if new_tax_id:
                    col_inserts.append((seq_id, new_tax_id))

            # Exécuter les mises à jour
            if col_inserts:
                cursor.executemany("""
                    INSERT OR REPLACE INTO seq_tax (seq_id, tax_id, backbone, status)
                    VALUES (?, ?, 'CoL', 'Correct')
                """, col_inserts)

            # Marquer comme correct
            cursor.executemany("""
                UPDATE seq_tax SET status = 'Correct'
                WHERE tax_id = ? AND backbone = 'Raw'
            """, [(tid,) for tid in original_ids])

            self.conn.commit()
            return len(successful_matches)
            
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Database update failed: {e}")
            raise

    def mark_failed_taxonomies(self):
        """Marque les taxonomies non matchées comme Failed."""
        cursor = self.conn.execute("""
            UPDATE seq_tax
            SET status = 'Failed'
            WHERE backbone = 'Raw' AND status = 'Pending'
        """)
        failed_count = cursor.rowcount
        logger.info(f"Marked {failed_count} taxonomies as Failed")
        return failed_count

    def remove_duplicate_taxonomies(self):
        """Supprime les doublons dans la table taxonomy."""
        cursor = self.conn.execute("""
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
            logger.info("No duplicate taxonomies found")
            return 0

        removed_count = 0
        for dup in duplicates:
            keep_id = dup['keep_id']
            all_ids = [int(x) for x in dup['all_ids'].split(',')]
            remove_ids = [x for x in all_ids if x != keep_id]
            
            if remove_ids:
                placeholders = ','.join('?' for _ in remove_ids)
                self.conn.execute(f"UPDATE seq_tax SET tax_id = ? WHERE tax_id IN ({placeholders})", 
                            (keep_id, *remove_ids))
                self.conn.execute(f"DELETE FROM taxonomy WHERE tax_id IN ({placeholders})", remove_ids)
                removed_count += len(remove_ids)
        
        logger.info(f"Removed {removed_count} duplicate taxonomies")
        return removed_count

    def run(self, batch_size: int = 100):
        """Lance le traitement principal."""
        start_time = time.time()
        self.report["run_summary"]["start_time"] = datetime.now().isoformat()
        
        logger.info("Starting CoL taxonomic matching...")
        
        try:
            # Récupérer les taxonomies
            taxonomies = self.get_pending_taxonomies()
            total_taxonomies = len(taxonomies)
            
            if not taxonomies:
                logger.info("No taxonomies found for processing")
                self._finalize_report(start_time, 0, 0, 0)
                return
            
            # Traitement par lots
            total_matched = 0
            for i in range(0, len(taxonomies), batch_size):
                batch = taxonomies[i:i + batch_size]
                results = self.process_batch(batch)
                matched_in_batch = self.update_database(results)
                total_matched += matched_in_batch
                
                progress = (i + len(batch)) / len(taxonomies) * 100
                logger.info(f"Progress: {progress:.1f}% - Batch matched: {matched_in_batch}/{len(batch)}")
            
            # Nettoyage final
            duplicates_removed = self.remove_duplicate_taxonomies()
            failed_marked = self.mark_failed_taxonomies()
            
            self._finalize_report(start_time, total_taxonomies, total_matched, failed_marked)
            
            logger.info(f"CoL matching completed successfully")
            logger.info(f"Matched: {total_matched}, Failed: {failed_marked}")
            
        except Exception as e:
            logger.error(f"Error in CoL processing: {e}", exc_info=True)
            self.report["run_summary"]["status"] = "Error"
            self.report["run_summary"]["error_message"] = str(e)
            raise

    def _finalize_report(self, start_time: float, total: int, matched: int, failed: int):
        """Finalise le rapport."""
        end_time = time.time()
        processing_time = end_time - start_time
        
        self.report["run_summary"].update({
            "end_time": datetime.now().isoformat(),
            "total_duration_seconds": round(processing_time, 2),
            "status": "Success"
        })
        
        self.report["statistics"].update({
            "total_taxonomies": total,
            "matched_taxonomies": matched,
            "failed_taxonomies": failed,
            "api_calls": self.api_client.api_calls,
            "failed_api_calls": self.api_client.failed_calls,
            "success_rate": (matched / total * 100) if total > 0 else 0,
            "avg_time_per_call": processing_time / max(self.api_client.api_calls, 1)
        })

    def write_report(self, output_path: str):
        """Écrit le rapport JSON."""
        logger.info(f"Writing report to: {output_path}")
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with output_file.open('w', encoding='utf-8') as f:
            json.dump(self.report, f, indent=2, ensure_ascii=False)

def main():
    """Point d'entrée principal."""
    try:
        # Configuration Snakemake ou standalone
        if 'snakemake' in globals():
            db_path = snakemake.params.db
            batch_size = snakemake.params.batch_size
            config = snakemake.config
            output_report = snakemake.output[0]
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            
            # Configuration du logging
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            logger.warning("Running in standalone test mode")
            db_path = "test_output.sqlite"
            config = {
                'CoL_api': {
                    'max_requests_per_second': 10,
                    'max_workers': 4
                }
            }
            batch_size=100
            output_report = "CoL_matching_report.json"
        
        # Traitement principal
        with CoLProcessor(db_path, config) as processor:
            processor.run(batch_size=batch_size)
            processor.write_report(output_report)
            
            # Affichage des résultats
            stats = processor.report["statistics"]
            logger.info("\n" + "="*50)
            logger.info("CoL MATCHING SUMMARY")
            logger.info("="*50)
            logger.info(f"Total taxonomies processed: {stats['total_taxonomies']}")
            logger.info(f"Successfully matched: {stats['matched_taxonomies']}")
            logger.info(f"Failed to match: {stats['failed_taxonomies']}")
            logger.info(f"Success rate: {stats['success_rate']:.1f}%")
            logger.info(f"API calls made: {stats['api_calls']}")
            logger.info("="*50)
        
    except Exception as e:
        logger.critical(f"Fatal error in CoL matching: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()