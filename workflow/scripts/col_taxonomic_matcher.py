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
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from sqlalchemy import text

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

    def __init__(self, db_manager, config: Dict):
        self.db_manager = db_manager
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
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.api_client.close()

    def _get_table(self, table_name: str) -> str:
        """Get schema-qualified table name for PostgreSQL, unqualified for SQLite."""
        if self.db_manager.db_type == "postgresql":
            schema = self.db_manager.connection_params.get('schema', 'public')
            return f"{schema}.{table_name}"
        return table_name

    def get_pending_taxonomies(self, backbone: str = 'Raw', status: str = 'Pending') -> List[TaxonomyRecord]:
        """Récupère les taxonomies à traiter."""
        taxonomy_table = self._get_table('taxonomy')
        seq_tax_table = self._get_table('seq_tax')

        query = text(f"""
            SELECT DISTINCT t.tax_id, t.kingdom, t.phylum, t.class_, t.order_,
                   t.family, t.genus, t.species
            FROM {taxonomy_table} t
            INNER JOIN {seq_tax_table} st ON t.tax_id = st.tax_id
            WHERE st.backbone = :backbone AND st.status = :status
            ORDER BY t.tax_id
        """)

        taxonomies = []

        with self.db_manager.get_session() as session:
            result = session.execute(query, {"backbone": backbone, "status": status})

            for row in result:
                taxonomy = TaxonomyRecord(
                    tax_id=row[0],
                    kingdom=row[1],
                    phylum=row[2],
                    class_=row[3],
                    order_=row[4],
                    family=row[5],
                    genus=row[6],
                    species=row[7]
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

        try:
            with self.db_manager.get_session() as session:
                # Mapper les IDs
                tax_id_map = {}
                for result in successful_matches:
                    col_taxonomy = result.col_taxonomy
                    tax_values = {
                        'kingdom': col_taxonomy.get('kingdom'),
                        'phylum': col_taxonomy.get('phylum'),
                        'class_': col_taxonomy.get('class_'),
                        'order_': col_taxonomy.get('order_'),
                        'family': col_taxonomy.get('family'),
                        'genus': col_taxonomy.get('genus'),
                        'species': col_taxonomy.get('species')
                    }

                    # Utiliser la méthode helper du db_manager
                    new_tax_id = self.db_manager.get_or_create_taxonomy(tax_values, result.original_tax_id)
                    tax_id_map[result.original_tax_id] = new_tax_id

                # Mise à jour des relations
                original_ids = list(tax_id_map.keys())

                # Préparer les mises à jour
                seq_tax_table = self._get_table('seq_tax')

                if self.db_manager.db_type == "postgresql":
                    result = session.execute(text(f"""
                        SELECT seq_id, tax_id FROM {seq_tax_table}
                        WHERE tax_id = ANY(:tax_ids) AND backbone = 'Raw'
                    """), {"tax_ids": original_ids})
                else:
                    placeholders = ','.join([f':id{i}' for i in range(len(original_ids))])
                    params = {f'id{i}': tax_id for i, tax_id in enumerate(original_ids)}
                    result = session.execute(text(f"""
                        SELECT seq_id, tax_id FROM {seq_tax_table}
                        WHERE tax_id IN ({placeholders}) AND backbone = 'Raw'
                    """), params)

                # Préparer les insertions
                col_inserts = []
                for seq_id, original_tax_id in result:
                    new_tax_id = tax_id_map.get(original_tax_id)
                    if new_tax_id:
                        col_inserts.append({
                            'seq_id': seq_id,
                            'tax_id': new_tax_id,
                            'backbone': 'CoL',
                            'status': 'Correct'
                        })

                # Bulk insert seq_tax relations
                if col_inserts:
                    self.db_manager.bulk_insert_seq_tax(col_inserts)

                # Marquer comme correct
                seq_tax_table = self._get_table('seq_tax')
                for tid in original_ids:
                    session.execute(text(f"""
                        UPDATE {seq_tax_table} SET status = 'Correct'
                        WHERE tax_id = :tax_id AND backbone = 'Raw'
                    """), {"tax_id": tid})

                session.commit()
                return len(successful_matches)

        except Exception as e:
            logger.error(f"Database update failed: {e}")
            raise

    def mark_failed_taxonomies(self):
        """Marque les taxonomies non matchées comme Failed."""
        seq_tax_table = self._get_table('seq_tax')
        with self.db_manager.get_session() as session:
            result = session.execute(text(f"""
                UPDATE {seq_tax_table}
                SET status = 'Failed'
                WHERE backbone = 'Raw' AND status = 'Pending'
            """))
            failed_count = result.rowcount
            session.commit()
            logger.info(f"Marked {failed_count} taxonomies as Failed")
            return failed_count

    def remove_duplicate_taxonomies(self):
        """Supprime les doublons dans la table taxonomy."""
        taxonomy_table = self._get_table('taxonomy')
        seq_tax_table = self._get_table('seq_tax')

        # Adapter la requête selon le type de base de données
        if self.db_manager.db_type == "postgresql":
            group_concat_expr = "STRING_AGG(CAST(tax_id AS TEXT), ',')"
        else:
            group_concat_expr = "GROUP_CONCAT(tax_id)"

        query = text(f"""
            SELECT kingdom, phylum, class_, order_, family, genus, species,
                    MIN(tax_id) as keep_id,
                    {group_concat_expr} as all_ids
            FROM {taxonomy_table}
            WHERE parent_tax_id IS NOT NULL
            GROUP BY kingdom, phylum, class_, order_, family, genus, species
            HAVING COUNT(*) > 1
        """)

        with self.db_manager.get_session() as session:
            result = session.execute(query)
            duplicates = result.fetchall()

            if not duplicates:
                logger.info("No duplicate taxonomies found")
                return 0

            removed_count = 0
            for dup in duplicates:
                keep_id = dup[7]  # keep_id (UUID string)
                all_ids = [x.strip() for x in dup[8].split(',')]  # all_ids (UUID strings)
                remove_ids = [x for x in all_ids if x != keep_id]

                if remove_ids:
                    # Mise à jour seq_tax
                    if self.db_manager.db_type == "postgresql":
                        session.execute(text(f"""
                            UPDATE {seq_tax_table} SET tax_id = :keep_id
                            WHERE tax_id = ANY(:remove_ids)
                        """), {"keep_id": keep_id, "remove_ids": remove_ids})

                        session.execute(text(f"""
                            DELETE FROM {taxonomy_table} WHERE tax_id = ANY(:remove_ids)
                        """), {"remove_ids": remove_ids})
                    else:
                        placeholders = ','.join([f':id{i}' for i in range(len(remove_ids))])
                        params = {f'id{i}': rid for i, rid in enumerate(remove_ids)}
                        params['keep_id'] = keep_id

                        session.execute(text(f"""
                            UPDATE {seq_tax_table} SET tax_id = :keep_id
                            WHERE tax_id IN ({placeholders})
                        """), params)

                        session.execute(text(f"""
                            DELETE FROM {taxonomy_table} WHERE tax_id IN ({placeholders})
                        """), {f'id{i}': rid for i, rid in enumerate(remove_ids)})

                    removed_count += len(remove_ids)

            session.commit()
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
                'database': {'path': db_path},
                'CoL_api': {
                    'max_requests_per_second': 10,
                    'max_workers': 4
                }
            }
            batch_size=100
            output_report = "CoL_matching_report.json"

        # Import et initialisation du database manager
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        # Traitement principal
        with CoLProcessor(db_manager, config) as processor:
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