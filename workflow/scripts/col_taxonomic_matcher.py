#!/usr/bin/env python3
"""
Improved CoL Taxonomic Matcher
workflow/scripts/CoL_taxonomic_matcher.py

Matches taxonomies against the Catalogue of Life (CoL) API using a modular,
responsibility-separated design.
"""

import sys
import logging
import json
import time
import re
import threading
import requests
from requests.adapters import HTTPAdapter
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Iterator
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from sqlalchemy import text

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)
logger = logging.getLogger(__name__)

@dataclass
class TaxonomyRecord:
    """Structured container for taxonomy data."""
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
    """Result of a taxonomy matching attempt."""
    original_tax_id: int
    col_taxonomy: Optional[Dict[str, str]] = None
    col_id: Optional[str] = None  # accepted CoL name-usage id (for reproducibility)
    success: bool = False
    error_message: Optional[str] = None

class CoLAPIClient:
    """ChecklistBank (Catalogue of Life) client using the `nameusage/search` endpoint.

    Rate limiting is THREAD-SAFE (shared lock) and polite: a realistic default rate,
    HTTP 429 handled with Retry-After/backoff, and a descriptive User-Agent. This
    replaces the previous `/match/nameusage` client whose limiter was not thread-safe
    and could not resolve synonyms (see README / SATIVA_modifications style notes).
    """

    HOST = "https://api.checklistbank.org"

    def __init__(self, max_requests_per_second: float = 8.0, timeout: int = 30,
                 contact: str = "https://github.com/Aaramis/CurateMake",
                 dataset: str = "3LR"):
        self.session = requests.Session()
        adapter = HTTPAdapter(pool_connections=20, pool_maxsize=20, max_retries=3, pool_block=True)
        self.session.mount('http://', adapter)
        self.session.mount('https://', adapter)
        self.session.headers.update({'User-Agent': f'CurateMake/1.0 (ITS reference DB curation; +{contact})'})

        # Thread-safe rate limiter. The cap is a safety ceiling against runaway loops;
        # effective throughput is bounded by worker concurrency, not this rate.
        rps = max(0.5, min(float(max_requests_per_second or 150), 500.0))
        self.request_interval = 1.0 / rps
        self._lock = threading.Lock()
        self._last_request = 0.0
        self.timeout = timeout

        # Metrics
        self.api_calls = 0
        self.failed_calls = 0

        # REPRODUCIBILITY: resolve and pin a specific, DOI'd CoL release. `dataset` may be
        # "3LR" (latest COL release, resolved & pinned at run start), a numeric release key,
        # or "3" (the live, moving COL project — NOT reproducible). All matching queries then
        # hit this fixed snapshot; its version + DOI are recorded in provenance.
        self.provenance = self._resolve_dataset(str(dataset))
        self.dataset_key = self.provenance["dataset_key"]

    def _resolve_dataset(self, dataset: str) -> Dict[str, Any]:
        meta = self._get(f"dataset/{dataset}") or {}
        key = meta.get("key")
        if not key:  # resolution failed -> fall back but flag it
            key = dataset if str(dataset).isdigit() else 3
            logger.warning(f"Could not resolve CoL dataset '{dataset}'; falling back to key {key} (provenance unknown)")
        prov = {
            "requested": dataset,
            "dataset_key": key,
            "alias": meta.get("alias"),
            "version": meta.get("version"),
            "issued": meta.get("issued"),
            "doi": meta.get("doi"),
            "title": meta.get("title"),
            "resolved_at": datetime.now().isoformat(),
            "api_host": self.HOST,
        }
        logger.info("CoL source pinned: %s (%s, issued %s) DOI=%s dataset_key=%s",
                    prov.get("alias"), prov.get("version"), prov.get("issued"), prov.get("doi"), key)
        return prov

    def _throttle(self):
        """Serialize + space out requests across all threads."""
        with self._lock:
            wait = self.request_interval - (time.time() - self._last_request)
            if wait > 0:
                time.sleep(wait)
            self._last_request = time.time()

    def _get(self, path: str, params: Optional[Dict[str, str]] = None) -> Optional[Dict]:
        for attempt in range(4):
            self._throttle()
            try:
                r = self.session.get(f"{self.HOST}/{path}", params=params or {}, timeout=self.timeout)
                with self._lock:
                    self.api_calls += 1
                if r.status_code == 200:
                    return r.json()
                if r.status_code == 429:  # rate limited -> respect Retry-After then retry
                    retry = r.headers.get("Retry-After")
                    time.sleep(min(float(retry) if retry else 2 ** attempt, 30) + 0.5)
                    continue
                logger.warning(f"CoL API error {r.status_code} for {path} params: {params}")
                with self._lock:
                    self.failed_calls += 1
                return None
            except requests.RequestException as e:
                logger.warning(f"CoL request failed for {path} {params}: {e}")
                time.sleep(2 ** attempt)
        with self._lock:
            self.failed_calls += 1
        return None

    def search(self, name: str, match_type: str = "EXACT", limit: int = 20) -> List[Dict]:
        """Search the pinned CoL release name index. type=EXACT avoids cross-kingdom epithet
        noise and returns the synonym's accepted usage in the SAME response."""
        d = self._get(f"dataset/{self.dataset_key}/nameusage/search",
                       {"q": name, "type": match_type, "content": "SCIENTIFIC_NAME", "limit": limit})
        return (d or {}).get("result", []) if d else []

    def close(self):
        """Close the HTTP session."""
        self.session.close()

_PROVISIONAL_WORDS = {
    "sp.", "sp", "cf.", "cf", "aff.", "aff", "spp", "spp.", "uncultured", "unidentified",
    "unknown", "unclassified", "environmental", "metagenome", "endophyte", "sample",
    "isolate", "clone", "voucher", "bacterium", "genomosp", "incertae",
}
# CoL rank name -> our TaxonomyRecord field
_RANK_KEY = {"domain": "kingdom", "kingdom": "kingdom", "phylum": "phylum",
             "class": "class_", "order": "order_", "family": "family",
             "genus": "genus", "species": "species"}


class TaxonomyMatcher:
    """Taxonomy matching against CoL using the EXACT-search strategy.

    `/nameusage/search` with `type=EXACT&content=SCIENTIFIC_NAME` (name only, no context
    in the query) returns the exact usage AND — for a synonym — its accepted name in ONE
    call. Homonyms are disambiguated by filtering the exact results on the expected
    kingdom. Misspellings fall back to a FUZZY search constrained to the same genus.
    Downstream ranks use the same class_/order_ keys the DB layer expects.
    """

    LEVELS = ["species", "genus", "family", "order_", "class_", "phylum", "kingdom"]

    def __init__(self, api_client: CoLAPIClient):
        self.api_client = api_client

    @staticmethod
    def _clean(value: Optional[str]) -> Optional[str]:
        if not value:
            return None
        v = str(value).strip()
        if v.lower() in ('', 'na', 'n/a', 'null', 'none', 'unknown', 'unclassified', 'unidentified'):
            return None
        v = re.sub(r"\s*\(.*", "", v)          # drop author citation, e.g. "(Schaeff.) Pers."
        v = v.replace('_', ' ').strip()
        return v or None

    @staticmethod
    def _is_provisional(name: Optional[str]) -> bool:
        if not name:
            return True
        words = [w.strip(".,;:") for w in name.lower().split()]
        if any(w in _PROVISIONAL_WORDS for w in words):
            return True
        if len(words) >= 2:                     # strain-code epithet, e.g. "sp. F123"
            ep = words[1]
            if ep and sum(c.isdigit() for c in ep) / len(ep) > 0.4:
                return True
        return False

    @staticmethod
    def _kingdom_of(classification: List[Dict]) -> Optional[str]:
        for c in classification or []:
            if c.get("rank") == "kingdom":
                return c.get("name")
        return None

    def _pick(self, results: List[Dict], expected_kingdom: Optional[str],
              require_genus: Optional[str] = None) -> Optional[Dict]:
        """Filter exact/fuzzy results by kingdom (homonyms) and optional genus (typo guard),
        preferring an accepted usage."""
        cand = results
        if expected_kingdom:
            kf = [r for r in cand
                  if (self._kingdom_of(r.get("classification")) or "").lower() == expected_kingdom.lower()]
            if kf:
                cand = kf
        if require_genus:
            cand = [r for r in cand if any(
                c.get("rank") == "genus" and (c.get("name", "").lower() == require_genus.lower())
                for c in r.get("classification", []))]
        if not cand:
            return None
        accepted = [r for r in cand if r.get("usage", {}).get("status") == "accepted"]
        return accepted[0] if accepted else cand[0]

    @staticmethod
    def _accepted_id(usage: Dict) -> Optional[str]:
        """The accepted name-usage id of a chosen result (for a synonym, the id of its
        accepted usage, not the synonym's own id)."""
        if usage.get("status") == "synonym":
            return (usage.get("accepted") or {}).get("id")
        return usage.get("id")

    def _resolve(self, result: Dict) -> Dict[str, str]:
        """Build a {kingdom..species} dict from a chosen result, resolving synonym->accepted."""
        usage = result.get("usage", {})
        tax: Dict[str, str] = {}
        # higher lineage (for a synonym, CoL's classification is the ACCEPTED placement)
        for c in result.get("classification", []):
            rank = (c.get("rank") or "").lower()
            key = _RANK_KEY.get(rank)
            if key and c.get("name"):
                tax[key] = c["name"]
        # deepest accepted name
        if usage.get("status") == "synonym":
            acc_name = usage.get("accepted", {}).get("name", {})
        else:
            acc_name = usage.get("name", {})
        if isinstance(acc_name, dict):
            key = _RANK_KEY.get((acc_name.get("rank") or "").lower())
            if key and acc_name.get("scientificName"):
                tax[key] = acc_name["scientificName"]
        return tax

    def match_taxonomy(self, taxonomy: TaxonomyRecord) -> MatchResult:
        """Match a single taxonomy against CoL with the EXACT-first strategy."""
        try:
            cleaned = {lvl: self._clean(getattr(taxonomy, lvl, None)) for lvl in self.LEVELS}
            exp_kingdom = cleaned.get("kingdom")
            # a provisional species (sp., cf., strain code...) must not be matched as a species
            if cleaned.get("species") and self._is_provisional(cleaned["species"]):
                cleaned["species"] = None

            for lvl in self.LEVELS:                       # deepest -> shallowest
                q = cleaned.get(lvl)
                if not q or self._is_provisional(q):
                    continue
                hit = self._pick(self.api_client.search(q, "EXACT"), exp_kingdom)
                if hit is None and lvl in ("species", "genus"):
                    # typo fallback: fuzzy, but constrained to the same genus to avoid noise
                    hit = self._pick(self.api_client.search(q, "FUZZY"), exp_kingdom,
                                     require_genus=cleaned.get("genus"))
                if hit is not None:
                    tax = self._resolve(hit)
                    # keep original, more-specific ranks CoL did not return
                    for r in self.LEVELS:
                        if cleaned.get(r) and r not in tax:
                            tax[r] = cleaned[r]
                    return MatchResult(original_tax_id=taxonomy.tax_id, col_taxonomy=tax,
                                       col_id=self._accepted_id(hit.get("usage", {})),
                                       success=bool(tax))

            return MatchResult(original_tax_id=taxonomy.tax_id, success=False)

        except Exception as e:
            logger.error(f"Error matching taxonomy {taxonomy.tax_id}: {e}")
            return MatchResult(original_tax_id=taxonomy.tax_id, success=False, error_message=str(e))

class CoLProcessor:
    """High-level processor for CoL matching."""

    def __init__(self, db_manager, config: Dict):
        self.db_manager = db_manager
        self.config = config

        # API configuration (polite defaults for the public ChecklistBank API)
        api_config = config.get('CoL_api', {})
        max_requests = api_config.get('max_requests_per_second', 150)
        # Concurrent CoL requests. ChecklistBank handles high concurrency well (~470 req/s,
        # 0 errors), so allow a generous cap; the shared thread-safe limiter is the safety net.
        max_workers = max(1, min(int(api_config.get('max_workers', 20)), 32))
        # CoL release to pin for reproducibility ("3LR" = latest release, resolved at run start)
        dataset = str(api_config.get('dataset', '3LR'))

        self.api_client = CoLAPIClient(max_requests_per_second=max_requests, dataset=dataset)
        self.matcher = TaxonomyMatcher(self.api_client)
        self.max_workers = max_workers

        # Report metadata — records the exact CoL release (version + DOI) used, so the run
        # is reproducible and citable (addresses reviewer note on CoL provenance).
        self.report = {
            "run_summary": {},
            "parameters": {
                "max_requests_per_second": max_requests,
                "max_workers": max_workers,
                "col_source": self.api_client.provenance,
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
        """Retrieve pending taxonomies that must be processed."""
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

        def fetch_rows():
            with self.db_manager.get_session() as session:
                result = session.execute(query, {"backbone": backbone, "status": status})
                return list(result)

        rows = self.db_manager._execute_with_retry(
            "col_match_get_pending_taxonomies",
            fetch_rows
        )

        taxonomies = []
        for row in rows:
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
        """Process a batch of taxonomies using a thread pool."""
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
        """Persist successful matches to the database."""
        successful_matches = [r for r in results if r.success and r.col_taxonomy]

        if not successful_matches:
            return 0

        def update_op():
            # Build taxonomy payload once to avoid per-row get_or_create calls
            taxonomy_data = []
            taxonomy_tuple_by_original = {}
            for result in successful_matches:
                col_taxonomy = result.col_taxonomy or {}
                tax_values = {
                    'kingdom': col_taxonomy.get('kingdom'),
                    'phylum': col_taxonomy.get('phylum'),
                    'class_': col_taxonomy.get('class_'),
                    'order_': col_taxonomy.get('order_'),
                    'family': col_taxonomy.get('family'),
                    'genus': col_taxonomy.get('genus'),
                    'species': col_taxonomy.get('species')
                }
                taxonomy_data.append(tax_values)
                taxonomy_tuple_by_original[result.original_tax_id] = (
                    tax_values['kingdom'],
                    tax_values['phylum'],
                    tax_values['class_'],
                    tax_values['order_'],
                    tax_values['family'],
                    tax_values['genus'],
                    tax_values['species']
                )

            tax_id_map_by_tuple = self.db_manager.bulk_insert_taxonomies(taxonomy_data)
            tax_id_map = {
                original_tax_id: tax_id_map_by_tuple.get(tax_tuple)
                for original_tax_id, tax_tuple in taxonomy_tuple_by_original.items()
            }

            # Update existing relations
            original_ids = [tid for tid in tax_id_map.keys() if tid]
            if not original_ids:
                return 0

            seq_tax_table = self._get_table('seq_tax')

            with self.db_manager.get_session() as session:
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

                # Prepare CoL backbone insertions
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

                # Record CoL provenance (accepted usage id + pinned release) on each matched
                # taxonomy row. This also stamps rows that already existed from the Raw
                # ingestion (shared lineage): a plain UPDATE, so it covers both cases.
                taxonomy_table = self._get_table('taxonomy')
                dataset_key = self.api_client.dataset_key
                seen_tid = set()
                col_prov_updates = []
                for result in successful_matches:
                    new_tid = tax_id_map.get(result.original_tax_id)
                    if new_tid and result.col_id and new_tid not in seen_tid:
                        seen_tid.add(new_tid)
                        col_prov_updates.append(
                            {"tid": new_tid, "cid": result.col_id, "dk": dataset_key})
                if col_prov_updates:
                    session.execute(text(f"""
                        UPDATE {taxonomy_table}
                        SET col_id = :cid, col_dataset_key = :dk
                        WHERE tax_id = :tid
                    """), col_prov_updates)

                # Mark legacy relations as correct in one update
                if self.db_manager.db_type == "postgresql":
                    session.execute(text(f"""
                        UPDATE {seq_tax_table} SET status = 'Correct'
                        WHERE tax_id = ANY(:tax_ids) AND backbone = 'Raw'
                    """), {"tax_ids": original_ids})
                else:
                    placeholders = ','.join([f':id{i}' for i in range(len(original_ids))])
                    params = {f'id{i}': tax_id for i, tax_id in enumerate(original_ids)}
                    session.execute(text(f"""
                        UPDATE {seq_tax_table} SET status = 'Correct'
                        WHERE tax_id IN ({placeholders}) AND backbone = 'Raw'
                    """), params)

                session.commit()

            return len(successful_matches)

        try:
            return self.db_manager._execute_with_retry(
                "col_match_update_database",
                update_op
            )
        except Exception as e:
            logger.error(f"Database update failed: {e}")
            raise

    def mark_failed_taxonomies(self):
        """Mark taxonomies that never matched as Failed."""
        seq_tax_table = self._get_table('seq_tax')
        def update_op():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    UPDATE {seq_tax_table}
                    SET status = 'Failed'
                    WHERE backbone = 'Raw' AND status = 'Pending'
                """))
                failed_count = result.rowcount
                session.commit()
                return failed_count

        failed_count = self.db_manager._execute_with_retry(
            "col_match_mark_failed",
            update_op
        )
        logger.info(f"Marked {failed_count} taxonomies as Failed")
        return failed_count

    def remove_duplicate_taxonomies(self):
        """Remove duplicate taxonomies from the taxonomy table."""
        taxonomy_table = self._get_table('taxonomy')
        seq_tax_table = self._get_table('seq_tax')

        # Adjust expressions per database type
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

        def dedupe_op():
            with self.db_manager.get_session() as session:
                result = session.execute(query)
                duplicates = result.fetchall()

                if not duplicates:
                    return 0

                removed_count = 0
                for dup in duplicates:
                    keep_id = dup[7]  # keep_id (UUID string)
                    all_ids = [x.strip() for x in dup[8].split(',')]  # all_ids (UUID strings)
                    remove_ids = [x for x in all_ids if x != keep_id]

                    if remove_ids:
                        # Update seq_tax entries
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
                return removed_count

        removed_count = self.db_manager._execute_with_retry(
            "col_match_remove_duplicates",
            dedupe_op
        )
        if removed_count == 0:
            logger.info("No duplicate taxonomies found")
        else:
            logger.info(f"Removed {removed_count} duplicate taxonomies")
        return removed_count

    def run(self, batch_size: int = 100):
        """Execute the main processing loop."""
        start_time = time.time()
        self.report["run_summary"]["start_time"] = datetime.now().isoformat()
        
        logger.info("Starting CoL taxonomic matching...")
        
        try:
            # Retrieve taxonomies
            taxonomies = self.get_pending_taxonomies()
            total_taxonomies = len(taxonomies)
            
            if not taxonomies:
                logger.info("No taxonomies found for processing")
                self._finalize_report(start_time, 0, 0, 0)
                return
            
            # Batch processing
            total_matched = 0
            for i in range(0, len(taxonomies), batch_size):
                batch = taxonomies[i:i + batch_size]
                results = self.process_batch(batch)
                matched_in_batch = self.update_database(results)
                total_matched += matched_in_batch
                
                progress = (i + len(batch)) / len(taxonomies) * 100
                logger.info(f"Progress: {progress:.1f}% - Batch matched: {matched_in_batch}/{len(batch)}")
            
            # Final cleanup
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
        """Finalize and persist statistics in the in-memory report."""
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
        """Write the JSON report to disk."""
        logger.info(f"Writing report to: {output_path}")
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with output_file.open('w', encoding='utf-8') as f:
            json.dump(self.report, f, indent=2, ensure_ascii=False)

def main():
    """Script entry point."""
    try:
        # Snakemake or standalone configuration
        if 'snakemake' in globals():
            db_path = snakemake.params.db
            batch_size = snakemake.params.batch_size
            config = snakemake.config
            output_report = snakemake.output[0]
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            # Logging configuration for Snakemake
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

        # Import and initialize the database manager
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        # Main processing
        with CoLProcessor(db_manager, config) as processor:
            processor.run(batch_size=batch_size)
            processor.write_report(output_report)

            # Display aggregated statistics
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
