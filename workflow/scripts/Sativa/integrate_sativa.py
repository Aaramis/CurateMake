#!/usr/bin/env python3
"""
Integrate SATIVA corrections into the database - Version
workflow/scripts/Sativa/integrate_sativa.py

Processes SATIVA .mis/.tax files and integrates taxonomic corrections with CoL validation.

Key features:
- CoL validation for all modifications (except simple deletions)
- Multi-level cache for CoL API calls
- Parallel processing of .mis files
- Bulk database inserts for performance
- Memory-efficient streaming of large files

Compatible: PostgreSQL + SQLite via DatabaseManager
"""

import sys
import json
import logging
import time
import requests
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Set, Iterator
from collections import defaultdict
from dataclasses import dataclass, asdict
from enum import Enum
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
import pickle

# Add utils + sibling scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
sys.path.insert(0, str(Path(__file__).parent.parent))
from db_manager import get_database_manager
# Reuse the single corrected CoL client/matcher (EXACT-search strategy, synonym resolution,
# thread-safe polite rate limiting, pinned DOI'd release + provenance). One source of truth
# for CoL logic across the pipeline (see col_taxonomic_matcher.py).
from col_taxonomic_matcher import CoLAPIClient, TaxonomyMatcher, TaxonomyRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

DISPLAY_TAXONOMY_LEVELS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
DB_TAXONOMY_LEVEL_KEYS = ["kingdom", "phylum", "class_", "order_", "family", "genus", "species"]


# ============================================================================
# Data Classes and Enums
# ============================================================================

class ModificationType(Enum):
    """Types of SATIVA modifications"""
    VALIDATED = "VALIDATED"              # In .tax but not in .mis
    DELETION_SIMPLE = "DELETION_SIMPLE"  # 1 level removed
    DELETION_MULTIPLE = "DELETION_MULTIPLE"  # 2+ levels removed
    SUBSTITUTION = "SUBSTITUTION"        # Name change at any level
    COMPLETION = "COMPLETION"            # Levels added


class CoLMatchQuality(Enum):
    """Quality of CoL match"""
    EXACT = "exact"          # CoL returned exactly what SATIVA proposed
    MODIFIED = "modified"    # CoL returned different taxonomy
    NO_MATCH = "no_match"    # CoL found nothing


@dataclass
class ModificationRecord:
    """Record from SATIVA .mis file"""
    seq_id: str
    seq_original_id: str
    mislabeled_level: str
    original_label: str
    proposed_label: str
    confidence: float
    original_taxonomy: str  # Semicolon-separated
    proposed_taxonomy: str  # Semicolon-separated
    per_rank_confidence: str
    barcode: str
    tax_level: str


@dataclass
class CoLValidationResult:
    """Result from CoL validation"""
    success: bool
    match_quality: Optional[CoLMatchQuality]
    col_taxonomy: Optional[Dict[str, str]]  # {kingdom, phylum, ...}
    error_message: Optional[str] = None


@dataclass
class IntegrationAction:
    """Actions to take for integrating a modification"""
    create_seq_tax: bool
    update_seq_tax_backbone: bool  # True for VALIDATED: update existing backbone to 'Sativa'
    create_warning: bool
    taxonomy_to_use: Dict[str, str]  # Taxonomy dict to insert/use
    seq_tax_status: str  # 'Correct' or 'Check'
    warning_type: Optional[str]
    col_validation_result: Optional[CoLValidationResult]


# ============================================================================
# CoL Validator with Caching
# ============================================================================

class CoLValidator:
    """
    CoL API client with multi-level caching and rate limiting.

    Cache hierarchy:
    1. In-memory dict (fast, limited size)
    2. File-based pickle cache (persistent)
    """

    def __init__(self, max_requests_per_second: int = 100,
                 timeout: int = 30, cache_enabled: bool = True,
                 cache_file: Optional[Path] = None, max_memory_cache: int = 100000,
                 dataset: str = "3LR"):
        # Shared CoL client + matcher (EXACT-search, synonym resolution, pinned DOI'd
        # release, thread-safe rate limiting).
        self._api = CoLAPIClient(max_requests_per_second=max_requests_per_second,
                                 timeout=timeout, dataset=dataset)
        self._matcher = TaxonomyMatcher(self._api)
        self.provenance = self._api.provenance
        self.base_url = "https://api.checklistbank.org/dataset/3LR/match/nameusage"
        self.session = requests.Session()

        # Session configuration
        adapter = HTTPAdapter(
            pool_connections=20,
            pool_maxsize=20,
            max_retries=3,
            pool_block=False
        )
        self.session.mount('http://', adapter)
        self.session.mount('https://', adapter)
        self.session.headers.update({'User-Agent': 'SATIVAIntegration/2.0'})

        # Rate limiting
        self.max_requests_per_second = max_requests_per_second
        self.request_interval = 1.0 / max_requests_per_second if max_requests_per_second > 0 else 0
        self.last_request_time = 0
        self.timeout = timeout

        # Cache configuration
        self.cache_enabled = cache_enabled
        self.memory_cache: Dict[str, CoLValidationResult] = {}
        self.max_memory_cache = max_memory_cache
        self.cache_file = cache_file or Path("cache/col_validation_cache.pkl")

        # Load file cache if exists
        if self.cache_enabled:
            self._load_file_cache()

        # Statistics
        self.stats = {
            'api_calls': 0,
            'cache_hits': 0,
            'failed_calls': 0
        }

    def _load_file_cache(self):
        """Load cache from pickle file"""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'rb') as f:
                    self.memory_cache = pickle.load(f)
                logger.info(f"Loaded {len(self.memory_cache)} entries from cache file")
            except Exception as e:
                logger.warning(f"Failed to load cache file: {e}")
                self.memory_cache = {}

    def save_cache(self):
        """Save cache to pickle file"""
        if not self.cache_enabled:
            return

        try:
            self.cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.cache_file, 'wb') as f:
                pickle.dump(self.memory_cache, f)
            logger.info(f"Saved {len(self.memory_cache)} entries to cache file")
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")

    @staticmethod
    def _taxonomy_to_cache_key(taxonomy: Dict[str, str]) -> str:
        """Convert taxonomy dict to cache key"""
        return ";".join([
            taxonomy.get('kingdom') or '',
            taxonomy.get('phylum') or '',
            taxonomy.get('class_') or '',
            taxonomy.get('order_') or '',
            taxonomy.get('family') or '',
            taxonomy.get('genus') or '',
            taxonomy.get('species') or ''
        ])

    def _get_from_cache(self, cache_key: str) -> Optional[CoLValidationResult]:
        """Get result from cache"""
        if not self.cache_enabled:
            return None

        if cache_key in self.memory_cache:
            self.stats['cache_hits'] += 1
            return self.memory_cache[cache_key]

        return None

    def _add_to_cache(self, cache_key: str, result: CoLValidationResult):
        """Add result to cache"""
        if not self.cache_enabled:
            return

        self.memory_cache[cache_key] = result

        # Evict oldest entries if cache is too large
        if len(self.memory_cache) > self.max_memory_cache:
            # Remove first 10% of entries
            to_remove = len(self.memory_cache) // 10
            for key in list(self.memory_cache.keys())[:to_remove]:
                del self.memory_cache[key]

    def _rate_limit(self):
        """Implement rate limiting for API calls"""
        if self.request_interval <= 0:
            return

        current_time = time.time()
        time_since_last = current_time - self.last_request_time

        if time_since_last < self.request_interval:
            time.sleep(self.request_interval - time_since_last)

        self.last_request_time = time.time()

    def _query_col_api(self, taxonomy: Dict[str, str]) -> Optional[Dict]:
        """Query CoL API with taxonomy"""
        try:
            self._rate_limit()

            # Build query parameters
            params = {"fuzzy": "true"}

            # Try from most specific to least specific
            taxonomy_levels = ['species', 'genus', 'family', 'order_', 'class_', 'phylum', 'kingdom']

            for level in taxonomy_levels:
                if taxonomy.get(level):
                    # Primary search term
                    params['q'] = taxonomy[level]

                    # Add context from higher levels
                    for context_level in taxonomy_levels[taxonomy_levels.index(level)+1:]:
                        if taxonomy.get(context_level):
                            param_name = context_level.replace('_', '')
                            params[param_name] = taxonomy[context_level]

                    break

            if 'q' not in params:
                return None

            response = self.session.get(self.base_url, params=params, timeout=self.timeout)
            self.stats['api_calls'] += 1

            if response.status_code != 200:
                logger.warning(f"CoL API error {response.status_code}")
                self.stats['failed_calls'] += 1
                return None

            data = response.json()

            if data.get("match") and "usage" in data:
                return data["usage"]

            return None

        except Exception as e:
            logger.warning(f"CoL API request failed: {e}")
            self.stats['failed_calls'] += 1
            return None

    @staticmethod
    def _extract_classification_from_usage(usage: Dict) -> Dict[str, str]:
        """Extract taxonomy classification from CoL usage response"""
        classification = {}

        # Extract from classification array
        if "classification" in usage:
            for taxon in usage["classification"]:
                if taxon.get('rank') and taxon.get('name'):
                    rank = taxon['rank'].lower()
                    if rank == 'class':
                        rank = 'class_'
                    elif rank == 'order':
                        rank = 'order_'
                    classification[rank] = taxon['name']

        # Add the main taxon
        if usage.get("rank") and usage.get("name"):
            rank = usage["rank"].lower()
            if rank == 'class':
                rank = 'class_'
            elif rank == 'order':
                rank = 'order_'
            classification[rank] = usage["name"]

        return classification

    @staticmethod
    def _taxonomies_match(tax1: Dict[str, str], tax2: Dict[str, str]) -> bool:
        """Check if two taxonomies are identical"""
        levels = ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species']

        for level in levels:
            val1 = tax1.get(level, '').strip().lower() if tax1.get(level) else ''
            val2 = tax2.get(level, '').strip().lower() if tax2.get(level) else ''

            if val1 != val2:
                return False

        return True

    def validate_taxonomy(self, taxonomy: Dict[str, str]) -> CoLValidationResult:
        """
        Validate a taxonomy against CoL.

        Args:
            taxonomy: Dict with keys kingdom, phylum, class_, order_, family, genus, species

        Returns:
            CoLValidationResult with match quality and returned taxonomy
        """
        # Check cache first
        cache_key = self._taxonomy_to_cache_key(taxonomy)
        cached = self._get_from_cache(cache_key)
        if cached:
            return cached

        # Query CoL via the shared corrected matcher (EXACT search + synonym resolution)
        rec = TaxonomyRecord(
            tax_id=0,
            kingdom=taxonomy.get('kingdom'), phylum=taxonomy.get('phylum'),
            class_=taxonomy.get('class_'), order_=taxonomy.get('order_'),
            family=taxonomy.get('family'), genus=taxonomy.get('genus'),
            species=taxonomy.get('species'))
        match = self._matcher.match_taxonomy(rec)
        self.stats['api_calls'] = self._api.api_calls
        self.stats['failed_calls'] = self._api.failed_calls

        if not match.success or not match.col_taxonomy:
            # No match found
            result = CoLValidationResult(
                success=False,
                match_quality=CoLMatchQuality.NO_MATCH,
                col_taxonomy=None
            )
        else:
            col_taxonomy = match.col_taxonomy

            # Compare with input taxonomy
            if self._taxonomies_match(taxonomy, col_taxonomy):
                # Exact match
                result = CoLValidationResult(
                    success=True,
                    match_quality=CoLMatchQuality.EXACT,
                    col_taxonomy=col_taxonomy
                )
            else:
                # Modified match
                result = CoLValidationResult(
                    success=True,
                    match_quality=CoLMatchQuality.MODIFIED,
                    col_taxonomy=col_taxonomy
                )

        # Add to cache
        self._add_to_cache(cache_key, result)

        return result

    def close(self):
        """Close sessions and save cache"""
        self.save_cache()
        self.session.close()
        self._api.close()


# ============================================================================
# SATIVA .mis/.tax Parser
# ============================================================================

class SativaMisParser:
    """Parse SATIVA .mis and .tax files"""

    MISSING_TAX_TOKENS = {"", "na", "n/a", "none", "null", "unknown", "unidentified"}

    @staticmethod
    def _normalize_tax_token(value: str | None) -> str | None:
        if value is None:
            return None
        cleaned = value.strip()
        if not cleaned:
            return None
        if cleaned.lower() in SativaMisParser.MISSING_TAX_TOKENS:
            return None
        return cleaned

    @staticmethod
    def parse_taxonomy_path_to_dict(tax_path: str) -> Dict[str, str]:
        """
        Convert semicolon-separated taxonomy path to dict.

        Args:
            tax_path: "Kingdom;Phylum;Class;Order;Family;Genus;Species"

        Returns:
            Dict with keys: kingdom, phylum, class_, order_, family, genus, species
        """
        levels = tax_path.split(';')

        # Normalize and remove empty/NA tokens
        levels = [SativaMisParser._normalize_tax_token(level) for level in levels]

        # Pad with None if fewer than 7 levels
        while len(levels) < 7:
            levels.append(None)

        return {
            'kingdom': levels[0],
            'phylum': levels[1],
            'class_': levels[2],
            'order_': levels[3],
            'family': levels[4],
            'genus': levels[5],
            'species': levels[6]
        }

    @staticmethod
    def parse_mis_file(mis_path: Path, barcode: str, tax_level: str) -> List[ModificationRecord]:
        """
        Parse a .mis file into ModificationRecords.

        Args:
            mis_path: Path to .mis file
            barcode: Barcode name (ITS, ITS1, ITS2)
            tax_level: Taxonomic level (family, genus, etc.)

        Returns:
            List of ModificationRecord objects
        """
        records = []

        try:
            with open(mis_path, 'r') as f:
                for line in f:
                    # Skip comment lines and empty lines
                    if line.startswith(';') or not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 8:
                        logger.warning(f"Invalid line in {mis_path.name}: {line[:50]}")
                        continue

                    try:
                        record = ModificationRecord(
                            seq_id=parts[0],
                            seq_original_id=parts[0],  # Will be updated from DB
                            mislabeled_level=parts[1],
                            original_label=parts[2],
                            proposed_label=parts[3],
                            confidence=float(parts[4]),
                            original_taxonomy=parts[5],
                            proposed_taxonomy=parts[6],
                            per_rank_confidence=parts[7],
                            barcode=barcode,
                            tax_level=tax_level
                        )
                        records.append(record)
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Failed to parse line in {mis_path.name}: {e}")
                        continue

            logger.debug(f"Parsed {len(records)} modifications from {mis_path.name}")

        except Exception as e:
            logger.error(f"Failed to read {mis_path}: {e}")

        return records

    @staticmethod
    def parse_tax_file(tax_path: Path) -> Set[str]:
        """
        Parse a .tax file to get all seq_ids.

        Args:
            tax_path: Path to .tax file

        Returns:
            Set of seq_ids present in the .tax file
        """
        seq_ids = set()

        try:
            with open(tax_path, 'r') as f:
                for line in f:
                    if not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) >= 1:
                        seq_ids.add(parts[0])

            logger.debug(f"Found {len(seq_ids)} sequences in {tax_path.name}")

        except Exception as e:
            logger.error(f"Failed to read {tax_path}: {e}")

        return seq_ids

    @staticmethod
    def identify_validated_sequences(tax_path: Path, modifications: List[ModificationRecord]) -> Set[str]:
        """
        Identify sequences that are in .tax but NOT in .mis (VALIDATED).

        Args:
            tax_path: Path to .tax file
            modifications: List of modifications from .mis file

        Returns:
            Set of seq_ids that are validated (in .tax but not in .mis)
        """
        tax_seq_ids = SativaMisParser.parse_tax_file(tax_path)
        mis_seq_ids = {mod.seq_id for mod in modifications}

        validated = tax_seq_ids - mis_seq_ids
        logger.debug(f"Identified {len(validated)} validated sequences")

        return validated


# ============================================================================
# Modification Classifier
# ============================================================================

class ModificationClassifier:
    """Classify SATIVA modifications into types"""

    @staticmethod
    def classify_modification(original_path: str, proposed_path: str) -> ModificationType:
        """
        Classify a SATIVA modification.

        Args:
            original_path: Original taxonomy path (semicolon-separated)
            proposed_path: Proposed taxonomy path (semicolon-separated)

        Returns:
            ModificationType enum
        """
        # Parse into lists, treating NA/empty as missing
        original_levels = [
            token
            for token in (SativaMisParser._normalize_tax_token(x) for x in original_path.split(';'))
            if token
        ]
        proposed_levels = [
            token
            for token in (SativaMisParser._normalize_tax_token(x) for x in proposed_path.split(';'))
            if token
        ]

        # Find common length and check if common levels are identical
        min_length = min(len(original_levels), len(proposed_levels))
        common_original = original_levels[:min_length]
        common_proposed = proposed_levels[:min_length]

        if common_original != common_proposed:
            # At least one common level changed -> SUBSTITUTION
            # This includes cases like: Fungi -> Plantae;Tracheophyta;...
            return ModificationType.SUBSTITUTION

        # Common levels are identical, check length difference
        if len(proposed_levels) > len(original_levels):
            # More levels � COMPLETION
            return ModificationType.COMPLETION

        elif len(proposed_levels) < len(original_levels):
            # Fewer levels and common levels identical -> Pure DELETION
            deletion_count = len(original_levels) - len(proposed_levels)
            if deletion_count == 1:
                return ModificationType.DELETION_SIMPLE
            else:
                return ModificationType.DELETION_MULTIPLE

        else:
            # Same length and common levels identical (already checked above)
            # This shouldn't happen in .mis file
            logger.warning(f"Identical taxonomies found in .mis: {original_path}")
            return ModificationType.VALIDATED

    @staticmethod
    def count_deleted_levels(original_path: str, proposed_path: str) -> List[str]:
        """
        Identify taxonomy levels that were deleted between the original and proposed paths.
        Returns the list of level names (kingdom → species) that disappeared.
        """
        orig = SativaMisParser.parse_taxonomy_path_to_dict(original_path)
        prop = SativaMisParser.parse_taxonomy_path_to_dict(proposed_path)

        deleted = []
        for display_name, key in zip(DISPLAY_TAXONOMY_LEVELS, DB_TAXONOMY_LEVEL_KEYS):
            orig_value = orig.get(key)
            prop_value = prop.get(key)

            if orig_value and not prop_value:
                deleted.append(display_name)

        return deleted


# ============================================================================
# Database Integrator
# ============================================================================

class DatabaseIntegrator:
    """Handle database operations for SATIVA integration"""

    def __init__(self, db_manager, col_validator: CoLValidator):
        self.db_manager = db_manager
        self.col_validator = col_validator
        self.classifier = ModificationClassifier()

    def _get_table(self, table_name: str) -> str:
        """Get schema-qualified table name"""
        if self.db_manager.db_type == "postgresql":
            schema = self.db_manager.connection_params.get('schema', 'public')
            return f"{schema}.{table_name}"
        return table_name

    def _run_with_retry(self, operation_name: str, func):
        """Execute a DB operation with retry."""
        return self.db_manager._execute_with_retry(operation_name, func)

    def get_original_tax_id_for_sequence(self, seq_id: str) -> Optional[str]:
        """
        Get the original (Raw) taxonomy ID for a sequence.

        Args:
            seq_id: Sequence UUID

        Returns:
            Original tax_id or None
        """
        from sqlalchemy import text

        seq_tax_table = self._get_table('seq_tax')

        def fetch_tax_id():
            with self.db_manager.get_session() as session:
                # Priority: Raw taxonomy
                result = session.execute(text(f"""
                    SELECT tax_id FROM {seq_tax_table}
                    WHERE seq_id = :seq_id AND backbone = 'Raw'
                    LIMIT 1
                """), {"seq_id": seq_id}).fetchone()

                if result:
                    return result[0]

                # Fallback: any taxonomy
                result = session.execute(text(f"""
                    SELECT tax_id FROM {seq_tax_table}
                    WHERE seq_id = :seq_id
                    LIMIT 1
                """), {"seq_id": seq_id}).fetchone()

                return result[0] if result else None

        return self.db_manager._execute_with_retry("sativa_get_original_tax_id", fetch_tax_id)

    def get_col_taxonomy_for_sequence(self, seq_id: str) -> Optional[Dict[str, str]]:
        """
        Get the CoL taxonomy for a sequence if it exists.

        Args:
            seq_id: Sequence UUID

        Returns:
            Taxonomy dict if CoL taxonomy exists, None otherwise
        """
        from sqlalchemy import text

        seq_tax_table = self._get_table('seq_tax')
        taxonomy_table = self._get_table('taxonomy')

        def fetch_col_taxonomy():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT t.kingdom, t.phylum, t.class_, t.order_, t.family, t.genus, t.species
                    FROM {seq_tax_table} st
                    JOIN {taxonomy_table} t ON st.tax_id = t.tax_id
                    WHERE st.seq_id = :seq_id AND st.backbone = 'CoL'
                    LIMIT 1
                """), {"seq_id": seq_id}).fetchone()

                if result:
                    return {
                        'kingdom': result[0] or '',
                        'phylum': result[1] or '',
                        'class_': result[2] or '',
                        'order_': result[3] or '',
                        'family': result[4] or '',
                        'genus': result[5] or '',
                        'species': result[6] or ''
                    }

                return None

        return self.db_manager._execute_with_retry("sativa_get_col_taxonomy", fetch_col_taxonomy)

    def check_taxonomy_exists(self, taxonomy: Dict[str, str]) -> Optional[str]:
        """
        Check if a taxonomy already exists in the database.

        Args:
            taxonomy: Dict with kingdom, phylum, class_, order_, family, genus, species

        Returns:
            tax_id if exists, None otherwise
        """
        from sqlalchemy import text

        taxonomy_table = self._get_table('taxonomy')

        def fetch_existing_taxonomy():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT tax_id FROM {taxonomy_table}
                    WHERE
                        COALESCE(kingdom, '') = :kingdom AND
                        COALESCE(phylum, '') = :phylum AND
                        COALESCE(class_, '') = :class_ AND
                        COALESCE(order_, '') = :order_ AND
                        COALESCE(family, '') = :family AND
                        COALESCE(genus, '') = :genus AND
                        COALESCE(species, '') = :species
                    LIMIT 1
                """), {
                    'kingdom': taxonomy.get('kingdom') or '',
                    'phylum': taxonomy.get('phylum') or '',
                    'class_': taxonomy.get('class_') or '',
                    'order_': taxonomy.get('order_') or '',
                    'family': taxonomy.get('family') or '',
                    'genus': taxonomy.get('genus') or '',
                    'species': taxonomy.get('species') or ''
                }).fetchone()

                return result[0] if result else None

        return self.db_manager._execute_with_retry("sativa_check_taxonomy", fetch_existing_taxonomy)

    def get_taxonomy_by_id(self, tax_id: str) -> Optional[Dict[str, str]]:
        """Retrieve a taxonomy dict by its ID."""
        from sqlalchemy import text

        taxonomy_table = self._get_table('taxonomy')

        def fetch_taxonomy():
            with self.db_manager.get_session() as session:
                row = session.execute(text(f"""
                    SELECT kingdom, phylum, class_, order_, family, genus, species
                    FROM {taxonomy_table}
                    WHERE tax_id = :tax_id
                    LIMIT 1
                """), {"tax_id": tax_id}).fetchone()

                if not row:
                    return None

                return {
                    'kingdom': row[0],
                    'phylum': row[1],
                    'class_': row[2],
                    'order_': row[3],
                    'family': row[4],
                    'genus': row[5],
                    'species': row[6],
                }

        return self.db_manager._execute_with_retry("sativa_get_taxonomy_by_id", fetch_taxonomy)

    def create_taxonomy(self, taxonomy: Dict[str, str], parent_tax_id: Optional[str]) -> str:
        """
        Create a new taxonomy entry.

        Args:
            taxonomy: Dict with taxonomy levels
            parent_tax_id: Parent taxonomy ID (for traceability)

        Returns:
            New tax_id
        """
        # Convert None values to empty strings for UniqueConstraint compatibility
        return self.db_manager.add_taxonomy(
            kingdom=taxonomy.get('kingdom') or '',
            phylum=taxonomy.get('phylum') or '',
            class_=taxonomy.get('class_') or '',
            order_=taxonomy.get('order_') or '',
            family=taxonomy.get('family') or '',
            genus=taxonomy.get('genus') or '',
            species=taxonomy.get('species') or '',
            parent_tax_id=parent_tax_id
        )

    def get_or_create_taxonomy(self, taxonomy: Dict[str, str], parent_tax_id: Optional[str]) -> str:
        """Get existing taxonomy or create new one"""
        existing_tax_id = self.check_taxonomy_exists(taxonomy)

        if existing_tax_id:
            return existing_tax_id

        return self.create_taxonomy(taxonomy, parent_tax_id)

    def check_seq_tax_exists(self, seq_id: str, tax_id: str, backbone: str) -> bool:
        """Check if a seq_tax link already exists"""
        from sqlalchemy import text

        seq_tax_table = self._get_table('seq_tax')

        def fetch_exists():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    SELECT 1 FROM {seq_tax_table}
                    WHERE seq_id = :seq_id AND tax_id = :tax_id AND backbone = :backbone
                """), {"seq_id": seq_id, "tax_id": tax_id, "backbone": backbone}).fetchone()

                return result is not None

        return self.db_manager._execute_with_retry("sativa_check_seq_tax_exists", fetch_exists)

    # ------------------------------------------------------------------ #
    # BULK helpers (avoid the per-sequence N+1: one chunked query instead #
    # of one query per sequence). Backend-aware: `= ANY(:ids)` on         #
    # PostgreSQL/YugabyteDB (single param, no variable limit), `IN (...)`  #
    # on SQLite. Chunked to keep queries under YugabyteDB RPC timeouts.   #
    # ------------------------------------------------------------------ #
    def _id_chunk_size(self) -> int:
        return 900 if self.db_manager.db_type == "sqlite" else 2000

    def _in_clause(self, column: str, chunk: List[str], prefix: str):
        """Return (sql_fragment, params) for `column IN <chunk>`, backend-aware."""
        if self.db_manager.db_type == "postgresql":
            return f"{column} = ANY(:{prefix}ids)", {f"{prefix}ids": list(chunk)}
        placeholders = ", ".join(f":{prefix}{i}" for i in range(len(chunk)))
        return f"{column} IN ({placeholders})", {f"{prefix}{i}": v for i, v in enumerate(chunk)}

    def _chunked_rows(self, seq_ids, sql_template: str, column: str, op_name: str,
                      extra_params: Optional[Dict] = None):
        """Run `sql_template` (with a `{clause}` placeholder) over chunks of seq_ids."""
        from sqlalchemy import text
        ids = list(seq_ids)
        cs = self._id_chunk_size()
        for i in range(0, len(ids), cs):
            chunk = ids[i:i + cs]
            clause, params = self._in_clause(column, chunk, "s")
            if extra_params:
                params = {**params, **extra_params}
            sql = sql_template.format(clause=clause)

            def run(sql=sql, params=params):
                with self.db_manager.get_session() as session:
                    return list(session.execute(text(sql), params))

            for row in self.db_manager._execute_with_retry(f"{op_name}_{i}", run):
                yield row

    def bulk_get_original_tax_ids(self, seq_ids) -> Dict[str, str]:
        """seq_id -> original tax_id (prefer backbone 'Raw', else any)."""
        st = self._get_table('seq_tax')
        raw, other = {}, {}
        sql = f"SELECT seq_id, tax_id, backbone FROM {st} WHERE {{clause}}"
        for seq_id, tax_id, backbone in self._chunked_rows(seq_ids, sql, "seq_id", "sativa_bulk_orig_tax"):
            if backbone == 'Raw':
                raw[seq_id] = tax_id
            else:
                other.setdefault(seq_id, tax_id)
        for sid, tid in other.items():
            raw.setdefault(sid, tid)
        return raw

    def bulk_get_col_taxonomy(self, seq_ids) -> Dict[str, Tuple[str, Dict[str, str]]]:
        """seq_id -> (col_tax_id, taxonomy_dict) for sequences with a 'CoL' backbone."""
        st = self._get_table('seq_tax')
        tt = self._get_table('taxonomy')
        sql = (f"SELECT st.seq_id, st.tax_id, t.kingdom, t.phylum, t.class_, t.order_, "
               f"t.family, t.genus, t.species FROM {st} st JOIN {tt} t ON st.tax_id = t.tax_id "
               f"WHERE {{clause}} AND st.backbone = 'CoL'")
        out = {}
        for r in self._chunked_rows(seq_ids, sql, "st.seq_id", "sativa_bulk_col_tax"):
            out[r[0]] = (r[1], {
                'kingdom': r[2] or '', 'phylum': r[3] or '', 'class_': r[4] or '',
                'order_': r[5] or '', 'family': r[6] or '', 'genus': r[7] or '', 'species': r[8] or ''
            })
        return out

    def bulk_get_existing_seq_tax(self, seq_ids, backbone: str) -> Set[Tuple[str, str]]:
        """Set of (seq_id, tax_id) already linked under the given backbone."""
        st = self._get_table('seq_tax')
        sql = f"SELECT seq_id, tax_id FROM {st} WHERE {{clause}} AND backbone = :bb"
        return {(r[0], r[1]) for r in self._chunked_rows(
            seq_ids, sql, "seq_id", "sativa_bulk_existing", extra_params={"bb": backbone})}

    def decide_integration_action(self, modification: ModificationRecord,
                                   mod_type: ModificationType,
                                   existing_col_taxonomy: Optional[Dict[str, str]]) -> IntegrationAction:
        """
        Decide what actions to take for a modification.

        New logic:
        - For each SATIVA suggestion, validate the proposed taxonomy via CoL.
        - Always create a Sativa seq_tax entry using the SATIVA proposed taxonomy.
        - Always create a warning entry to record CoL outcome (COL_MATCHED,
          COL_MODIFIED_SATIVA, or COL_NO_MATCH).

        Args:
            modification: The modification record
            mod_type: The classified modification type
            existing_col_taxonomy: Existing CoL taxonomy for this sequence (if any)

        Returns:
            IntegrationAction with all necessary information
        """
        proposed_tax_dict = SativaMisParser.parse_taxonomy_path_to_dict(
            modification.proposed_taxonomy
        )
        col_result = self.col_validator.validate_taxonomy(proposed_tax_dict)

        warning_type = 'COL_NO_MATCH'

        if col_result.success and col_result.match_quality == CoLMatchQuality.EXACT:
            warning_type = 'COL_MATCHED'
        elif col_result.match_quality == CoLMatchQuality.MODIFIED:
            warning_type = 'COL_MODIFIED_SATIVA'

        seq_tax_status = 'Correct' if (
            col_result.success and col_result.match_quality == CoLMatchQuality.EXACT
        ) else 'Check'

        taxonomy_to_use = proposed_tax_dict
        if col_result.success and col_result.col_taxonomy:
            taxonomy_to_use = self._merge_col_with_sativa(
                proposed_tax_dict,
                col_result.col_taxonomy
            )

        return IntegrationAction(
            create_seq_tax=True,
            update_seq_tax_backbone=False,
            create_warning=True,
            taxonomy_to_use=taxonomy_to_use,
            seq_tax_status=seq_tax_status,
            warning_type=warning_type,
            col_validation_result=col_result
        )

    @staticmethod
    def _merge_col_with_sativa(sativa_tax: Dict[str, str],
                               col_tax: Dict[str, str]) -> Dict[str, str]:
        """
        Use CoL taxonomy where available, but keep SATIVA values for deletions.

        A deletion is when CoL returns a missing/empty rank while SATIVA had one.
        """
        merged = {}
        for key in ('kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species'):
            col_value = col_tax.get(key)
            sativa_value = sativa_tax.get(key)
            merged[key] = col_value if col_value else sativa_value
        return merged

    def prepare_seq_tax_entry(self, seq_id: str, tax_id: str, status: str,
                             confidence: float) -> Dict:
        """Prepare a seq_tax entry dict for bulk insert"""
        return {
            'seq_id': seq_id,
            'tax_id': tax_id,
            'backbone': 'Sativa',
            'status': status,
            'confidence': confidence
        }

    def update_seq_tax_backbone(self, seq_id: str, confidence: float) -> bool:
        """
        Update existing seq_tax backbone to 'Sativa' for validated sequences.
        Priority: CoL > Raw (updates the most recent/best taxonomy).

        Args:
            seq_id: Sequence ID
            confidence: SATIVA confidence score

        Returns:
            bool: True if updated, False if no seq_tax found
        """
        from sqlalchemy import text

        seq_tax_table = self._get_table('seq_tax')

        # Define backbone priority: CoL preferred over Raw
        priority_order = ['CoL', 'Raw']

        def update_op():
            with self.db_manager.get_session() as session:
                # Fetch the backbones already linked to this sequence
                all_backbones = session.execute(text(f"""
                    SELECT backbone FROM {seq_tax_table}
                    WHERE seq_id = :seq_id
                """), {"seq_id": seq_id}).fetchall()

                if not all_backbones:
                    logger.warning(f"No seq_tax entries at all for {seq_id}")
                    return False

                existing_backbones = [row[0] for row in all_backbones]
                logger.debug(f"Seq {seq_id} has backbones: {existing_backbones}")

                # Find existing seq_tax with priority: CoL > Raw
                for backbone in priority_order:
                    if backbone not in existing_backbones:
                        continue

                    # Update this backbone to 'Sativa'
                    result = session.execute(text(f"""
                        UPDATE {seq_tax_table}
                        SET backbone = 'Sativa',
                            status = 'Correct',
                            confidence = :confidence
                        WHERE seq_id = :seq_id AND backbone = :backbone
                    """), {"seq_id": seq_id, "backbone": backbone, "confidence": confidence})

                    if result.rowcount > 0:
                        logger.debug(f"Updated seq_tax for {seq_id}: {backbone} → Sativa")
                        return True

                logger.warning(f"No CoL or Raw seq_tax found for {seq_id} (has: {existing_backbones})")
                return False

        return self.db_manager._execute_with_retry("sativa_update_seq_tax_backbone", update_op)

    def prepare_warning_entry(self, modification: ModificationRecord,
                             warning_type: str, original_tax_id: str,
                             col_result: Optional[CoLValidationResult]) -> Dict:
        """Prepare a warning entry dict for bulk insert"""
        # Map ModificationType to database-allowed values
        mod_type = self.classifier.classify_modification(
            modification.original_taxonomy,
            modification.proposed_taxonomy
        )
        db_mod_type = mod_type.value

        entry = {
            'seq_id': modification.seq_id,
            'seq_original_id': modification.seq_original_id,
            'original_tax_id': original_tax_id,
            'warning_type': warning_type,
            'sativa_modification_type': db_mod_type,
            'sativa_confidence': modification.confidence,
            'sativa_original_taxonomy': modification.original_taxonomy,
            'sativa_proposed_taxonomy': modification.proposed_taxonomy,
            'barcode': modification.barcode,
            'tax_level': modification.tax_level,
            'mislabeled_level': modification.mislabeled_level,
            'per_rank_confidence': modification.per_rank_confidence,
            'col_validation_attempted': col_result is not None
        }

        # Add deleted levels for deletions
        if 'DELETION' in warning_type:
            deleted_levels = self.classifier.count_deleted_levels(
                modification.original_taxonomy,
                modification.proposed_taxonomy
            )
            entry['deleted_levels'] = json.dumps(deleted_levels) if deleted_levels else None
        else:
            entry['deleted_levels'] = None

        # Add CoL validation info
        if col_result:
            if col_result.col_taxonomy:
                col_tax_path = ";".join([
                    col_result.col_taxonomy.get('kingdom') or '',
                    col_result.col_taxonomy.get('phylum') or '',
                    col_result.col_taxonomy.get('class_') or '',
                    col_result.col_taxonomy.get('order_') or '',
                    col_result.col_taxonomy.get('family') or '',
                    col_result.col_taxonomy.get('genus') or '',
                    col_result.col_taxonomy.get('species') or ''
                ])
                entry['col_returned_taxonomy'] = col_tax_path
            else:
                entry['col_returned_taxonomy'] = None

            entry['col_match_quality'] = col_result.match_quality.value if col_result.match_quality else None
        else:
            entry['col_returned_taxonomy'] = None
            entry['col_match_quality'] = None

        return entry

    def bulk_insert_seq_tax(self, entries: List[Dict], batch_size: int = 1000):
        """Bulk insert seq_tax entries in batches.

        Args:
            entries: List of entry dicts to insert
            batch_size: Number of entries per batch (default 1000)
        """
        if not entries:
            return

        from sqlalchemy import text

        seq_tax_table = self._get_table('seq_tax')

        # Process in batches
        for i in range(0, len(entries), batch_size):
            batch = entries[i:i + batch_size]

            def insert_batch(batch_entries=batch):
                with self.db_manager.get_session() as session:
                    if self.db_manager.db_type == "postgresql":
                        # Use PostgreSQL-specific INSERT ... ON CONFLICT with executemany
                        stmt = text(f"""
                            INSERT INTO {seq_tax_table} (seq_id, tax_id, backbone, status, confidence)
                            VALUES (:seq_id, :tax_id, :backbone, :status, :confidence)
                            ON CONFLICT (seq_id, tax_id, backbone) DO NOTHING
                        """)
                        for entry in batch_entries:
                            session.execute(stmt, entry)
                    else:
                        # SQLite: Use INSERT OR IGNORE
                        stmt = text(f"""
                            INSERT OR IGNORE INTO {seq_tax_table} (seq_id, tax_id, backbone, status, confidence)
                            VALUES (:seq_id, :tax_id, :backbone, :status, :confidence)
                        """)
                        for entry in batch_entries:
                            session.execute(stmt, entry)

            self.db_manager._execute_with_retry(f"sativa_bulk_insert_seq_tax_batch_{i}", insert_batch)

    def bulk_insert_warnings(self, entries: List[Dict], batch_size: int = 1000):
        """Bulk insert sativa_warnings entries in batches.

        Args:
            entries: List of entry dicts to insert
            batch_size: Number of entries per batch (default 1000)
        """
        if not entries:
            return

        from sqlalchemy import text
        import uuid

        warnings_table = self._get_table('sativa_warnings')

        # Add warning_id to all entries upfront
        for entry in entries:
            if 'warning_id' not in entry:
                entry['warning_id'] = str(uuid.uuid4())

        # Process in batches
        for i in range(0, len(entries), batch_size):
            batch = entries[i:i + batch_size]

            def insert_batch(batch_entries=batch):
                with self.db_manager.get_session() as session:
                    stmt = text(f"""
                        INSERT INTO {warnings_table} (
                            warning_id, seq_id, seq_original_id, original_tax_id,
                            warning_type, sativa_modification_type, sativa_confidence,
                            sativa_original_taxonomy, sativa_proposed_taxonomy,
                            deleted_levels, col_validation_attempted, col_returned_taxonomy,
                            col_match_quality, barcode, tax_level, mislabeled_level,
                            per_rank_confidence
                        ) VALUES (
                            :warning_id, :seq_id, :seq_original_id, :original_tax_id,
                            :warning_type, :sativa_modification_type, :sativa_confidence,
                            :sativa_original_taxonomy, :sativa_proposed_taxonomy,
                            :deleted_levels, :col_validation_attempted, :col_returned_taxonomy,
                            :col_match_quality, :barcode, :tax_level, :mislabeled_level,
                            :per_rank_confidence
                        )
                    """)
                    for entry in batch_entries:
                        session.execute(stmt, entry)

            self.db_manager._execute_with_retry(f"sativa_bulk_insert_warnings_batch_{i}", insert_batch)


# ============================================================================
# Main Processing Logic
# ============================================================================

def create_sativa_warnings_table(db_manager):
    """Create sativa_warnings table if it doesn't exist.

    Table creation and index creation are done in separate transactions
    to avoid YugabyteDB timeout issues with CREATE INDEX.
    """
    from sqlalchemy import text

    integrator = DatabaseIntegrator(db_manager, None)
    warnings_table = integrator._get_table('sativa_warnings')

    # Step 1: Create table
    def create_table_op():
        with db_manager.get_session() as session:
            if db_manager.db_type == "postgresql":
                session.execute(text(f"""
                    CREATE TABLE IF NOT EXISTS {warnings_table} (
                        warning_id UUID PRIMARY KEY,
                        seq_id UUID NOT NULL,
                        seq_original_id TEXT NOT NULL,
                        original_tax_id UUID NOT NULL,
                    warning_type TEXT NOT NULL CHECK (warning_type IN (
                        'DELETION_SIMPLE', 'DELETION_MULTIPLE',
                        'COL_MATCHED', 'COL_MODIFIED_SATIVA', 'COL_NO_MATCH'
                    )),
                    sativa_modification_type TEXT CHECK (sativa_modification_type IN (
                        'DELETION_SIMPLE', 'DELETION_MULTIPLE', 'SUBSTITUTION', 'COMPLETION'
                    )),
                        sativa_confidence FLOAT,
                        sativa_original_taxonomy TEXT,
                        sativa_proposed_taxonomy TEXT,
                        deleted_levels TEXT,
                        col_validation_attempted BOOLEAN DEFAULT FALSE,
                        col_returned_taxonomy TEXT,
                        col_match_quality TEXT,
                        barcode TEXT NOT NULL,
                        tax_level TEXT NOT NULL,
                        mislabeled_level TEXT,
                        per_rank_confidence TEXT,
                        created_at TIMESTAMP DEFAULT NOW()
                    )
                """))
            else:  # SQLite
                session.execute(text(f"""
                    CREATE TABLE IF NOT EXISTS {warnings_table} (
                        warning_id TEXT PRIMARY KEY,
                        seq_id TEXT NOT NULL,
                        seq_original_id TEXT NOT NULL,
                        original_tax_id TEXT NOT NULL,
                    warning_type TEXT NOT NULL CHECK (warning_type IN (
                        'DELETION_SIMPLE', 'DELETION_MULTIPLE',
                        'COL_MATCHED', 'COL_MODIFIED_SATIVA', 'COL_NO_MATCH'
                    )),
                    sativa_modification_type TEXT CHECK (sativa_modification_type IN (
                        'DELETION_SIMPLE', 'DELETION_MULTIPLE', 'SUBSTITUTION', 'COMPLETION'
                    )),
                        sativa_confidence REAL,
                        sativa_original_taxonomy TEXT,
                        sativa_proposed_taxonomy TEXT,
                        deleted_levels TEXT,
                        col_validation_attempted INTEGER DEFAULT 0,
                        col_returned_taxonomy TEXT,
                        col_match_quality TEXT,
                        barcode TEXT NOT NULL,
                        tax_level TEXT NOT NULL,
                        mislabeled_level TEXT,
                        per_rank_confidence TEXT,
                        created_at TEXT DEFAULT CURRENT_TIMESTAMP
                    )
                """))

    db_manager._execute_with_retry("sativa_create_warnings_table", create_table_op)

    # Step 2: Create indexes separately (each in its own transaction with retry)
    index_statements = [
        ("idx_sativa_warnings_seq", f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_seq ON {warnings_table}(seq_id)"),
        ("idx_sativa_warnings_type", f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_type ON {warnings_table}(warning_type)"),
        ("idx_sativa_warnings_barcode", f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_barcode ON {warnings_table}(barcode)"),
    ]

    for idx_name, idx_sql in index_statements:
        def create_index_op(sql=idx_sql):
            with db_manager.get_session() as session:
                session.execute(text(sql))

        try:
            db_manager._execute_with_retry(f"sativa_create_index_{idx_name}", create_index_op)
        except Exception as e:
            # Index creation failure should not stop the pipeline
            logger.warning(f"Could not create index {idx_name}: {e}")

    logger.info("sativa_warnings table created/verified")


def process_single_mis_file(mis_file: Path, barcode: str, tax_level: str,
                            col_validator: CoLValidator, db_manager,
                            config: Dict) -> Dict:
    """
    Process a single .mis file.

    Args:
        mis_file: Path to .mis file
        barcode: Barcode name
        tax_level: Taxonomic level
        col_validator: CoL validator instance
        db_manager: Database manager
        config: Configuration dict

    Returns:
        Statistics dict
    """
    stats = defaultdict(int)
    parser = SativaMisParser()
    integrator = DatabaseIntegrator(db_manager, col_validator)

    # Parse .mis file
    modifications = parser.parse_mis_file(mis_file, barcode, tax_level)
    stats['total_modifications'] = len(modifications)

    tax_dir = mis_file.parent
    candidate_by_stem = tax_dir / f"{mis_file.stem}.tax"
    candidate_by_taxon = tax_dir / f"{tax_dir.name}.tax"
    tax_file = None
    if candidate_by_stem.exists():
        tax_file = candidate_by_stem
    elif candidate_by_taxon.exists():
        tax_file = candidate_by_taxon

    if not modifications and not tax_file:
        logger.warning(f"No modifications and no .tax file for {mis_file.name}")
        return dict(stats)

    # Initialize batches
    seq_tax_batch = []
    warnings_batch = []

    # Parse .tax file and identify validated sequences
    validated_seq_ids = []
    if tax_file and tax_file.exists():
        validated_seq_ids = parser.identify_validated_sequences(tax_file, modifications)
        stats['validated_sequences'] = len(validated_seq_ids)
        stats['tax_file_sequences'] = len(parser.parse_tax_file(tax_file))
        logger.info(
            "[%s] %d sequences in .tax, %d modifications in .mis, %d classified as VALIDATED",
            mis_file.name,
            stats['tax_file_sequences'],
            len(modifications),
            len(validated_seq_ids)
        )

    # BULK PRE-FETCH — resolve everything the per-sequence loops need in a few chunked,
    # backend-aware queries instead of ~4 queries per sequence. This is the key win for
    # YugabyteDB (each query = a network round-trip + distributed transaction).
    all_seq_ids = set(validated_seq_ids) | {mod.seq_id for mod in modifications}
    orig_map = integrator.bulk_get_original_tax_ids(all_seq_ids)                 # seq_id -> raw tax_id
    col_map = integrator.bulk_get_col_taxonomy(all_seq_ids)                      # seq_id -> (tax_id, dict)
    existing_sativa = integrator.bulk_get_existing_seq_tax(all_seq_ids, 'Sativa')  # {(seq_id, tax_id)}

    def _flush_seq_tax():
        nonlocal seq_tax_batch
        if len(seq_tax_batch) >= config['performance']['bulk_insert_size']:
            integrator.bulk_insert_seq_tax(seq_tax_batch)
            stats['seq_tax_created'] += len(seq_tax_batch)
            seq_tax_batch = []

    # Process validated sequences: reuse the existing CoL (else original) tax_id — no
    # per-sequence DB, and no taxonomy creation (the CoL/original taxonomy already exists).
    for seq_id in validated_seq_ids:
        try:
            col_entry = col_map.get(seq_id)
            if col_entry:
                tax_id = col_entry[0]
                if (seq_id, tax_id) not in existing_sativa:
                    seq_tax_batch.append(integrator.prepare_seq_tax_entry(seq_id, tax_id, 'Correct', 1.0))
                    existing_sativa.add((seq_id, tax_id))
                    stats['validated_with_col'] += 1
                else:
                    stats['duplicates_skipped'] += 1
            else:
                original_tax_id = orig_map.get(seq_id)
                if not original_tax_id:
                    logger.warning(f"No original taxonomy found for validated sequence {seq_id}")
                    stats['validated_without_col'] += 1
                    continue
                if (seq_id, original_tax_id) not in existing_sativa:
                    seq_tax_batch.append(integrator.prepare_seq_tax_entry(seq_id, original_tax_id, 'Correct', 1.0))
                    existing_sativa.add((seq_id, original_tax_id))
                    stats['validated_without_col'] += 1
                else:
                    stats['duplicates_skipped'] += 1
            _flush_seq_tax()
        except Exception as e:
            logger.error(f"Error processing validated sequence {seq_id}: {e}")
            stats['errors'] += 1

    # Process modifications from .mis

    for mod in modifications:
        try:
            # Original + existing CoL taxonomy from the pre-fetched maps (no per-seq query)
            original_tax_id = orig_map.get(mod.seq_id)
            if not original_tax_id:
                logger.warning(f"No original taxonomy for {mod.seq_id}")
                stats['errors'] += 1
                continue

            _col_entry = col_map.get(mod.seq_id)
            col_taxonomy = _col_entry[1] if _col_entry else None

            # Classify modification
            mod_type = integrator.classifier.classify_modification(
                mod.original_taxonomy,
                mod.proposed_taxonomy
            )
            stats[f'type_{mod_type.value}'] += 1

            # Decide action with new logic
            action = integrator.decide_integration_action(mod, mod_type, col_taxonomy)

            # Skip if no action needed
            if not action.create_seq_tax and not action.create_warning:
                stats['skipped_no_action'] += 1
                continue

            # Create seq_tax if needed
            if action.create_seq_tax:
                # Log taxonomy before processing
                logger.debug(f"Processing taxonomy for {mod.seq_id}: {action.taxonomy_to_use}")

                # Check for None values in taxonomy
                none_keys = [k for k, v in action.taxonomy_to_use.items() if v is None]
                if none_keys:
                    logger.debug(f"Found None values in taxonomy for {mod.seq_id} at keys: {none_keys}")

                # Get or create taxonomy (the SATIVA-proposed one may be new; modifications
                # are the minority, so a per-item get_or_create here is acceptable).
                tax_id = integrator.get_or_create_taxonomy(
                    action.taxonomy_to_use,
                    parent_tax_id=original_tax_id
                )

                # Duplicate check against the pre-fetched Sativa links (no per-item query)
                if (mod.seq_id, tax_id) in existing_sativa:
                    stats['duplicates_skipped'] += 1
                    continue

                # Prepare entry
                seq_tax_batch.append(integrator.prepare_seq_tax_entry(
                    mod.seq_id, tax_id, action.seq_tax_status, mod.confidence
                ))
                existing_sativa.add((mod.seq_id, tax_id))

                stats['seq_tax_created_pending'] += 1

            # Create warning if needed
            if action.create_warning:
                warnings_batch.append(integrator.prepare_warning_entry(
                    mod, action.warning_type, original_tax_id, action.col_validation_result
                ))
                stats['warnings_created_pending'] += 1

            # Bulk insert when batches are full
            if len(seq_tax_batch) >= config['performance']['bulk_insert_size']:
                integrator.bulk_insert_seq_tax(seq_tax_batch)
                stats['seq_tax_created'] += len(seq_tax_batch)
                seq_tax_batch = []

            if len(warnings_batch) >= config['performance']['bulk_insert_size']:
                integrator.bulk_insert_warnings(warnings_batch)
                stats['warnings_created'] += len(warnings_batch)
                warnings_batch = []

        except Exception as e:
            logger.error(f"Error processing modification {mod.seq_id}: {e}")
            logger.error(f"Modification details - Type: {mod_type}, Barcode: {mod.barcode}, Tax level: {mod.tax_level}")
            logger.error(f"Original taxonomy: {mod.original_taxonomy}")
            logger.error(f"Proposed taxonomy: {mod.proposed_taxonomy}")
            if 'action' in locals():
                logger.error(f"Taxonomy to use: {action.taxonomy_to_use}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            stats['errors'] += 1

    # Insert remaining
    if seq_tax_batch:
        integrator.bulk_insert_seq_tax(seq_tax_batch)
        stats['seq_tax_created'] += len(seq_tax_batch)

    if warnings_batch:
        integrator.bulk_insert_warnings(warnings_batch)
        stats['warnings_created'] += len(warnings_batch)

    return dict(stats)


def scan_sativa_output_directory(sativa_dir: Path, barcodes: List[str],
                                 tax_levels: List[str]) -> List[Tuple[Path, str, str]]:
    """
    Scan SATIVA output directory for .mis files.

    Returns:
        List of (mis_file_path, barcode, tax_level) tuples
    """
    mis_files = []

    for barcode in barcodes:
        for tax_level in tax_levels:
            alignment_dir = sativa_dir / barcode / f"alignments_{tax_level}"

            if not alignment_dir.exists():
                logger.debug("Directory not found: %s", alignment_dir)
                continue

            for taxon_dir in alignment_dir.iterdir():
                if not taxon_dir.is_dir():
                    continue

                # SATIVA outputs are prefixed with "{barcode}_{tax_level}_{taxon}.mis"
                mis_candidates = sorted(taxon_dir.glob("*.mis"))
                if not mis_candidates:
                    logger.debug("No .mis file found in %s", taxon_dir)
                    continue

                mis_file = mis_candidates[0]
                mis_files.append((mis_file, barcode, tax_level))

    return mis_files


def merge_stats(global_stats: Dict, file_stats: Dict):
    """Merge file statistics into global statistics"""
    for key, value in file_stats.items():
        if isinstance(value, (int, float)):
            global_stats[key] = global_stats.get(key, 0) + value


def generate_reports(stats: Dict, config: Dict, output_dir: Path):
    """Generate JSON and CSV reports"""
    output_dir.mkdir(parents=True, exist_ok=True)

    # JSON report
    if config['reporting']['export_json']:
        json_path = output_dir / "sativa_integration_report.json"
        with open(json_path, 'w') as f:
            json.dump(stats, f, indent=2, default=str)
        logger.info(f"JSON report written to {json_path}")

    # CSV reports will be added in future version
    logger.info("Report generation completed")


def main(snakemake):
    """Main entry point"""
    start_time = time.time()

    logger.info("=" * 70)
    logger.info("SATIVA INTEGRATION v2.0")
    logger.info("=" * 70)

    # Extract configuration
    config = snakemake.config
    sativa_dir = Path(snakemake.params.sativa_output_dir)
    barcodes = snakemake.params.barcodes
    tax_levels = snakemake.params.tax_levels
    output_report = Path(snakemake.output.report)

    # Initialize database
    db_manager = get_database_manager(config)

    # Create sativa_warnings table
    if config['sativa_integration']['database']['create_warnings_table']:
        create_sativa_warnings_table(db_manager)

    # Initialize CoL validator
    col_config = config['sativa_integration']['col_api']
    cache_config = config['sativa_integration']['cache']

    col_validator = CoLValidator(
        max_requests_per_second=col_config['max_requests_per_second'],
        timeout=col_config['timeout_per_request'],
        cache_enabled=cache_config['enabled'],
        cache_file=Path(cache_config['file_cache_path']),
        max_memory_cache=cache_config['max_memory_entries'],
        dataset=col_config.get('dataset', '3LR')
    )
    logger.info("CoL source (integration): %s", col_validator.provenance)

    # Scan for .mis files
    mis_files = scan_sativa_output_directory(sativa_dir, barcodes, tax_levels)
    logger.info(f"Found {len(mis_files)} .mis files to process")

    if not mis_files:
        logger.warning("No .mis files found")
        generate_reports({}, config['sativa_integration'], output_report.parent)
        return

    # Process files in parallel
    global_stats = defaultdict(int)
    global_stats['metadata'] = {
        'start_time': datetime.now().isoformat(),
        'total_mis_files': len(mis_files)
    }

    perf_config = config['sativa_integration']['performance']

    with ThreadPoolExecutor(max_workers=perf_config['max_workers_files']) as executor:
        futures = {
            executor.submit(
                process_single_mis_file,
                mis_file, barcode, tax_level,
                col_validator, db_manager, config['sativa_integration']
            ): (mis_file, barcode, tax_level)
            for mis_file, barcode, tax_level in mis_files
        }

        completed = 0
        for future in as_completed(futures):
            completed += 1
            mis_file, barcode, tax_level = futures[future]

            try:
                file_stats = future.result()
                merge_stats(global_stats, file_stats)
                logger.info(f"[{completed}/{len(mis_files)}]  {mis_file.name}: "
                          f"{file_stats.get('total_modifications', 0)} modifications")
            except Exception as e:
                logger.error(f"[{completed}/{len(mis_files)}]  {mis_file.name}: {e}")
                global_stats['failed_files'] += 1

    # Save CoL cache
    col_validator.close()

    # Finalize statistics
    global_stats['metadata']['end_time'] = datetime.now().isoformat()
    global_stats['metadata']['runtime_seconds'] = time.time() - start_time
    global_stats['col_api_stats'] = col_validator.stats

    # Generate reports
    generate_reports(dict(global_stats), config['sativa_integration'], output_report.parent)

    logger.info("=" * 70)
    logger.info("INTEGRATION COMPLETED")
    logger.info(f"Runtime: {global_stats['metadata']['runtime_seconds']:.2f}s")
    logger.info(f"seq_tax created: {global_stats.get('seq_tax_created', 0)}")
    logger.info(f"Warnings created: {global_stats.get('warnings_created', 0)}")
    logger.info(f"CoL API calls: {col_validator.stats['api_calls']}")
    logger.info(f"Cache hits: {col_validator.stats['cache_hits']}")
    logger.info("=" * 70)


if __name__ == "__main__":
    if 'snakemake' in globals():
        main(snakemake)
    else:
        logger.error("This script must be run through Snakemake")
        sys.exit(1)
