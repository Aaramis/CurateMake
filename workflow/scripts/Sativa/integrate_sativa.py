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

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
from db_manager import get_database_manager

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
                 cache_file: Optional[Path] = None, max_memory_cache: int = 100000):
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

        # Query CoL API
        usage = self._query_col_api(taxonomy)

        if not usage:
            # No match found
            result = CoLValidationResult(
                success=False,
                match_quality=CoLMatchQuality.NO_MATCH,
                col_taxonomy=None
            )
        else:
            # Extract taxonomy from response
            col_taxonomy = self._extract_classification_from_usage(usage)

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
        """Close session and save cache"""
        self.save_cache()
        self.session.close()


# ============================================================================
# SATIVA .mis/.tax Parser
# ============================================================================

class SativaMisParser:
    """Parse SATIVA .mis and .tax files"""

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

        # Remove empty strings and strip whitespace
        levels = [level.strip() if level and level.strip() else None for level in levels]

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
        # Parse into lists
        original_levels = [x.strip() for x in original_path.split(';') if x.strip()]
        proposed_levels = [x.strip() for x in proposed_path.split(';') if x.strip()]

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

    def get_taxonomy_by_id(self, tax_id: str) -> Optional[Dict[str, str]]:
        """Retrieve a taxonomy dict by its ID."""
        from sqlalchemy import text

        taxonomy_table = self._get_table('taxonomy')

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

        with self.db_manager.get_session() as session:
            result = session.execute(text(f"""
                SELECT 1 FROM {seq_tax_table}
                WHERE seq_id = :seq_id AND tax_id = :tax_id AND backbone = :backbone
            """), {"seq_id": seq_id, "tax_id": tax_id, "backbone": backbone}).fetchone()

            return result is not None

    def decide_integration_action(self, modification: ModificationRecord,
                                   mod_type: ModificationType,
                                   existing_col_taxonomy: Optional[Dict[str, str]]) -> IntegrationAction:
        """
        Decide what actions to take for a modification.

        New logic:
        - No modification (VALIDATED): If CoL exists -> create new entry with same taxonomy + backbone=sativa
                                       If no CoL -> do nothing
        - Modification (deletion/completion/substitution):
            If CoL exists:
                - Validate new taxonomy via CoL API
                - If validation OK -> create new entry with new taxonomy + backbone=sativa
                - If validation fails -> create new entry with old CoL taxonomy + backbone=sativa
            If no CoL exists:
                - Validate new taxonomy via CoL API
                - If validation OK -> create new entry with new taxonomy + backbone=sativa
                - If validation fails -> do nothing

        Args:
            modification: The modification record
            mod_type: The classified modification type
            existing_col_taxonomy: Existing CoL taxonomy for this sequence (if any)

        Returns:
            IntegrationAction with all necessary information
        """
        if mod_type == ModificationType.VALIDATED:
            # No modification from SATIVA
            if existing_col_taxonomy:
                # CoL exists -> create new entry with same taxonomy + backbone=sativa
                return IntegrationAction(
                    create_seq_tax=True,
                    update_seq_tax_backbone=False,
                    create_warning=False,
                    taxonomy_to_use=existing_col_taxonomy,
                    seq_tax_status='Correct',
                    warning_type=None,
                    col_validation_result=None
                )
            else:
                # No CoL -> do nothing
                return IntegrationAction(
                    create_seq_tax=False,
                    update_seq_tax_backbone=False,
                    create_warning=False,
                    taxonomy_to_use={},
                    seq_tax_status='',
                    warning_type=None,
                    col_validation_result=None
                )

        elif mod_type == ModificationType.DELETION_SIMPLE:
            # Simple deletion: only create warning, no seq_tax
            return IntegrationAction(
                create_seq_tax=False,
                update_seq_tax_backbone=False,
                create_warning=True,
                taxonomy_to_use={},
                seq_tax_status='',
                warning_type='DELETION_SIMPLE',
                col_validation_result=None
            )

        else:
            # DELETION_MULTIPLE, SUBSTITUTION, COMPLETION
            # -> Validate new SATIVA taxonomy with CoL
            proposed_tax_dict = SativaMisParser.parse_taxonomy_path_to_dict(
                modification.proposed_taxonomy
            )

            col_result = self.col_validator.validate_taxonomy(proposed_tax_dict)

            if col_result.success and col_result.match_quality in [CoLMatchQuality.EXACT, CoLMatchQuality.MODIFIED]:
                # CoL validation successful -> use CoL returned taxonomy
                return IntegrationAction(
                    create_seq_tax=True,
                    update_seq_tax_backbone=False,
                    create_warning=col_result.match_quality == CoLMatchQuality.MODIFIED,
                    taxonomy_to_use=col_result.col_taxonomy,
                    seq_tax_status='Correct' if col_result.match_quality == CoLMatchQuality.EXACT else 'Check',
                    warning_type='COL_MODIFIED_SATIVA' if col_result.match_quality == CoLMatchQuality.MODIFIED else None,
                    col_validation_result=col_result
                )
            else:
                # CoL validation failed (NO_MATCH)
                if existing_col_taxonomy:
                    # CoL exists -> use old CoL taxonomy + backbone=sativa
                    return IntegrationAction(
                        create_seq_tax=True,
                        update_seq_tax_backbone=False,
                        create_warning=True,
                        taxonomy_to_use=existing_col_taxonomy,
                        seq_tax_status='Check',
                        warning_type='COL_NO_MATCH',
                        col_validation_result=col_result
                    )
                else:
                    # No CoL exists -> do nothing
                    return IntegrationAction(
                        create_seq_tax=False,
                        update_seq_tax_backbone=False,
                        create_warning=True,
                        taxonomy_to_use={},
                        seq_tax_status='',
                        warning_type='COL_NO_MATCH',
                        col_validation_result=col_result
                    )

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

        with self.db_manager.get_session() as session:
            # Debug: check what backbones exist for this sequence
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

    def prepare_warning_entry(self, modification: ModificationRecord,
                             warning_type: str, original_tax_id: str,
                             col_result: Optional[CoLValidationResult]) -> Dict:
        """Prepare a warning entry dict for bulk insert"""
        # Map ModificationType to database-allowed values
        mod_type = self.classifier.classify_modification(
            modification.original_taxonomy,
            modification.proposed_taxonomy
        )
        # Map DELETION_SIMPLE and DELETION_MULTIPLE to just DELETION for DB constraint
        db_mod_type = 'DELETION' if 'DELETION' in mod_type.value else mod_type.value

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

    def bulk_insert_seq_tax(self, entries: List[Dict]):
        """Bulk insert seq_tax entries"""
        if not entries:
            return

        from sqlalchemy import text
        from sqlalchemy.dialects.postgresql import insert as pg_insert

        seq_tax_table = self._get_table('seq_tax')

        with self.db_manager.get_session() as session:
            if self.db_manager.db_type == "postgresql":
                # Use PostgreSQL-specific INSERT ... ON CONFLICT
                for entry in entries:
                    session.execute(text(f"""
                        INSERT INTO {seq_tax_table} (seq_id, tax_id, backbone, status, confidence)
                        VALUES (:seq_id, :tax_id, :backbone, :status, :confidence)
                        ON CONFLICT (seq_id, tax_id, backbone) DO NOTHING
                    """), entry)
            else:
                # SQLite: Use INSERT OR IGNORE
                for entry in entries:
                    session.execute(text(f"""
                        INSERT OR IGNORE INTO {seq_tax_table} (seq_id, tax_id, backbone, status, confidence)
                        VALUES (:seq_id, :tax_id, :backbone, :status, :confidence)
                    """), entry)

            session.commit()

    def bulk_insert_warnings(self, entries: List[Dict]):
        """Bulk insert sativa_warnings entries"""
        if not entries:
            return

        from sqlalchemy import text
        import uuid

        warnings_table = self._get_table('sativa_warnings')

        with self.db_manager.get_session() as session:
            for entry in entries:
                # Add warning_id
                entry['warning_id'] = str(uuid.uuid4())

                session.execute(text(f"""
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
                """), entry)

            session.commit()


# ============================================================================
# Main Processing Logic
# ============================================================================

def create_sativa_warnings_table(db_manager):
    """Create sativa_warnings table if it doesn't exist"""
    from sqlalchemy import text

    integrator = DatabaseIntegrator(db_manager, None)
    warnings_table = integrator._get_table('sativa_warnings')

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
                        'COL_MODIFIED_SATIVA', 'COL_NO_MATCH'
                    )),
                    sativa_modification_type TEXT CHECK (sativa_modification_type IN (
                        'DELETION', 'SUBSTITUTION', 'COMPLETION'
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

            # Create indices
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_seq ON {warnings_table}(seq_id)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_type ON {warnings_table}(warning_type)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_barcode ON {warnings_table}(barcode)"))

        else:  # SQLite
            session.execute(text(f"""
                CREATE TABLE IF NOT EXISTS {warnings_table} (
                    warning_id TEXT PRIMARY KEY,
                    seq_id TEXT NOT NULL,
                    seq_original_id TEXT NOT NULL,
                    original_tax_id TEXT NOT NULL,
                    warning_type TEXT NOT NULL CHECK (warning_type IN (
                        'DELETION_SIMPLE', 'DELETION_MULTIPLE',
                        'COL_MODIFIED_SATIVA', 'COL_NO_MATCH'
                    )),
                    sativa_modification_type TEXT CHECK (sativa_modification_type IN (
                        'DELETION', 'SUBSTITUTION', 'COMPLETION'
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

            # Create indices
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_seq ON {warnings_table}(seq_id)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_type ON {warnings_table}(warning_type)"))
            session.execute(text(f"CREATE INDEX IF NOT EXISTS idx_sativa_warnings_barcode ON {warnings_table}(barcode)"))

        session.commit()

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

        # Process validated sequences with new logic
        for seq_id in validated_seq_ids:
            try:
                original_tax_id = integrator.get_original_tax_id_for_sequence(seq_id)
                col_taxonomy = integrator.get_col_taxonomy_for_sequence(seq_id)

                if col_taxonomy:
                    # CoL exists -> create new entry with same taxonomy + backbone=sativa
                    tax_id = integrator.get_or_create_taxonomy(col_taxonomy, parent_tax_id=original_tax_id)

                    if not integrator.check_seq_tax_exists(seq_id, tax_id, 'Sativa'):
                        seq_tax_batch.append(integrator.prepare_seq_tax_entry(
                            seq_id, tax_id, 'Correct', 1.0
                        ))
                        stats['validated_with_col'] += 1
                    else:
                        stats['duplicates_skipped'] += 1
                else:
                    # Pas de CoL : réutiliser la taxonomie d'origine (Raw) pour créer un lien Sativa
                    if not original_tax_id:
                        logger.warning(f"No original taxonomy found for validated sequence {seq_id}")
                        stats['validated_without_col'] += 1
                        continue

                    original_taxonomy = integrator.get_taxonomy_by_id(original_tax_id)
                    if not original_taxonomy:
                        logger.warning(
                            "Unable to load original taxonomy record %s for %s",
                            original_tax_id, seq_id
                        )
                        stats['validated_without_col'] += 1
                        continue

                    if not integrator.check_seq_tax_exists(seq_id, original_tax_id, 'Sativa'):
                        seq_tax_batch.append(integrator.prepare_seq_tax_entry(
                            seq_id, original_tax_id, 'Correct', 1.0
                        ))
                        stats['validated_without_col'] += 1
                    else:
                        stats['duplicates_skipped'] += 1
            except Exception as e:
                logger.error(f"Error processing validated sequence {seq_id}: {e}")
                stats['errors'] += 1

    # Process modifications from .mis

    for mod in modifications:
        try:
            # Get original taxonomy
            original_tax_id = integrator.get_original_tax_id_for_sequence(mod.seq_id)
            if not original_tax_id:
                logger.warning(f"No original taxonomy for {mod.seq_id}")
                stats['errors'] += 1
                continue

            # Get existing CoL taxonomy
            col_taxonomy = integrator.get_col_taxonomy_for_sequence(mod.seq_id)

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

                # Get or create taxonomy
                tax_id = integrator.get_or_create_taxonomy(
                    action.taxonomy_to_use,
                    parent_tax_id=original_tax_id
                )

                # Check for duplicates
                if integrator.check_seq_tax_exists(mod.seq_id, tax_id, 'Sativa'):
                    stats['duplicates_skipped'] += 1
                    continue

                # Prepare entry
                seq_tax_batch.append(integrator.prepare_seq_tax_entry(
                    mod.seq_id, tax_id, action.seq_tax_status, mod.confidence
                ))

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
        max_memory_cache=cache_config['max_memory_entries']
    )

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
