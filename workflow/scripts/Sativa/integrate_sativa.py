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
import re
import time
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple, Set, Iterator
from collections import defaultdict
from dataclasses import dataclass, asdict
from enum import Enum
from concurrent.futures import ThreadPoolExecutor, as_completed
import pickle

# Add utils + sibling scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
sys.path.insert(0, str(Path(__file__).parent.parent))
from db_manager import get_database_manager, pg_hint
# Reuse the single corrected CoL client/matcher (EXACT-search strategy, synonym resolution,
# thread-safe polite rate limiting, pinned DOI'd release + provenance). One source of truth
# for CoL logic across the pipeline (see col_taxonomic_matcher.py).
from sativa_tokens import (normalize_tax_token as _normalize_tax_token,
                           MISSING_TAX_TOKENS as _MISSING_TAX_TOKENS,
                           TAXCLUSTER_RE as _TAXCLUSTER_RE)
from col_taxonomic_matcher import (CoLAPIClient, TaxonomyMatcher, TaxonomyRecord,
                                   match_taxonomies_parallel)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def _attach_snakemake_log(snakemake) -> None:
    """Send this script's log to the rule's `log:` file.

    basicConfig alone writes to stderr, which Snakemake does not capture for
    `script:` directives, so a failing clade would otherwise leave no diagnostic.
    """
    try:
        log_attr = getattr(snakemake, 'log', None)
        if not log_attr:
            return
        log_file = log_attr[0] if isinstance(log_attr, (list, tuple)) else log_attr
        if any(getattr(h, 'baseFilename', None) == str(Path(log_file).absolute())
               for h in logger.handlers):
            return
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(log_file, mode="w")
        handler.setFormatter(logging.Formatter(
            '%(asctime)s | %(levelname)s | %(name)s | %(message)s'))
        logger.addHandler(handler)
    except Exception as exc:  # noqa: BLE001 - logging must never break the job
        logger.warning("could not attach file logging: %s", exc)


import uuid as _uuid

# IDENTITY columns of sativa_warnings. FROZEN — the primary key is a uuid5 over exactly
# these, so appending one would give every existing row a different id and orphan the
# 1.5 M rows already in the table. New descriptive columns go in _WARNING_EXTRA_COLUMNS
# below, which the key deliberately ignores.
_WARNING_ID_COLUMNS = (
    'warning_id', 'seq_id', 'seq_original_id', 'original_tax_id',
    'warning_type', 'sativa_modification_type', 'sativa_confidence',
    'sativa_original_taxonomy', 'sativa_proposed_taxonomy',
    'deleted_levels', 'col_validation_attempted', 'col_returned_taxonomy',
    'col_match_quality', 'barcode', 'tax_level', 'mislabeled_level',
    'per_rank_confidence',
)

# Descriptive only: how many independent splits agreed, out of how many. NULL when
# sativa.replicates == 1 (no vote was taken). Because these sit outside the key, a row
# re-integrated with a different replicate count keeps its original id and the
# ON CONFLICT DO NOTHING below leaves the stored counts alone — clear the table if you
# want a different `replicates` setting reflected in rows that already exist.
_WARNING_EXTRA_COLUMNS = ('sativa_votes', 'sativa_replicates')

# Columns of sativa_warnings, in the order used by the bulk insert.
_WARNING_COLUMNS = _WARNING_ID_COLUMNS + _WARNING_EXTRA_COLUMNS

# Fixed namespace -> the same row always yields the same warning_id, across
# processes and across runs. Do not change it: it would orphan existing rows.
_WARNING_ID_NAMESPACE = _uuid.UUID('b7f0c3a2-5d41-5e88-9c2b-1f6a4e0d7c35')


def _warning_id(entry: Dict[str, Any]) -> str:
    """Deterministic, content-addressed primary key for a sativa_warnings row."""
    key = '|'.join(str(entry.get(c) or '') for c in _WARNING_ID_COLUMNS[1:])
    return str(_uuid.uuid5(_WARNING_ID_NAMESPACE, key))


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
    # Set only when the record comes from a consensus over several splits (sativa_vote).
    votes: Optional[int] = None
    n_replicates: Optional[int] = None


@dataclass
class CoLValidationResult:
    """Result from CoL validation"""
    success: bool
    match_quality: Optional[CoLMatchQuality]
    col_taxonomy: Optional[Dict[str, str]]  # {kingdom, phylum, ...}
    col_id: Optional[str] = None  # accepted CoL name-usage id (for reproducibility)
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
        # NB: every HTTP call goes through self._api (CoLAPIClient), which owns
        # its own session, throttle and 429 handling. The validator used to keep
        # a second requests.Session plus a _rate_limit() of its own; neither was
        # ever used, so both are gone.
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
        match = self._matcher.match_taxonomy(self._record_from_taxonomy(taxonomy))
        self._refresh_stats()
        result = self._result_from_match(taxonomy, match)
        self._add_to_cache(cache_key, result)
        return result

    @staticmethod
    def _record_from_taxonomy(taxonomy: Dict[str, str], tax_id: int = 0) -> TaxonomyRecord:
        """Build the TaxonomyRecord the shared CoL matcher expects."""
        return TaxonomyRecord(
            tax_id=tax_id,
            kingdom=taxonomy.get('kingdom'), phylum=taxonomy.get('phylum'),
            class_=taxonomy.get('class_'), order_=taxonomy.get('order_'),
            family=taxonomy.get('family'), genus=taxonomy.get('genus'),
            species=taxonomy.get('species'))

    def _refresh_stats(self) -> None:
        self.stats['api_calls'] = self._api.api_calls
        self.stats['failed_calls'] = self._api.failed_calls

    def _result_from_match(self, taxonomy: Dict[str, str], match) -> CoLValidationResult:
        """Turn a MatchResult into a CoLValidationResult.

        Shared by the serial path (validate_taxonomy) and the parallel one
        (prefetch), so both classify a CoL answer identically.
        """
        if not match.success or not match.col_taxonomy:
            return CoLValidationResult(
                success=False,
                match_quality=CoLMatchQuality.NO_MATCH,
                col_taxonomy=None
            )
        quality = (CoLMatchQuality.EXACT
                   if self._taxonomies_match(taxonomy, match.col_taxonomy)
                   else CoLMatchQuality.MODIFIED)
        return CoLValidationResult(
            success=True,
            match_quality=quality,
            col_taxonomy=match.col_taxonomy,
            col_id=match.col_id
        )

    def prefetch(self, taxonomies: List[Dict[str, str]], max_workers: int = 4) -> int:
        """Resolve every distinct lineage against CoL up front, in parallel.

        Delegates to match_taxonomies_parallel -- the same concurrent CoL query
        the early pipeline uses via CoLProcessor.process_batch -- and drops the
        answers into the cache, so the per-modification loop below becomes a
        cache lookup instead of a serial request per lineage.
        """
        pending: Dict[str, Dict[str, str]] = {}
        for tax in taxonomies:
            if not tax:
                continue
            key = self._taxonomy_to_cache_key(tax)
            if key in pending or self._get_from_cache(key) is not None:
                continue
            pending[key] = tax
        if not pending:
            return 0

        ordered = list(pending.items())
        records = [self._record_from_taxonomy(tax, tax_id=i)
                   for i, (_key, tax) in enumerate(ordered)]
        for match in match_taxonomies_parallel(self._matcher, records, max_workers):
            key, tax = ordered[match.original_tax_id]
            self._add_to_cache(key, self._result_from_match(tax, match))
        self._refresh_stats()
        return len(pending)

    def close(self):
        """Close sessions and save cache"""
        self.save_cache()
        self._api.close()


# ============================================================================
# SATIVA .mis/.tax Parser
# ============================================================================

class SativaMisParser:
    """Parse SATIVA .mis and .tax files"""

    MISSING_TAX_TOKENS = _MISSING_TAX_TOKENS

    _TAXCLUSTER_RE = _TAXCLUSTER_RE

    _normalize_tax_token = staticmethod(_normalize_tax_token)

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
    def parse_mis_file(mis_path: Path, barcode: str, tax_level: str,
                       only_native_taxon: Optional[str] = None,
                       min_confidence: float = 0.0) -> List[ModificationRecord]:
        """
        Parse a .mis file into ModificationRecords.

        The same parser reads a raw SATIVA .mis and the consensus written by
        `vote_sativa.py`: the consensus uses the identical 8-column layout and merely
        appends Votes / Replicates / NativeTaxon. Columns 0..7 are read the same way in
        both cases, so there is exactly one path from a SATIVA verdict to the database.

        Args:
            mis_path: Path to .mis file
            barcode: Barcode name (ITS, ITS1, ITS2)
            tax_level: Taxonomic level (family, genus, etc.)
            only_native_taxon: consensus only — keep the sequences this clade owns.
            min_confidence: drop flags below this confidence. 0.0 keeps everything, which
                is the historical behaviour. See config sativa_integration.min_confidence
                for why 0.511 is the interesting value.
                The consensus spans a whole barcode while integration is scattered per
                clade, and column 10 says which native clade each sequence belongs to.

        Returns:
            List of ModificationRecord objects
        """
        records = []
        below_cutoff = 0

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

                    if only_native_taxon is not None:
                        if len(parts) < 11 or parts[10] != only_native_taxon:
                            continue

                    def _int_col(idx):
                        if len(parts) > idx and parts[idx].strip().isdigit():
                            return int(parts[idx])
                        return None

                    if min_confidence > 0.0:
                        try:
                            if float(parts[4]) < min_confidence:
                                below_cutoff += 1
                                continue
                        except ValueError:
                            pass        # malformed: let the parse below report it

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
                            tax_level=tax_level,
                            votes=_int_col(8),
                            n_replicates=_int_col(9)
                        )
                        records.append(record)
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Failed to parse line in {mis_path.name}: {e}")
                        continue

            if below_cutoff:
                logger.info("%s: %d flag(s) dropped below confidence %.3f, %d kept",
                            mis_path.name, below_cutoff, min_confidence, len(records))

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
        Classify a SATIVA modification, comparing the two lineages RANK BY RANK.

        The ranks must be compared in place, never after compacting away the empty ones.
        This taxonomy is full of mid-lineage holes -- `class` is absent from most of it --
        and squeezing them out shifts every deeper rank one slot left, so a proposal that
        merely FILLS a hole looks like it renamed everything below it. Measured on the real
        shapes: `Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA` ->
        `Fungi;Ascomycota;Sordariomycetes;Diaporthales;Valsaceae;Valsa;NA` was reported as a
        SUBSTITUTION when it adds two ranks and contradicts none, and the reverse case
        (dropping `class` alone) was reported as SUBSTITUTION instead of DELETION_SIMPLE.
        `count_deleted_levels` below was already positional, so the two disagreed.

        Rules:
          * any rank filled on BOTH sides with different names -> SUBSTITUTION
          * ranks only gained                                  -> COMPLETION
          * ranks only lost                                    -> DELETION_SIMPLE/MULTIPLE
          * some gained AND some lost, nothing contradicted    -> SUBSTITUTION, on purpose:
            the shape moved in both directions, so it takes the cautious path (CoL
            validation) rather than being filed as a clean addition or removal.

        Args:
            original_path: Original taxonomy path (semicolon-separated)
            proposed_path: Proposed taxonomy path (semicolon-separated)

        Returns:
            ModificationType enum
        """
        def slots(path: str):
            """(7 rank slots, depth) where depth = 1 + index of the deepest populated rank.

            A blank BELOW the depth asserts "this lineage has no such rank" -- common here,
            most lineages have no `class`. A blank at or beyond the depth is an abstention:
            the lineage simply stopped. The two must not be conflated.
            """
            vals = [SativaMisParser._normalize_tax_token(x) or ""
                    for x in (path or "").split(';')]
            n = len(DISPLAY_TAXONOMY_LEVELS)
            vals = (vals + [""] * n)[:n]
            depth = max((i for i, v in enumerate(vals) if v), default=-1) + 1
            return vals, depth

        orig, o_depth = slots(original_path)
        prop, p_depth = slots(proposed_path)

        contradicted = added = removed = 0
        for i, (o, p) in enumerate(zip(orig, prop)):
            o_said, p_said = i < o_depth, i < p_depth
            if o_said and p_said:
                if o != p:
                    # name vs different name, or name vs asserted absence: a real conflict
                    contradicted += 1
            elif p_said and p:
                added += 1          # the proposal reaches a rank the original never claimed
            elif o_said and o:
                removed += 1        # the proposal stops short of a rank the original claimed

        if contradicted:
            return ModificationType.SUBSTITUTION
        if added and removed:
            return ModificationType.SUBSTITUTION
        if added:
            return ModificationType.COMPLETION
        if removed:
            return (ModificationType.DELETION_SIMPLE if removed == 1
                    else ModificationType.DELETION_MULTIPLE)

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

    # Ranks that define a lineage, in the order used as the taxonomy cache key.
    _LINEAGE_RANKS = ('kingdom', 'phylum', 'class_', 'order_', 'family', 'genus', 'species')

    def __init__(self, db_manager, col_validator: CoLValidator):
        self.db_manager = db_manager
        self.col_validator = col_validator
        self.classifier = ModificationClassifier()
        # Per-clade memo of lineage -> tax_id. check_taxonomy_exists is not sargable
        # (COALESCE around every rank defeats uq_taxonomy_lineage), so each miss is a
        # full scan of taxonomy; the .mis files repeat each proposed lineage ~5x on
        # average, so caching removes ~80% of those scans. Safe: a tax_id never changes
        # once known, and a lineage created here is cached under the id just created.
        self._tax_id_cache: Dict[tuple, str] = {}
        self._tax_cache_preloaded = False

    def _get_table(self, table_name: str) -> str:
        """Get schema-qualified table name"""
        if self.db_manager.db_type == "postgresql":
            schema = self.db_manager.connection_params.get('schema', 'public')
            return f"{schema}.{table_name}"
        return table_name


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


    def create_taxonomy(self, taxonomy: Dict[str, str], parent_tax_id: Optional[str],
                        col_id: Optional[str] = None, col_dataset_key: Optional[str] = None) -> str:
        """
        Create a new taxonomy entry.

        Args:
            taxonomy: Dict with taxonomy levels
            parent_tax_id: Parent taxonomy ID (for traceability)
            col_id, col_dataset_key: Optional CoL provenance (accepted usage id + release).

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
            parent_tax_id=parent_tax_id,
            col_id=col_id,
            col_dataset_key=col_dataset_key
        )

    def _stamp_col_provenance(self, tax_id: str, col_id: str, col_dataset_key: Optional[str]) -> None:
        """Record CoL provenance on an existing taxonomy row (only if not already set)."""
        from sqlalchemy import text
        taxonomy_table = self._get_table('taxonomy')

        def op():
            with self.db_manager.get_session() as session:
                session.execute(text(f"""
                    UPDATE {taxonomy_table}
                    SET col_id = :cid, col_dataset_key = :dk
                    WHERE tax_id = :tid AND col_id IS NULL
                """), {"cid": col_id, "dk": col_dataset_key, "tid": tax_id})
                session.commit()

        self.db_manager._execute_with_retry("sativa_stamp_col_provenance", op)

    def preload_taxonomy_cache(self) -> int:
        """Load every known lineage into the in-memory cache in one read.

        check_taxonomy_exists wraps each rank in COALESCE, which defeats
        uq_taxonomy_lineage and turns every lookup into a full scan of taxonomy
        (measured on this cluster: 0.29 s idle, 0.72 s under load, i.e. ~20 min
        of scans for the ~1700 distinct lineages of one clade). Reading the whole
        table costs 2.3 s for 468k rows, so one read replaces the lot.

        The key mirrors the COALESCE semantics exactly -- NULL and '' collapse --
        because both spellings exist in the table (add_taxonomy stores '' while
        older rows hold NULL). Duplicates resolve to the smallest tax_id so the
        choice is stable across jobs, where the old LIMIT 1 without ORDER BY was
        arbitrary.
        """
        from sqlalchemy import text

        taxonomy_table = self._get_table('taxonomy')
        ranks = ", ".join(self._LINEAGE_RANKS)

        def fetch():
            with self.db_manager.get_session() as session:
                return list(session.execute(
                    text(f"SELECT tax_id, {ranks} FROM {taxonomy_table}")))

        rows = self.db_manager._execute_with_retry("sativa_preload_taxonomy", fetch)
        for row in rows:
            key = tuple((v or '') for v in row[1:])
            current = self._tax_id_cache.get(key)
            if current is None or str(row[0]) < current:
                self._tax_id_cache[key] = str(row[0])
        self._tax_cache_preloaded = True
        return len(rows)

    def get_or_create_taxonomy(self, taxonomy: Dict[str, str], parent_tax_id: Optional[str],
                              col_id: Optional[str] = None, col_dataset_key: Optional[str] = None) -> str:
        """Get existing taxonomy or create new one (stamping CoL provenance when given)."""
        cache_key = tuple((taxonomy.get(rank) or '') for rank in self._LINEAGE_RANKS)

        existing_tax_id = self._tax_id_cache.get(cache_key)
        if existing_tax_id is None:
            # With the cache preloaded from the whole table, a miss means the lineage
            # was absent when the snapshot was taken, so the non-sargable
            # check_taxonomy_exists would only repeat that answer at the price of a
            # full scan. Go straight to create_taxonomy: add_taxonomy does its own
            # existence check with an indexable filter_by and settles the race with
            # another writer through IntegrityError.
            existing_tax_id = (None if self._tax_cache_preloaded
                               else self.check_taxonomy_exists(taxonomy))

            if existing_tax_id is None:
                new_tax_id = self.create_taxonomy(taxonomy, parent_tax_id, col_id, col_dataset_key)
                self._tax_id_cache[cache_key] = new_tax_id
                return new_tax_id

            self._tax_id_cache[cache_key] = existing_tax_id

        # Stamping stays outside the cache: a later call may carry a col_id the first
        # one did not, and _stamp_col_provenance only writes when col_id IS NULL.
        if col_id:
            self._stamp_col_provenance(existing_tax_id, col_id, col_dataset_key)
        return existing_tax_id


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
        # Force the pkey (seq_id HASH) batched lookup; see pg_hint() for the measurements.
        sql = (f"{pg_hint([('st', 'seq_tax_pkey')])}"
               f"SELECT st.seq_id, st.tax_id, t.kingdom, t.phylum, t.class_, t.order_, "
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
        # Force the pkey (seq_id HASH) batched lookup; see pg_hint() for the measurements.
        sql = (f"{pg_hint([('st', 'seq_tax_pkey')])}"
               f"SELECT st.seq_id, st.tax_id FROM {st} st WHERE {{clause}} AND st.backbone = :bb")
        return {(r[0], r[1]) for r in self._chunked_rows(
            seq_ids, sql, "st.seq_id", "sativa_bulk_existing", extra_params={"bb": backbone})}

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
            'col_validation_attempted': col_result is not None,
            # NULL unless the row came from a consensus over several splits
            'sativa_votes': modification.votes,
            'sativa_replicates': modification.n_replicates
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

    _vote_columns_checked = False

    def _ensure_vote_columns(self, warnings_table: str) -> None:
        """Idempotently add sativa_votes / sativa_replicates to an existing table.

        Once per process, not per batch. Both backends tolerate being told twice:
        PostgreSQL/YugabyteDB via IF NOT EXISTS, SQLite by raising on a duplicate column,
        which is caught here.
        """
        if DatabaseIntegrator._vote_columns_checked:
            return
        from sqlalchemy import text
        for col in _WARNING_EXTRA_COLUMNS:
            def add_col(column=col):
                with self.db_manager.get_session() as session:
                    if self.db_manager.db_type == 'postgresql':
                        session.execute(text(
                            f"ALTER TABLE {warnings_table} "
                            f"ADD COLUMN IF NOT EXISTS {column} INTEGER"))
                    else:
                        try:
                            session.execute(text(
                                f"ALTER TABLE {warnings_table} ADD COLUMN {column} INTEGER"))
                        except Exception as exc:          # already there
                            if 'duplicate column' not in str(exc).lower():
                                raise
            try:
                self.db_manager._execute_with_retry(f"sativa_add_column_{col}", add_col)
            except Exception as exc:
                logger.warning("could not add %s to %s: %s", col, warnings_table, exc)
        DatabaseIntegrator._vote_columns_checked = True

    def bulk_insert_warnings(self, entries: List[Dict], batch_size: int = 1000):
        """Bulk insert sativa_warnings entries in batches.

        Args:
            entries: List of entry dicts to insert
            batch_size: Number of entries per batch (default 1000)
        """
        if not entries:
            return

        from sqlalchemy import text

        warnings_table = self._get_table('sativa_warnings')

        # Write the historical column set unless a vote actually took place. With
        # sativa.replicates == 1 the statement is byte-for-byte the one this pipeline has
        # always issued, so an existing database needs no migration at all. When a vote IS
        # present, make sure the two columns exist first — CREATE TABLE IF NOT EXISTS never
        # adds a column to a table that is already there.
        has_votes = any(e.get('sativa_votes') is not None for e in entries)
        if has_votes:
            self._ensure_vote_columns(warnings_table)
        cols = _WARNING_COLUMNS if has_votes else _WARNING_ID_COLUMNS

        # Content-addressed primary key. The previous random uuid4 made every
        # replay collide on sativa_warnings_pkey: an ambiguous commit (the batch
        # lands server-side but the client sees an error) replayed the SAME ids,
        # and that violation is permanent, so the retry loop slept 48 minutes
        # before failing the whole clade. Deriving the id from the row content
        # makes both an internal retry and a full task re-run converge instead.
        for entry in entries:
            if 'warning_id' not in entry:
                entry['warning_id'] = _warning_id(entry)

        for i in range(0, len(entries), batch_size):
            batch = entries[i:i + batch_size]

            def insert_batch(batch_entries=batch):
                with self.db_manager.get_session() as session:
                    if self.db_manager.db_type == 'postgresql':
                        # One statement per page instead of one round-trip per
                        # row (measured on this cluster: 3.9k -> 54k rows/s).
                        from psycopg2.extras import execute_values
                        raw = session.connection().connection
                        dbapi = getattr(raw, 'dbapi_connection', raw)
                        values = [tuple(e.get(c) for c in cols) for e in batch_entries]
                        with dbapi.cursor() as cur:
                            execute_values(
                                cur,
                                f"INSERT INTO {warnings_table} ({', '.join(cols)}) VALUES %s "
                                f"ON CONFLICT (warning_id) DO NOTHING",
                                values,
                                page_size=500,
                            )
                    else:
                        stmt = text(
                            f"INSERT OR IGNORE INTO {warnings_table} ({', '.join(cols)}) "
                            f"VALUES ({', '.join(':' + c for c in cols)})"
                        )
                        # Bind every declared column explicitly: a named parameter raises
                        # if the key is simply absent from the dict, which happens as soon
                        # as a batch mixes voted and unvoted rows.
                        for entry in batch_entries:
                            session.execute(stmt, {c: entry.get(c) for c in cols})

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
                        'VALIDATED', 'DELETION_SIMPLE', 'DELETION_MULTIPLE', 'SUBSTITUTION', 'COMPLETION'
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
                        sativa_votes INTEGER,
                        sativa_replicates INTEGER,
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
                        'VALIDATED', 'DELETION_SIMPLE', 'DELETION_MULTIPLE', 'SUBSTITUTION', 'COMPLETION'
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
                        sativa_votes INTEGER,
                        sativa_replicates INTEGER,
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
                            config: Dict,
                            records: Optional[List[ModificationRecord]] = None) -> Dict:
    """
    Process a single .mis file.

    Args:
        mis_file: Path to .mis file. Still used to locate the sibling .tax even when
            `records` is supplied, since that is what marks the VALIDATED sequences.
        barcode: Barcode name
        tax_level: Taxonomic level
        col_validator: CoL validator instance
        db_manager: Database manager
        config: Configuration dict
        records: pre-parsed modifications. Supplied by the consensus path (sativa_vote):
            the flags that survived the vote for this clade. Everything downstream —
            CoL prefetch, classification, seq_tax, warnings — is the code below,
            unchanged, so the vote adds no second route into the database.

    Returns:
        Statistics dict
    """
    stats = defaultdict(int)
    parser = SativaMisParser()
    integrator = DatabaseIntegrator(db_manager, col_validator)

    # One read of taxonomy instead of a full scan per distinct lineage.
    n_tax = integrator.preload_taxonomy_cache()
    logger.info("[%s] taxonomy cache preloaded: %d rows, %d distinct lineages",
                mis_file.name, n_tax, len(integrator._tax_id_cache))

    # Parse .mis file (or take the consensus records the vote already selected)
    modifications = (parser.parse_mis_file(mis_file, barcode, tax_level,
                                           min_confidence=float(config.get('min_confidence', 0.0)))
                     if records is None else records)
    stats['total_modifications'] = len(modifications)
    if records is not None:
        stats['from_consensus'] = 1
        voted = [m.votes for m in modifications if m.votes is not None]
        if voted:
            stats['consensus_min_votes'] = min(voted)
            stats['consensus_replicates'] = modifications[0].n_replicates or 0

    # Resolve every distinct proposed lineage against CoL up front and in
    # parallel, the same way the early pipeline does it. Without this the
    # per-modification loop issues one serial request per lineage.
    # NB: this function receives config['sativa_integration'], not the root config.
    col_workers = int(config.get('col_api', {}).get('max_workers', 8))
    prefetched = col_validator.prefetch(
        [SativaMisParser.parse_taxonomy_path_to_dict(m.proposed_taxonomy)
         for m in modifications if m.proposed_taxonomy],
        max_workers=col_workers)
    if prefetched:
        logger.info("[%s] CoL prefetch: %d distinct lineages resolved with %d workers",
                    mis_file.name, prefetched, min(col_workers, 32))

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
                # When CoL validated the proposal, stamp the accepted CoL usage id +
                # pinned release for reproducibility (col_id is None on a CoL no-match).
                _col_val = action.col_validation_result
                _col_id = _col_val.col_id if _col_val else None
                tax_id = integrator.get_or_create_taxonomy(
                    action.taxonomy_to_use,
                    parent_tax_id=original_tax_id,
                    col_id=_col_id,
                    col_dataset_key=(integrator.col_validator.provenance.get("dataset_key")
                                     if _col_id else None)
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


def integrate_one_clade(snakemake):
    """Integrate a single SATIVA clade (scattered rule sativa_integrate_task).

    Runs as soon as this clade's SATIVA job finishes, so the DB/CoL I/O (which is
    the bulk of integration wall-clock, ~0 CPU) overlaps with the remaining SATIVA
    compute of other clades. Concurrency across clades is gated by the db_lock
    resource (serial on SQLite, parallel on PostgreSQL/YugabyteDB); concurrent
    taxonomy inserts are made safe by add_taxonomy's IntegrityError handling.
    The sativa_warnings table and all ORM tables are created earlier by the schema
    step, so nothing is created here.
    """
    _attach_snakemake_log(snakemake)
    start_time = time.time()
    config = snakemake.config
    barcode = snakemake.wildcards.barcode
    tax_level = snakemake.wildcards.level
    taxon = snakemake.wildcards.taxon
    clade_dir = (Path(snakemake.params.sativa_output_dir) / barcode
                 / f"alignments_{tax_level}" / taxon)
    out_report = Path(snakemake.output.report)

    db_manager = get_database_manager(config)

    col_config = config['sativa_integration']['col_api']
    cache_config = config['sativa_integration']['cache']
    # Per-clade cache file: avoids concurrent writes to a single shared pickle when
    # clade jobs run in parallel (the in-memory cache still covers repeated lineages
    # within the clade, which is where most cache hits come from).
    clade_cache = (clade_dir / "col_integration_cache.pkl") if cache_config['enabled'] else None
    col_validator = CoLValidator(
        max_requests_per_second=col_config['max_requests_per_second'],
        timeout=col_config['timeout_per_request'],
        cache_enabled=cache_config['enabled'],
        cache_file=clade_cache,
        max_memory_cache=cache_config['max_memory_entries'],
        dataset=col_config.get('dataset', '3LR')
    )

    stats = defaultdict(int)
    stats['metadata'] = {
        'barcode': barcode, 'tax_level': tax_level, 'taxon': taxon,
        'start_time': datetime.now().isoformat()
    }

    # With sativa.replicates > 1 the clade integrates the barcode-wide consensus restricted
    # to the sequences it owns, not its own raw .mis: a single-split flag is ~63% artefact
    # (2026-08-26 measurement on three independent ITS2 splits). Sequences that did not
    # reach the threshold are absent from the consensus, so the .tax below marks them
    # VALIDATED — no reproducible objection to their current label.
    consensus_in = getattr(snakemake.input, "consensus", None)
    if isinstance(consensus_in, (list, tuple)):
        consensus_in = consensus_in[0] if consensus_in else None

    # Only clades with a .mis file are integrated (same rule as the previous full-scan
    # behaviour: scan_sativa_output_directory required a .mis).
    mis_candidates = [p for p in sorted(clade_dir.glob("*.mis")) if p.name != "consensus.mis"]

    records = None
    if consensus_in:
        records = SativaMisParser.parse_mis_file(
            Path(consensus_in), barcode, tax_level, only_native_taxon=taxon,
            min_confidence=float(config['sativa_integration'].get('min_confidence', 0.0)))
        stats['metadata']['consensus'] = str(consensus_in)
        logger.info("[%s] consensus mode: %d flag(s) survived the vote for this clade",
                    taxon, len(records))
        if not mis_candidates:
            # SATIVA always writes a .mis, but keep a usable path for the .tax lookup.
            mis_candidates = [clade_dir / f"{barcode}_{tax_level}_{taxon}.mis"]

    if mis_candidates:
        try:
            file_stats = process_single_mis_file(
                mis_candidates[0], barcode, tax_level,
                col_validator, db_manager, config['sativa_integration'],
                records=records)
            merge_stats(stats, file_stats)
        except Exception as e:
            # Do NOT swallow this: writing out_report below would mark the clade as
            # integrated, and Snakemake would report success for a clade whose rows
            # never reached the database.
            logger.exception("[%s] integration failed: %s", taxon, e)
            stats['failed_files'] += 1
            raise
    else:
        logger.info(f"[{taxon}] no .mis file — nothing to integrate")

    col_validator.close()
    stats['metadata']['end_time'] = datetime.now().isoformat()
    stats['metadata']['runtime_seconds'] = time.time() - start_time
    stats['col_api_stats'] = dict(col_validator.stats)

    out_report.parent.mkdir(parents=True, exist_ok=True)
    with open(out_report, 'w') as f:
        json.dump(dict(stats), f, indent=2, default=str)
    logger.info(f"[{taxon}] integrated in {stats['metadata']['runtime_seconds']:.1f}s "
                f"({stats.get('total_modifications', 0)} modifications) -> {out_report}")


def main(snakemake):
    """Main entry point (legacy full-scan mode; the DAG now scatters per clade)."""
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
        # Per-clade scatter (sativa_integrate_task) vs. legacy full-scan (main).
        if getattr(snakemake, 'rule', None) == 'sativa_integrate_task':
            integrate_one_clade(snakemake)
        else:
            main(snakemake)
    else:
        logger.error("This script must be run through Snakemake")
        sys.exit(1)
