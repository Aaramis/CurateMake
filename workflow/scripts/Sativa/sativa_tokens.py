"""Normalisation of the taxon tokens SATIVA writes in its .mis files.

Deliberately DEPENDENCY-FREE: the standard library only.

This lived on `SativaMisParser` in integrate_sativa.py, and vote_sativa.py imported it
from there to keep one single normalisation. That import chain was
    vote_sativa -> integrate_sativa -> col_taxonomic_matcher -> requests
and `requests` is not declared in workflow/envs/sativa.yaml, so the vote job died with
ModuleNotFoundError before it could even open its log file -- a 0-byte log, no consensus,
and every SATIVA correction silently dropped downstream (a sequence missing from the
consensus is read as VALIDATED). Pulling an ORM and an HTTP client in to normalise a
string was the actual mistake; this module keeps the shared behaviour without them.
"""
import re

# What SATIVA/UNITE write when a rank carries no usable name.
MISSING_TAX_TOKENS = {"", "na", "n/a", "none", "null", "unknown", "unidentified"}

# SATIVA merges taxa whose sequences it cannot tell apart into one pseudo-taxon named
# `__TAXCLUSTER<n>__<name>` (epa_trainer.check_identical_ranks), where <name> is just the
# first member of the merged set. Downstream normalisation turns the underscores into
# spaces, hence both forms. It is an internal marker, not a taxonomic name: the only thing
# it asserts is that the rank could NOT be resolved, so it is treated as unclassified.
# The raw proposal is still recorded in sativa_warnings.
TAXCLUSTER_RE = re.compile(r"^_*TAXCLUSTER\d*_*\s*", re.IGNORECASE)


def normalize_tax_token(value):
    """Return a clean taxon name, or None when the rank carries no usable name."""
    if value is None:
        return None
    cleaned = value.strip()
    if not cleaned:
        return None
    if TAXCLUSTER_RE.match(cleaned):
        return None
    if cleaned.lower() in MISSING_TAX_TOKENS:
        return None
    return cleaned
