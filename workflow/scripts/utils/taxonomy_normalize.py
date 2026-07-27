"""Shared taxonomy normalisation for all source formatters.

Every source (UNITE, BOLD, PLANiTS, CALeDNA, PERSONAL) previously reimplemented its
own rank-cleaning with a different, hard-coded placeholder set, so the same
uninformative label was dropped by one source and kept as a real taxon by another.
This module is the single source of truth:

  * clean_rank(value)     -> cleaned rank string, or None if it is a placeholder
  * clean_species(value)  -> normalised binomial 'Genus epithet', or None if the
                             species is not actually identified
  * is_placeholder(value) -> True when the value carries no real taxonomic information

Design goals:
  - one placeholder vocabulary covering NA/none/null/unknown/unidentified/unclassified/
    unassigned/undetermined, "Incertae sedis", "environmental sample", "uncultured",
    "metagenome", etc.;
  - species not identified to species ("Genus sp.", "sp.", "unknown Genus",
    "Genus cf./aff./nr. x", OTU/strain codes) collapse to None so the genus is kept
    without inventing a pseudo-species;
  - binomial normalisation: drop infraspecific qualifiers (subsp./ssp./var./f./forma)
    and author citations/strain codes, keep the hybrid marker.
"""

import re

# Exact tokens (after lower-casing and underscore->space) that mean "no taxon here".
_PLACEHOLDER_EXACT = frozenset({
    "", "na", "n/a", "n.a.", "none", "null", "nan", "undefined",
    "unknown", "unidentified", "unclassified", "unassigned", "undetermined",
    "unidentifiable", "unresolved", "incertae sedis", "environmental sample",
    "environmental", "uncultured", "uncultured organism", "metagenome",
    "unclassified organism", "root", "cellular organisms",
})

# Substrings that make the whole value a placeholder wherever they appear.
_PLACEHOLDER_SUBSTR = (
    "incertae sedis", "environmental sample", "uncultured", "unclassified",
    "metagenome", "endophyte of", " unclassified",
)

# Leading words that mark a species as not identified (e.g. "unknown Asterochloris").
_UNIDENTIFIED_LEAD = frozenset({
    "unknown", "unidentified", "unclassified", "uncultured", "environmental",
})

# Epithet tokens that indicate the species is not resolved to species rank.
_SP_EPITHET = frozenset({
    "sp", "spp", "cf", "aff", "nr", "indet", "gen", "genosp", "sect", "var",
})

# Domain / superkingdom labels that some sources put in the "kingdom" field (e.g.
# CALeDNA stores "Eukaryota"). A domain is not a kingdom, so callers treat these as
# "kingdom unknown" and let downstream harmonisation fill the real kingdom.
_DOMAIN_TOKENS = frozenset({
    "eukaryota", "eukarya", "bacteria", "archaea", "prokaryota", "prokaryotae",
    "cellular organisms",
})

# Infraspecific qualifiers: everything from here on is dropped from the binomial.
_INFRASPECIFIC = frozenset({
    "subsp", "ssp", "var", "f", "forma", "subvar", "subf", "cv", "morph", "race",
})

_MULTISPACE = re.compile(r"\s+")


def _basic_clean(value):
    """Underscore->space, collapse whitespace, strip. None-safe."""
    if value is None:
        return None
    v = str(value).replace("_", " ")
    v = _MULTISPACE.sub(" ", v).strip()
    return v


def is_placeholder(value):
    """True when the value carries no real taxonomic information."""
    v = _basic_clean(value)
    if not v:
        return True
    low = v.lower()
    if low in _PLACEHOLDER_EXACT:
        return True
    return any(sub in low for sub in _PLACEHOLDER_SUBSTR)


def is_domain(value):
    """True if the value is a domain / superkingdom rather than a kingdom."""
    v = _basic_clean(value)
    return bool(v) and v.lower() in _DOMAIN_TOKENS


def clean_rank(value):
    """Clean a non-species rank. Returns the cleaned string or None."""
    if is_placeholder(value):
        return None
    return _basic_clean(value)


def _strip_token(tok):
    """Lower-case a token with trailing punctuation removed (for classification)."""
    return tok.strip(".,;:()[]").lower()


def clean_species(value):
    """Normalise a species value.

    Returns a 'Genus epithet' binomial, or None when the species is not actually
    identified (placeholders, 'Genus sp.', 'unknown Genus', qualifiers, etc.).
    """
    v = _basic_clean(value)
    if not v or is_placeholder(v):
        return None

    # Normalise the hybrid multiplication sign to a plain 'x' token.
    v = v.replace("×", " x ")
    v = _MULTISPACE.sub(" ", v).strip()

    tokens = v.split(" ")
    if not tokens:
        return None

    # "unknown Genus" / "unidentified Genus" ... -> not identified.
    if _strip_token(tokens[0]) in _UNIDENTIFIED_LEAD:
        return None

    # A lone token that is an 'sp.'-like marker -> not identified.
    if len(tokens) == 1:
        return None if _strip_token(tokens[0]) in _SP_EPITHET else v

    genus = tokens[0]

    # Hybrid: "Genus x epithet" -> keep "Genus x epithet".
    if _strip_token(tokens[1]) == "x" and len(tokens) >= 3:
        return f"{genus} x {tokens[2]}"

    epithet = tokens[1]
    if _strip_token(epithet) in _SP_EPITHET:
        # "Genus sp.", "Genus cf. ...", "Genus aff. ..." -> species not resolved.
        return None

    # Keep the binomial only; drop authors, strain/OTU codes, infraspecific parts.
    return f"{genus} {epithet}"
