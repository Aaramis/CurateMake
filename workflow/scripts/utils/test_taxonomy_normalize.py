"""Test cases for taxonomy_normalize, covering taxonomic problems seen across sources.

Run: python workflow/scripts/utils/test_taxonomy_normalize.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from taxonomy_normalize import (clean_rank, clean_species, is_placeholder, is_domain,
                                clean_kingdom, is_legacy_kingdom, EUKARYOTE_KINGDOMS)

# (input, expected) for is_domain
DOMAIN_CASES = [
    ("Eukaryota", True), ("eukaryota", True), ("Eukarya", True),
    ("Bacteria", True), ("Archaea", True), ("cellular organisms", True),
    ("Fungi", False), ("Plantae", False), ("Animalia", False),
    ("", False), (None, False), ("Chromista", False),
]

# (input, expected) for non-species ranks
RANK_CASES = [
    # --- real placeholders that MUST become None ---
    ("Incertae_sedis", None),        # UNITE: 519x in mini set
    ("Incertae sedis", None),
    ("unidentified", None),
    ("Unidentified", None),          # case-insensitive (UNITE was case-sensitive)
    ("unclassified", None),
    ("unclassified Fungi", None),
    ("environmental sample", None),
    ("uncultured", None),
    ("uncultured fungus", None),
    ("metagenome", None),
    ("NA", None), ("N/A", None), ("none", None), ("null", None),
    ("unknown", None), ("", None), ("   ", None), (None, None),
    # --- real taxa that MUST be kept ---
    ("Fungi", "Fungi"),
    ("Basidiomycota", "Basidiomycota"),
    ("GS2", "GS2"),                  # UNITE early-diverging clade code: NOT a placeholder
    ("Agaricomycetes", "Agaricomycetes"),
    ("Russula", "Russula"),
    ("Eukaryota", "Eukaryota"),
    ("  Aspergillus_ ", "Aspergillus"),   # underscore + whitespace
]

# (input, expected) for species
SPECIES_CASES = [
    # --- "not identified to species" -> None (genus kept elsewhere) ---
    ("Aspergillus sp", None),        # UNITE: 978x "X_sp"
    ("Aspergillus_sp", None),
    ("Aspergillus sp.", None),
    ("Bipolaris sp", None),
    ("Trebouxia sp. OTU I01", None), # BOLD strain/OTU code
    ("Asterochloris sp. LV-2018", None),
    ("Aspergillus sp. NRRL 145", None),
    ("unknown Asterochloris", None), # PLANiTS: 61x "unknown Genus"
    ("unknown Scenedesmus", None),
    ("unknown Chlamydomonadales", None),
    ("Cortinarius cf. balteatus", None),   # cf. = uncertain -> not resolved
    ("Russula aff. nigricans", None),
    ("sp", None), ("sp.", None),
    ("NA", None), ("unidentified", None), ("", None), (None, None),
    # --- genuine binomials -> kept, normalised ---
    ("Russula_abietina", "Russula abietina"),
    ("Amanita muscaria", "Amanita muscaria"),
    ("Amanita muscaria (L.) Lam.", "Amanita muscaria"),   # author citation dropped
    ("Primula veris subsp. macrocalyx", "Primula veris"), # infraspecific dropped
    ("Primula nutans ssp. finmarchia", "Primula nutans"),
    ("Scorzonera purpurea subsp. rosea", "Scorzonera purpurea"),
    ("Salix × fragilis", "Salix x fragilis"),             # hybrid preserved
    ("Salix x fragilis", "Salix x fragilis"),
]


# (input, expected) for clean_kingdom -- see _LEGACY_KINGDOM_TOKENS for the rationale.
KINGDOM_CASES = [
    # kingdoms the pipeline keeps
    ("Fungi", "Fungi"), ("Plantae", "Plantae"), ("Animalia", "Animalia"),
    ("Chromista", "Chromista"), ("Protozoa", "Protozoa"),
    # BOLD junk: 6493 sequences on the 2026-08-25 run carried "Mixture"
    ("Mixture", None), ("mixture", None), ("Unspecified", None), ("unspecified", None),
    # legacy grouping: 5813 BOLD sequences, members belong under Chromista/Plantae
    ("Protista", None), ("protista", None), ("Protoctista", None),
    # domains are not kingdoms (already covered by is_domain, checked again here)
    ("Eukaryota", None), ("Bacteria", None),
    # ordinary placeholders still rejected
    ("unidentified", None), ("NA", None), ("", None), (None, None),
]

LEGACY_CASES = [
    ("Protista", True), ("protista", True), ("Protoctista", True),
    ("Chromista", False), ("Protozoa", False), ("Fungi", False),
    ("", False), (None, False),
]


def run():
    failures = []
    for value, expected in RANK_CASES:
        got = clean_rank(value)
        if got != expected:
            failures.append(("clean_rank", value, expected, got))
    for value, expected in SPECIES_CASES:
        got = clean_species(value)
        if got != expected:
            failures.append(("clean_species", value, expected, got))
    for value, expected in DOMAIN_CASES:
        got = is_domain(value)
        if got != expected:
            failures.append(("is_domain", value, expected, got))
    for value, expected in KINGDOM_CASES:
        got = clean_kingdom(value)
        if got != expected:
            failures.append(("clean_kingdom", value, expected, got))
    for value, expected in LEGACY_CASES:
        got = is_legacy_kingdom(value)
        if got != expected:
            failures.append(("is_legacy_kingdom", value, expected, got))
    # every kingdom clean_kingdom lets through must be one we keep
    for k in ("Fungi", "Plantae", "Animalia", "Chromista", "Protozoa"):
        if clean_kingdom(k) not in EUKARYOTE_KINGDOMS:
            failures.append(("whitelist", k, "in EUKARYOTE_KINGDOMS", clean_kingdom(k)))

    total = (len(RANK_CASES) + len(SPECIES_CASES) + len(DOMAIN_CASES)
             + len(KINGDOM_CASES) + len(LEGACY_CASES) + 5)
    if failures:
        print(f"FAILED {len(failures)}/{total}")
        for fn, value, expected, got in failures:
            print(f"  {fn}({value!r}) -> {got!r}  (expected {expected!r})")
        return 1
    print(f"PASSED {total}/{total}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
