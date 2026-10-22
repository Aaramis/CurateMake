"""Test cases for ModificationClassifier, focused on mid-lineage holes.

Run: python workflow/scripts/Sativa/test_classify_modification.py

Most lineages in this database have no `class` rank at all. A blank rank is therefore
NOT always missing information:
  * blank WITH populated ranks after it  -> asserts "this lineage has no such rank"
  * blank with nothing after it          -> abstention, the lineage just stops there
The classifier used to compact the blanks away before comparing, which shifted every
deeper rank one slot left, so merely FILLING a hole was reported as a SUBSTITUTION of
everything below it, and dropping `class` alone was reported as a SUBSTITUTION instead
of a DELETION. count_deleted_levels was already positional, so the two disagreed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from integrate_sativa import ModificationClassifier  # noqa: E402

C = ModificationClassifier.classify_modification
D = ModificationClassifier.count_deleted_levels

CASES = [
    # --- the hole is asserted: changing that assertion is a substitution ---------------
    ("trou comble au class",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA",
     "Fungi;Ascomycota;Sordariomycetes;Diaporthales;Valsaceae;Valsa;NA", "SUBSTITUTION"),
    ("class retire, reste identique",
     "Fungi;Ascomycota;Sordariomycetes;Diaporthales;Valsaceae;NA;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA", "SUBSTITUTION"),
    # --- trailing blanks are abstentions: filling them completes ----------------------
    ("genre ajoute apres un NA final",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;NA", "COMPLETION"),
    ("lignee courte prolongee",
     "Fungi;Ascomycota;NA", "Fungi;Ascomycota;NA;Diaporthales", "COMPLETION"),
    ("espece ajoutee",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;Valsa mali", "COMPLETION"),
    # --- real renames -----------------------------------------------------------------
    ("changement de famille",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Diaporthaceae;NA;NA", "SUBSTITUTION"),
    ("changement de regne",
     "Plantae;Tracheophyta;NA;Asterales;NA;NA;NA",
     "Fungi;Basidiomycota;NA", "SUBSTITUTION"),
    # --- removals ---------------------------------------------------------------------
    ("genre supprime",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA", "DELETION_SIMPLE"),
    ("famille et genre supprimes",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;NA",
     "Fungi;Ascomycota;NA;Diaporthales;NA;NA;NA", "DELETION_MULTIPLE"),
]

passed = failed = 0
for name, orig, prop, expected in CASES:
    got = C(orig, prop).value
    if got == expected:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL  {name}: got {got}, expected {expected}")

# count_deleted_levels is positional and must stay consistent with the classifier
CONSISTENCY = [
    ("genre supprime",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA", ["genus"]),
    ("famille et genre supprimes",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;NA",
     "Fungi;Ascomycota;NA;Diaporthales;NA;NA;NA", ["family", "genus"]),
    ("trou au class des deux cotes: rien de supprime",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;Valsa;NA", []),
]
for name, orig, prop, expected in CONSISTENCY:
    got = D(orig, prop)
    if got == expected:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL  deleted_levels/{name}: got {got}, expected {expected}")

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
