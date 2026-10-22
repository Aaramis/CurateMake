"""Test cases for the CoL matcher's rank handling.

Run: python workflow/scripts/test_col_taxonomic_matcher.py

Covers the three guards added after Edi Prifti's review of this branch, run against fish
data with CoL release COL26.10:

  * A rank CoL leaves unresolved must stay EMPTY. For an incertae sedis family CoL returns
    no order, and the matcher used to reinstate the source's order, so the same species
    came out as Perciformes from a FishBase-derived source and Labriformes from an
    NCBI-derived one: 40 of the 102 species seen by both methods disagreed inside the
    layer whose whole purpose is to make them comparable. Only ranks MORE SPECIFIC than
    the matched one may be inherited.
  * A result must sit at the rank that was queried. The 2026-07-27 probe (15.CoL_logic_test)
    showed a misspelling plus a kingdom hint resolving to the kingdom itself and being
    accepted as a match.
  * A FUZZY hit must be a near miss, not any suggestion CoL ranks first: without a bound,
    `niger` returns birds and `grisescens` returns moths.

Deliberately offline: the fake client answers from a table, so the test never touches
ChecklistBank and cannot drift with a CoL release.
"""
import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "col_taxonomic_matcher", Path(__file__).with_name("col_taxonomic_matcher.py"))
col = importlib.util.module_from_spec(_SPEC)
sys.modules["col_taxonomic_matcher"] = col
_SPEC.loader.exec_module(col)

passed = failed = 0


def check(name, ok):
    global passed, failed
    passed += bool(ok)
    failed += (not ok)
    print(f"  {name:<58} {'OK' if ok else 'ECHEC'}")


def lineage(*pairs):
    """CoL returns the lineage ABOVE the usage as a list of {rank, name}."""
    return [{"rank": r, "name": n} for r, n in pairs]


def hit(name, rank, classification, status="accepted", accepted_name=None):
    usage = {"status": status, "name": {"scientificName": name, "rank": rank}, "id": "u1"}
    if status == "synonym":
        usage["accepted"] = {"id": "u2",
                             "name": {"scientificName": accepted_name, "rank": rank}}
    return {"usage": usage, "classification": classification}


class FakeClient:
    """Answers EXACT/FUZZY from a table keyed by the queried name."""

    def __init__(self, exact=None, fuzzy=None):
        self.exact, self.fuzzy = exact or {}, fuzzy or {}

    def search(self, name, match_type="EXACT", limit=20):
        return (self.exact if match_type == "EXACT" else self.fuzzy).get(name, [])


def record(**kw):
    return col.TaxonomyRecord(tax_id=1, **kw)


print("\nRangs non resolus par CoL")

# Labridae is incertae sedis in CoL: the lineage stops at class, no order.
FISH = lineage(("kingdom", "Animalia"), ("phylum", "Chordata"),
               ("class", "Actinopterygii"), ("family", "Labridae"), ("genus", "Thalassoma"))
FISH_HIT = {"Thalassoma hardwicke": [hit("Thalassoma hardwicke", "species", FISH)]}
COMMON = dict(kingdom="Animalia", phylum="Chordata", class_="Actinopterygii",
              family="Labridae", genus="Thalassoma", species="Thalassoma hardwicke")

m = col.TaxonomyMatcher(FakeClient(exact=FISH_HIT))
out = m.match_taxonomy(record(order_="Perciformes", **COMMON))
check("l ordre absent de CoL n est pas repris de la source", "order_" not in out.col_taxonomy)
check("ce que CoL renvoie est conserve",
      out.col_taxonomy.get("family") == "Labridae"
      and out.col_taxonomy.get("species") == "Thalassoma hardwicke")

video = m.match_taxonomy(record(order_="Perciformes", **COMMON))
edna = m.match_taxonomy(record(order_="Labriformes", **COMMON))
check("deux sources aux ordres divergents convergent", video.col_taxonomy == edna.col_taxonomy)

# A match at family must still keep the source's genus and species.
FAM = lineage(("kingdom", "Fungi"), ("phylum", "Basidiomycota"),
              ("order", "Russulales"), ("family", "Russulaceae"))
m = col.TaxonomyMatcher(FakeClient(exact={"Russulaceae": [hit("Russulaceae", "family", FAM)]}))
out = m.match_taxonomy(record(kingdom="Fungi", family="Russulaceae",
                              genus="Russula", species="Russula azurea"))
check("les rangs plus specifiques restent herites",
      out.col_taxonomy.get("genus") == "Russula"
      and out.col_taxonomy.get("species") == "Russula azurea")
check("l ordre vient de CoL et non de la source", out.col_taxonomy.get("order_") == "Russulales")


print("\nGarde de rang")

m = col.TaxonomyMatcher(FakeClient(exact={"Aspergilus niger": [hit("Fungi", "kingdom", [])]}))
out = m.match_taxonomy(record(kingdom="Fungi", species="Aspergilus niger"))
check("un resultat au rang du regne est refuse pour une espece", not out.success)

SUB = lineage(("kingdom", "Fungi"), ("phylum", "Basidiomycota"),
              ("family", "Russulaceae"), ("genus", "Russula"))
m = col.TaxonomyMatcher(FakeClient(exact={"Russula azurea": [
    hit("Russula azurea var. minor", "subspecies", SUB)]}))
out = m.match_taxonomy(record(kingdom="Fungi", species="Russula azurea"))
check("un rang non cartographie est garde pour sa lignee",
      out.success and out.col_taxonomy.get("genus") == "Russula")


print("\nGarde de distance sur le repli fuzzy")

RUS = lineage(("kingdom", "Fungi"), ("family", "Russulaceae"), ("genus", "Russula"))
m = col.TaxonomyMatcher(FakeClient(fuzzy={"Russula azurae": [hit("Russula azurea", "species", RUS)]}))
out = m.match_taxonomy(record(kingdom="Fungi", genus="Russula", species="Russula azurae"))
check("une faute de frappe a 1 edition est corrigee",
      out.success and out.col_taxonomy.get("species") == "Russula azurea")

m = col.TaxonomyMatcher(FakeClient(fuzzy={"Russula azurae": [hit("Russula brevipes", "species", RUS)]}))
out = m.match_taxonomy(record(kingdom="Fungi", genus="Russula", species="Russula azurae"))
check("une suggestion lointaine est refusee",
      out.col_taxonomy is None or out.col_taxonomy.get("species") != "Russula brevipes")

m = col.TaxonomyMatcher(FakeClient(fuzzy={"Russula olivascen": [
    hit("Russula olivascens", "species", RUS, status="synonym", accepted_name="Russula alutacea")]}))
out = m.match_taxonomy(record(kingdom="Fungi", genus="Russula", species="Russula olivascen"))
check("la distance se mesure sur le nom matche, pas sur l accepte",
      out.success and out.col_taxonomy.get("species") == "Russula alutacea")


print("\nFonction de distance")

for a, b, expected in [("Russula", "Russula", 0), ("Russula", "Rusula", 1),
                       ("Russula", "Rusulla", 2)]:
    check(f"distance({a!r}, {b!r}) == {expected}", col._edit_distance(a, b, cap=3) == expected)
check("un nom tres different depasse le plafond", col._edit_distance("Cortinarius", "Pinaceae", cap=3) > 3)
check("insensible a la casse, tolerante aux valeurs nulles",
      col._edit_distance("Russula", "russula") == 0
      and col._edit_distance(None, "Russula", cap=2) > 2
      and col._edit_distance(None, None) == 0)

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
