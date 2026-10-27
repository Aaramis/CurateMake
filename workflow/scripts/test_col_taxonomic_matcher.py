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


def hit(name, rank, classification, status="accepted", accepted_name=None,
        accepted_rank=None):
    usage = {"status": status, "name": {"scientificName": name, "rank": rank}, "id": "u1"}
    if "synonym" in status:
        # An unresolved synonym is modelled as an `accepted` block without a usable name,
        # which is one of the two shapes CoL returns for one.
        usage["accepted"] = ({"id": "u2", "name": {"scientificName": accepted_name,
                                                   "rank": accepted_rank or rank}}
                             if accepted_name else {"id": "u2"})
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

print("\nFaute dans le genre lui-meme")

RUSACEAE = lineage(("kingdom", "Fungi"), ("family", "Russulaceae"))
WITH_GENUS = RUSACEAE + lineage(("genus", "Russula"))
SRC = dict(kingdom="Fungi", family="Russulaceae", genus="Rusula", species="Rusula azurea")

# the genus is misspelled AND CoL knows the species under the corrected genus
m = col.TaxonomyMatcher(FakeClient(fuzzy={
    "Rusula": [hit("Russula", "genus", RUSACEAE)],
    "Rusula azurea": [hit("Russula azurea", "species", WITH_GENUS)]}))
out = m.match_taxonomy(record(**SRC))
check("un genre mal orthographie est corrige", out.col_taxonomy.get("genus") == "Russula")
check("le binome entier est recupere", out.col_taxonomy.get("species") == "Russula azurea")

# the genus is misspelled and CoL does not know the species: the stale epithet is dropped
m = col.TaxonomyMatcher(FakeClient(fuzzy={"Rusula": [hit("Russula", "genus", RUSACEAE)]}))
out = m.match_taxonomy(record(**SRC))
check("le genre est corrige meme sans l espece", out.col_taxonomy.get("genus") == "Russula")
check("l epithete perime n est pas recopie sous le nouveau genre",
      "species" not in out.col_taxonomy)

# a near genus in ANOTHER family must not be accepted
m = col.TaxonomyMatcher(FakeClient(fuzzy={
    "Rusula": [hit("Russula", "genus", lineage(("kingdom", "Fungi"), ("family", "Amanitaceae")))]}))
out = m.match_taxonomy(record(**SRC))
check("un genre proche d une autre famille est refuse", not out.success)

# the ordinary case must not regress
m = col.TaxonomyMatcher(FakeClient(exact={"Russula": [hit("Russula", "genus", RUSACEAE)]}))
out = m.match_taxonomy(record(kingdom="Fungi", family="Russulaceae",
                              genus="Russula", species="Russula azurea"))
check("un binome coherent reste herite", out.col_taxonomy.get("species") == "Russula azurea")

check("la garde accepte un binome coherent",
      col.TaxonomyMatcher._binomial_agrees("Russula azurea", "Russula"))
check("la garde refuse un binome contredit",
      not col.TaxonomyMatcher._binomial_agrees("Carangoides ferdau", "Ferdauia"))
check("la garde laisse passer ce qu elle ne peut pas comparer",
      col.TaxonomyMatcher._binomial_agrees("Russula azurea", None)
      and col.TaxonomyMatcher._binomial_agrees(None, "Russula"))


print("\nMarqueurs d absence: le piege de la sous-chaine")

# Historical bug (fixed in the August rewrite, locked in here): the placeholder test was
# a SUBSTRING test against 'na', so every name containing those two letters was nulled.
# For ITS, Cortinarius alone is a large share of the data.
for name in ("Cortinarius", "Pinaceae", "Gentianaceae", "Solanaceae", "Ranunculaceae",
             "Anaptychia", "Nannizzia", "Hypomyces chrysospermus"):
    check(f"{name!r} n est pas pris pour une valeur absente",
          col.TaxonomyMatcher._clean(name) == name)

for blank in ("", "   ", "NA", "na", "N/A", "null", "None", "unknown", "UNCLASSIFIED",
              "unidentified"):
    check(f"{blank!r} est bien traite comme absent",
          col.TaxonomyMatcher._clean(blank) is None)

check("une citation d auteur est retiree",
      col.TaxonomyMatcher._clean("Russula (Schaeff.) Pers.") == "Russula")
check("les tirets bas redeviennent des espaces",
      col.TaxonomyMatcher._clean("Russula_azurea") == "Russula azurea")


print("\nNoms provisoires")

for prov in ("Russula sp.", "Russula sp", "Russula cf. azurea", "Russula aff. azurea",
             "uncultured fungus", "unidentified Basidiomycota", "Russula F123",
             "environmental sample"):
    check(f"{prov!r} est refuse comme nom exploitable",
          col.TaxonomyMatcher._is_provisional(prov))

for real in ("Russula azurea", "Cortinarius", "Thalassoma hardwicke"):
    check(f"{real!r} reste exploitable", not col.TaxonomyMatcher._is_provisional(real))


print("\nSynonymes resolus vers le nom accepte")

# CoL's EXACT search returns a synonym's accepted usage in the same response. The
# January matcher recorded the synonym instead, which is the opposite of harmonising.
m = col.TaxonomyMatcher(FakeClient(exact={"Russula olivascens": [
    hit("Russula olivascens", "species", RUS, status="synonym",
        accepted_name="Russula alutacea")]}))
out = m.match_taxonomy(record(kingdom="Fungi", genus="Russula", species="Russula olivascens"))
check("un synonyme est remplace par le nom accepte",
      out.col_taxonomy.get("species") == "Russula alutacea")
check("l identifiant conserve est celui de l usage accepte", out.col_id == "u2")


print("\nHomonymes entre regnes")

# Morus is a plant (mulberry) and a bird (gannet). The 2026-07-27 probe showed a bare
# search returning either one; the kingdom of the source has to decide.
MORUS = [hit("Morus", "genus", lineage(("kingdom", "Animalia"), ("family", "Sulidae"))),
         hit("Morus", "genus", lineage(("kingdom", "Plantae"), ("family", "Moraceae")))]
m = col.TaxonomyMatcher(FakeClient(exact={"Morus": MORUS}))
plant = m.match_taxonomy(record(kingdom="Plantae", genus="Morus"))
bird = m.match_taxonomy(record(kingdom="Animalia", genus="Morus"))
check("le Morus vegetal est choisi pour une source Plantae",
      plant.col_taxonomy.get("family") == "Moraceae")
check("le Morus oiseau est choisi pour une source Animalia",
      bird.col_taxonomy.get("family") == "Sulidae")


print("\nHomonymes a l interieur d un meme regne")

# Chromis is a damselfish (Pomacentridae, accepted) and a cichlid (Cichlidae, only
# provisionally accepted). The kingdom cannot separate them, so the family the source
# declares has to. Reported by Edi Prifti on fish data.
CHROMIS = [
    hit("Chromis", "genus", lineage(("kingdom", "Animalia"), ("order", "Perciformes"),
                                    ("family", "Pomacentridae"))),
    hit("Chromis", "genus", lineage(("kingdom", "Animalia"), ("order", "Cichliformes"),
                                    ("family", "Cichlidae")), status="provisionally accepted"),
]
m = col.TaxonomyMatcher(FakeClient(exact={"Chromis": CHROMIS}))
out = m.match_taxonomy(record(kingdom="Animalia", family="Cichlidae", genus="Chromis"))
check("la famille brute depart age l homonyme", out.col_taxonomy.get("family") == "Cichlidae")
out = m.match_taxonomy(record(kingdom="Animalia", family="Pomacentridae", genus="Chromis"))
check("et dans l autre sens aussi", out.col_taxonomy.get("family") == "Pomacentridae")

out = m.match_taxonomy(record(kingdom="Animalia", genus="Chromis"))
check("sans famille, l ambiguite est signalee", bool(out.review))
check("le signalement nomme les deux lignees",
      out.review and "Pomacentridae" in out.review and "Cichlidae" in out.review)

# the order alone must be enough when the family is missing
out = m.match_taxonomy(record(kingdom="Animalia", order_="Cichliformes", genus="Chromis"))
check("l ordre brut depart age aussi", out.col_taxonomy.get("family") == "Cichlidae")

# a preference must never turn a usable match into no match
out = m.match_taxonomy(record(kingdom="Animalia", family="Famille_inexistante", genus="Chromis"))
check("une famille absente de CoL ne fait pas perdre le match", out.success)


print("\nSynonymes ambigus")

# A synonym CoL does not resolve: the name has to be kept, but the record must say so.
UNRES = [hit("Moolgarda", "genus", lineage(("kingdom", "Animalia"), ("family", "Mugilidae")),
             status="synonym", accepted_name=None)]
m = col.TaxonomyMatcher(FakeClient(exact={"Moolgarda": UNRES}))
out = m.match_taxonomy(record(kingdom="Animalia", family="Mugilidae", genus="Moolgarda"))
check("un synonyme non resolu est signale", bool(out.review))
check("le motif dit pourquoi", out.review and "synonym" in out.review)

# a clean, resolved match carries no flag at all
m = col.TaxonomyMatcher(FakeClient(exact={"Russula": [hit("Russula", "genus", RUSACEAE)]}))
out = m.match_taxonomy(record(kingdom="Fungi", family="Russulaceae", genus="Russula"))
check("un match propre n est pas signale", out.review is None)


print("\nHomonyme vivant contre synonyme")

# CoL has several synonym statuses, and only the plain one may be followed. All of them
# have to be RECOGNISED as synonyms, though, or an ambiguous synonym looks like a rival
# taxon of the same name. Reported by Edi Prifti on the fish test set.
for st in ("synonym", "ambiguous synonym", "misapplied name"):
    check(f"{st!r} compte comme un synonyme", col.TaxonomyMatcher._is_any_synonym(st))
for st in ("accepted", "provisionally accepted", None):
    check(f"{st!r} ne compte pas comme un synonyme",
          not col.TaxonomyMatcher._is_any_synonym(st))

# One accepted usage surrounded by synonyms of the same name is settled: nomenclature has
# decided. CoL answers "Pentacerotidae" with the accepted fish family and a synonym of the
# same name pointing at a sea-star family.
PENTA = [hit("Pentacerotidae", "family", lineage(("kingdom", "Animalia"),
                                                 ("order", "Perciformes"))),
         hit("Pentacerotidae", "family", lineage(("kingdom", "Animalia")),
             status="synonym", accepted_name="Stauranderasteridae")]
m = col.TaxonomyMatcher(FakeClient(exact={"Pentacerotidae": PENTA}))
out = m.match_taxonomy(record(kingdom="Animalia", family="Pentacerotidae"))
check("un usage accepte entoure de synonymes n est pas signale", out.review is None)
check("et c est bien lui qui est retenu",
      out.col_taxonomy.get("family") == "Pentacerotidae")

# A merely provisional usage against an ambiguous synonym of the same name, with the name
# accepted NOWHERE: CoL has no settled opinion, so the pick must be flagged. This is
# COL26.10 on "Moolgarda", and the rule Edi Prifti asked for.
MOOL = [hit("Moolgarda", "genus", lineage(("kingdom", "Animalia"), ("order", "Mugiliformes"),
                                          ("family", "Mugilidae")),
            status="provisionally accepted"),
        hit("Moolgarda", "genus", lineage(("kingdom", "Animalia"), ("order", "Mugiliformes"),
                                          ("family", "Mugilidae")),
            status="ambiguous synonym", accepted_name="Crenimugil")]
m = col.TaxonomyMatcher(FakeClient(exact={"Moolgarda": MOOL}))
out = m.match_taxonomy(record(kingdom="Animalia", family="Mugilidae", genus="Moolgarda"))
check("provisoire + synonyme ambigu, nom accepte nulle part : signale", bool(out.review))
check("le motif dit qu aucun usage n est accepte",
      out.review and "no accepted usage" in out.review)
check("le nom provisoire est quand meme retenu",
      out.col_taxonomy.get("genus") == "Moolgarda")

# Same shape, but the name IS accepted elsewhere at the same rank: it is a homonym, and
# the family of the source resolves it. COL26.10 on "Chromis", given family Cichlidae.
CHROM = [hit("Chromis", "genus", lineage(("kingdom", "Animalia"), ("order", "Perciformes"),
                                         ("family", "Pomacentridae"))),
         hit("Chromis", "genus", lineage(("kingdom", "Animalia"), ("order", "Cichliformes"),
                                         ("family", "Cichlidae")),
             status="provisionally accepted"),
         hit("Chromis", "genus", lineage(("kingdom", "Animalia"), ("order", "Cichliformes"),
                                         ("family", "Cichlidae")),
             status="ambiguous synonym", accepted_name="Tilapia")]
m = col.TaxonomyMatcher(FakeClient(exact={"Chromis": CHROM}))
out = m.match_taxonomy(record(kingdom="Animalia", family="Cichlidae", genus="Chromis"))
check("le meme cas, mais nom accepte ailleurs : pas signale", out.review is None)
check("et la famille brute a bien tranche", out.col_taxonomy.get("family") == "Cichlidae")
out = m.match_taxonomy(record(kingdom="Animalia", family="Pomacentridae", genus="Chromis"))
check("et l autre branche reste accessible",
      out.col_taxonomy.get("family") == "Pomacentridae" and out.review is None)


print("\nBinome coherent avec le genre")

# CoL often resolves a species synonym to an INFRASPECIFIC accepted name. We have no slot
# below species, so the name must be filed under species, and _result_rank_key has to say
# the same or _pick drops the candidate from a species query. COL26.10 answers
# "Anemone patens" with a plain synonym of "Pulsatilla halleri subsp. halleri".
SUB = [hit("Anemone patens", "species",
           lineage(("kingdom", "Plantae"), ("family", "Ranunculaceae"),
                   ("genus", "Pulsatilla")),
           status="synonym", accepted_name="Pulsatilla halleri subsp. halleri",
           accepted_rank="subspecies")]
check("un nom accepte infraspecifique compte comme une espece",
      col.TaxonomyMatcher._result_rank_key(SUB[0]) == "species")
m = col.TaxonomyMatcher(FakeClient(exact={"Anemone patens": SUB}))
out = m.match_taxonomy(record(kingdom="Plantae", family="Ranunculaceae",
                              genus="Anemone", species="Anemone patens"))
check("l espece est tiree du nom infraspecifique",
      out.col_taxonomy.get("species") == "Pulsatilla halleri")
check("et le genre reste celui de la lignee",
      out.col_taxonomy.get("genus") == "Pulsatilla")

# A plain synonym can be followed; an ambiguous one cannot. With nothing accepted, the
# plain one must win, or the only candidate able to give a species is thrown away.
BOTH = [hit("Anemone patens", "species",
            lineage(("kingdom", "Plantae"), ("family", "Ranunculaceae"),
                    ("genus", "Pulsatilla")),
            status="ambiguous synonym", accepted_name="Pulsatilla patens subsp. patens",
            accepted_rank="subspecies"),
        hit("Anemone patens", "species",
            lineage(("kingdom", "Plantae"), ("family", "Ranunculaceae"),
                    ("genus", "Pulsatilla")),
            status="synonym", accepted_name="Pulsatilla halleri subsp. halleri",
            accepted_rank="subspecies")]
m = col.TaxonomyMatcher(FakeClient(exact={"Anemone patens": BOTH}))
out = m.match_taxonomy(record(kingdom="Plantae", family="Ranunculaceae",
                              genus="Anemone", species="Anemone patens"))
check("le synonyme simple prime sur le synonyme ambigu",
      out.col_taxonomy.get("species") == "Pulsatilla halleri")
check("et le desaccord est signale", bool(out.review))

# An ambiguous synonym we refuse to follow, with the accepted lineage around it: the row
# would otherwise carry genus Synotis and species "Senecio aureus".
AMB = [hit("Senecio aureus", "species",
           lineage(("kingdom", "Plantae"), ("family", "Asteraceae"), ("genus", "Synotis")),
           status="ambiguous synonym", accepted_name="Synotis cappa")]
m = col.TaxonomyMatcher(FakeClient(exact={"Senecio aureus": AMB}))
out = m.match_taxonomy(record(kingdom="Plantae", family="Asteraceae",
                              genus="Senecio", species="Senecio aureus"))
check("un binome contredisant le genre n est jamais ecrit",
      out.col_taxonomy.get("species") in (None, ""))
check("le genre accepte est conserve", out.col_taxonomy.get("genus") == "Synotis")

# the ordinary case must be untouched: species synonym -> accepted binomial, both coherent
OK = [hit("Lavatera trimestris", "species",
          lineage(("kingdom", "Plantae"), ("family", "Malvaceae"), ("genus", "Malva")),
          status="synonym", accepted_name="Malva trimestris")]
m = col.TaxonomyMatcher(FakeClient(exact={"Lavatera trimestris": OK}))
out = m.match_taxonomy(record(kingdom="Plantae", family="Malvaceae",
                              genus="Lavatera", species="Lavatera trimestris"))
check("un synonyme resolu normalement garde son binome",
      out.col_taxonomy.get("species") == "Malva trimestris")
check("sans signalement", out.review is None)


print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
