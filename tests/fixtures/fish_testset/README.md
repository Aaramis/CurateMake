# Fish test set for the CurateMake CoL matcher

Raw identifications from a paired study of fish surveyed with baited video (BRUVS) and eDNA metabarcoding
(New Caledonia), with the expected harmonisation against **Catalogue of Life COL26.10**
(ChecklistBank key `316595`, doi:10.48580/dgzh4). Always pin this release: other releases can give
other answers.

## Files

- `raw_records.tsv`: 1,439 records, one per distinct raw identification. These are the input to the matcher.
  - `source = BRUVS` (472 records): video IDs. Lineage from FishBase, e.g. "Perciformes" in the broad sense, Scaridae, Serranidae s.l.
  - `source = eDNA` (967 records): MOTU assignments (ecotag, teleo 12S). NCBI lineage, e.g. class "Actinopteri",
    orders such as Labriformes or Spariformes, an NCBI placeholder `Scopelarchus sp. ausc01` and a trinomial.
  - Columns: `kingdom` (always Animalia), `class`, `order`, `family`, `genus`, `species`, plus `verbatim` and
    `verbatim_rank` (the rank of the identification). Ranks finer than `verbatim_rank` are empty.
- `expected.tsv`: the expected output per `record_id`: `expected_class` … `expected_species` and
  `expected_status` (`Correct` / `Check`), with a `note` on the 469 records that exercise a rule.
- `edge_cases.tsv`: 10 synthetic cases (typos, homonyms, synonym, provisional name).

## Rules the expected output follows

1. **Synonyms**: replaced by the accepted name, with the accepted lineage (38 records, e.g. *Carangoides ferdau* →
   *Ferdauia ferdau*, *Chromis atripes* → *Pycnochromis atripes*, *Plectroglyphidodon lacrymatus* → *Stegastes lacrymatus*).
2. **Ranks CoL leaves unplaced stay empty.** Never fill them from the source: 364 records whose family is
   *incertae sedis* in CoL (Labridae, Pomacentridae, Lutjanidae…) must have an **empty order**. Filling them
   gives the same species two orders across sources (40 of the 102 species shared by both sources).
3. **The unresolved rank itself keeps its raw name.** 37 eDNA records assigned to an NCBI-only order (Labriformes,
   Spariformes…) match only higher up; their order keeps the raw value.
4. **Rank variants**: NCBI class "Actinopteri" is a superclass in CoL, so `expected_class` is empty for those records
   rather than Actinopteri.
5. **Provisional names** (`sp.`, strain codes) are resolved at genus, never at species.
6. **Ambiguous names → Check**: *Moolgarda* (3 records) is both a provisionally accepted genus and an ambiguous
   synonym of *Crenimugil* in COL26.10. Either name is acceptable (`Moolgarda|Crenimugil`), but the status
   must be `Check`.

## Status of CurateMake `restructure-preview` @ ef149e3 on this set (22 Oct 2026)

- 1,439/1,439 matched. Species, genus, family, order and class are identical to `expected.tsv`, except *Moolgarda*:
  it gives `Moolgarda` (fine) but not `Check`, because the provisionally accepted usage is picked and the ambiguous synonym ignored.
- `review = ambiguous` on 6 records, 4 of them noise:
  - *Nectamia* ×2 and *Notoscopelus*: the same taxon returned at genus and subgenus rank.
  - *Centrophorus* (vs a katydid synonym) and *Pentacerotidae* (vs a sea-star family): the source family/order already settled the choice.

  Suggestion: ignore candidates whose rank differs from the query (subgenus), and don't flag when the family/order
  preference leaves a single candidate.
- `edge_cases.tsv`: all as expected, except case 9 (*Moolgarda*, `Check` missing).

## Provenance

Expected output = the CurateMake-method R implementation used in the paper (`harmonise_taxonomy_col.R`),
cross-checked with CurateMake. Matching rules are hierarchical (species → class), kingdom + source family as
context, fuzzy ≤ 2 edits at the same rank, synonyms → accepted, and unresolved finer ranks keep raw names.
