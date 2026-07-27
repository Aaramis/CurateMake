# SATIVA modifications in CurateMake — leave-one-out via EPA-ng

**Status: integrated as the default (2026-07-27).**

This document justifies the modification made to the SATIVA copy bundled with the pipeline
(`workflow/scripts/Sativa/sativa/`).

## 1. Summary

SATIVA's mislabel-detection step (the **leave-one-out** test) is now performed with
**EPA-ng** (Barbera et al. 2019, *Syst. Biol.*) instead of the historical **RAxML `-f O`**
implementation. Everything else in SATIVA is unchanged: reference-tree construction
(`mfresolv`), taxonomic labelling of the branches, the `assign_taxonomy_maxsum`
aggregation, the `check_seq_tax_labels` decision rule and the confidence cutoff all remain
**native SATIVA code**. Only the *placement engine* of the leave-one-out changes.

## 2. Motivation

RAxML's leave-one-out (`raxmlHPC -f O`) **does not scale**: beyond ~1500–2000 sequences per
taxonomic group it hits the 6 h timeout **without producing any result** (observed on
Trichocomaceae 2575 seqs, Calloriaceae 9422, Acacia 1881, etc.). Since SATIVA accounts for
~80 % of the pipeline's compute time, this is the main bottleneck.

EPA-ng re-implements phylogenetic placement (SIMD + multithreaded) and is 1–2 orders of
magnitude faster, while emitting a standard `.jplace`.

## 3. How it works

Files:
- **`epac/epang_l1o.py`** (new) — the `run_epang_l1o()` helper.
- **`sativa.py`** — `run_leave_seq_out_test()` calls EPA-ng by default (see §5).
- **`../../envs/sativa.yaml`** — adds the `epa-ng=0.3.8` dependency.

Helper algorithm (leave-one-out reformulated as **k-fold**, K=5 by default):
1. Partition the leaves into K folds. For each fold: prune its leaves from the reference
   tree and alignment (`ete3.prune`, branch lengths preserved), then place the fold's
   sequences with EPA-ng.
2. **Keep ALL placements** (`--filter-acc-lwr 0.99999 --filter-max 100000`). EPA-ng's
   default (`--filter-max=7`) truncated ~80 % of the LWR mass and corrupted the confidence —
   this setting is essential.
3. **Remap EPA-ng edge numbers** to SATIVA's `refjson` `B=` numbering via leaf-bipartition
   matching (exact bijection). The placements are then directly consumable by `classify_seq`
   through `bid_taxonomy_map`.
4. Model: EPA-ng reads the `RAxML_info.mfresolv_*` produced by the tree-building step
   (GTRGAMMA); falls back to `GTR+G` if absent.

`run_epang_l1o()` returns a placement list in the exact format of
`EpaJsonParser.get_placement()`; the remainder of `run_leave_seq_out_test()` is unchanged.

## 4. Validation

- **Strict concordance (Acaulium group, 106 seqs)**: per-sequence leave-one-out, running the
  real SATIVA code on the EPA-ng placements → **100 % recall, perfect 20/20 rank agreement**
  vs SATIVA-RAxML (EPA-ng is a strict superset).
- **Full-pipeline benchmark (mini dataset `data/raw`, SQLite, 66 groups, T=1)**:

  | | Native SATIVA (RAxML `-f O`) | SATIVA + EPA-ng (default) | speed-up |
  |---|---|---|---|
  | SATIVA time (Σ runtime) | 11,127 s (3.1 h) | **1,294 s (21.6 min)** | **~8.6×** |
  | Total pipeline time | ~3.9 h | **~1.2 h** | **~3.3×** |
  | Mislabels detected | 1,412 | 1,486 | — |
  | Concordance (k-fold) | — | recall 84 % · precision 80 % · **rank 98 %** | — |

  Example: Aspergillaceae ITS (1047 seqs) 2094 s → 207 s (T=1) / 65 s (T=8).

The k-fold scheme is an **approximation** of the strict leave-one-out: aggregate concordance
is 84 % (recall) / 80 % (precision) with 98 % rank agreement on the shared mislabels. EPA-ng
flags slightly more mislabels than RAxML (1486 vs 1412).

## 5. Configuration

- **Default**: EPA-ng (no variable to set).
- **Historical RAxML fallback**: `export CURATEMAKE_SATIVA_L1O=raxml` (discouraged — does not
  scale).
- **Number of folds**: `export CURATEMAKE_EPANG_FOLDS=5` (default). More folds ⇒ closer to the
  strict leave-one-out, but slower.
- **EPA-ng binary**: resolved from `PATH` (provided by `sativa.yaml`); overridable via
  `CURATEMAKE_EPANG_BIN`.

## 6. Known limitations / future work

- **Small groups** (< ~300 seqs): EPA-ng can be slightly slower than RAxML (process start-up
  cost), but the absolute time there is negligible; the net gain comes from the large groups.
- SATIVA's **final EPA confirmation** step (`run_epa_once`, run only on the flagged suspects)
  still uses RAxML — could be switched to EPA-ng later.
- Model GTRGAMMA (EPA-ng, from `RAxML_info`) vs GTRCAT (RAxML `-f O`): a minor difference is
  expected, with no notable effect on the decisions (98 % rank agreement).

---

*Refs: Barbera et al. 2019, EPA-ng, Syst. Biol. 68(2):365–369. Kozlov et al. 2016, SATIVA.*
