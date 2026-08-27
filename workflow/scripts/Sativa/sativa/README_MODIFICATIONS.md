# SATIVA modifications in CurateMake — leave-one-out via EPA-ng

**Status: integrated as the default (2026-07-27).**

This document justifies the modification made to the SATIVA copy bundled with the pipeline
(`workflow/scripts/Sativa/sativa/`).

## 1. Summary

**Both placement steps** of SATIVA now run on **EPA-ng** (Barbera et al. 2019, *Syst.
Biol.*) instead of RAxML: the **leave-one-out** test (`-f O`, pass 1) and the **final
confirmation** of the flagged suspects (`-f v`, pass 2). Everything else in SATIVA is
unchanged: reference-tree construction (`mfresolv`), taxonomic labelling of the branches,
the `assign_taxonomy_maxsum` aggregation, the `check_seq_tax_labels` decision rule and the
confidence cutoff all remain **native SATIVA code** — the files that hold them
(`epac/classify_util.py`, `taxonomy_util.py`, `json_util.py`, `raxml_util.py`, `msa.py`)
are byte-identical to upstream v0.9.3. The whole modification is two hunks in `sativa.py`,
one hunk in `epac/config.py` and the new `epac/epang_l1o.py`.

## 2. Motivation

RAxML's leave-one-out (`raxmlHPC -f O`) **does not scale**: beyond ~1500–2000 sequences per
taxonomic group it hits the 6 h timeout **without producing any result** (observed on
Trichocomaceae 2575 seqs, Calloriaceae 9422, Acacia 1881, etc.). Since SATIVA accounts for
~80 % of the pipeline's compute time, this is the main bottleneck.

EPA-ng re-implements phylogenetic placement (SIMD + multithreaded) and is 1–2 orders of
magnitude faster, while emitting a standard `.jplace`.

## 3. How it works

Files:
- **`epac/epang_l1o.py`** (new) — the `run_epang_l1o()` (pass 1) and `run_epang_final()`
  (pass 2) helpers.
- **`sativa.py`** — `run_leave_seq_out_test()` and `run_epa_once()` both call EPA-ng by
  default (see §5).
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

**Pass 2 — final confirmation (`run_epang_final()`).** SATIVA re-places the flagged
suspects on the mislabel-free tree to confirm them; that step used RAxML `-f v`
(`run_epa_once`). It now runs EPA-ng on the same inputs: reference = the leaves of the
pruned tree, queries = the alignment minus those leaves (i.e. the suspects), same model
(`RAxML_info.mfresolv_*`) and same "keep all placements" filter as pass 1. The `.jplace`
it writes is self-consistent (its own tree + placements), so `EpaJsonParser` consumes it
directly and SATIVA rebuilds `bid_tax_map` from that tree — no edge remapping needed here.
The log line is `[epang-final] N suspects re-placed`.

## 4. Validation

Measured in `18.sativa_speedup/` (see its `RESULTS.md` for the figure, the tables and the
method). Reference: unmodified SATIVA v0.9.3. Subsets of one clustered ITS alignment, 100
to 5402 sequences, on 32 cores. Every timed run on local disk and alone on the machine:
on the project filesystem the I/O dominates and adds several-fold noise, which is what the
earlier numbers in this file suffered from.

**Speed** (total SATIVA runtime, K=5):

| n | upstream, 1 thread | upstream, 8 threads | EPA-ng, 1 thread | EPA-ng, 8 threads | EPA-ng, 16 threads |
|---|---|---|---|---|---|
| 400 | 35.9 s | 26.9 s | 5.5 s | 3.3 s | 3.3 s |
| 1600 | 257.4 s | 214.1 s | 34.8 s | 12.7 s | 11.7 s |
| 5402 | 1536 s | 1157 s | 723 s | 129 s | 87 s |

Threads buy RAxML little here (242 alignment columns leave nothing to split); EPA-ng
parallelises over queries. Above 500 taxa SATIVA switches its own RAxML to GTRCAT, and
above 1000 taxa to a partial-branch heuristic, so the single-thread ratio narrows with size.

**Agreement**, per sequence, cutoff 0.4 both sides:

| n | K=5 recall / precision | K=25 recall / precision |
|---|---|---|
| 800 | 0.833 / 0.851 | 0.958 / 0.939 |
| 1600 | 0.839 / 0.756 | 0.898 / 0.855 |
| 5402 | 0.783 / 0.753 | 0.793 / 0.800 |

To be read against how much SATIVA moves against itself: at n=800, changing its own model
gives 0.917 / 0.898, and its own python-2 build gives 0.875 / 0.857.

**Where the k-fold wins.** With 33 mislabels injected into three clades (k = 1, 2, 3, 5
sequences relabelled with a distant genus), unmodified SATIVA finds 28, the EPA-ng version
33. Every miss happens when at least two sequences share the wrong label: the strict
leave-one-out removes one and the remaining copies confirm the error, while a fold removes
several at once.

## 5. Configuration

- **Default**: EPA-ng for both passes (no variable to set).
- **Historical RAxML fallback**: `export CURATEMAKE_SATIVA_L1O=raxml` — the same switch
  reverts *both* passes to RAxML (`-f O` and `-f v`). Discouraged: does not scale.
- **Number of folds**: `export CURATEMAKE_EPANG_FOLDS=25` (default since 2026-08-18, was 5).
  More folds ⇒ closer to the strict leave-one-out, but slower. K=25 costs about 2.7× the
  placement time of K=5 and lifts agreement at n=800 from 0.875/0.840 to 0.979/0.904, still
  6× to 7× faster than RAxML. It is also the only setting whose calls came out stable across
  every configuration measured. K equal to the number of leaves is the strict leave-one-out; it
  costs 128 s at n=400 against 36 s for the RAxML it replaces, for the same calls as K=25.
- **Reference-tree model**: `sativa.cfg` is upstream's, `raxml_model=auto`. It was briefly
  pinned to `GTRGAMMA` (2026-08-18) so that EPA-ng would read a fitted gamma shape, and that
  was wrong at scale: RAxML aborts on an assertion in `optimizeModel.c` above roughly 10 000
  taxa, which is SATIVA's own `GAMMA_UPPER_THRES`, and every clustered MSA at or above 10 000
  sequences died. Reverted on 2026-08-19.
  The problem it was solving is real and is now handled in `epac/epang_l1o.py` instead: under
  `auto` a tree above 500 taxa is built with GTRCAT, whose `RAxML_info` carries the
  placeholder `alpha: 1.000000` because CAT estimates no shape, and EPA-ng, which always
  places under GAMMA, took it for an estimate. The module now checks whether that file holds
  an actual estimate and, if not, hands EPA-ng nothing and lets it fit its own model. At 5402
  sequences, K=25: 576 calls with the placeholder, 608 with EPA-ng fitting, 634 with a shape
  fitted by RAxML under GTRGAMMA. Fitting costs nothing measurable (140 s against 150 s).
- **EPA-ng binary**: resolved from `PATH` (provided by `sativa.yaml`); overridable via
  `CURATEMAKE_EPANG_BIN`.
- **Speed and determinism knobs**: see section 6. All of them default to the fast, ordered
  behaviour, so there is nothing to set for normal use. They keep the `SATIVA_EPANG_*`
  naming of the published module rather than the pipeline's `CURATEMAKE_*`, so that
  `epac/epang_l1o.py` here stays byte-identical to the one in
  `18.sativa_speedup/sativa_epang/` and updating it is a copy.

## 6. Speed of the leave-one-out (2026-08-19)

`epac/epang_l1o.py` was reworked for speed. **The calls do not change**: at both sizes
measured the `.mis` file is byte-identical to the one the previous module produced, twice
in a row.

Wall clock of the leave-one-out plus the classification, K=25, 8 threads, with the
reference tree built once and reused so that only this part varies:

| | n=1600 | n=5402 |
|---|---|---|
| before | 33.2 s | 233 s |
| after | 13.0 s | 150 s |
| | 2.6× | 1.6× |

Three changes carry it, and none of them touches what SATIVA decides:

1. **Pruning the reference tree per fold.** ete3's `prune()` calls `get_common_ancestor` for
   every kept node, once per fold. An O(n) pass over flattened arrays replaces it and
   reproduces its output exactly, topology and branch lengths both. 28.0 s → 0.8 s at 5402
   taxa, K=25.
2. **Mapping EPA-ng edges back to the `B=` numbering.** This built a leaf-name set per node
   plus its complement, O(n²), once per fold. Each clade is now summarised by the XOR of a
   random 128-bit value per leaf, which makes it O(n). 20.5 s → 7.6 s under the same
   conditions. `SATIVA_EPANG_MAP_CHECK=1` runs the old map alongside and reports any
   disagreement; there is none.
3. **Running several folds at once**, with the thread budget split between them. EPA-ng
   parallelises over queries, so this only pays when a fold has few of them: at 5402 taxa,
   K=25 (216 queries a fold) four at once is 1.31×, while K=5 (1080 a fold) four at once is
   *slower* than one. The module picks the number itself from that ratio and caps it at four,
   because each concurrent EPA-ng holds its own copy of the reference (~3 GB at 5402 taxa).
   `SATIVA_EPANG_FOLD_JOBS` overrides it; raise it on a large node.

**What to expect here.** The gain shrinks with alignment size, because EPA-ng itself is 77 %
to 88 % of the leave-one-out and the work removed above is the python around it. The
pipeline's clustered alignments sit at the top of the measured range, so budget about 1.6×,
not 2.6×.

**Reproducibility, which was not there before.** Two runs of the same command on the same
data did not necessarily give the same file. EPA-ng above one thread returns the queries in
whatever order its threads finish, and when two candidate labels tie,
`assign_taxonomy_maxsum` keeps whichever it saw first. The placements are now sorted before
SATIVA sees them, which costs nothing measurable; `SATIVA_EPANG_SORT=0` opts out. Three runs
at n=1600, K=5 gave one file with it and two different files without. It is intermittent
rather than systematic: the same three runs at n=5402 agreed even without it.

One cause is *not* fixed and is upstream behaviour: the reference tree step. RAxML PTHREADS
does not fix the order it reduces alignment columns in, so `-p` does not pin the tree the
way it looks like it should. Three full runs of the same command at 1600 sequences, one
after another, produced three different reference trees.

That does not automatically move the calls, though. Those same three runs wrote byte
identical `.mis` files. Where the trees differ is not where the flagged sequences sit, at
least on this alignment, so the tree varying is not on its own a reason to expect a
different answer.

Detail and the full matrix, including the approximations that were measured and rejected:
`18.sativa_speedup/OPTIMISATIONS.md`.

## 7. Known limitations / future work

- **Small groups** (< ~300 seqs): EPA-ng can be slightly slower than RAxML (process start-up
  cost), but the absolute time there is negligible; the net gain comes from the large groups.
- SATIVA's **final EPA confirmation** (`run_epa_once`, run only on the flagged suspects)
  also runs on EPA-ng since 2026-08-18; RAxML `-f v` was fragile on the small pruned trees.
  It is a plain placement (no k-fold), so it is not an approximation — but it *is* a second
  point where the engine differs from upstream, and any comparison against unmodified
  SATIVA should say so.
- **Model**: EPA-ng always runs GTRGAMMA. Since `sativa.cfg` pins the reference tree to
  GTRGAMMA, the shape parameter it reads is now fitted; before that it silently used 1.0 on
  any tree above 500 taxa. The RAxML fallback still switches to GTRCAT at scale, so the two
  engines do not run the same model when that fallback is used.
- **Duplicated wrong labels**: with 33 mislabels injected into three clades, unmodified
  SATIVA flags 28 and this version 33, and all five misses are cases where at least two
  sequences share the wrong label. What recovers them is EPA-ng, not the k-fold: holding out
  one sequence at a time also finds all 33 once EPA-ng does the placing. Why RAxML misses
  them is not established. This is a behavioural difference from upstream, in this project's
  favour, and it should be presented as such rather than as noise.
- **Confidence cutoff**: pooled over all sizes, calls the reference makes between 0.40 and
  0.50 confidence are reproduced about half the time, against ~0.8 above 0.5, and above 0.5
  confidence stops predicting reproducibility. The pipeline runs at `-C 0.4`; 0.5 would drop
  the unstable band at the cost of 39 calls out of 788.

---

*Refs: Barbera et al. 2019, EPA-ng, Syst. Biol. 68(2):365–369. Kozlov et al. 2016, SATIVA.*
