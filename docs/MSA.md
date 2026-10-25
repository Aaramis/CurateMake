# MSA stage — bottom-up construction (SUPERSEDED)

> **This page describes a strategy the pipeline no longer runs.** Bottom-up construction
> was replaced by the eHMM shared frame, in `rules/ehmm_msa.smk` and
> `scripts/MSA/ehmm/build_ehmm_msa.py`, whose own docstring is the reference for what runs
> today. The rule file this page refers to, `rules/bottomup_msa.smk`, and the design note
> it cites, `INTEGRATION_DESIGN.md`, were both deleted with it.
>
> It is kept because it records why bottom-up replaced the top-down path, which still
> matters: that reasoning is one of the four steps in `alignment-strategies.md`, the
> plain-language history of every strategy tried here. Read that page first.

The MSA stage turns the post-ITSx sequence set into reference alignments and one HMM
profile per emergent taxonomic unit, then hands SATIVA the alignments worth validating.

Two rule files were involved:

- `rules/msa.smk` — reduced to `prepare_msa` (DB → `results/full.csv`) and the report rule.
  This one still exists, and `rules/msa_groups.smk` now merges what eHMM emits.
- `rules/bottomup_msa.smk` — the construction itself (six rules). Deleted.

`msa.bottomup.enabled` selected it. The older top-down path (align a whole family, then
descend when the alignment is poor) had already been removed at that point.

## Why bottom-up

Top-down aligned each family in one shot, dropped sequences that didn't overlap
(`trimAl -seqoverlap`), then re-classified the ~3M discarded sequences against the HMM
library — that re-classification was the dominant cost. Bottom-up changes two things:

- it builds from the **species level up** and merges upward only while the alignment
  stays coherent, so references are cleaner and SATIVA gets real contrast to work with;
- it trims with `trimAl -gappyout`, which drops sparse columns but **keeps every
  sequence**, so the re-classification residual shrinks by ~1000× (only true singletons
  are left to place).

## Unit of work

Scatter is per `(barcode, family)`; the whole intra-family tree (species → genus →
family) runs inside a single job. Families are independent, so this parallelises across
the node/cluster without a deep checkpoint DAG. Barcodes are never mixed — ITS, ITS1 and
ITS2 are different regions and get separate trees.

`checkpoint bottomup_split` (`split_families.py`) streams `results/full.csv` in chunks
and writes one CSV per `(barcode, family)` plus `families_manifest.json`.

## The floor (per species)

For each species, `build_family_msa.py`:

1. clusters raw sequences by identity — `vsearch --cluster_fast`, `cluster_id = 0.97`;
2. aligns each cluster (MAFFT) and trims it (`trimAl -gappyout`).

Clustering first is deliberate: a genome can carry several divergent ITS copies, so one
species can legitimately hold two very different barcodes. Each copy becomes its own leaf
MSA instead of being forced into one bad alignment.

Groups larger than `max_group` (300) are sub-sampled before alignment — past a few
hundred sequences the alignment no longer improves — while every membership is still
recorded for downstream bookkeeping.

## The climb and the gate

Leaves are merged upward with `mafft --merge`, which preserves the sub-alignments instead
of re-aligning them. At each rank a quality gate decides whether to promote the merge:

- `occ ≥ gate_occ` (default 40) — percentage of columns filled in ≥50% of sequences;
- `ratio ≤ gate_ratio` (default 2.5) — alignment length ÷ median ungapped length.

Pass → the group becomes a single reference at that rank. Fail → the sub-groups stay
separate. The **stop level is emergent per lineage**: species for divergent fungi, often
genus or family for plants. An MSA below `min_ref` (2 sequences) can't be validated on
its own and goes to the residual.

## Residual, contamination, cross-family climb

- **`bottomup_residual`** (`residual_classify_integrate.py`) — the residual (singletons,
  unplaced) is classified against the HMM profile library with `hmmsearch`, and each
  orphan is added to its best MSA with `hmmalign` (mafft fallback). Small by construction.
- **`bottomup_contamination`** (`coarse_contamination.py`) — a coarse cross-kingdom /
  phylum check: each reference is classified against the full profile library and its
  homology compared to its label. This catches gross contamination (e.g. a fungal
  sequence labelled as a plant) that SATIVA structurally cannot see — SATIVA only detects
  fine mislabels *inside* an already-coherent clade. Flagged MSAs are re-homed to their
  true clade.
- **`climb_lineage`** (`climb_lineage.py`) — an unbounded cross-family climb
  (family → order → class → phylum) that merges family references while the same gate
  passes. Runs after the residual step, on already-enriched MSAs.

## Collection and SATIVA hand-off

`bottomup_collect` (`collect_bottomup.py`) aggregates the per-family manifests and the
climbed manifest into `bottomup_summary.json` and emits the SATIVA manifest. Only MSAs
with real validation power are sent: at least `sativa_min_seqs` (4) sequences and, when
`sativa_min_taxa > 1`, at least that many distinct taxa — a single-species MSA has
nothing to contrast. Re-homed / emptied MSAs are excluded.

## Rules and scripts

| Rule | Script | Does |
|---|---|---|
| `bottomup_split` (checkpoint) | `split_families.py` | `full.csv` → one CSV per (barcode, family) + manifest |
| `build_family` | `build_family_msa.py` | floor + gated climb → family MSAs + HMM profiles |
| `bottomup_contamination` | `coarse_contamination.py` | cross-kingdom detection + re-homing |
| `bottomup_residual` | `residual_classify_integrate.py` | `hmmsearch` reclassification + `hmmalign` integration |
| `climb_lineage` | `climb_lineage.py` | unbounded family → … → phylum climb |
| `bottomup_collect` | `collect_bottomup.py` | summary + SATIVA manifest |

Dependency order (via inputs): `build_family` → `contamination` → `residual` →
`climb_lineage` → `collect`. Outputs live under `results/MSA/bottomup/`.

## Configuration

```yaml
msa:
  bottomup:
    enabled: true
    cluster_id: 0.97        # vsearch identity at the floor
    gate_ratio: 2.5         # climb if (length / median length) <= this
    gate_occ: 40.0          # ... and column occupancy (>=50% filled) >= this %
    max_group: 300          # sub-sample cap per cluster/species
    min_ref: 2              # smaller MSAs go to the residual
    threads_per_family: 2   # few threads/family -> many families in parallel
    climb_ranks: ["order","class","phylum"]   # cross-family climb; [] stops at family
    sativa_min_seqs: 4
    sativa_min_taxa: 1      # >1 skips single-species MSAs (can cut coverage sharply)
    contamination_flag_ranks: ["kingdom","phylum"]
```

`gate_ratio`, `gate_occ` and `cluster_id` are the calibration knobs: a tighter gate
lowers the emergent stop level, giving more but smaller MSAs.

## Determinism

Reproducible by construction: inputs are sorted before `vsearch` (single thread), MAFFT
runs as FFT-NS-2, single thread, input order preserved, `-gappyout` is deterministic, and
merge order is sorted.

## Tips

```bash
# Rebuild one family
snakemake --use-conda results/MSA/bottomup/out/ITS__Russulaceae/manifest.json

# Stop-level distribution and per-barcode counts
jq '.stop_levels, .per_barcode' results/MSA/bottomup/bottomup_summary.json

# Residual reinsertion (should be small; only unassigned is a true loss)
jq '{residual:.n_residual, reinserted:.n_integrated, unassigned:.n_unassigned}' \
   results/MSA/bottomup/residual_report.json
```
