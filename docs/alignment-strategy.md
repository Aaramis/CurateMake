# Alignment strategy for a hypervariable marker

**Status: REJECTED at Gate 1 on 2026-10-22. Nothing here was implemented, and on
this evidence nothing should be.** The result is kept because it answers the obvious
reviewer question, "why not align per clade?", with a measurement rather than an opinion. The current eHMM path
remains the only one. This document exists so the reasoning, the measurements and the
conditions for abandoning the idea are written down before any code is touched, because
the alignment step is the most consequential one in the pipeline.

## The question

ITS is a spacer. Outside the conserved 5.8S and the SSU/LSU flanks it is under no
structural constraint, so homology between distant lineages is not merely hard to detect,
it is largely absent. 16S is the opposite: constrained over its whole length, which is why
a single global alignment and a single profile work there. Methods that are sound for 16S
are not automatically sound here, and the strategy should be judged against ITS rather
than inherited.

## What the pipeline's own numbers already say

These are from the comments in `config/config.yaml`, measured on real ITS data during the
eHMM work:

| configuration | reads assigned |
|---|---:|
| 4 clade-level HMMs | 20.4 % |
| 108 UPP-style profiles | 72.8 % |
| backbone subsets of 25 | 60.1 % |
| backbone subsets of 10 | 66.1 % |
| backbone 9x larger, 1 264 profiles | 69.1 % |

The last row is the important one. Nine times more backbone buys three points, and about
30 % of reads remain unassigned by any profile. A plateau like that is not a tuning
problem. It says the barcode is not describable by profile HMMs at any granularity, and
that the method is already at its ceiling.

A measurement taken for this document says the same thing from the other side. Per-sequence
occupancy inside the group alignments the pipeline actually produced in the August 2026
full run:

| barcode | median occupancy | distribution |
|---|---:|---|
| ITS1 | ~0.37 | smooth, unimodal |
| ITS | ~0.51 | smooth, unimodal |

In an alignment of genuinely homologous sequences occupancy is high. At 0.37 most columns
of the frame belong to someone else's lineage. And the absence of bimodality matters
separately: there is no coverage threshold that separates complete sequences from
fragments, because the width is set by the diversity of the group, not by the completeness
of the sequence. That is why no coverage filter is proposed here.

## The structural problem

One alignment is being asked to do two incompatible jobs.

Detecting cross-kingdom contamination requires the contaminant **and its true relatives**
in the same MSA; that is what kingdom-interleaved groups and `scaffold_every_group` buy.
Placing a sequence correctly at genus or species requires columns that are **positionally
homologous**, which for ITS is true only inside a clade.

A frame that spans kingdoms is therefore a shared coordinate system, not an alignment in
the sense placement needs. The simulation results match that reading precisely: the
entropy gain runs from -53 % at phylum to -2 % at species, MISIDENTIFICATION stays at 0 %
in the controlled benchmark, and cross-kingdom contamination is where the arm performs
best. The design is good at the coarse job and weak at the fine one, which is exactly what
it trades.

## The proposal: two independent paths, no reassembly

To be explicit, because it is easy to misread: **nothing is stitched back together.** The
5.8S is never spliced onto the spacers. Two paths run on two different materials and their
verdicts meet only at integration, where they are already treated as different error
classes.

**Coarse path, on 5.8S.** One global alignment, one profile, about 158 columns. It answers
one question: is this sequence in the right kingdom, phylum, perhaps order. This replaces
what the kingdom interleaving buys today, at the cost of a small alignment instead of 122
groups of 20 000.

**Fine path, on the spacers.** Groups by clade: family by default, split at genus when too
large, merged upwards when too small. Each group aligned on its own, with its own columns.
No global frame, no interleaving.

### Is there enough 5.8S? Measured, yes

ITSx already runs with `--save_regions all`, so 5.8S is extracted and then discarded: only
`["ITS", "ITS1", "ITS2"]` travel downstream. The August run left 1 650 806 distinct
sequences carrying a 5.8S on disk, 140 to 168 nt long with a median of 158 — a 28 nt spread,
the signature of a conserved region and the opposite of the spacer.

Coverage of the 5 088 644 non-duplicate entries, counting a sequence as coverable when it
carries a 5.8S itself or when its ITSx parent does:

| | entries | share |
|---|---:|---:|
| own 5.8S | 1 246 991 | 24.51 % |
| through the ITSx parent | 2 984 000 | 58.64 % |
| **coverable** | **4 230 991** | **83.15 %** |
| not coverable | 857 653 | 16.85 % |

By barcode, full-length ITS is 99 % coverable, ITS1 74 %, ITS2 66 %. The uncoverable
remainder is overwhelmingly records that arrived already trimmed, with no parent holding a
full ITS. They need an alignment-free fallback or an explicit "not screened" status; they
must not be silently treated as clean.

## What is not known

- Whether a 158-column 5.8S profile actually separates kingdoms on this data. Believed,
  not demonstrated.
- Whether per-clade alignment improves genus and species verdicts at all. This is the
  hypothesis, not a result.
- The group-size distribution per clade, and how many families fall below a usable size.

## Gates

Nothing is replaced. If this is built, it is added behind a configuration switch with the
current path as the default, and each gate must pass before the next is attempted.

**Gate 1 — does per-clade alignment change the fine verdicts?** Three to five well-populated
families across kingdoms. Each aligned alone, SATIVA run on it, verdicts at genus and
species compared with the global-frame verdicts for exactly the same sequences. If the
verdicts barely move, the global frame is acceptable, this proposal is wrong, and the work
stops here. About one day.

**Gate 2 — is 5.8S discriminative?** Only if Gate 1 passes. Coverage is already measured at
83.15 %. What remains is whether a 5.8S profile separates kingdoms well enough to replace
what interleaving detects today.

**Gate 3 — does it hold end to end?** Implemented as an alternative strategy, run on the
mini dataset, then on one palier of the simulation grid, with arm scores compared against
the current ones. A full-database run only after that.

## Gate 1 result: the proposal fails

Three families were aligned on their own and placed with SATIVA, then their verdicts
compared with the global frame's verdicts for exactly the same sequences.

| family | kingdom | sequences | columns alone | flags per-clade | flags global frame | shared | same rank |
|---|---|---:|---:|---:|---:|---:|---:|
| Enchytraeidae | Animalia | 1 359 | 3 331 | 23 | 65 | 16 | 0 % |
| Meruliaceae | Fungi | 3 685 | 7 340 | 549 | 784 | 139 | 58 % |
| Caryophyllaceae | Plantae | 2 303 | 1 805 | 425 | 70 | 39 | 49 % |

**Per-clade alignment is not denser.** ITS sequences of roughly 600 nt spread over 3 331
columns inside a single family of annelids and 7 340 inside a single family of fungi,
against 472 for the global frame. The hypervariability is not a between-kingdom effect
that clade grouping removes; it is present between genera of one family. Caryophyllaceae
at 1 805 columns shows the width is strongly clade-dependent, which makes the cost
unpredictable.

**It costs more, not less.** At roughly 330 bytes per taxon-by-column cell, a global-frame
task of 20 000 sequences over 472 columns needs about 3.1 GB per fold, while one family of
3 685 sequences over 7 340 columns needs about 8.9 GB -- three times the memory for five
times fewer sequences. The narrowness of the global frame is doing real work: it is what
bounds the placement cost. Observed during the test: four concurrent epa-ng processes at
6.5 GB each for one family.

**It cannot validate the rank it groups on.** In Enchytraeidae 22 of 23 flags land at
Family rank, in a tree that contains only one family. SATIVA has no reference for "not
Enchytraeidae" there, so those flags are artefacts of the missing alternative rather than
findings. The global frame flagged the same sequences at Species and Genus. This is the
mirror of the singleton-taxon trap: a clade-restricted reference cannot detect an error at
or above the rank that defines the clade -- which is precisely what the kingdom
interleaving exists to buy.

**A confound to state.** The global-frame verdicts are post-vote, ten replicates with a
two-thirds consensus; the per-clade run is a single pass. A single pass is known to carry
about 38 % artefacts in its unstable population against 4 % after the vote, so a large part
of the per-clade-only flags are likely artefacts rather than discoveries. "Different" must
not be read as "better". Removing the confound would mean rerunning the per-clade arm with
the vote, which is worth doing only if the first two findings are ever overturned.

**Scope.** Three families, one per kingdom, no Chromista family in the usable size band.
This is enough to stop the proposal, not enough to characterise the behaviour.

## What should not change

The backbone selection by cluster representatives is a good, measured decision inside the
current design: it assigns 81.1 % of reads against 66.1 % for taxonomy-stratified
selection. ITSx extraction and the barcode reannotation are sound and are not in question
here.
