# Alignment strategies, in plain language

Four strategies have been tried on this pipeline. Three were abandoned and one is running
today. This page says what each one did, what went wrong, and what the measurements were,
without assuming you know the tools. The technical page for the current one is
`MSA.md`; the reasoning behind rejecting the fourth is in `alignment-strategy.md`.

## Why this step is the hard one

Before the pipeline can tell whether a sequence is labelled correctly, it has to line up
sequences against each other so that the same position means the same thing in every one
of them. That is the alignment.

For ITS this is unusually hard. ITS is hypervariable: two sequences from different genera
may share almost nothing, and even sequences of about 600 letters can end up spread across
thousands of alignment columns, each mostly empty. Every strategy below is an attempt to
get around that, and the differences between them come down to one question: **what do you
align together, and what do you align against?**

## 1. Top-down: align a whole family, then descend

The first approach aligned an entire family in one go. Sequences that did not overlap
enough with the rest were dropped, and the pipeline then went back and tried to re-attach
those discarded sequences by matching them against a library of profiles.

**Why it was abandoned.** That re-attachment step became the single most expensive part of
the pipeline, with roughly three million sequences to place. The alignment was throwing
away most of its input and then paying to get it back.

## 2. Bottom-up: start at species, merge upwards while it holds

The replacement went the other way. It aligned each species on its own, then merged
species into genera, genera into families, and stopped climbing as soon as the merged
alignment stopped being coherent. It also changed how columns are trimmed: instead of
dropping sequences that fit poorly, it drops sparse columns and keeps every sequence.

**What it bought.** The leftover pile of unplaced sequences shrank by about a thousandfold.
Only genuine singletons were left over.

**Why it was abandoned.** It produces one alignment per taxonomic unit. That is good for
placing a sequence precisely, but it makes one kind of error invisible. If a plant sequence
is mislabelled as a fungus, you can only notice by seeing it sit next to its true
relatives. In a fungus-only alignment, there are no plants to sit next to.

## 3. eHMM: one shared frame per barcode (what runs today)

The current strategy builds a single set of columns per barcode and makes every sequence
fit into it. Nothing is ever realigned afterwards. It runs in four steps.

**Build the frame.** Cluster the sequences, take one representative per large cluster,
and align only those. That alignment defines the columns, once and for all.

The choice of representatives matters more than it looks. A sample balanced by kingdom and
phylum lets the pipeline assign 66.1 % of reads. Cluster representatives reach 81.1 %,
because what has to be covered is the variety of *sequences*, not the variety of *names*.

Columns that almost nobody occupies are then thrown away. Without that, a cross-kingdom
ITS alignment spreads roughly 470-letter reads over about 6 700 columns, so any given
column is filled 7 % of the time, and a profile built on it stops recognising even the
sequences it was trained on: 5 out of 17 in the measurement. Keeping only columns filled
at least 10 % of the time leaves about 1 000 columns, and recognition rises to about 77 %.

**Project everyone onto it.** Sequences are bucketed by clade, one profile is built per
bucket over the whole frame, and each read is matched to its best bucket, with the correct
orientation detected at the same time. The read is then placed into the frame. Every
sequence ends up as a row of exactly the same width.

**Rescue the misses.** The reads no profile recognises are mostly short fragments, around
166 letters against a typical 457. A fast similarity search finds them against reads that
are already placed, and they inherit that read's bucket and orientation.

**Cut into working groups.** The frame is sliced into groups small enough to validate, and
kingdoms are dealt out in rotation rather than in blocks, with the backbone added to every
group. Otherwise the small kingdoms get used up first and the last groups contain only two
kingdoms out of five, which is exactly the situation where a contaminant becomes
undetectable.

## 4. Per-clade alignment: proposed, measured, rejected

The obvious objection to the current strategy is that a frame spanning all kingdoms cannot
be positionally meaningful, so why not align each clade separately and keep a small
separate alignment for the conserved region in the middle? That was written up as a formal
proposal with pass-or-fail gates, and it failed the first one. Three families were aligned
on their own and compared with the shared frame on the same sequences.

| family | sequences | columns alone | columns in the shared frame |
|---|---:|---:|---:|
| Enchytraeidae | 1 359 | 3 331 | 472 |
| Meruliaceae | 3 685 | 7 340 | 472 |
| Caryophyllaceae | 2 303 | 1 805 | 472 |

**It is not denser.** Sequences of about 600 letters still spread over thousands of columns
*inside a single family*. The hypervariability is not something that clade grouping
removes, because it is already there between genera of the same family.

**It costs more.** One family of 3 685 sequences over 7 340 columns needs about three times
the memory of a shared-frame task holding 20 000 sequences, which is five times more
sequences. The narrowness of the shared frame is what keeps the cost bounded.

**It cannot check the rank it groups on.** In the annelid family, 22 of 23 flags landed at
family rank, inside a tree containing one single family. There is no alternative family to
compare against, so those flags are artefacts of what is missing rather than findings.

## Where that leaves things

| strategy | status | the thing that killed it, or kept it |
|---|---|---|
| Top-down | removed | re-attaching discarded sequences dominated the cost |
| Bottom-up | removed | cannot see cross-kingdom contamination |
| eHMM shared frame | **current** | bounded cost, and contamination stays visible |
| Per-clade | rejected before implementation | wider, more expensive, and blind at its own rank |

The honest summary is that the current design is good at the coarse job, spotting a
sequence that sits in the wrong kingdom or the wrong family, and weak at the fine one,
separating close species. That is a trade-off, not an oversight: the measurements above say
the alternative is worse on both cost and coarse detection, and nothing yet measured says
it would be better on the fine one.

One experiment is not written up here because it was measured but never committed: giving
each kingdom its own frame instead of one shared frame. It should be written up or redone
before anyone concludes anything from it.
