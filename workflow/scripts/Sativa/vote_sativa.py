#!/usr/bin/env python3
"""Consensus vote over several independent SATIVA splits of the SAME barcode.

Why
---
The group boundaries handed to SATIVA are an artefact of the traversal order in
`slice_groups`, not of biology, and the neighbourhood drives the verdict hard. Measured on
ITS2 (2026-08-26, three independent splits of the same 1 277 566 sequences):

  * only 43% of the flags raised by any one split are reproduced by another;
  * a SATIVA confidence of exactly 1.000 still leaves 29% unreproduced -- the confidence
    scores the placement WITHIN the given neighbourhood, so by construction it cannot see
    the arbitrariness of that neighbourhood;
  * fitting the observed 0/1/2/3-vote counts identifies three populations: 41.8% never
    flagged, 47.3% flagged ~38% of the time (the artefact zone), 11.0% flagged by every
    split. The 29.2% flag rate the pipeline reports is 11% real + 18% coin tosses;
  * agreement collapses with depth: among sequences flagged by >= 2 splits the proposed
    kingdom agrees 99.7% of the time, the order 91.4%, the genus 66.0%, the species 29.1%.
    That last line is the whole reason the vote is taken RANK BY RANK and the lineage
    truncated at the first rank without a majority, instead of being accepted or rejected
    whole.

Output contract
---------------
The consensus is written in the SATIVA `.mis` layout, so `integrate_sativa.py` consumes it
through the exact same `SativaMisParser` path as a raw SATIVA output -- no second code path
into the database. Two extra columns (Votes, Replicates) are appended; `parse_mis_file`
reads columns 0..7 and ignores the rest, so an old reader still works.

A sequence below the threshold is simply absent from the consensus. Integration then finds
it in the clade's `.tax` but not in the `.mis` and classifies it VALIDATED -- which is
precisely the intended meaning: no reproducible objection to its current label.

Each consensus line also carries the NATIVE clade that owns the sequence, taken from the
seed-0 `group_members.tsv`. Integration is scattered per clade, and a backbone sequence sits
in every clade's `.tax` as scaffold, so filtering the consensus on the `.tax` alone would
hand the same sequence to all 122 ITS clades. `group_members.tsv` lists members only, which
is exactly the one-sequence-one-clade attribution the scatter needs.
"""
import argparse
import json
import logging
import math
import sys
from fractions import Fraction
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sativa_tokens import normalize_tax_token  # noqa: E402
# NB: deliberately NOT `from integrate_sativa import SativaMisParser`. That chain pulls
# col_taxonomic_matcher -> requests, which workflow/envs/sativa.yaml does not declare,
# and the job then died at import -- before opening its log, so the failure was a
# 0-byte log and an absent consensus rather than an error.

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
log = logging.getLogger("vote_sativa")

N_RANKS = 7
_norm = normalize_tax_token

# "this rank does not exist for this lineage" -- a positive claim, unlike an abstention.
# A sentinel object so it can never collide with a taxon name.
ABSENT = object()


def parse_threshold(value: str) -> Fraction:
    """Accept `2/3`, `1/2`, or a decimal.

    Fractions are taken EXACTLY. A decimal is taken literally, which is rarely what is
    meant for small replicate counts: 0.667 x 3 = 2.001, so ceil() asks for 3 votes out of
    3, not 2. Write `2/3` and the threshold is 2 at N=3 and 7 at N=10, as intended.
    """
    text = str(value).strip()
    try:
        return Fraction(text)                    # handles "2/3" and "0.667" alike
    except (ValueError, ZeroDivisionError):
        raise SystemExit(f"--vote-threshold: cannot read {value!r} (try 2/3 or 0.5)")


def _path_to_list(path: str) -> Tuple[List[str], int]:
    """Semicolon path -> (7 normalised slots, number of ranks AS WRITTEN).

    The written length matters: SATIVA emits as many ranks as it resolved -- three for
    `Fungi;Basidiomycota;NA`, seven for `Plantae;Tracheophyta;NA;Asterales;NA;NA;NA` --
    and the consensus must emit the same number, otherwise N identical replicates stop
    reproducing their own source. Padding to 7 and trimming trailing NA are both wrong.
    """
    raw = (path or "").split(";")
    vals = [_norm(v) or "" for v in raw]
    return (vals + [""] * N_RANKS)[:N_RANKS], min(len(raw), N_RANKS)


def _conf_to_list(per_rank: str, n: int) -> List[float]:
    out = []
    for tok in (per_rank or "").split(";"):
        try:
            out.append(float(tok))
        except ValueError:
            out.append(0.0)
    return (out + [0.0] * n)[:n]


class Vote:
    """All the verdicts a single sequence collected, one entry per replicate."""

    __slots__ = ("by_rep",)

    def __init__(self):
        # replicate -> list of
        #   (proposed[7], original[7], per_rank[7], confidence, prop_len, orig_len)
        self.by_rep: Dict[int, List[Tuple[List[str], List[str], List[float],
                                          float, int, int]]] = defaultdict(list)

    def representative(self, rep: int):
        """One verdict per replicate.

        A backbone sequence is scaffolded into EVERY group, so within one replicate it can
        be flagged several times. That must still count as a single vote, otherwise a
        scaffold read alone could clear the threshold. The modal proposed lineage wins,
        ties broken by confidence.
        """
        items = self.by_rep[rep]
        if len(items) == 1:
            return items[0]
        modal = Counter(tuple(it[0]) for it in items).most_common(1)[0][0]
        cands = [it for it in items if tuple(it[0]) == modal]
        return max(cands, key=lambda it: it[3])


def read_mis(path: Path) -> List[List[str]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(";") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            rows.append(p)
    return rows


def read_members(path: Path) -> Dict[str, str]:
    """seq_id -> native clade tag (members only; scaffold copies are not listed)."""
    owner = {}
    with open(path) as fh:
        header = fh.readline()
        if not header.lower().startswith("tag"):      # no header line: rewind
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                owner[parts[1]] = parts[0]
    return owner


def collect(tasks: List[Dict], sativa_dir: Path) -> Tuple[Dict[str, Vote], Counter]:
    votes: Dict[str, Vote] = defaultdict(Vote)
    seen = Counter()
    for t in tasks:
        rep = int(t.get("replicate", 0))
        clade = sativa_dir / t["barcode"] / f"alignments_{t['tax_level']}" / t["taxon"]
        # A missing summary means SATIVA has not run for this clade -- which, if it happens
        # for ALL of them, means this job was scheduled BEFORE its inputs existed. That is
        # only possible if the rule's input function handed Snakemake an empty list, and it
        # must be an error, never an empty consensus (see the check in main()).
        if not (clade / "sativa.summary.json").exists():
            seen["clades_without_summary"] += 1
        mis = [p for p in sorted(clade.glob("*.mis")) if p.name != "consensus.mis"]
        if not mis:
            seen["clades_without_mis"] += 1
            continue
        seen["clades_read"] += 1
        for p in read_mis(mis[0]):
            prop, prop_len = _path_to_list(p[6])
            orig, orig_len = _path_to_list(p[5])
            try:
                conf = float(p[4])
            except ValueError:
                conf = 0.0
            votes[p[0]].by_rep[rep].append(
                (prop, orig, _conf_to_list(p[7], N_RANKS), conf, prop_len, orig_len))
            seen["records_read"] += 1
    return votes, seen


def consensus_line(seq_id: str, vote: Vote, n_rep: int, threshold: int,
                   native_taxon: str, stats: Counter) -> Optional[str]:
    reps = sorted(vote.by_rep)
    n_votes = len(reps)
    stats[f"votes_{n_votes}"] += 1
    if n_votes < threshold:
        stats["below_threshold"] += 1
        return None

    picked = [vote.representative(r) for r in reps]

    # Rank-by-rank majority. A rank left empty by EVERY voter is legitimately absent (the
    # pipeline's taxonomy has no class rank at all: 0 of 224 501 double-flagged ITS2
    # sequences carry one) and must not stop the walk -- only an actual disagreement does.
    # A rank left blank is NOT always silence. In this taxonomy some lineages genuinely
    # have no `class`, and a blank followed by populated ranks ASSERTS that absence. Only a
    # blank with nothing after it is an abstention -- the lineage simply stopped there.
    # So each voter has an opinion at rank i only while i is above its deepest populated
    # rank, and that opinion is either a name or ABSENT. Both are votes; abstentions are not.
    depths = [max((j for j, v in enumerate(pk[0]) if v), default=-1) + 1 for pk in picked]

    cons: List[str] = []
    cons_conf: List[float] = []
    truncated_at = None
    for i in range(N_RANKS):
        opinions = [(pk[0][i] or ABSENT) for pk, d in zip(picked, depths) if i < d]
        if not opinions:
            # past every voter's lineage: trailing padding, never a conflict
            cons.append("")
            cons_conf.append(sum(pk[2][i] for pk in picked) / n_votes)
            continue
        counts = Counter(opinions)
        top, cnt = counts.most_common(1)[0]
        distinct_names = {o for o in opinions if o is not ABSENT}

        if cnt * 2 > n_votes:                      # strict majority of ALL voters
            chosen = top
        elif len(distinct_names) >= 2:
            # two or more DIFFERENT names, none with a majority: a real taxonomic conflict.
            # Everything below is unreliable too (a genus belongs to a family), stop here.
            truncated_at = i
            break
        else:
            # at most one name in play, the rest asserting the rank does not exist:
            # never a conflict, and a name without majority is not asserted
            chosen = ABSENT

        agree = [pk[2][i] for pk in picked if (pk[0][i] or ABSENT) == chosen]
        cons.append("" if chosen is ABSENT else chosen)
        cons_conf.append(sum(agree) / len(agree) if agree
                         else sum(pk[2][i] for pk in picked) / n_votes)
    if truncated_at is not None:
        stats["truncated"] += 1
        stats[f"truncated_at_rank_{truncated_at}"] += 1
    # Trailing empty ranks are KEPT: SATIVA itself writes them (`Plantae;Tracheophyta;NA;
    # Asterales;NA;NA;NA` with seven confidences), and trimming them made the consensus
    # differ from its own source when every replicate agreed -- which broke the invariant
    # that N identical splits must reproduce the original .mis line for line.
    if not any(cons):
        stats["dropped_no_consensus"] += 1
        return None
    # emit exactly as many ranks as the voters wrote, capped by where the vote stopped
    keep = min(len(cons), max(pk[4] for pk in picked))
    cons, cons_conf = cons[:keep], cons_conf[:keep]

    # The original label is the database's, identical across replicates; take the modal
    # anyway rather than trusting that.
    orig = list(Counter(tuple(pk[1]) for pk in picked).most_common(1)[0][0])
    orig_keep = max(pk[5] for pk in picked)

    # Shallowest rank where the consensus departs from the original -- covers a renamed
    # rank AND a truncation (the DELETION case, where the consensus is simply shorter).
    level_idx = None
    for i in range(N_RANKS):
        o = orig[i] if i < len(orig) else ""
        c = cons[i] if i < len(cons) else ""
        if o != c:
            level_idx = i
            break
    if level_idx is None:
        stats["dropped_agrees_with_original"] += 1
        return None

    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    fields = [
        seq_id,
        ranks[level_idx],
        orig[level_idx] or "NA",
        (cons[level_idx] if level_idx < len(cons) else "") or "NA",
        f"{sum(pk[3] for pk in picked) / n_votes:.3f}",
        ";".join(v or "NA" for v in orig[:orig_keep]) or "NA",
        ";".join(v or "NA" for v in cons),
        ";".join(f"{c:.3f}" for c in cons_conf),
        str(n_votes),
        str(n_rep),
        native_taxon,
    ]
    stats["kept"] += 1
    return "\t".join(fields)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="merged SATIVA task manifest")
    ap.add_argument("--sativa-dir", required=True)
    ap.add_argument("--barcode", required=True)
    ap.add_argument("--level", required=True)
    ap.add_argument("--replicates", type=int, required=True)
    ap.add_argument("--native-members", required=True,
                    help="seed-0 group_members.tsv: which native clade owns each sequence")
    ap.add_argument("--vote-threshold", default="2/3",
                    help="fraction of the replicates that must agree, e.g. 2/3 or 1/2")
    ap.add_argument("--output", required=True, help="consensus .mis")
    ap.add_argument("--report", required=True)
    ap.add_argument("--log-file")
    a = ap.parse_args()

    if a.log_file:
        Path(a.log_file).parent.mkdir(parents=True, exist_ok=True)
        h = logging.FileHandler(a.log_file)
        h.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
        logging.getLogger().addHandler(h)

    manifest = json.load(open(a.manifest))
    tasks = [t for t in manifest.get("tasks", [])
             if t["barcode"] == a.barcode and t["tax_level"] == a.level]
    reps_present = sorted({int(t.get("replicate", 0)) for t in tasks})
    if len(reps_present) != a.replicates:
        log.warning("manifest carries %d replicate(s) %s but --replicates is %d",
                    len(reps_present), reps_present, a.replicates)
    n_rep = len(reps_present)
    frac = parse_threshold(a.vote_threshold)
    threshold = max(1, math.ceil(frac * n_rep))
    log.info("%s/%s: %d task(s) over %d replicate(s); threshold = %d/%d votes",
             a.barcode, a.level, len(tasks), n_rep, threshold, n_rep)

    owner = read_members(Path(a.native_members))
    log.info("native membership: %d sequences over %d clade(s)",
             len(owner), len(set(owner.values())))

    votes, seen = collect(tasks, Path(a.sativa_dir))
    log.info("read %d records from %d clade(s) (%d without .mis, %d without summary) "
             "-> %d distinct sequences",
             seen["records_read"], seen["clades_read"], seen["clades_without_mis"],
             seen["clades_without_summary"], len(votes))

    # ---- fail LOUD ---------------------------------------------------------------------
    # An empty consensus is not a verdict. By design a sequence absent from it is taken as
    # VALIDATED, so writing an empty file silently discards every SATIVA correction and the
    # pipeline degrades to CoL-only output with zero errors reported -- exactly the
    # regression seen on 2026-08-26 (39.3% correction, identical to CoL, against 79.4% for
    # the single-split path). Refuse instead.
    if tasks and seen["clades_without_summary"] == len(tasks):
        raise SystemExit(
            f"{a.barcode}/{a.level}: none of the {len(tasks)} clade(s) has a "
            f"sativa.summary.json. This job ran BEFORE SATIVA -- its input function must "
            f"have returned an empty list, so Snakemake had no ordering constraint.")
    if tasks and seen["clades_read"] == 0:
        raise SystemExit(
            f"{a.barcode}/{a.level}: {len(tasks)} clade(s) declared, not one .mis found "
            f"under {a.sativa_dir}. Refusing to write an empty consensus, which would be "
            f"read downstream as 'nothing is mislabelled'.")

    stats = Counter()
    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write(f";Consensus of {n_rep} independent SATIVA splits of {a.barcode}, "
                 f"threshold {threshold}/{n_rep} votes.\n")
        fh.write(";Ranks are voted independently and the lineage is truncated at the first "
                 "rank without a strict majority.\n")
        fh.write(";SeqID\tMislabeledLevel\tOriginalLabel\tProposedLabel\tConfidence"
                 "\tOriginalTaxonomyPath\tProposedTaxonomyPath\tPerRankConfidence"
                 "\tVotes\tReplicates\tNativeTaxon\n")
        for seq_id in sorted(votes):
            native = owner.get(seq_id)
            if native is None:
                # Every projected sequence is a member of exactly one native group, so a
                # miss means consensus and native split disagree about the universe.
                stats["dropped_no_native_clade"] += 1
                continue
            line = consensus_line(seq_id, votes[seq_id], n_rep, threshold, native, stats)
            if line:
                fh.write(line + "\n")

    report = {
        "barcode": a.barcode, "tax_level": a.level,
        "replicates": n_rep, "replicate_ids": reps_present,
        "vote_threshold": str(frac), "min_votes": threshold,
        "clades_read": seen["clades_read"],
        "clades_without_mis": seen["clades_without_mis"],
        "records_read": seen["records_read"],
        "clades_without_summary": seen["clades_without_summary"],
        "sequences_flagged_at_least_once": len(votes),
        "consensus_written": stats["kept"],
        "below_threshold": stats["below_threshold"],
        "dropped_no_consensus": stats["dropped_no_consensus"],
        "dropped_agrees_with_original": stats["dropped_agrees_with_original"],
        "dropped_no_native_clade": stats["dropped_no_native_clade"],
        "lineages_truncated": stats["truncated"],
        "votes_histogram": {str(k): stats[f"votes_{k}"] for k in range(1, n_rep + 1)},
        "truncation_by_rank": {str(k): stats[f"truncated_at_rank_{k}"]
                               for k in range(N_RANKS) if stats[f"truncated_at_rank_{k}"]},
    }
    Path(a.report).parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(a.report, "w"), indent=2)
    log.info("consensus: %d of %d flagged sequences kept (%.1f%%), %d truncated -> %s",
             stats["kept"], len(votes),
             100.0 * stats["kept"] / len(votes) if votes else 0.0,
             stats["truncated"], out)


if __name__ == "__main__":
    main()
