#!/usr/bin/env python3
"""
COARSE contamination-detection tier (cross-kingdom / cross-clade).

SATIVA only detects FINE mislabels within an alignable clade. A coarse contaminant
(e.g. a Fungi labelled as a Plant) is isolated by the homology clustering -> its own
cluster/MSA -> SATIVA finds it consistent -> missed. But its HOMOLOGY contradicts its
LABEL, at a high rank: that is THE signal.

Method: for each reference MSA, classify one representative against the FULL profile
library (all kingdoms). The best profile from ANOTHER family gives the "true" taxonomy
by homology. If it differs from the MSA's label at a high rank (kingdom/phylum...),
with a strong score -> CONTAMINATION flagged.

Requires NO cross-kingdom alignment (just a classification) -> the "ITS not alignable
across kingdoms" problem does not arise.

SCALES via Snakemake: the classification is O(M^2) (one rep per MSA against every
profile), so it is split into three steps for parallelism:
  prepare : build the profile library + one rep per MSA, split into batches, + an index.
  search  : one hmmsearch per rep-batch against the library (Snakemake scatter).
  report  : parse the batch hits, flag label-vs-homology disagreements, re-home.
A legacy all-in-one `run` subcommand is kept for small-scale / standalone use.
"""
import argparse
import json
import logging
import os
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from external_tools import run_hmmalign_add, run_mafft_add, existing_msa_ids

logger = logging.getLogger("coarse_contamination")
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def read_all_taxed(fasta):
    """[(raw_seq_id, tax_dict, degapped_seq)] for ALL sequences (taxonomized headers)."""
    out, header, s = [], None, []

    def flush():
        if header is None:
            return
        sid = header.split("|")[0]
        tax = {r: "" for r in RANKS}
        for part in header.split("|")[1:]:
            if "=" in part:
                k, v = part.split("=", 1)
                if k in tax:
                    tax[k] = v.replace("_", " ")
        out.append((sid, tax, "".join(s).replace("-", "").upper()))

    for line in open(fasta):
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            header, s = line[1:].split()[0], []
        else:
            s.append(line)
    flush()
    return out


def read_first_record(fasta):
    """(seq_id, tax_list, degapped_seq) of the 1st sequence, taxonomy from the header."""
    header, seq = None, []
    with open(fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:]
            elif header is not None:
                seq.append(line)
    if header is None:
        return None
    tax = {r: "" for r in RANKS}
    for part in header.split("|")[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            if k in tax:
                tax[k] = v.replace("_", " ")
    return header.split("|")[0], [tax[r] for r in RANKS], "".join(seq).replace("-", "").upper()


# --------------------------------------------------------------------- shared core
def build_index(family_manifests):
    """Read all family manifests -> (msa_tax, msa_info, reps, profiles).
    msa_tax: msa_id -> tax_list ; msa_info: msa_id -> {msa_path, profile_path}
    reps: [(msa_id, rep_seq)] ; profiles: [profile_path]."""
    msa_tax, msa_info, reps, profiles = {}, {}, [], []
    for mpath in family_manifests:
        m = json.load(open(mpath, encoding="utf-8"))
        for e in m.get("msas", []):
            rec = read_first_record(e["msa_path"]) if os.path.exists(e.get("msa_path", "")) else None
            if not rec:
                continue
            _sid, tax, seq = rec
            msa_tax[e["msa_id"]] = tax
            msa_info[e["msa_id"]] = {"msa_path": e["msa_path"], "profile_path": e.get("profile_path")}
            reps.append((e["msa_id"], seq))            # rep header = msa_id (self-exclusion)
            if e.get("profile_path"):
                profiles.append(e["profile_path"])
    return msa_tax, msa_info, reps, profiles


def parse_hits(tbl_files):
    """hmmsearch tblout(s) -> {rep_msa_id: [(score, profile_msa_id), ...]} (self excluded)."""
    hits = defaultdict(list)
    for tbl in tbl_files:
        if not tbl or not os.path.exists(tbl):
            continue
        for line in open(tbl):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            rep, prof, score = f[0], f[2], float(f[5])   # target=rep, query=profile, score
            if prof == rep:
                continue
            hits[rep].append((score, prof))
    return hits


def flag_contaminants(hits, msa_tax, flag_ranks, min_score, top_k):
    """Compare each rep's label to the homology CONSENSUS of its top-K hits."""
    idx = {r: RANKS.index(r) for r in flag_ranks}
    flagged = []
    for rep, hlist in hits.items():
        hlist.sort(reverse=True)
        if not hlist or hlist[0][0] < min_score:
            continue
        top = hlist[:top_k]
        lab = msa_tax.get(rep, [""] * 7)
        for rank in flag_ranks:
            i = idx[rank]
            vals = [msa_tax.get(p, [""] * 7)[i].strip() for _, p in top]
            vals = [v for v in vals if v]
            if not vals:
                continue
            hom = Counter(vals).most_common(1)[0][0]     # majority kingdom/phylum of the top-K
            lv = lab[i].strip()
            if lv and hom and lv != hom:
                flagged.append({
                    "msa_id": rep, "rank": rank, "label": lv, "homology": hom,
                    "consensus": f"{Counter(vals).most_common(1)[0][1]}/{len(vals)}",
                    "best_hit_profile": top[0][1], "score": round(top[0][0], 1),
                    "label_lineage": ";".join(lab),
                })
                break
    return flagged


def rehome(flagged, msa_info, threads):
    """Re-route flagged sequences to their TRUE clade (best-hit), then empty the source
    MSA + its profile. Idempotent. Returns (emptied_msa_list, n_rehomed)."""
    emptied, n_rehomed = [], 0
    for f in flagged:
        src, tgt = f["msa_id"], f["best_hit_profile"]
        src_info, tgt_info = msa_info.get(src), msa_info.get(tgt)
        if not src_info or not tgt_info:
            continue
        src_msa, tgt_msa = src_info["msa_path"], tgt_info["msa_path"]
        if not os.path.exists(src_msa) or not os.path.exists(tgt_msa):
            continue
        recs = read_all_taxed(src_msa)
        if not recs:
            if os.path.getsize(src_msa) == 0:          # already emptied (partial re-run)
                emptied.append(src)
                f["rehomed_to"] = tgt
            continue
        df = pd.DataFrame([{"seq_id": sid, "sequence": seq, **tax} for sid, tax, seq in recs])
        already = existing_msa_ids(tgt_msa)            # idempotence
        new_df = df[~df["seq_id"].isin(already)]
        ok = new_df.empty
        if not ok:
            if tgt_info.get("profile_path"):
                ok = run_hmmalign_add(new_df, tgt_msa, tgt_msa,
                                      tgt_info["profile_path"], "sequence", RANKS, "seq_id")
            if not ok:
                ok = run_mafft_add(new_df, tgt_msa, tgt_msa,
                                   "sequence", RANKS, threads, "seq_id")
        if ok:
            open(src_msa, "w").close()                 # empty the source MSA (contaminant)
            prof = src_info.get("profile_path")
            if prof and os.path.exists(prof):
                open(prof, "w").close()                # empty the profile TOO (unclassifiable downstream)
            emptied.append(src)
            n_rehomed += len(new_df)
            f["rehomed_to"] = tgt
    return sorted(set(emptied)), n_rehomed


def write_report(path, msa_tax, flagged, emptied, n_rehomed, n_profiles, flag_ranks):
    report = {"n_msa": len(msa_tax), "n_profiles": n_profiles, "flag_ranks": flag_ranks,
              "n_flagged": len(flagged), "n_rehomed": n_rehomed,
              "emptied_msa": emptied, "n_rehomed_msa": len(emptied), "flagged": flagged}
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(path, "w"), indent=2)
    by_rank = defaultdict(int)
    for f in flagged:
        by_rank[f["rank"]] += 1
    logger.info("Coarse contamination: %d MSAs flagged (%s) out of %d, %d re-homed.",
                len(flagged), dict(by_rank), len(msa_tax), n_rehomed)


# --------------------------------------------------------------------- subcommands
def cmd_prepare(args):
    """Build the reps target DB + PROFILE-LIBRARY shards + index (for the scattered search).

    hmmsearch iterates over the query PROFILES, so parallelism must split the profiles
    (not the reps): each shard searches its slice of the library against ALL reps, and
    every profile is processed exactly once across the shards (wall-clock ~ /n_shards).
    """
    msa_tax, msa_info, reps, profiles = build_index(args.family_manifests)
    profiles = sorted(set(profiles))
    # Only worth searching when there are >=2 reps AND at least one profile.
    reps = sorted(reps, key=lambda r: r[0]) if (len(reps) >= 2 and profiles) else []

    # reps.fa = the single shared target DB (all reps).
    Path(args.out_reps).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_reps, "w") as fh:
        for mid, seq in reps:
            fh.write(f">{mid}\n{seq}\n")

    index = {"n_profiles": len(profiles),
             "msas": {mid: {"tax": msa_tax[mid], **msa_info[mid]} for mid in msa_tax}}
    Path(args.out_index).parent.mkdir(parents=True, exist_ok=True)
    json.dump(index, open(args.out_index, "w"))

    shards_dir = Path(args.shards_dir)
    shards_dir.mkdir(parents=True, exist_ok=True)
    for old in shards_dir.glob("shard_*.hmm"):
        old.unlink()
    bs = max(1, args.batch_size)
    n_shards = 0
    if reps:
        for i in range(0, len(profiles), bs):
            with open(shards_dir / f"shard_{i // bs:04d}.hmm", "w") as out:
                for p in profiles[i:i + bs]:
                    out.write(open(p).read())
            n_shards += 1
    logger.info("Prepare: %d MSAs, %d reps, %d profiles -> %d shards (size %d).",
                len(msa_tax), len(reps), len(profiles), n_shards, bs)


def cmd_report(args):
    """Parse the scattered hits, flag, re-home, and write the report."""
    index = json.load(open(args.index, encoding="utf-8"))
    msa_tax = {mid: e["tax"] for mid, e in index.get("msas", {}).items()}
    msa_info = {mid: {"msa_path": e.get("msa_path"), "profile_path": e.get("profile_path")}
                for mid, e in index.get("msas", {}).items()}
    n_profiles = index.get("n_profiles", 0)

    hits = parse_hits(args.hits)
    flagged = flag_contaminants(hits, msa_tax, args.flag_ranks, args.min_score, args.top_k)
    emptied, n_rehomed = ([], 0) if args.no_rehome else rehome(flagged, msa_info, args.threads)
    write_report(args.output_report, msa_tax, flagged, emptied, n_rehomed, n_profiles, args.flag_ranks)


def cmd_run(args):
    """Legacy all-in-one (small scale / standalone): build index, one hmmsearch, flag, re-home."""
    msa_tax, msa_info, reps, profiles = build_index(args.family_manifests)
    if len(reps) < 2 or not profiles:
        write_report(args.output_report, msa_tax, [], [], 0, len(set(profiles)), args.flag_ranks)
        logger.info("Nothing to check (msa=%d, profiles=%d).", len(reps), len(profiles))
        return
    with tempfile.TemporaryDirectory(prefix="coarse_") as tmp:
        lib, repfa, tbl = f"{tmp}/lib.hmm", f"{tmp}/reps.fa", f"{tmp}/hits.tbl"
        with open(lib, "w") as out:
            for p in sorted(set(profiles)):
                out.write(open(p).read())
        with open(repfa, "w") as fh:
            for mid, seq in reps:
                fh.write(f">{mid}\n{seq}\n")
        subprocess.run(["hmmsearch", "--tblout", tbl, "-E", str(args.max_evalue),
                        "--cpu", str(args.threads), "--noali", lib, repfa],
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, check=True)
        hits = parse_hits([tbl])
    flagged = flag_contaminants(hits, msa_tax, args.flag_ranks, args.min_score, args.top_k)
    emptied, n_rehomed = ([], 0) if args.no_rehome else rehome(flagged, msa_info, args.threads)
    write_report(args.output_report, msa_tax, flagged, emptied, n_rehomed, len(set(profiles)), args.flag_ranks)


def parse_args():
    p = argparse.ArgumentParser(description="Coarse cross-kingdom contamination check (scatterable).")
    sub = p.add_subparsers(dest="cmd", required=True)

    pp = sub.add_parser("prepare", help="Build reps target DB + profile-library shards + index.")
    pp.add_argument("--family-manifests", nargs="*", default=[])
    pp.add_argument("--out-reps", required=True, help="reps.fa (all MSA representatives, shared target DB).")
    pp.add_argument("--out-index", required=True)
    pp.add_argument("--shards-dir", required=True, help="Dir of shard_*.hmm (profile-library shards).")
    pp.add_argument("--batch-size", type=int, default=200, help="Profiles per shard (fewer -> more parallel).")
    pp.add_argument("--log-file")
    pp.set_defaults(func=cmd_prepare)

    pr = sub.add_parser("report", help="Parse scattered hits, flag, re-home, write report.")
    pr.add_argument("--hits", nargs="*", default=[])
    pr.add_argument("--index", required=True)
    pr.add_argument("--output-report", required=True)
    pr.add_argument("--flag-ranks", nargs="+", default=["kingdom", "phylum"])
    pr.add_argument("--min-score", type=float, default=80.0)
    pr.add_argument("--top-k", type=int, default=5)
    pr.add_argument("--threads", type=int, default=8)
    pr.add_argument("--no-rehome", action="store_true")
    pr.add_argument("--log-file")
    pr.set_defaults(func=cmd_report)

    rn = sub.add_parser("run", help="Legacy all-in-one (small scale / standalone).")
    rn.add_argument("--family-manifests", nargs="*", default=[])
    rn.add_argument("--output-report", required=True)
    rn.add_argument("--flag-ranks", nargs="+", default=["kingdom", "phylum"])
    rn.add_argument("--min-score", type=float, default=80.0)
    rn.add_argument("--top-k", type=int, default=5)
    rn.add_argument("--max-evalue", type=float, default=1e-5)
    rn.add_argument("--threads", type=int, default=8)
    rn.add_argument("--no-rehome", action="store_true")
    rn.add_argument("--log-file")
    rn.set_defaults(func=cmd_run)
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if getattr(args, "log_file", None):
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        logger.addHandler(logging.FileHandler(args.log_file, mode="w"))
    args.func(args)


if __name__ == "__main__":
    main()
