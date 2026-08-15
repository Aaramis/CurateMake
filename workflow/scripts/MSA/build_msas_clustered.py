#!/usr/bin/env python3
"""
Clustered MSA construction for ONE barcode — the sequence-first replacement for the
(barcode, family) scatter + gated climb.

Rationale (validated in 10.DOI_CurateMake/simulation/SIMULATION_REPORT.md): grouping sequences by
their declared FAMILY silos a contaminant inside the clade of its (possibly wrong) label, so SATIVA
never sees it next to its true relatives and cannot correct it. Grouping by SEQUENCE instead — one
clustering over the whole barcode — puts a mislabelled read with its true kin, so cross-kingdom /
cross-clade contamination becomes visible. Taxonomy is then used only to MERGE clusters (not to
partition), up to a size cap that keeps each SATIVA MSA feasible (SATIVA tolerates very gappy
alignments, so size — not occupancy — is the real bound).

Pipeline (per barcode):
  1. mmseqs2 cluster the whole barcode (linear, scales to millions).
  2. Absorb clusters with < min_ref sequences into the SEQUENCE-nearest larger cluster (mmseqs
     search of small reps vs large reps).
  3. MSA per cluster: MAFFT -> trimAl -gappyout (trims columns, never drops sequences).
  4. Merge clusters that SHARE a taxonomy label (union-find), packing to <= max_group sequences.
  5. Pad any under-sized group with the taxonomically-farthest spare clusters (cross-kingdom) up to
     ~max_group, so every SATIVA MSA has a broad reference frame to place outliers against.
  6. Global-align each merged group (MAFFT) + trimAl, emit <out>/msa/<msa_id>.fasta and a
     SATIVA-compatible manifest {tasks:[{barcode, tax_level, taxon, taxon_label, msa_path, n}]}.

Reads full.csv (columns: seq_id, sequence, kingdom..species, barcode) in chunks.
"""
import argparse, json, logging, os, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
GROUP_RANKS = ["kingdom", "phylum", "class"]
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("clustered_msa")


def run(cmd, **kw):
    kw.setdefault("check", True)
    kw.setdefault("stdout", subprocess.DEVNULL)
    kw.setdefault("stderr", subprocess.DEVNULL)
    return subprocess.run(cmd, **kw)


def load_barcode(full_csv, barcode, chunksize=200000):
    """Yield (seq_id, seq, tax_dict) for one barcode from full.csv (pandas chunks)."""
    import pandas as pd
    keep = ["seq_id", "sequence", *RANKS, "barcode"]
    for chunk in pd.read_csv(full_csv, chunksize=chunksize, dtype=str, keep_default_na=False):
        cols = [c for c in keep if c in chunk.columns]
        chunk = chunk[cols]
        if "barcode" in chunk.columns:
            chunk["barcode"] = chunk["barcode"].replace("", "ITS")
        else:
            chunk["barcode"] = "ITS"
        chunk = chunk[chunk["barcode"] == barcode]
        for r in chunk.itertuples(index=False):
            d = r._asdict()
            if d.get("sequence"):
                yield d["seq_id"], d["sequence"], {rk: d.get(rk, "") or "NA" for rk in RANKS}


def write_fasta(ids, seqs, path):
    with open(path, "w") as f:
        for sid in ids:
            f.write(f">{sid}\n{seqs[sid]}\n")


def taxed_header(sid, tax):
    """Header in the pipeline-canonical format SATIVA/auto_sativa expects:
    >seq_id|kingdom=..|phylum=..|...|species=Genus_species

    The '|' separator + spaces->'_' escaping are load-bearing: MAFFT/trimAl truncate
    FASTA names at the first whitespace (a raw 'species=Genus species' would lose the
    epithet), and auto_sativa splits the header on '|' (a ';'-separated header parses to
    an all-'NA' .tax, so SATIVA would receive no reference taxonomy)."""
    parts = [sid]
    for r in RANKS:
        v = (tax.get(r) or "").strip().replace(" ", "_")
        if v:
            parts.append(f"{r}={v}")
    return "|".join(parts)


def mmseqs_cluster(fasta, out_dir, ident, cov, threads, env):
    """mmseqs easy-cluster -> {rep_id: [member_ids]}."""
    pref = out_dir / "mm"
    tmp = out_dir / "mmtmp"; tmp.mkdir(exist_ok=True)
    run(["conda", "run", "-n", env, "mmseqs", "easy-cluster", str(fasta), str(pref), str(tmp),
         "--min-seq-id", str(ident), "-c", str(cov), "--threads", str(threads), "-v", "1"])
    clusters = defaultdict(list)
    with open(f"{pref}_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            clusters[rep].append(mem)
    return clusters


def absorb_small(clusters, seqs, min_ref, out_dir, ident, threads, env):
    """Reassign members of clusters with < min_ref members to the sequence-nearest large cluster."""
    big = {r: ms for r, ms in clusters.items() if len(ms) >= min_ref}
    small = {r: ms for r, ms in clusters.items() if len(ms) < min_ref}
    if not big or not small:
        return clusters
    big_reps = out_dir / "big_reps.fasta"; small_reps = out_dir / "small_reps.fasta"
    write_fasta(list(big), seqs, big_reps)
    write_fasta(list(small), seqs, small_reps)
    res = out_dir / "absorb.m8"; tmp = out_dir / "abtmp"; tmp.mkdir(exist_ok=True)
    try:
        run(["conda", "run", "-n", env, "mmseqs", "easy-search", str(small_reps), str(big_reps),
             str(res), str(tmp), "--min-seq-id", str(max(0.5, ident - 0.2)), "--threads",
             str(threads), "-v", "1"])
    except subprocess.CalledProcessError:
        res = None
    nearest, best = {}, {}
    if res and res.exists():
        for line in open(res):                  # m8: query, target, pident, ...
            f = line.rstrip("\n").split("\t")
            q, t, pid = f[0], f[1], float(f[2])
            if pid > best.get(q, 0.0):
                best[q] = pid; nearest[q] = t
    merged = {r: list(ms) for r, ms in big.items()}
    for r, ms in small.items():
        tgt = nearest.get(r)
        if tgt in merged:
            merged[tgt].extend(ms)
        else:
            merged[r] = ms                      # no hit -> keep as its own (tiny) cluster
    return merged


def _subsample(members, k):
    """Deterministic representative subsample of a cluster: sorted, then evenly spaced, <= k.
    Clusters are near-identical (>= cluster_id, e.g. 97%), so any k reps span their variation and a
    giant same-species cluster need not blow one MSA past the feasibility cap. ALL members are still
    carried for accounting / verdict propagation (the non-sampled ones inherit the cluster's verdict).
    k <= 0 disables the cap."""
    ms = sorted(members)
    if k <= 0 or len(ms) <= k:
        return ms
    step = len(ms) / k
    return [ms[int(i * step)] for i in range(k)]


def group_clusters(clusters, tax, max_group, max_cluster_reps):
    """Pack CLUSTERS into groups whose MSA stays <= max_group sequences, INTERLEAVING kingdoms so
    every bounded group spans kingdoms. Each cluster contributes at most `max_cluster_reps`
    representatives to its MSA; ALL its members are carried for accounting.

    Kingdom-interleaved packing is load-bearing (validated at scale, cf. SIMULATION_REPORT.md /
    compare_corrupt/large): a contaminant grouped by its WRONG label still meets its true-kingdom kin
    in the same bounded MSA, so SATIVA can move it across kingdoms — clade-coherent packing separates
    them and misses cross-kingdom contamination (11.8% vs 54.4% at 41k reads).

    We pack clusters directly and ENFORCE the cap on the REP count. A union-find that merges every
    pair of clusters sharing a label does NOT bound size — contaminations are exactly the cross-links
    that chain all clades into one giant component (empirically: 331 cross-kingdom contaminants
    collapsed 41 k reads into a single 40 k 'group'). Cutting a new group at the cap guarantees every
    SATIVA MSA stays feasible.

    Returns a list of groups; each group is a list of (reps, all_members) cluster units. Deterministic
    given a deterministic cluster set (sorted members, tie-broken orderings)."""
    def kkey(ms):
        return Counter(tuple(tax[m].get(rk, "NA") for rk in GROUP_RANKS) for m in ms).most_common(1)[0][0]
    # unit = (reps_for_msa, all_members); reps subsampled + sorted -> deterministic
    units = []
    for members in clusters.values():
        allm = sorted(members)
        units.append((_subsample(allm, max_cluster_reps), allm))
    # interleave kingdoms so every bounded group spans kingdoms (broad reference frame -> a
    # mis-grouped contaminant still meets its true-kingdom kin)
    by_k = defaultdict(list)
    for u in units:
        by_k[kkey(u[0])[0]].append(u)
    for k in by_k:                                  # deterministic order within a kingdom
        by_k[k].sort(key=lambda u: (-len(u[0]), u[0][0]))
    order, ks = [], sorted(by_k)
    while any(by_k[k] for k in ks):
        for k in ks:
            if by_k[k]:
                order.append(by_k[k].pop(0))
    groups, cur, cur_n = [], [], 0
    for u in order:
        r = len(u[0])
        if cur and cur_n + r > max_group:           # cap enforced (on MSA size) between clusters
            groups.append(cur); cur, cur_n = [], 0
        cur.append(u); cur_n += r
    if cur:
        groups.append(cur)
    return groups


def mafft_trim(ids, seqs, tax, tag, out_dir, threads, env):
    gdir = out_dir / "work" / tag; gdir.mkdir(parents=True, exist_ok=True)
    raw = gdir / "raw.fasta"
    with open(raw, "w") as f:
        for sid in ids:
            f.write(f">{taxed_header(sid, tax[sid])}\n{seqs[sid]}\n")
    aln, trimmed = gdir / "aln.fasta", gdir / "trim.fasta"
    with open(aln, "w") as fo:
        run(["conda", "run", "-n", env, "mafft", "--retree", "2", "--maxiterate", "0",
             "--quiet", "--thread", str(threads), "--inputorder", str(raw)], stdout=fo)
    # trimAl -gappyout : trim columns, keep every sequence
    run(["conda", "run", "-n", env, "trimal", "-in", str(aln), "-out", str(trimmed), "-gappyout"])
    return trimmed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-csv", required=True)
    ap.add_argument("--barcode", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cluster-id", type=float, default=0.90)
    ap.add_argument("--cluster-cov", type=float, default=0.80)
    ap.add_argument("--min-ref", type=int, default=4)
    ap.add_argument("--max-group", type=int, default=10000)
    ap.add_argument("--max-cluster-reps", type=int, default=0,
                    help="cap the reps ONE cluster contributes to its MSA (0 = no cap). Bounds the MSA "
                         "when a giant near-identical cluster would otherwise blow past --max-group.")
    ap.add_argument("--sativa-min-seqs", type=int, default=4)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--mmseqs-env", default="mmseqs")
    ap.add_argument("--tool-env", default="its-pipeline-dev",
                    help="conda env providing mafft + trimal")
    args = ap.parse_args()
    out = Path(args.out_dir); (out / "msa").mkdir(parents=True, exist_ok=True)

    seqs, tax = {}, {}
    for sid, s, t in load_barcode(args.full_csv, args.barcode):
        seqs[sid] = s; tax[sid] = t
    log.info("barcode %s: %d sequences", args.barcode, len(seqs))
    if not seqs:
        json.dump({"tasks": [], "stats": {"barcode": args.barcode, "n_sequences": 0, "n_clusters": 0,
                                          "n_groups": 0, "n_msa": 0, "ref_seqs": 0, "residual": 0}},
                  open(out / "sativa_manifest.json", "w"))
        return

    # sort the input -> mmseqs/grouping determinism (reproducibility, cf. INTEGRATION_DESIGN §5.1)
    allfa = out / "all.fasta"; write_fasta(sorted(seqs), seqs, allfa)
    clusters = mmseqs_cluster(allfa, out, args.cluster_id, args.cluster_cov, args.threads, args.mmseqs_env)
    log.info("clusters: %d", len(clusters))
    clusters = absorb_small(clusters, seqs, args.min_ref, out, args.cluster_id, args.threads, args.mmseqs_env)
    log.info("clusters after absorbing <%d: %d", args.min_ref, len(clusters))
    groups = group_clusters(clusters, tax, args.max_group, args.max_cluster_reps)
    log.info("merged groups (cap %d, <=%d reps/cluster): %d", args.max_group, args.max_cluster_reps, len(groups))

    tasks = []
    ref_seqs = residual = 0
    members_out = open(out / "group_members.tsv", "w")   # tag<TAB>member -> verdict propagation
    members_out.write("tag\tmember_seq_id\n")
    for i, units in enumerate(groups):
        msa_members = [m for reps, _ in units for m in reps]
        all_members = [m for _, allm in units for m in allm]
        if len(msa_members) < args.sativa_min_seqs:
            residual += len(all_members)                 # too few refs for SATIVA -> unprocessed
            continue
        tag = f"{args.barcode}__clustered__{i:04d}"
        trimmed = mafft_trim(msa_members, seqs, tax, tag, out, args.threads, args.tool_env)
        msa_path = out / "msa" / f"{tag}.fasta"
        os.replace(trimmed, msa_path)
        label = Counter(f"{tax[m].get('kingdom','NA')};{tax[m].get('phylum','NA')}"
                        for m in msa_members).most_common(1)[0][0]
        tasks.append({"barcode": args.barcode, "tax_level": "clustered", "taxon": tag,
                      "taxon_label": label, "msa_path": str(msa_path),
                      "n": len(all_members), "n_msa_seqs": len(msa_members)})
        for m in all_members:
            members_out.write(f"{tag}\t{m}\n")
        ref_seqs += len(all_members)
        log.info("  %s: %d reps / %d members -> %s", tag, len(msa_members), len(all_members), msa_path.name)
    members_out.close()
    stats = {"barcode": args.barcode, "n_sequences": len(seqs), "n_clusters": len(clusters),
             "n_groups": len(groups), "n_msa": len(tasks), "ref_seqs": ref_seqs, "residual": residual}
    json.dump({"tasks": tasks, "stats": stats}, open(out / "sativa_manifest.json", "w"), indent=2)
    log.info("wrote %d SATIVA tasks (ref_seqs=%d, residual=%d)", len(tasks), ref_seqs, residual)


if __name__ == "__main__":
    main()
