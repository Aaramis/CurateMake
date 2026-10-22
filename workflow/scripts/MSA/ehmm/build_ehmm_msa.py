#!/usr/bin/env python3
"""Backbone + ensemble-HMM MSA for ONE barcode — the only MSA strategy.

Reads full.csv, and emits the manifest contract msa_groups.smk merges:
  <out>/msa/<msa_id>.fasta, <out>/group_members.tsv, <out>/sativa_manifest.json

Four stages (prototype + measurements in 20.eHMM_MSA_prototype/README.md):
  1. BACKBONE — representatives of the largest mmseqs clusters (--backbone-source cluster, the
     default) aligned once; a (kingdom,phylum) stratified sample is the legacy alternative. A profile
     HMM models SEQUENCE space, so the backbone must cover sequence space: on the real ITS2 the
     stratified sample assigns 66.1% of reads and cluster representatives 81.1%. This alignment IS the single
     column frame; nothing is ever realigned afterwards. Only columns occupied by >= --min-col-occ
     of the backbone are kept: a raw cross-kingdom ITS alignment spreads ~470 bp reads over ~6700
     columns (7% occupancy) and an HMM built on that no longer recognises even its own training
     rows (measured: 5/17). At >=10% (~1000 columns) nhmmer recognises ~77% of reads.
  2. PROJECT — one HMM per clade bucket, built over the WHOLE frame (Stockholm + `#=GC RF` +
     `hmmbuild --hand`, because hmmbuild silently drops all-gap columns even with --symfrac 0.0,
     which would make bucket HMMs narrower than the frame). nhmmer picks the best bucket AND the
     strand; hmmalign then places the read and only match columns are kept, so every read becomes a
     row of exactly `frame` columns in one shared coordinate system.
  3. RECOVER (optional) — reads nhmmer misses are short fragments (median 166 bp vs 457). mmseqs
     finds them against the already-projected reads; they inherit that read's bucket and strand.
  4. SLICE — kingdom round-robin packing into groups <= max_group, optionally with the backbone
     injected into every group as a taxonomic scaffold, since a contaminant is only detectable when
     its TRUE clade sits in the same MSA (a plain round-robin exhausts the small kingdoms first and
     leaves the last groups with 2 kingdoms out of 5).
"""
import argparse
import json
import logging
import os
import random
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("ehmm_msa")

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def run(cmd, stdout=None):
    r = subprocess.run([str(c) for c in cmd], stdout=stdout,
                       stderr=subprocess.PIPE, text=(stdout is None))
    if r.returncode != 0:
        sys.stderr.write((r.stderr or "")[-2000:] if isinstance(r.stderr, str) else "")
        raise RuntimeError("command failed: " + " ".join(str(c) for c in cmd[:4]))
    return r


def tool(binary, hmmer_bin):
    return str(Path(hmmer_bin) / binary) if hmmer_bin else binary


def mmseqs_cmd(mmseqs_env):
    """How to invoke mmseqs: directly when it is on PATH, through `conda run` otherwise.

    Under Snakemake --use-conda the rule's own environment provides it (envs/msa.yaml
    declares mmseqs2), so the plain binary is the normal case. The `conda run` fallback
    keeps a deployment that relies on a hand-made named environment working. Before this,
    the named environment was the ONLY path: on any machine without an env literally
    called "mmseqs" the step died with "Not a conda environment", which is not something
    a reader of the repository could have guessed, since no env file declared it.
    """
    direct = shutil.which("mmseqs")
    if direct:
        return [direct]
    if mmseqs_env:
        log.warning("mmseqs absent du PATH, repli sur `conda run -n %s`", mmseqs_env)
        return ["conda", "run", "-n", mmseqs_env, "mmseqs"]
    raise RuntimeError("mmseqs introuvable sur le PATH et aucun --mmseqs-env fourni")


def read_fasta(path):
    h, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(buf)
            h, buf = line[1:].rstrip("\n"), []
        else:
            buf.append(line.strip())
    if h is not None:
        yield h, "".join(buf)


def load_barcode(full_csv, barcode, chunksize=200000):
    """(seq_id, sequence, taxonomy) for one barcode."""
    import pandas as pd
    keep = ["seq_id", "sequence", *RANKS, "barcode"]
    for chunk in pd.read_csv(full_csv, chunksize=chunksize, dtype=str, keep_default_na=False):
        cols = [c for c in keep if c in chunk.columns]
        chunk = chunk[cols]
        chunk["barcode"] = chunk["barcode"].replace("", "ITS") if "barcode" in chunk.columns else "ITS"
        chunk = chunk[chunk["barcode"] == barcode]
        for r in chunk.itertuples(index=False):
            d = r._asdict()
            if d.get("sequence"):
                yield d["seq_id"], d["sequence"], {rk: (d.get(rk) or "NA") for rk in RANKS}


def taxed_header(sid, tax):
    """`id|kingdom=..|..` — the format auto_sativa parses."""
    parts = [sid]
    for r in RANKS:
        v = (tax.get(r) or "").strip().replace(" ", "_")
        if v and v != "NA":
            parts.append(f"{r}={v}")
    return "|".join(parts)


# ----------------------------------------------------------------------------- 1. backbone
def cluster_backbone(seqs, work, target, ident, threads, mmseqs_env):
    """Backbone rows = representatives of the largest mmseqs clusters.

    A profile HMM models SEQUENCE space, so the backbone has to cover sequence space. Stratifying on
    the declared (kingdom, phylum) instead covers LABEL space -- the very thing this pipeline exists
    to distrust -- and with only ~60 strata the proportional quota handed 81% of the rows to Fungi.
    Measured on the real ITS2 (1.4M seqs, 1300-row budget, subsets of 10, same 10k queries):
    taxonomy-stratified assigns 66.1% of reads, cluster representatives 81.1%. Taking the LARGEST
    clusters also beats spreading the picks over the size-ranked list (81.1% vs 70.2%): covering the
    bulk of the data pays more than covering rare sequence space. The frame comes out wider too
    (512 vs 468 columns), because coherent representatives align to each other more tightly.
    """
    fa = work / "all_for_clust.fasta"
    with open(fa, "w") as f:
        for sid, seq in seqs.items():
            f.write(f">{sid}\n{seq}\n")
    pref = work / "lc"
    run(mmseqs_cmd(mmseqs_env) + ["easy-linclust", fa, pref, work / "lctmp",
         "--min-seq-id", ident, "-c", "0.8", "--kmer-per-seq", "300", "--threads", threads, "-v", "1"])
    size = defaultdict(int)
    for line in open(f"{pref}_cluster.tsv"):
        size[line.split("\t")[0]] += 1
    ranked = sorted(size.items(), key=lambda kv: (-kv[1], kv[0]))
    log.info("linclust: %d clusters (largest %d seqs, %d singletons)",
             len(ranked), ranked[0][1] if ranked else 0, sum(1 for _, n in ranked if n == 1))
    return [r for r, _ in ranked[:target] if r in seqs]


def taxonomy_backbone(seqs, tax, target, min_per_stratum):
    """Legacy backbone: (kingdom, phylum) stratified sample. Kept for comparison only."""
    import statistics
    strata = defaultdict(list)
    for sid in seqs:
        strata[(tax[sid].get("kingdom", "NA"), tax[sid].get("phylum", "NA"))].append(sid)
    med = statistics.median(len(s) for s in seqs.values())
    total = sum(len(v) for v in strata.values())
    picked = []
    for key in sorted(strata):
        pool = sorted(strata[key], key=lambda sid: (abs(len(seqs[sid]) - med), sid))
        picked.extend(pool[:max(min_per_stratum, round(target * len(pool) / total))])
    return picked


def build_backbone(seqs, tax, work, target, min_per_stratum, aligner, threads,
                   source, cluster_id, mmseqs_env):
    picked = (cluster_backbone(seqs, work, target, cluster_id, threads, mmseqs_env)
              if source == "cluster" else
              taxonomy_backbone(seqs, tax, target, min_per_stratum))
    kings = Counter(tax[s].get("kingdom", "NA") for s in picked)
    log.info("backbone (%s): %d rows, kingdoms %s", source, len(picked), dict(kings.most_common(4)))
    raw = work / "backbone_raw.fasta"
    with open(raw, "w") as f:
        for sid in picked:
            f.write(f">{sid}\n{seqs[sid]}\n")
    aln = work / "backbone.aln.fasta"
    with open(aln, "w") as fo:
        run([aligner, "--adjustdirection", "--thread", threads, "--quiet", raw], stdout=fo)
    return aln, set(picked)


def frame_from_backbone(aln_path, min_col_occ):
    rows = {h.split()[0].removeprefix("_R_").split(";")[0]: s.upper()
            for h, s in read_fasta(aln_path)}
    raw_w = len(next(iter(rows.values())))
    occ = [0] * raw_w
    for row in rows.values():
        for j, c in enumerate(row):
            if c != "-":
                occ[j] += 1
    n = len(rows)
    frame = [j for j, o in enumerate(occ) if o / n >= min_col_occ]
    rows = {sid: "".join(row[j] for j in frame) for sid, row in rows.items()}
    log.info("frame: %d of %d backbone columns kept (occupancy >= %.0f%%)",
             len(frame), raw_w, 100 * min_col_occ)
    return rows, len(frame)


# ----------------------------------------------------------------------------- 2. projection
def make_buckets(bb_rows, tax, bucket_mode, bucket_rank, min_bucket, subset_size):
    """Split the backbone into the profiles of the ensemble.

    `upp` (default) is the granularity UPP actually uses: many small, tight profiles rather than one
    per clade. Measured on the real ITS1 (1.4 M seqs, same backbone, same 10 000 queries): 4
    kingdom-level HMMs assign 20.4% of reads, 108 UPP-style HMMs assign 72.8%. A single HMM asked to
    model 6 781 fungal genera is too diluted to recognise anything -- the same dilution that already
    bit us one level up, on the column frame. Subsets are cut INSIDE (kingdom, phylum, order) so each
    profile stays taxonomically tight; a GLOBAL profile is always kept as a fallback.
    """
    if bucket_mode == "clade":
        b = defaultdict(list)
        for sid, row in bb_rows.items():
            b[tax.get(sid, {}).get(bucket_rank, "NA")].append((sid, row))
        out = {k: v for k, v in sorted(b.items()) if len(v) >= min_bucket}
    else:
        by_clade = defaultdict(list)
        for sid, row in bb_rows.items():
            t = tax.get(sid, {})
            by_clade[(t.get("kingdom", "NA"), t.get("phylum", "NA"), t.get("order", "NA"))].append(
                (sid, row))
        out = {}
        for clade, rows in sorted(by_clade.items()):
            rows.sort()
            for i in range(0, len(rows), subset_size):
                chunk = rows[i:i + subset_size]
                if len(chunk) < 3 and out:            # a stray tail joins the previous profile
                    out[list(out)[-1]].extend(chunk)
                else:
                    name = f"{clade[0]}_{clade[1]}_{clade[2]}_{i // subset_size}"
                    out[name.replace(" ", "_").replace("/", "_")[:60]] = chunk
    out["GLOBAL"] = [(sid, row) for sid, row in bb_rows.items()]
    return out


def project(bb_rows, width, tax, seqs, work, hmmer_bin, bucket_mode, bucket_rank, min_bucket,
            subset_size, threads, nhmmer_jobs=4, nhmmer_threads=None):
    keep = make_buckets(bb_rows, tax, bucket_mode, bucket_rank, min_bucket, subset_size)
    log.info("ensemble: %d profiles (mode=%s, subset_size=%d)", len(keep), bucket_mode, subset_size)

    hmms = {}
    for k, rows in keep.items():
        sto = work / f"bucket_{k}.sto"
        with open(sto, "w") as f:
            f.write("# STOCKHOLM 1.0\n")
            for sid, row in rows:
                f.write(f"{sid}  {row}\n")
            f.write(f"#=GC RF  {'x' * width}\n//\n")
        hmm = work / f"bucket_{k}.hmm"
        run([tool("hmmbuild", hmmer_bin), "--dna", "--hand", "--cpu", threads,
             "-n", f"b_{k}"[:60], hmm, sto])
        hmms[k] = hmm

    queries = work / "queries.fasta"
    with open(queries, "w") as f:
        for sid, s in seqs.items():
            f.write(f">{sid}\n{s}\n")

    # nhmmer scales well but NOT linearly: measured 2026-08-25 on one profile x 100k reads,
    # --cpu 4/8/16/32 gave 238.8/116.8/63.3/38.7 s, i.e. 107 reads/s/core at 8 threads against
    # 80.8 at 32. Within the 150-core cgroup quota, many small processes therefore beat a few
    # wide ones: 18 x 8 threads moves ~15400 reads/s where 4 x 32 moved ~10300.
    n_threads = int(nhmmer_threads or threads)

    # One nhmmer per profile, run CONCURRENTLY and RESUMABLY.
    #
    # Sequential, this loop is the whole cost of the step: on the ITS barcode a single profile
    # takes ~14 min (2.29M reads x 624 bp against a 940-column model), so 138 profiles is well
    # over a day while one nhmmer at --cpu 32 leaves ~120 of the 150 cgroup cores idle.
    # The profiles are independent -- each writes its own table and nothing is shared -- so a
    # small pool of subprocesses is all that is needed; parsing still happens afterwards, in
    # this thread, exactly as before.
    #
    # Resumable: HMMER terminates a complete table with a "# [ok]" line, which is a reliable
    # marker that the profile finished. A run interrupted mid-profile leaves that table
    # truncated, so it is recomputed while every finished profile is reused.
    def _table_complete(path):
        try:
            size = path.stat().st_size
        except OSError:
            return False
        with open(path, "rb") as fh:
            fh.seek(max(0, size - 64))
            return b"# [ok]" in fh.read()

    def _scan(item):
        k, hmm = item
        tbl = work / f"nhmmer_{k}.tbl"
        if _table_complete(tbl):
            return k, tbl, True
        run([tool("nhmmer", hmmer_bin), "--cpu", n_threads, "--tblout", tbl, "-o", os.devnull,
             hmm, queries])
        return k, tbl, False

    workers = max(1, min(int(nhmmer_jobs), len(hmms)))
    log.info("nhmmer: %d profiles, %d concurrent x %s threads", len(hmms), workers, n_threads)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        scanned = list(pool.map(_scan, list(hmms.items())))
    reused = sum(1 for _, _, r in scanned if r)
    if reused:
        log.info("nhmmer: %d/%d profile(s) reused from a previous run", reused, len(scanned))

    best = {}
    for k, tbl, _ in scanned:
        for line in open(tbl):
            if line.startswith("#"):
                continue
            f = line.split()
            qid, strand, score = f[0].split(";")[0], f[11], float(f[13])
            if qid not in best or score > best[qid][0]:
                best[qid] = (score, k, strand)
    log.info("nhmmer assigned %d/%d reads", len(best), len(seqs))
    return hmms, best


def hmmalign_bucket(hmm, items, work, tag, width, hmmer_bin, out_rows):
    if not items:
        return 0
    qf = work / f"q_{tag}.fasta"
    with open(qf, "w") as f:
        for sid, s in items:
            f.write(f">{sid}\n{s}\n")
    af = work / f"a_{tag}.afa"
    run([tool("hmmalign", hmmer_bin), "--dna", "--outformat", "afa", "-o", af, hmm, qf])
    n = 0
    for h, s in read_fasta(af):
        row = "".join(c for c in s if c == "-" or c.isupper())
        if len(row) == width:
            out_rows[h.split()[0]] = row
            n += 1
    return n


def hmmalign_buckets(per_bucket, hmms, work, width, hmmer_bin, out_rows, jobs, tag_prefix=""):
    """Project every bucket into the shared column frame, buckets running in parallel.

    hmmalign has no --cpu in HMMER 3.4 (checked on the shipped binary; nhmmer does), so one
    bucket occupies exactly one core. Measured on the 2026-08-25 ITS build: 138 buckets took
    17.0 h of the 33.2 h total, on a single core, while the node quota is 150. Buckets are
    independent -- own profile, own temp files, disjoint seq_ids -- so they parallelise the
    same way the nhmmer scan above already does. Each worker fills a local dict; the merge
    happens here, so nothing relies on dict writes being atomic.
    """
    if not per_bucket:
        return 0
    workers = max(1, min(int(jobs), len(per_bucket)))

    def one(item):
        bucket, items = item
        local = {}
        n = hmmalign_bucket(hmms[bucket], items, work, f"{tag_prefix}{bucket}",
                            width, hmmer_bin, local)
        return local, n

    total = 0
    log.info("hmmalign: %d bucket(s), %d concurrent (hmmalign is single-threaded)",
             len(per_bucket), workers)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for local, n in pool.map(one, list(per_bucket.items())):
            out_rows.update(local)
            total += n
    return total


# ----------------------------------------------------------------------------- 3. rescue
def map_bucket(name, hmms):
    """Translate a bucket name coming from ANOTHER barcode into a local one.

    Buckets are named from the lineage (bucket_<Kingdom>_<Phylum>_<Order>_<n>), so the
    three barcodes share most of them: measured on the 2026-08-25 run, 100 of ITS1's 136
    buckets and 104 of ITS2's 138 exist under the same name in ITS. Fall back on the same
    lineage with another index, then on the kingdom, then give up.
    """
    if name in hmms:
        return name
    base = name.rsplit("_", 1)[0]
    for k in hmms:
        if k.rsplit("_", 1)[0] == base:
            return k
    parts = name.split("_")
    if len(parts) > 2:
        prefix = "_".join(parts[:2])
        for k in sorted(hmms):
            if k.startswith(prefix):
                return k
    return None


def read_cross_refs(paths):
    """Load >seq_id|bucket records emitted by another barcode. Returns {seq_id: (seq, bucket)}."""
    out = {}
    for p in paths or []:
        if not p or not os.path.exists(p):
            continue
        sid = bucket = None
        chunks = []
        with open(p) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if sid:
                        out[sid] = ("".join(chunks), bucket)
                    head = line[1:]
                    sid, _, bucket = head.partition("|")
                    chunks = []
                else:
                    chunks.append(line)
        if sid:
            out[sid] = ("".join(chunks), bucket)
    return out


def recover(seqs, rows, best, hmms, work, width, hmmer_bin, mmseqs_env, min_id, min_cov, threads,
            cross=None, bucket_of=None, align_jobs=1):
    missing = [sid for sid in seqs if sid not in rows]
    if not missing:
        # An empty dict, NOT 0: the caller reads .get("rescued_cross_barcode") on whatever
        # comes back, and both are falsy where the result is only tested for truthiness.
        # Returning 0 here crashed the run whenever every read was projected, which never
        # happens on the full database and always happens on a small clean dataset.
        return {}
    q = work / "unassigned.fasta"
    with open(q, "w") as f:
        for sid in missing:
            f.write(f">{sid}\n{seqs[sid]}\n")
    # Targets: this barcode's own projected reads, plus -- when given -- the projected reads
    # of another barcode. ITS1 and ITS2 are exact substrings of their ITS parent (checked:
    # 200/200 both ways), and ITS projects ~100% of its reads against ~90% for ITS1/ITS2, so
    # the ITS frame is where the missing 267k reads can still be placed. --cov-mode 2 scores
    # coverage of the QUERY, so a 177nt ITS1 inside a 570nt ITS still passes -c 0.60.
    # ITS1 <-> ITS2 is deliberately NOT wired: they sit either side of 5.8S and share no
    # sequence at all (0 common 25-mers over 100 same-parent pairs), so any hit is noise.
    cross = cross or {}
    tgt = work / "targets.fasta"
    with open(tgt, "w") as f:
        for sid, row in rows.items():
            f.write(f">{sid}\n{row.replace('-', '')}\n")
        for sid, (seq, _bucket) in cross.items():
            f.write(f">x|{sid}\n{seq}\n")
    hits = work / "recover.m8"
    run(mmseqs_cmd(mmseqs_env) + ["easy-search", q, tgt, hits,
         work / "rec_tmp", "--search-type", "3", "--min-seq-id", min_id, "-c", min_cov,
         "--cov-mode", "2", "--max-seqs", "3", "--threads", threads, "-v", "1"])
    pick = {}
    from_cross = {}
    dropped = 0
    for line in open(hits):
        f = line.rstrip("\n").split("\t")
        qid, tid, bits = f[0].split(";")[0], f[1], float(f[11])
        is_cross = tid.startswith("x|")
        if is_cross:
            foreign = cross.get(tid[2:], (None, None))[1]
            bucket = map_bucket(foreign, hmms) if foreign else None
            if bucket is None:
                dropped += 1
                continue
        else:
            bucket = best.get(tid, (0, "GLOBAL", "+"))[1]
        if qid not in pick or bits > pick[qid][0]:
            pick[qid] = (bits, bucket, int(f[8]) > int(f[9]))
            from_cross[qid] = is_cross
    per = defaultdict(list)
    for sid, (_b, bucket, minus) in pick.items():
        per[bucket].append((sid, seqs[sid].translate(COMP)[::-1] if minus else seqs[sid]))
    before = set(rows)
    if bucket_of is not None:
        for bucket, items in per.items():
            for sid, _s in items:
                bucket_of[sid] = bucket
    n = hmmalign_buckets(per, hmms, work, width, hmmer_bin, rows, align_jobs, tag_prefix="rec_")
    # Split the rescue by where the winning target came from, so the manifest can report it the
    # way the other inter-barcode steps do (cf. barcode_reannotation.json's by_change).
    rescued = set(rows) - before
    n_cross = sum(1 for sid in rescued if from_cross.get(sid))
    stats = {"unassigned": len(missing), "rescued": n,
             "rescued_same_barcode": n - n_cross, "rescued_cross_barcode": n_cross,
             "cross_targets": len(cross), "cross_hits_no_local_bucket": dropped}
    if cross:
        log.info("cross-barcode targets: %d sequences, %d hit(s) dropped (no local bucket)",
                 len(cross), dropped)
        log.info("rescued %d reads via another barcode", n_cross)
    log.info("rescued %d of %d unassigned reads", n, len(missing))
    return stats


# ----------------------------------------------------------------------------- 4. slicing
def slice_groups(rows, tax, backbone_ids, out, barcode, max_group, min_group, scaffold,
                 n_input, unprojected, rescue=None, split_seed=0):
    by_k = defaultdict(list)
    for sid in sorted(rows):
        by_k[tax.get(sid, {}).get("kingdom", "NA")].append(sid)
    # split_seed != 0 shuffles within each kingdom before the round-robin, which produces a
    # DIFFERENT but equally legitimate packing: the group boundaries are an artefact of the
    # traversal order, not of biology. Measured on the 2026-08-25 run, the neighbourhood
    # drives the verdict hard -- of the 1200 scaffold sequences present in all 122 ITS
    # groups, 984 (89.4%) got a different proposed lineage depending on the group, with a
    # median of 6 distinct proposals. Re-slicing with several seeds and voting is the way to
    # separate a real correction from an artefact of who happened to be in the same task.
    # seed 0 keeps the historical deterministic order.
    if split_seed:
        rnd = random.Random(split_seed)
        for k in by_k:
            rnd.shuffle(by_k[k])
    order, ks = [], sorted(by_k)
    while any(by_k[k] for k in ks):
        for k in ks:
            if by_k[k]:
                order.append(by_k[k].pop(0))
    scaf = sorted(sid for sid in backbone_ids if sid in rows) if scaffold else []
    budget = max(min_group, max_group - len(scaf))
    groups, cur = [], []
    for sid in order:
        if len(cur) >= budget:
            groups.append(cur); cur = []
        cur.append(sid)
    if cur:
        groups.append(cur)
    if scaffold:
        log.info("scaffold: %d backbone refs added to each of the %d groups", len(scaf), len(groups))

    (out / "msa").mkdir(parents=True, exist_ok=True)
    tasks, ref_seqs, residual = [], 0, 0
    members = open(out / "group_members.tsv", "w")
    members.write("tag\tmember_seq_id\n")
    for i, g in enumerate(groups):
        if len(g) < min_group:
            residual += len(g)
            continue
        ids = g + [s for s in scaf if s not in set(g)]
        width = len(next(iter(rows.values())))
        keep_col = [False] * width
        for sid in ids:
            for j, c in enumerate(rows[sid]):
                if c != "-":
                    keep_col[j] = True
        cols = [j for j, k in enumerate(keep_col) if k]
        tag = f"{barcode}__groups__{i:04d}"
        path = out / "msa" / f"{tag}.fasta"
        with open(path, "w") as f:
            for sid in ids:
                f.write(">" + taxed_header(sid, tax.get(sid, {})) + "\n"
                        + "".join(rows[sid][j] for j in cols) + "\n")
        label = Counter(f"{tax[s].get('kingdom','NA')};{tax[s].get('phylum','NA')}"
                        for s in g).most_common(1)[0][0]
        tasks.append({"barcode": barcode, "tax_level": "groups", "taxon": tag,
                      "taxon_label": label, "msa_path": str(path),
                      "n": len(g), "n_msa_seqs": len(ids)})
        for sid in g:                                   # scaffold refs belong to their own group
            members.write(f"{tag}\t{sid}\n")
        ref_seqs += len(g)
        log.info("  %s: %d reads (+%d scaffold) x %d cols", tag, len(g), len(ids) - len(g), len(cols))
    members.close()
    # n_sequences MUST be the barcode's input count, and every read that reaches no MSA MUST land in
    # `residual` -- reporting only the projected reads made a third of the data vanish silently and
    # the downstream report then claimed 100% coverage.
    residual += len(unprojected)
    stats = {"barcode": barcode, "n_sequences": n_input, "n_clusters": len(rows),
             "n_groups": len(groups), "n_msa": len(tasks), "ref_seqs": ref_seqs,
             "residual": residual, "n_projected": len(rows), "n_unprojected": len(unprojected)}
    if rescue:
        stats["rescue"] = rescue
    if unprojected:
        with open(out / "unprojected_seq_ids.txt", "w") as f:
            f.write("\n".join(sorted(unprojected)) + "\n")
        pct = 100 * len(unprojected) / n_input
        (log.warning if pct >= 5 else log.info)(
            "%d of %d reads (%.1f%%) reached NO MSA and are reported as residual "
            "(ids in unprojected_seq_ids.txt)", len(unprojected), n_input, pct)
    json.dump({"tasks": tasks, "stats": stats}, open(out / "sativa_manifest.json", "w"), indent=2)
    log.info("wrote %d SATIVA tasks (ref_seqs=%d, residual=%d)", len(tasks), ref_seqs, residual)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-csv", required=True)
    ap.add_argument("--barcode", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--backbone-target", type=int, default=1200)
    ap.add_argument("--backbone-source", default="cluster", choices=["cluster", "taxonomy"])
    ap.add_argument("--backbone-cluster-id", type=float, default=0.90)
    ap.add_argument("--backbone-min-per-stratum", type=int, default=3)
    ap.add_argument("--backbone-aligner", default="mafft-linsi")
    ap.add_argument("--min-col-occ", type=float, default=0.10)
    ap.add_argument("--bucket-mode", default="upp", choices=["upp", "clade"])
    ap.add_argument("--subset-size", type=int, default=25)
    ap.add_argument("--bucket-rank", default="kingdom")
    ap.add_argument("--min-bucket", type=int, default=20)
    ap.add_argument("--recover", action="store_true")
    ap.add_argument("--split-seed", type=int, default=0,
                    help="shuffle sequences within each kingdom before packing groups; 0 keeps "
                         "the historical deterministic order. Different seeds give different "
                         "neighbourhoods for the same reads -- the basis for a consensus vote.")
    ap.add_argument("--hmmalign-jobs", type=int, default=16,
                    help="buckets aligned concurrently; hmmalign is single-threaded so this is "
                         "the only way to use more than one core for the projection")
    ap.add_argument("--cross-ref", action="append", default=[],
                    help="projected_refs.fasta of another barcode, used as extra rescue targets")
    ap.add_argument("--recover-min-id", type=float, default=0.75)
    ap.add_argument("--recover-min-qcov", type=float, default=0.60)
    ap.add_argument("--scaffold", action="store_true")
    ap.add_argument("--max-group", type=int, default=20000)
    ap.add_argument("--min-group", type=int, default=4)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--nhmmer-threads", type=int, default=0,
                    help="--cpu given to each nhmmer; 0 falls back to --threads. Prefer many "
                         "small processes: per-core throughput peaks around 8 (see scan()).")
    ap.add_argument("--nhmmer-jobs", type=int, default=4,
                    help="profiles scanned concurrently, each with --threads. One nhmmer at "
                         "--cpu 32 leaves ~120 of the 150 cgroup cores idle; 4 x 32 fills them.")
    ap.add_argument("--hmmer-bin", default="")
    ap.add_argument("--mmseqs-env", default="mmseqs")
    ap.add_argument("--tool-env", default="its-pipeline-dev")
    a = ap.parse_args()

    work = Path(a.work_dir); work.mkdir(parents=True, exist_ok=True)
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)

    # mafft/mmseqs live in a conda env; resolve it onto PATH once (a per-call `conda run` costs 1-2 s)
    # Optional convenience: if these named environments exist, put them on PATH once, so
    # a per-call `conda run` does not cost 1-2 s. Absent ones are skipped and logged, not
    # fatal: under --use-conda the rule's own environment already provides the binaries.
    for env in {a.tool_env, a.mmseqs_env}:
        if not env:
            continue
        p = subprocess.run(["conda", "run", "-n", env, "printenv", "CONDA_PREFIX"],
                           capture_output=True, text=True)
        if p.stdout.strip():
            os.environ["PATH"] = f"{p.stdout.strip()}/bin{os.pathsep}" + os.environ.get("PATH", "")
        else:
            log.info("environnement nomme '%s' absent, on s en passe", env)

    seqs, tax = {}, {}
    for sid, s, t in load_barcode(a.full_csv, a.barcode):
        seqs[sid] = s; tax[sid] = t
    log.info("barcode %s: %d sequences", a.barcode, len(seqs))
    if not seqs:
        json.dump({"tasks": [], "stats": {"barcode": a.barcode, "n_sequences": 0, "n_clusters": 0,
                                          "n_groups": 0, "n_msa": 0, "ref_seqs": 0, "residual": 0}},
                  open(out / "sativa_manifest.json", "w"))
        return

    aln, bb_ids = build_backbone(seqs, tax, work, a.backbone_target, a.backbone_min_per_stratum,
                                 a.backbone_aligner, a.threads, a.backbone_source,
                                 a.backbone_cluster_id, a.mmseqs_env)
    bb_rows, width = frame_from_backbone(aln, a.min_col_occ)
    hmms, best = project(bb_rows, width, tax, seqs, work, a.hmmer_bin, a.bucket_mode,
                         a.bucket_rank, a.min_bucket, a.subset_size, a.threads,
                         a.nhmmer_jobs, a.nhmmer_threads)

    per_bucket = defaultdict(list)
    for sid, (_score, bucket, strand) in best.items():
        if sid in seqs:
            per_bucket[bucket].append((sid, seqs[sid].translate(COMP)[::-1] if strand == "-"
                                       else seqs[sid]))
    rows = {}
    bucket_of = {}
    for bucket, items in per_bucket.items():
        for sid, _s in items:
            bucket_of[sid] = bucket
    hmmalign_buckets(per_bucket, hmms, work, width, a.hmmer_bin, rows, a.hmmalign_jobs)
    log.info("projected %d reads into the %d-column frame", len(rows), width)

    rescue_stats = {}
    if a.recover:
        cross = read_cross_refs(a.cross_ref)
        if cross:
            log.info("cross-barcode rescue enabled: %d reference reads loaded", len(cross))
        rescue_stats = recover(seqs, rows, best, hmms, work, width, a.hmmer_bin, a.mmseqs_env,
                               a.recover_min_id, a.recover_min_qcov, a.threads,
                               cross=cross, bucket_of=bucket_of, align_jobs=a.hmmalign_jobs)
        if rescue_stats.get("rescued_cross_barcode"):
            # Same shape as barcode_reannotation.json's by_change: "<source>-><target>": n
            rescue_stats["by_source"] = {
                f"{src}->{a.barcode}": rescue_stats["rescued_cross_barcode"]
                for src in sorted({os.path.basename(os.path.dirname(p)) for p in a.cross_ref})
            }
    log.info("total projected: %d/%d (%.1f%%)", len(rows), len(seqs), 100 * len(rows) / len(seqs))

    # Side-car other barcodes rescue against: every projected read, ungapped, tagged with the
    # bucket that accepted it. Written for every barcode; only ITS is consumed today.
    refs_path = out / "projected_refs.fasta"
    out.mkdir(parents=True, exist_ok=True)
    with open(refs_path, "w") as fh:
        for sid, row in rows.items():
            fh.write(f">{sid}|{bucket_of.get(sid, 'GLOBAL')}\n{row.replace('-', '')}\n")
    log.info("wrote %s (%d reads)", refs_path, len(rows))

    unprojected = {sid for sid in seqs if sid not in rows}
    slice_groups(rows, tax, bb_ids, out, a.barcode, a.max_group, a.min_group, a.scaffold,
                 len(seqs), unprojected, rescue_stats, a.split_seed)


if __name__ == "__main__":
    main()
