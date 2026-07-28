#!/usr/bin/env python3
"""
Classification + integration of the RESIDUAL orphans from the bottom-up construction,
plus the sequences that had NO family (they cannot enter a family tree, so they are
recovered here by homology instead of being dropped).

Sources of orphans:
  - per-family residual.fasta (singletons + unplaced sequences);
  - the _no_family.csv produced by split_families.

Steps:
  1. concatenate the HMM profiles of all reference MSAs -> library
  2. hmmsearch the orphans against the library -> best MSA per sequence (by homology;
     no barcode constraint — placement follows real homology, not the annotated barcode,
     which may itself be wrong; a true ITS1 will not score against an ITS2 profile anyway)
  3. add each orphan to its assigned MSA (hmmalign --mapali, fallback mafft --add)
Emits a JSON report. The assigned MSAs are mutated in place.
"""

import argparse
import csv
import json
import logging
import os
import tempfile
from pathlib import Path

import pandas as pd

from external_tools import run_hmmsearch, run_hmmalign_add, run_mafft_add, existing_msa_ids

logger = logging.getLogger("residual_classify_integrate")
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def read_fasta(path):
    """Returns [(raw_seq_id, tax_dict, sequence)] by parsing the taxonomized headers
    (>seq_id|kingdom=..|..). The raw id = before the first '|'."""
    recs, header, s = [], None, []

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
        recs.append((sid, tax, "".join(s)))

    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            header, s = line[1:].split()[0], []
        else:
            s.append(line)
    flush()
    return recs


def read_no_family_csv(path):
    """Returns [(seq_id, tax_dict, sequence)] from the _no_family.csv."""
    recs = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("seq_id") or "").strip()
            seq = (row.get("sequence") or "").replace("-", "").replace(".", "").upper()
            if not sid or not seq:
                continue
            tax = {r: (row.get(r) or "").strip() for r in RANKS}
            recs.append((sid, tax, seq))
    return recs


def parse_args():
    p = argparse.ArgumentParser(description="Classify+integrate residual orphans + no-family seqs (bottom-up).")
    p.add_argument("--family-manifests", nargs="*", default=[])
    p.add_argument("--no-family-csv", default=None,
                   help="split_families _no_family.csv: seqs without a family, recovered by homology.")
    p.add_argument("--output-report", required=True)
    p.add_argument("--max-evalue", type=float, default=1e-3)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--log-file")
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        h = logging.FileHandler(args.log_file, mode="w")
        h.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(h)

    # 1. map msa_id -> entry + collect the orphans (per-family residual.fasta + no-family CSV).
    msa_map = {}
    profile_paths = []
    orphan_recs = []            # (seq_id, tax_dict, seq)
    for mpath in args.family_manifests:
        m = json.load(open(mpath, encoding="utf-8"))
        for e in m.get("msas", []):
            msa_map[e["msa_id"]] = e
            if e.get("profile_path"):
                profile_paths.append(e["profile_path"])
        res = Path(mpath).parent / "residual.fasta"
        if res.exists():
            orphan_recs.extend(read_fasta(str(res)))

    n_no_family = 0
    if args.no_family_csv and os.path.exists(args.no_family_csv):
        nf = read_no_family_csv(args.no_family_csv)
        orphan_recs.extend(nf)
        n_no_family = len(nf)

    report = {"n_residual": len(orphan_recs), "n_no_family": n_no_family,
              "n_profiles": len(profile_paths), "n_classified": 0, "n_integrated": 0,
              "n_unassigned": 0, "per_msa": {}}

    if not orphan_recs or not profile_paths:
        logger.info("Nothing to do (orphans=%d, profiles=%d).", len(orphan_recs), len(profile_paths))
        Path(args.output_report).parent.mkdir(parents=True, exist_ok=True)
        json.dump(report, open(args.output_report, "w"), indent=2)
        return

    with tempfile.TemporaryDirectory(prefix="resid_") as tmp:
        # 2. profile library (flat concat, no hmmpress)
        lib = f"{tmp}/library.hmm"
        with open(lib, "w") as out:
            for p in sorted(profile_paths):
                with open(p) as fh:
                    out.write(fh.read())

        # df with the taxonomy (columns kingdom..species) to propagate it to the integration
        df = pd.DataFrame(
            [{"seq_id": sid, "sequence": seq, **tax} for sid, tax, seq in orphan_recs]
        ).drop_duplicates("seq_id")
        # hmmsearch: bare headers (tax_columns=[]) -> raw seq_id for the reclassification.
        # Placement is by homology only (no barcode filter): a mislabelled barcode should
        # not block a real homology hit, and a true ITS1 won't score against an ITS2 profile.
        assignments = run_hmmsearch(
            sequences_df=df, hmm_db_path=lib, seq_col="sequence",
            max_thread_workers=args.threads, tax_columns=[], id_col="seq_id",
            max_evalue=args.max_evalue,
        )
        report["n_classified"] = int(len(assignments))
        report["n_unassigned"] = len(df) - int(len(assignments))

        # 3. integration per assigned MSA (WITH taxonomy -> SATIVA can validate the orphans)
        for msa_id, grp in assignments.groupby("assigned_msa"):
            info = msa_map.get(msa_id)
            if not info:
                logger.warning("Profile %s with no known MSA; %d sequences ignored.", msa_id, len(grp))
                continue
            mpath = info.get("msa_path")
            # MSA emptied by the contamination re-homing (contaminant): its profile still lingers
            # in the lib, so an orphan can be assigned to it -> do NOT integrate there.
            if not mpath or not os.path.exists(mpath) or os.path.getsize(mpath) == 0:
                logger.warning("MSA %s empty/absent (re-homed); %d sequences not integrated.",
                               msa_id, len(grp))
                report["n_unassigned"] += len(grp)
                continue
            # Idempotence: do not re-add sequences already present (partial re-run:
            # build_family is cached on its manifest, not on the in-place-mutated FASTA).
            already = existing_msa_ids(mpath)
            subset = df[df["seq_id"].isin(grp["seq_id"]) & ~df["seq_id"].isin(already)]
            if subset.empty:
                continue
            ok = False
            if info.get("profile_path"):
                ok = run_hmmalign_add(subset, mpath, mpath,
                                      info["profile_path"], "sequence", RANKS, "seq_id")
            if not ok:
                ok = run_mafft_add(subset, mpath, mpath,
                                   "sequence", RANKS, args.threads, "seq_id")
            if ok:
                report["n_integrated"] += len(subset)
                report["per_msa"][msa_id] = int(len(subset))

    Path(args.output_report).parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(args.output_report, "w"), indent=2)
    logger.info("Residual+no-family: %d orphans (%d no-family), %d classified, %d integrated, %d unassigned.",
                report["n_residual"], report["n_no_family"], report["n_classified"],
                report["n_integrated"], report["n_unassigned"])


if __name__ == "__main__":
    main()
