#!/usr/bin/env python3
"""
Split full.csv into one CSV per (barcode, family) for the bottom-up scatter.

full.csv is large (~2.5 GB) and contains quoted fields (internal commas)
=> read in CHUNKS with pandas (robust to quoting, bounded memory). Each
(barcode, family) is written to <families_dir>/<fam_id>.csv, and a manifest lists
the families for Snakemake (checkpoint).

Sequences WITHOUT a family cannot enter a family tree -> written to
<families_dir>/_no_family.csv (handled as residual / downstream hmmsearch classif).
"""

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Dict

import pandas as pd

logger = logging.getLogger("split_families")

KEEP_COLS = ["seq_id", "sequence", "kingdom", "phylum", "class", "order",
             "family", "genus", "species", "barcode"]


def sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "NA"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Split full.csv into per-(barcode,family) CSVs.")
    p.add_argument("--input-csv", required=True)
    p.add_argument("--families-dir", required=True)
    p.add_argument("--output-manifest", required=True)
    p.add_argument("--barcodes", nargs="+", default=["ITS", "ITS1", "ITS2"])
    p.add_argument("--chunksize", type=int, default=200000)
    p.add_argument("--log-file")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        h = logging.FileHandler(args.log_file, mode="w")
        h.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(h)

    fam_dir = Path(args.families_dir)
    if fam_dir.exists():
        for f in fam_dir.glob("*.csv"):
            f.unlink()
    fam_dir.mkdir(parents=True, exist_ok=True)

    barcodes = set(args.barcodes)
    seen: Dict[str, Dict] = {}      # fam_id -> {barcode, family, n}
    written_header = set()
    no_family = 0

    def emit(fam_id: str, df: pd.DataFrame, header: bool):
        df.to_csv(fam_dir / f"{fam_id}.csv", mode="a", header=header, index=False)

    total = 0
    for chunk in pd.read_csv(args.input_csv, chunksize=args.chunksize, dtype=str, keep_default_na=False):
        chunk = chunk[[c for c in KEEP_COLS if c in chunk.columns]]
        # missing barcode -> 'ITS' (consistent with the extraction)
        if "barcode" in chunk.columns:
            chunk["barcode"] = chunk["barcode"].replace("", "ITS")
        chunk = chunk[chunk["barcode"].isin(barcodes)]
        total += len(chunk)

        nofam = chunk[chunk["family"].str.strip() == ""]
        if len(nofam):
            no_family += len(nofam)
            emit("_no_family", nofam, "_no_family" not in written_header)
            written_header.add("_no_family")

        withfam = chunk[chunk["family"].str.strip() != ""]
        for (bc, fam), grp in withfam.groupby(["barcode", "family"], sort=False):
            fam_id = f"{sanitize(bc)}__{sanitize(fam)}"
            emit(fam_id, grp, fam_id not in written_header)
            written_header.add(fam_id)
            e = seen.setdefault(fam_id, {"barcode": bc, "family": fam, "n": 0})
            e["n"] += len(grp)

    families = [{"fam_id": k, "barcode": v["barcode"], "family": v["family"],
                 "n_seq": v["n"], "csv_path": str(fam_dir / f"{k}.csv")}
                for k, v in sorted(seen.items())]
    manifest = {"n_families": len(families), "total_sequences": total,
                "n_no_family": no_family, "families": families}
    Path(args.output_manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_manifest, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("Split: %d families, %d sequences (%d without a family).",
                len(families), total, no_family)


if __name__ == "__main__":
    main()
