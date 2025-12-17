#!/usr/bin/env python3
"""
Standalone alignment runner for a single (barcode, clade) task.
Reads the global manifest, extracts the sequences for the requested task,
and invokes the existing mafft wrapper (external_tools.run_mafft).
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from external_tools import run_mafft


def setup_logger(log_path: Path) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"run_alignment_task.{log_path.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_path, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MAFFT alignment for a single manifest task.")
    parser.add_argument("--manifest", required=True, help="Path to manifest JSON from task_manifest.py")
    parser.add_argument("--barcode", required=True, help="Barcode identifier (e.g., ITS1)")
    parser.add_argument("--clade", required=True, help="Clade/task identifier (matches manifest entry)")
    parser.add_argument("--input-csv", required=True, help="Full CSV containing sequences")
    parser.add_argument("--output-folder", required=True, help="Base output directory for alignments")
    parser.add_argument("--seq-col", required=True, help="Sequence column in the CSV")
    parser.add_argument("--id-col", required=True, help="ID column in the CSV")
    parser.add_argument("--tax-columns", required=True, nargs='+', help="Taxonomy columns list")
    parser.add_argument("--log-file", required=True, help="Log file path")
    return parser.parse_args()


def load_task(manifest_path: Path, barcode: str, clade: str) -> dict:
    with open(manifest_path, "r", encoding="utf-8") as fh:
        manifest = json.load(fh)
    barcode_manifest = manifest["barcodes"].get(barcode)
    if not barcode_manifest:
        raise ValueError(f"Barcode {barcode} not found in manifest {manifest_path}")
    for task in barcode_manifest["tasks"] + barcode_manifest.get("large_tasks", []):
        if task["clade_name"] == clade:
            return task
    raise ValueError(f"Clade {clade} not found in manifest for barcode {barcode}")


def main():
    args = parse_args()
    logger = setup_logger(Path(args.log_file))

    logger.info("Loading manifest %s", args.manifest)
    task = load_task(Path(args.manifest), args.barcode, args.clade)

    logger.info("Loading CSV %s", args.input_csv)
    df = pd.read_csv(args.input_csv)

    seq_ids = task["sequence_ids"]
    subset = df[df[args.id_col].isin(seq_ids)]
    if subset.empty:
        logger.error("No sequences found for task %s/%s", args.barcode, args.clade)
        raise SystemExit(1)

    logger.info("Running MAFFT for %s/%s (%d sequences)", args.barcode, args.clade, len(subset))
    output_path = run_mafft(
        sequences_df=subset,
        barcode=args.barcode,
        clade_name=args.clade,
        level=task["level"],
        output_folder=args.output_folder,
        seq_col=args.seq_col,
        tax_columns=args.tax_columns,
        id_col=args.id_col
    )

    if not output_path:
        logger.error("MAFFT failed for %s/%s", args.barcode, args.clade)
        raise SystemExit(1)

    logger.info("Alignment written to %s", output_path)


if __name__ == "__main__":
    main()
