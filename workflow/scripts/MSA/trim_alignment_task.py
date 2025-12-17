#!/usr/bin/env python3
"""
Wrapper to run trimAl on a single alignment task.
"""

import argparse
import json
import logging
import shutil
from pathlib import Path

from external_tools import run_trimal_and_identify_removed


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"trim_alignment.{log_file.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Trim a single MSA alignment.")
    parser.add_argument("--input-msa", required=True, help="Path to raw alignment fasta")
    parser.add_argument("--trimmed-msa", required=True, help="Output path for trimmed alignment")
    parser.add_argument("--barcode", required=True)
    parser.add_argument("--clade", required=True)
    parser.add_argument("--min-overlap", type=float, required=True)
    parser.add_argument("--res-overlap", type=float, required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--removed-json", help="Optional JSON file to record removed sequence IDs")
    return parser.parse_args()


def main():
    args = parse_args()
    logger = setup_logger(Path(args.log_file))

    trimmed_path = Path(args.trimmed_msa)
    trimmed_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.input_msa, trimmed_path)

    logger.info("Trimming MSA %s for %s/%s", trimmed_path, args.barcode, args.clade)
    removed_ids = run_trimal_and_identify_removed(
        input_msa_path=str(trimmed_path),
        seq_overlap_threshold=int(args.min_overlap),
        seq_resoverlap_threshold=float(args.res_overlap)
    )

    if args.removed_json:
        removed_path = Path(args.removed_json)
        removed_path.parent.mkdir(parents=True, exist_ok=True)
        with open(removed_path, "w", encoding="utf-8") as fh:
            json.dump({"removed_ids": list(removed_ids)}, fh, indent=2)

    logger.info("Trimmed alignment saved to %s (%d seq removed)", trimmed_path, len(removed_ids))


if __name__ == "__main__":
    main()
