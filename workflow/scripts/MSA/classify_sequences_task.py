#!/usr/bin/env python3
"""
Classify a batch of sequences against a combined HMM database using hmmscan.
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from external_tools import run_hmmsearch


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"classify_sequences.{log_file.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify sequences batch with hmmscan.")
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--seq-ids-json", required=True, help="JSON with `sequence_ids` list")
    parser.add_argument("--combined-hmm", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--seq-col", required=True)
    parser.add_argument("--id-col", required=True)
    parser.add_argument("--tax-columns", required=True, nargs='+')
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--max-evalue", type=float, default=1e-3,
                        help="Maximum acceptable full-sequence E-value for hmmscan hits.")
    return parser.parse_args()


def main():
    args = parse_args()
    logger = setup_logger(Path(args.log_file))

    df = pd.read_csv(args.input_csv)
    with open(args.seq_ids_json, "r", encoding="utf-8") as fh:
        seq_ids = json.load(fh)["sequence_ids"]

    subset = df[df[args.id_col].isin(seq_ids)]
    if subset.empty:
        logger.warning("No sequences to classify for this batch.")
        Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[args.id_col, "assigned_msa", "score"]).to_csv(args.output_csv, index=False)
        return

    logger.info("Classifying %d sequences against %s", len(subset), args.combined_hmm)
    assignments = run_hmmsearch(
        sequences_df=subset,
        hmm_db_path=args.combined_hmm,
        seq_col=args.seq_col,
        max_thread_workers=args.threads,
        tax_columns=args.tax_columns,
        id_col=args.id_col,
        max_evalue=args.max_evalue
    )

    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(args.output_csv, index=False)
    logger.info("Assignments saved to %s", args.output_csv)


if __name__ == "__main__":
    main()
