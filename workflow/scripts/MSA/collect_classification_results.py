#!/usr/bin/env python3
"""
Aggregate hmmscan batch results into a single CSV and compute reassignment stats.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


logger = logging.getLogger("collect_classification_results")


def setup_logging(verbose: bool = False, log_file: Optional[Path] = None) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(log_file, mode="w")
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logging.getLogger().addHandler(handler)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate classification batches and compute statistics.")
    parser.add_argument("--batch-csvs", nargs="*", default=[], help="Hmmscan batch assignment CSVs")
    parser.add_argument("--unprocessed-csv", required=True, help="CSV with source metadata for sequences")
    parser.add_argument("--output-assignments", required=True, help="Combined assignments CSV")
    parser.add_argument("--output-json", required=True, help="Classification statistics JSON")
    parser.add_argument("--id-col", required=True)
    parser.add_argument("--barcode-col", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", help="Optional log file path")
    return parser.parse_args()


def generate_reassignment_statistics(assignments_df: pd.DataFrame, barcode_col: str) -> Dict[str, Dict[str, float]]:
    """Compute reassignment stats similar to HierarchicalAligner."""
    if assignments_df.empty:
        return {}

    stats: Dict[str, Dict[str, float]] = {}
    df = assignments_df.copy()
    df['target_barcode'] = df['assigned_msa'].str.extract(r'^(ITS\d?)_')[0]
    df['target_barcode'].fillna('Unknown', inplace=True)

    barcodes = ['ITS', 'ITS1', 'ITS2']
    for source in barcodes:
        source_total = (df[barcode_col] == source).sum()
        for target in barcodes:
            count = ((df[barcode_col] == source) & (df['target_barcode'] == target)).sum()
            percentage = (count / source_total * 100) if source_total > 0 else 0
            stats[f"{source}_to_{target}"] = {
                "count": int(count),
                "percentage": float(percentage)
            }
    return stats


def load_batches(paths: List[str]) -> List[pd.DataFrame]:
    dataframes: List[pd.DataFrame] = []
    for path in paths:
        if not Path(path).exists():
            continue
        df = pd.read_csv(path)
        if not df.empty:
            dataframes.append(df)
    return dataframes


def run(batch_csvs, unprocessed_csv, output_assignments, output_json, id_col, barcode_col, verbose=False, log_file=None):
    log_path = Path(log_file) if log_file else None
    setup_logging(verbose, log_path)

    batch_dfs = load_batches(batch_csvs)
    if batch_dfs:
        assignments = pd.concat(batch_dfs, ignore_index=True)
    else:
        assignments = pd.DataFrame(columns=[id_col, "assigned_msa", "score"])

    unprocessed_df = pd.read_csv(unprocessed_csv)[[id_col, barcode_col]].rename(
        columns={barcode_col: "source_barcode"}
    )
    merged = assignments.merge(unprocessed_df, on=id_col, how="left")

    Path(output_assignments).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_assignments, index=False)
    logger.info("Saved combined assignments to %s (%d rows)", output_assignments, len(merged))

    stats = generate_reassignment_statistics(merged, "source_barcode")
    summary = {
        "total_assignments": int(len(merged)),
        "reassignment_stats": stats
    }

    Path(output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    logger.info("Classification stats saved to %s", output_json)


def main():
    if "snakemake" in globals():
        input_batches = list(snakemake.input.get("batches", []))
        unprocessed_csv = snakemake.input.get("unprocessed")
        output_assignments = snakemake.output.get("assignments")
        output_json = snakemake.output.get("stats")
        id_col = snakemake.params.get("id_col")
        barcode_col = snakemake.params.get("barcode_col")
        log_file = snakemake.log[0] if snakemake.log else None
        run(input_batches, unprocessed_csv, output_assignments, output_json, id_col, barcode_col, False, log_file)
    else:
        args = parse_args()
        run(
            args.batch_csvs,
            args.unprocessed_csv,
            args.output_assignments,
            args.output_json,
            args.id_col,
            args.barcode_col,
            args.verbose,
            args.log_file
        )


if __name__ == "__main__":
    main()
