#!/usr/bin/env python3
"""
Scatter worker: add ONE clade's assigned sequences to its MSA via ``mafft --add``.

Invocation is identical to the legacy monolith (same run_mafft_add call, same
flags) — the only change is that each clade is now an independent Snakemake job,
so clades run in parallel instead of one serial loop. Writes a per-clade result
JSON that the gather step (collect_integration_results.py) aggregates.
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from external_tools import run_mafft_add, run_hmmalign_add
from add_assignments_to_msa import existing_sequence_ids


def setup_logger(log_file: str) -> logging.Logger:
    path = Path(log_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"integrate_clade.{path.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(path, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.handlers = [handler]
    return logger


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Add one clade's assignments to its MSA.")
    p.add_argument("--task-json", required=True, help="Per-clade task JSON.")
    p.add_argument("--unprocessed", required=True, help="CSV of sequences sent to classification.")
    p.add_argument("--output-json", required=True, help="Per-clade result JSON to write.")
    p.add_argument("--id-col", default="seq_id")
    p.add_argument("--seq-col", default="sequence")
    p.add_argument("--tax-columns", nargs="+", required=True)
    p.add_argument("--max-threads", type=int, default=2)
    p.add_argument("--method", choices=["mafft", "hmmalign"], default="mafft",
                   help="Integration method: mafft --add (validated) or hmmalign --mapali (faster).")
    p.add_argument("--log-file", required=True)
    return p.parse_args()


def write_result(path: str, obj: dict) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)


def main() -> None:
    args = parse_args()
    logger = setup_logger(args.log_file)

    with open(args.task_json, encoding="utf-8") as fh:
        task = json.load(fh)

    barcode, clade, level = task["barcode"], task["clade"], task["level"]
    target_path = Path(task["target_path"])
    seq_ids = task["seq_ids"]

    result = {
        "barcode": barcode,
        "clade": clade,
        "level": level,
        "msa_path": str(target_path),
        "added_sequences": 0,
    }

    if not target_path.exists():
        logger.warning("Target MSA missing: %s", target_path)
        write_result(args.output_json, result)
        return

    unprocessed = (pd.read_csv(args.unprocessed)
                     .drop_duplicates(subset=[args.id_col])
                     .set_index(args.id_col))

    present = [s for s in seq_ids if s in unprocessed.index]
    missing = [s for s in seq_ids if s not in unprocessed.index]
    if missing:
        logger.warning("%s/%s - %d sequences missing from unprocessed pool: %s",
                       barcode, clade, len(missing), ", ".join(map(str, missing[:5])))

    subset = unprocessed.loc[present].reset_index()

    # Dedup against sequences already present (idempotent re-runs), same as legacy.
    current_ids = existing_sequence_ids(target_path)
    subset = subset[~subset[args.id_col].isin(current_ids)]
    if subset.empty:
        logger.info("%s/%s - nothing to add (all present or missing).", barcode, clade)
        write_result(args.output_json, result)
        return

    profile_path = task.get("profile_path")
    method = args.method
    logger.info("%s/%s - adding %d sequences to %s via %s",
                barcode, clade, len(subset), target_path, method)

    success = False
    if method == "hmmalign":
        if profile_path and Path(profile_path).exists() and Path(profile_path).stat().st_size > 0:
            success = run_hmmalign_add(
                new_sequences_df=subset,
                existing_msa_path=str(target_path),
                output_msa_path=str(target_path),
                profile_path=profile_path,
                seq_col=args.seq_col,
                tax_columns=args.tax_columns,
                id_col=args.id_col,
            )
        if not success:
            # Fall back to MAFFT so one bad/missing profile never blocks the step.
            logger.warning("%s/%s - hmmalign failed/unavailable, falling back to mafft --add.",
                           barcode, clade)
            method = "mafft (fallback)"
    if not success:
        success = run_mafft_add(
            new_sequences_df=subset,
            existing_msa_path=str(target_path),
            output_msa_path=str(target_path),
            seq_col=args.seq_col,
            tax_columns=args.tax_columns,
            max_thread_workers=args.max_threads,
            id_col=args.id_col,
        )

    if not success:
        # Best-effort, like the legacy loop: log and record failure, do NOT
        # crash the whole step for one clade. The gather counts these.
        logger.error("%s/%s - integration failed (both methods).", barcode, clade)
        result["failed"] = True
        write_result(args.output_json, result)
        return

    result["method"] = method
    result["added_sequences"] = int(len(subset))
    write_result(args.output_json, result)


if __name__ == "__main__":
    main()
