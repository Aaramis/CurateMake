#!/usr/bin/env python3
"""
Append hmmscan classification assignments back into their corresponding MSA files.

This script mimics the legacy HierarchicalAligner integration step: for every
sequence that received an hmmscan assignment, append the sequence to the target
MSA using ``mafft --add`` so future runs can start from enriched alignments.
"""

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from external_tools import run_mafft_add
import utils


def setup_logger(log_file: Optional[Path] = None) -> logging.Logger:
    logger = logging.getLogger("msa_assignment_integration")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.handlers = [handler]

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(file_handler)

    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Append hmmscan assignments back into MSAs.")
    parser.add_argument("--assignments", required=True, help="CSV with hmmscan assignments.")
    parser.add_argument("--unprocessed", required=True, help="CSV of sequences sent to hmmscan.")
    parser.add_argument("--manifest", required=True, help="MSA task manifest JSON.")
    parser.add_argument("--msa-dir", required=True, help="Base directory for MSA outputs.")
    parser.add_argument("--id-col", default="seq_id")
    parser.add_argument("--seq-col", default="sequence")
    parser.add_argument("--tax-columns", nargs="+", required=True)
    parser.add_argument("--max-threads", type=int, default=4)
    parser.add_argument("--output-summary", required=True, help="JSON summary to write.")
    parser.add_argument("--log-file", help="Optional log file path.")
    return parser.parse_args()


def normalize_clade(label: str) -> str:
    """Remove barcode and suffix artefacts from hmmscan profile names."""
    if not label:
        return ""
    # Drop barcode prefix if still present.
    name = re.sub(r"^ITS\d?_?", "", label)
    # Remove legacy suffixes.
    for suffix in (".trimmed", ".msa", ".hmm"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.endswith("_msa"):
        name = name[:-4]
    return name


def read_manifest(manifest_path: Path, msa_dir: Path) -> Dict[Tuple[str, str], Dict[str, Path]]:
    """Build lookup from (barcode, clade) to alignment path and level."""
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    lookup: Dict[Tuple[str, str], Dict[str, Path]] = {}
    for barcode, payload in manifest.get("barcodes", {}).items():
        tasks: List[Dict] = []
        tasks.extend(payload.get("tasks", []))
        tasks.extend(payload.get("large_tasks", []))

        for task in tasks:
            clade = task["clade_name"]
            level = task["level"]
            base_dir = msa_dir / barcode / f"alignments_{level}"
            msa_path = base_dir / f"{clade}_msa.fasta"
            trimmed_path = base_dir / f"{clade}_msa.trimmed.fasta"
            lookup[(barcode, clade)] = {
                "level": level,
                "msa_path": msa_path,
                "trimmed_path": trimmed_path
            }
    return lookup


def existing_sequence_ids(msa_path: Path) -> set[str]:
    """Collect base seq_ids already present in an MSA to avoid duplicates."""
    if not msa_path.exists():
        return set()
    ids = set()
    for header in utils.get_fasta_ids(str(msa_path)):
        base_id = header.split("|", 1)[0]
        ids.add(base_id)
    return ids


def integrate_assignments(assignments_csv: Path,
                          unprocessed_csv: Path,
                          manifest_lookup: Dict[Tuple[str, str], Dict[str, Path]],
                          msa_dir: Path,
                          id_col: str,
                          seq_col: str,
                          tax_columns: List[str],
                          max_threads: int,
                          logger: logging.Logger) -> Dict:
    assignments = pd.read_csv(assignments_csv)
    if assignments.empty:
        logger.info("No hmmscan assignments detected; nothing to integrate.")
        return {"total_assignments": 0, "total_sequences_added": 0, "details": []}

    # Filter out rows without a target profile.
    assignments = assignments.dropna(subset=["assigned_msa"])
    if assignments.empty:
        logger.info("Assignments CSV contains no valid hmmscan hits; skipping integration.")
        return {"total_assignments": 0, "total_sequences_added": 0, "details": []}

    unprocessed = pd.read_csv(unprocessed_csv).drop_duplicates(subset=[id_col])
    unprocessed_index = unprocessed.set_index(id_col)

    assignments["target_barcode"] = assignments["assigned_msa"].str.extract(r"^(ITS\d?)_")[0]
    assignments["clade_name"] = assignments["assigned_msa"].apply(normalize_clade)
    assignments = assignments.dropna(subset=["target_barcode", "clade_name"])

    summary_entries = []
    total_added = 0

    for (barcode, clade_name), group in assignments.groupby(["target_barcode", "clade_name"]):
        lookup_key = (barcode, clade_name)
        info = manifest_lookup.get(lookup_key)
        if not info:
            logger.warning("No MSA path found for %s/%s; skipping %d sequences.",
                           barcode, clade_name, len(group))
            continue

        target_path = info["trimmed_path"] if info["trimmed_path"].exists() else info["msa_path"]
        if not target_path.exists():
            logger.warning("MSA file missing for %s/%s (%s); skipping.",
                           barcode, clade_name, target_path)
            continue

        seq_ids = group[id_col].unique().tolist()
        missing = [seq for seq in seq_ids if seq not in unprocessed_index.index]
        if missing:
            logger.warning("%s/%s - %d sequences missing from unprocessed pool: %s",
                           barcode, clade_name, len(missing), ", ".join(missing[:5]))

        subset = unprocessed_index.loc[[seq for seq in seq_ids if seq in unprocessed_index.index]].reset_index()
        if subset.empty:
            continue

        current_ids = existing_sequence_ids(target_path)
        subset = subset[~subset[id_col].isin(current_ids)]
        if subset.empty:
            logger.info("%s/%s - all %d sequences already present.", barcode, clade_name, len(seq_ids))
            continue

        logger.info("%s/%s - adding %d sequences to %s", barcode, clade_name, len(subset), target_path)
        success = run_mafft_add(
            new_sequences_df=subset,
            existing_msa_path=str(target_path),
            output_msa_path=str(target_path),
            seq_col=seq_col,
            tax_columns=tax_columns,
            max_thread_workers=max_threads,
            id_col=id_col
        )

        if not success:
            logger.error("%s/%s - failed to update MSA.", barcode, clade_name)
            continue

        total_added += len(subset)
        summary_entries.append({
            "barcode": barcode,
            "clade": clade_name,
            "level": info["level"],
            "msa_path": str(target_path),
            "added_sequences": int(len(subset))
        })

    result = {
        "total_assignments": int(len(assignments)),
        "total_sequences_added": int(total_added),
        "details": summary_entries
    }

    return result


def run_cli(args: argparse.Namespace) -> None:
    logger = setup_logger(Path(args.log_file) if args.log_file else None)
    summary = integrate_assignments(
        assignments_csv=Path(args.assignments),
        unprocessed_csv=Path(args.unprocessed),
        manifest_lookup=read_manifest(Path(args.manifest), Path(args.msa_dir)),
        msa_dir=Path(args.msa_dir),
        id_col=args.id_col,
        seq_col=args.seq_col,
        tax_columns=args.tax_columns,
        max_threads=args.max_threads,
        logger=logger
    )

    output_path = Path(args.output_summary)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    logger.info("Integration summary written to %s (added %d sequences)",
                output_path, summary["total_sequences_added"])


def run_from_snakemake(snakemake) -> None:
    params = snakemake.params
    logger = setup_logger(Path(snakemake.log[0]) if snakemake.log else None)
    summary = integrate_assignments(
        assignments_csv=Path(snakemake.input.assignments),
        unprocessed_csv=Path(snakemake.input.unprocessed),
        manifest_lookup=read_manifest(Path(snakemake.input.manifest), Path(params["msa_dir"])),
        msa_dir=Path(params["msa_dir"]),
        id_col=params["id_col"],
        seq_col=params["seq_col"],
        tax_columns=params["tax_columns"],
        max_threads=int(params["max_thread_workers"]),
        logger=logger
    )
    output_path = Path(snakemake.output.summary)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    logger.info("Integration summary written to %s (added %d sequences)",
                output_path, summary["total_sequences_added"])


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake(snakemake)  # type: ignore[name-defined]
    else:
        run_cli(parse_args())
