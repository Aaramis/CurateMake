#!/usr/bin/env python3
"""
Prepare the classification batches by collecting all sequences that still need
to be processed after the alignment/trim phases.

This script:
  - Reads the task manifest (per-barcode tasks + unprocessed IDs)
  - Scans the MSA output directory for trim removal JSON files
  - Extracts the corresponding sequences from the original CSV
  - Chunks the sequence IDs into hmmscan-friendly batches and writes one JSON per batch
  - Emits a manifest JSON summarizing batches + per-barcode counts
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


logger = logging.getLogger("prepare_classification_batches")


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
    parser = argparse.ArgumentParser(description="Prepare classification batches for MSA.")
    parser.add_argument("--input-csv", required=True, help="Full CSV with sequence data")
    parser.add_argument("--manifest-json", required=True, help="Task manifest JSON")
    parser.add_argument("--msa-dir", required=True, help="Base directory containing MSA outputs")
    parser.add_argument("--output-csv", required=True, help="CSV of all sequences to classify")
    parser.add_argument("--output-manifest", required=True, help="JSON manifest describing batches")
    parser.add_argument("--batches-dir", required=True, help="Directory to store per-batch JSON files")
    parser.add_argument("--batch-size", type=int, required=True, help="Number of sequences per batch")
    parser.add_argument("--id-col", required=True, help="Sequence identifier column")
    parser.add_argument("--barcode-col", required=True, help="Barcode column name")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.add_argument("--log-file", help="Optional log file path")
    return parser.parse_args()


def load_manifest(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def normalize_sequence_id(raw_id) -> str:
    """Normalize IDs from FASTA headers or JSON to match CSV seq_id values."""
    if raw_id is None:
        return ""
    raw_str = str(raw_id)
    if "|" in raw_str:
        return raw_str.split("|", 1)[0]
    return raw_str


def find_removed_sequences(msa_dir: Path, barcode: str) -> List[str]:
    """Collect IDs of sequences removed by trimAl for a specific barcode."""
    removed_ids: List[str] = []
    pattern = msa_dir / barcode
    if not pattern.exists():
        return removed_ids

    for removed_file in pattern.rglob("*_removed.json"):
        try:
            with open(removed_file, "r", encoding="utf-8") as fh:
                data = json.load(fh)
                for seq_id in data.get("removed_ids", []):
                    normalized = normalize_sequence_id(seq_id)
                    if normalized:
                        removed_ids.append(normalized)
        except (FileNotFoundError, json.JSONDecodeError) as exc:
            logger.warning("Failed to read %s: %s", removed_file, exc)
    return removed_ids


def chunk_list(items: List[str], chunk_size: int) -> List[List[str]]:
    if chunk_size <= 0:
        return [items]
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]


def ensure_clean_dir(path: Path) -> None:
    if path.exists():
        for file in path.glob("*"):
            if file.is_file():
                file.unlink()
    else:
        path.mkdir(parents=True, exist_ok=True)


def main():
    args = parse_args()
    log_path = Path(args.log_file) if args.log_file else None
    setup_logging(args.verbose, log_path)

    manifest = load_manifest(Path(args.manifest_json))
    df = pd.read_csv(args.input_csv)
    msa_dir = Path(args.msa_dir)

    barcode_info = manifest.get("barcodes", {})
    sequences_to_classify: Dict[str, List[str]] = {}

    for barcode, payload in barcode_info.items():
        unprocessed_ids = set(payload.get("unprocessed_sequence_ids", []))
        removed_ids = find_removed_sequences(msa_dir, barcode)
        unprocessed_ids.update(removed_ids)

        if unprocessed_ids:
            sequences_to_classify[barcode] = sorted(unprocessed_ids)
        logger.info(
            "%s - %d unprocessed + %d removed -> %d total sequences for classification",
            barcode,
            len(payload.get("unprocessed_sequence_ids", [])),
            len(removed_ids),
            len(unprocessed_ids)
        )

    all_ids = sorted({seq_id for ids in sequences_to_classify.values() for seq_id in ids})

    if not all_ids:
        logger.info("No sequences require classification; emitting empty files.")
        Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[args.id_col, args.barcode_col, "source_barcode"]).to_csv(args.output_csv, index=False)
        ensure_clean_dir(Path(args.batches_dir))
        batch_path = Path(args.batches_dir) / "batch_00001.json"
        with open(batch_path, "w", encoding="utf-8") as fh:
            json.dump({"batch_id": "batch_00001", "sequence_ids": []}, fh, indent=2)
        manifest_data = {
            "batch_size": args.batch_size,
            "total_sequences": 0,
            "barcode_summary": {},
            "batches": [{
                "batch_id": "batch_00001",
                "json_path": str(batch_path),
                "size": 0
            }]
        }
        Path(args.output_manifest).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output_manifest, "w", encoding="utf-8") as fh:
            json.dump(manifest_data, fh, indent=2)
        return

    sequences_df = df[df[args.id_col].isin(all_ids)].copy()
    sequences_df["source_barcode"] = sequences_df[args.barcode_col]
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    sequences_df.to_csv(args.output_csv, index=False)
    logger.info("Saved %d sequences to %s", len(sequences_df), args.output_csv)

    ensure_clean_dir(Path(args.batches_dir))

    batches = []
    seq_ids = sequences_df[args.id_col].tolist()
    for index, chunk in enumerate(chunk_list(seq_ids, args.batch_size), start=1):
        batch_id = f"batch_{index:05d}"
        batch_path = Path(args.batches_dir) / f"{batch_id}.json"
        with open(batch_path, "w", encoding="utf-8") as fh:
            json.dump({"batch_id": batch_id, "sequence_ids": chunk}, fh, indent=2)
        batches.append({
            "batch_id": batch_id,
            "json_path": str(batch_path),
            "size": len(chunk)
        })

    barcode_summary = {}
    for barcode, ids in sequences_to_classify.items():
        barcode_summary[barcode] = {
            "total_sequences": len(ids)
        }

    manifest_data = {
        "batch_size": args.batch_size,
        "total_sequences": len(seq_ids),
        "barcode_summary": barcode_summary,
        "batches": batches
    }

    Path(args.output_manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_manifest, "w", encoding="utf-8") as fh:
        json.dump(manifest_data, fh, indent=2)

    logger.info(
        "Classification manifest ready: %d batches (%d sequences)",
        len(batches),
        len(seq_ids)
    )


if __name__ == "__main__":
    main()
