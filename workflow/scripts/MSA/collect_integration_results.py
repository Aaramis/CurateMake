#!/usr/bin/env python3
"""
Gather step: merge the per-clade integration results into the single assignment
summary JSON (same path/schema the rest of the pipeline expects), replacing the
summary previously written by the monolithic add_assignments_to_msa.py.
"""

import argparse
import json
import logging
from pathlib import Path

logger = logging.getLogger("collect_integration_results")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate per-clade integration results.")
    p.add_argument("--result-jsons", nargs="*", default=[], help="Per-clade result JSONs.")
    p.add_argument("--manifest", required=True, help="Integration manifest JSON.")
    p.add_argument("--output-summary", required=True, help="Final assignment summary JSON.")
    p.add_argument("--log-file")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(args.log_file, mode="w")
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logging.getLogger().addHandler(handler)

    with open(args.manifest, encoding="utf-8") as fh:
        manifest = json.load(fh)
    total_assignments = int(manifest.get("total_assignments", 0))

    details = []
    total_added = 0
    failed = 0
    for rj in args.result_jsons:
        with open(rj, encoding="utf-8") as fh:
            r = json.load(fh)
        if r.get("failed"):
            failed += 1
        added = int(r.get("added_sequences", 0))
        total_added += added
        if added > 0:
            details.append({
                "barcode": r["barcode"],
                "clade": r["clade"],
                "level": r["level"],
                "msa_path": r["msa_path"],
                "added_sequences": added,
            })

    summary = {
        "total_assignments": total_assignments,
        "total_sequences_added": int(total_added),
        "failed_clades": failed,
        "details": details,
    }
    out = Path(args.output_summary)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    logger.info("Integration summary written to %s (added %d sequences, %d clade jobs, %d failed).",
                out, total_added, len(args.result_jsons), failed)


if __name__ == "__main__":
    main()
