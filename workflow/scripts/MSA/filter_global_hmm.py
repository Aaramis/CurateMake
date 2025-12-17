#!/usr/bin/env python3
"""
Filter low-quality HMM profiles using hmmstat metrics.

Profiles failing the configured thresholds are removed from the combined HMM
before hmmscan classification to avoid overly permissive matches.
"""

import argparse
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("filter_global_hmm")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def run_hmmstat(hmm_path: Path) -> str:
    command = ["hmmstat", str(hmm_path)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"hmmstat failed: {result.stderr.strip()}")
    return result.stdout


def parse_hmmstat_output(text: str) -> List[Dict]:
    rows: List[Dict] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        rows.append({
            "idx": int(parts[0]),
            "name": parts[1],
            "accession": parts[2],
            "nseq": float(parts[3]),
            "eff_nseq": float(parts[4]),
            "M": float(parts[5]),
            "relent": float(parts[6]),
            "info": float(parts[7]),
            "p_relE": float(parts[8]),
            "compKL": float(parts[9]),
        })
    return rows


def filter_profiles(rows: List[Dict], min_relent: float, min_prele: float) -> Dict[str, bool]:
    status: Dict[str, bool] = {}
    for row in rows:
        keep = (row["relent"] >= min_relent) and (row["p_relE"] >= min_prele)
        status[row["name"]] = keep
    return status


def write_filtered_hmm(input_hmm: Path, output_hmm: Path, keep: Dict[str, bool], logger: logging.Logger) -> Dict[str, int]:
    output_hmm.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    removed = 0

    with input_hmm.open("r", encoding="utf-8") as infile, output_hmm.open("w", encoding="utf-8") as outfile:
        current_block: List[str] = []
        current_name: str = ""

        for line in infile:
            current_block.append(line)
            if line.startswith("NAME"):
                current_name = line.split()[1]
            if line.startswith("//"):
                if keep.get(current_name, True):
                    outfile.writelines(current_block)
                    kept += 1
                else:
                    removed += 1
                current_block = []
                current_name = ""

    logger.info("Kept %d HMM profiles; removed %d.", kept, removed)
    return {"kept": kept, "removed": removed}


def save_report(report_path: Path, hmmstat_raw: str, rows: List[Dict], keep_map: Dict[str, bool], summary: Dict[str, int]) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "summary": summary,
        "profiles": [
            {
                **row,
                "kept": keep_map.get(row["name"], True)
            } for row in rows
        ],
        "hmmstat_raw": hmmstat_raw
    }
    report_path.write_text(json.dumps(payload, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter HMM profiles using hmmstat metrics.")
    parser.add_argument("--input-hmm", required=True)
    parser.add_argument("--output-hmm", required=True)
    parser.add_argument("--report-json", required=True)
    parser.add_argument("--min-relent", type=float, default=0.2)
    parser.add_argument("--min-prele", type=float, default=0.2)
    parser.add_argument("--log-file", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    logger = setup_logger(Path(args.log_file))

    input_hmm = Path(args.input_hmm)
    output_hmm = Path(args.output_hmm)
    report_path = Path(args.report_json)

    logger.info("Running hmmstat on %s", input_hmm)
    hmmstat_raw = run_hmmstat(input_hmm)
    rows = parse_hmmstat_output(hmmstat_raw)
    keep_map = filter_profiles(rows, args.min_relent, args.min_prele)

    logger.info("Filtering HMM profiles (min_relent=%s, min_p_relE=%s)", args.min_relent, args.min_prele)
    summary = write_filtered_hmm(input_hmm, output_hmm, keep_map, logger)

    save_report(report_path, hmmstat_raw, rows, keep_map, summary)
    logger.info("Filtering complete. Report written to %s", report_path)


if __name__ == "__main__":
    main()
