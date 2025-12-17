#!/usr/bin/env python3
"""
Combine all individual HMM profiles for a barcode into a single hmmpress-ready database.
"""

import argparse
import logging
from pathlib import Path

from external_tools import create_hmm_database


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"combine_hmm.{log_file.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.handlers = [handler]
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine barcode-specific HMM profiles.")
    parser.add_argument("--hmm-dir", required=True, help="Directory containing individual .hmm files")
    parser.add_argument("--output-hmm", required=True, help="Path to combined HMM database")
    parser.add_argument("--log-file", required=True, help="Path to log file")
    return parser.parse_args()


def run(hmm_dir: str, output_hmm: str, log_path: Path):
    logger = setup_logger(log_path)
    Path(output_hmm).parent.mkdir(parents=True, exist_ok=True)

    logger.info("Combining HMM profiles from %s into %s", hmm_dir, output_hmm)
    success = create_hmm_database(hmm_dir, output_hmm)
    if not success:
        logger.error("Failed to combine HMM profiles in %s", hmm_dir)
        raise SystemExit(1)

    logger.info("Combined HMM database ready: %s", output_hmm)


def main():
    if "snakemake" in globals():
        params = snakemake.params if hasattr(snakemake, "params") else {}
        input_dir = params.get("hmm_dir") or params.get("hmm_dir_path")
        output_hmm = snakemake.output[0] if len(snakemake.output) == 1 else snakemake.output.get("combined")
        log_file = snakemake.log[0] if snakemake.log else Path("logs/MSA/combine_hmm.log")
    else:
        args = parse_args()
        input_dir = args.hmm_dir
        output_hmm = args.output_hmm
        log_file = Path(args.log_file)

    run(input_dir, output_hmm, Path(log_file))


if __name__ == "__main__":
    main()
