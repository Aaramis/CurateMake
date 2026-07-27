#!/usr/bin/env python3
"""
Combine all barcode-level HMM databases into a single global HMM file
(plain text, consumed directly by hmmsearch — no hmmpress index needed).
"""

import argparse
import logging
from pathlib import Path

from external_tools import combine_hmm_databases


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"combine_global_hmm.{log_file.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.handlers = [handler]
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine barcode HMM databases into a single global HMM.")
    parser.add_argument("--input-hmms", nargs="+", required=True, help="List of per-barcode HMM database paths")
    parser.add_argument("--output-hmm", required=True, help="Combined HMM path")
    parser.add_argument("--log-file", required=True, help="Log file path")
    return parser.parse_args()


def run(input_hmms, output_hmm, log_file: Path):
    logger = setup_logger(log_file)
    output_path = Path(output_hmm)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Combining %d HMM databases into %s", len(input_hmms), output_hmm)
    success = combine_hmm_databases(input_hmms, output_hmm)
    if not success:
        logger.error("Failed to combine HMM databases.")
        raise SystemExit(1)

    logger.info("Global HMM ready: %s", output_hmm)


def main():
    if "snakemake" in globals():
        input_hmms = snakemake.input if isinstance(snakemake.input, list) else [snakemake.input]
        output_hmm = snakemake.output[0] if isinstance(snakemake.output, list) else snakemake.output
        log_file = Path(snakemake.log[0]) if snakemake.log else Path("logs/MSA/global_hmm.log")
    else:
        args = parse_args()
        input_hmms = args.input_hmms
        output_hmm = args.output_hmm
        log_file = Path(args.log_file)

    run(input_hmms, output_hmm, log_file)


if __name__ == "__main__":
    main()
