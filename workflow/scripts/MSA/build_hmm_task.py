#!/usr/bin/env python3
"""
Build an HMM from a trimmed alignment using hmmbuild.
"""

import argparse
import logging
from pathlib import Path

from external_tools import run_hmmbuild


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"build_hmm.{log_file.stem}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build HMM for a trimmed alignment.")
    parser.add_argument("--trimmed-msa", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--barcode", required=True)
    parser.add_argument("--clade", required=True)
    parser.add_argument("--log-file", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    logger = setup_logger(Path(args.log_file))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hmm_path = out_dir / f"{args.clade}.hmm"

    success = run_hmmbuild(
        msa_path=args.trimmed_msa,
        hmm_path=str(hmm_path),
        barcode=args.barcode
    )

    if not success:
        logger.error("hmmbuild failed for %s/%s", args.barcode, args.clade)
        raise SystemExit(1)

    logger.info("HMM saved to %s", hmm_path)


if __name__ == "__main__":
    main()
