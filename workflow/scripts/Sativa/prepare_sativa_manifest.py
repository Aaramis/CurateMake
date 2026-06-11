#!/usr/bin/env python3
"""
Scan the MSA directory and build the manifest of SATIVA tasks.
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
)
logger = logging.getLogger("prepare_sativa_manifest")


def infer_taxon_from_filename(msa_file: Path) -> str:
    """
    Derive a taxon name from an MSA filename that ends with
    `_msa.fasta` or `_msa.trimmed.fasta`.
    """
    name = msa_file.name
    if name.endswith("_msa.trimmed.fasta"):
        return name[: -len("_msa.trimmed.fasta")]
    if name.endswith("_msa.fasta"):
        return name[: -len("_msa.fasta")]
    return msa_file.stem.replace("_msa", "")


def is_fasta_empty(fasta_path: Path) -> bool:
    """Check if a FASTA file is empty (no sequences)."""
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        return True
    with fasta_path.open("r") as f:
        for line in f:
            if line.startswith(">"):
                return False
    return True


def list_msa_files(msa_dir: Path, barcodes: List[str], tax_levels: List[str]) -> List[Tuple[str, str, str, Path]]:
    """
    Enumerate every available MSA file for the requested barcodes/levels.
    Returns: List of tuples (barcode, tax_level, taxon, msa_path)
    """
    results: List[Tuple[str, str, str, Path]] = []
    skipped_empty = 0

    if not msa_dir.exists():
        logger.warning("MSA directory does not exist: %s", msa_dir)
        return results

    for barcode in barcodes:
        barcode_dir = msa_dir / barcode
        if not barcode_dir.exists():
            continue

        for level in tax_levels:
            align_dir = barcode_dir / f"alignments_{level}"
            if not align_dir.exists():
                continue

            trimmed = sorted(align_dir.glob("*_msa.trimmed.fasta"))
            candidates = trimmed if trimmed else sorted(align_dir.glob("*_msa.fasta"))

            for msa_file in candidates:
                if is_fasta_empty(msa_file):
                    skipped_empty += 1
                    continue
                taxon = infer_taxon_from_filename(msa_file)
                results.append((barcode, level, taxon, msa_file))

    if skipped_empty > 0:
        logger.info("Skipped %d empty MSA file(s)", skipped_empty)

    return results


def build_manifest(msa_dir: Path, barcodes: List[str], tax_levels: List[str]) -> Dict:
    """Build the JSON manifest structure."""
    tasks = []
    for barcode, level, taxon, msa_path in list_msa_files(msa_dir, barcodes, tax_levels):
        tasks.append(
            {
                "barcode": barcode,
                "tax_level": level,
                "taxon": taxon,
                "msa_path": str(msa_path),
            }
        )

    manifest = {
        "generated_at": datetime.now().isoformat(),
        "msa_dir": str(msa_dir),
        "total_tasks": len(tasks),
        "barcodes": sorted(set(barcodes)),
        "tax_levels": sorted(set(tax_levels)),
        "tasks": tasks,
    }
    return manifest


def main():
    try:
        if "snakemake" in globals():
            msa_dir = Path(snakemake.params.input_msa_dir)
            barcodes = snakemake.params.barcodes
            tax_levels = snakemake.params.tax_levels
            output_path = Path(snakemake.output.manifest)

            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
            logger.addHandler(file_handler)
        else:
            raise RuntimeError("This script must run within Snakemake.")

        logger.info("Scanning MSA directory: %s", msa_dir)
        manifest = build_manifest(msa_dir, barcodes, tax_levels)
        logger.info("Detected %d SATIVA tasks", manifest["total_tasks"])

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)
        logger.info("Manifest written to %s", output_path)

    except Exception as exc:
        logger.exception("Failed to build SATIVA manifest: %s", exc)
        raise


if __name__ == "__main__":
    main()
