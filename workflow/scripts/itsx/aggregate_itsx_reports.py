#!/usr/bin/env python3
"""
Aggregate the per-barcode ITSx reports into a single consolidated report.

Designed to run inside Snakemake after barcode-parallel ITSx jobs complete.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Any

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_barcode_report(report_path: Path) -> Dict[str, Any]:
    """
    Load a single barcode JSON report.

    Args:
        report_path: Path to the barcode JSON report.

    Returns:
        Dictionary containing the parsed report data.
    """
    try:
        with open(report_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to read {report_path}: {e}")
        raise


def calculate_totals(barcodes_data: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Compute global totals from every barcode report.

    Args:
        barcodes_data: Mapping of barcode name to its report data.

    Returns:
        Dictionary with aggregated totals.
    """
    totals = {
        "total_processed": 0,
        "total_subregions_found": 0,
        "total_subregions_inserted": 0,
        "total_reconstructed_its": 0,
        "total_no_detection": 0
    }

    for barcode, stats in barcodes_data.items():
        totals["total_processed"] += stats.get("processed", 0)
        totals["total_no_detection"] += stats.get("no_detection", 0)

        # Subregions
        subregions = stats.get("subregions", {})
        for region, region_stats in subregions.items():
            totals["total_subregions_found"] += region_stats.get("found", 0)
            totals["total_subregions_inserted"] += region_stats.get("inserted", 0)

        # Reconstructed ITS
        reconstructed = stats.get("reconstructed_its", {})
        totals["total_reconstructed_its"] += reconstructed.get("inserted", 0)

    return totals


def aggregate_reports(input_files: list, output_file: Path) -> None:
    """
    Aggregate every barcode report into a final JSON document.

    Args:
        input_files: List of paths to each barcode report.
        output_file: Destination path for the consolidated report.
    """
    logger.info(f"Aggregating {len(input_files)} ITSx reports...")

    # Load every individual report
    barcodes_data = {}
    for report_path in input_files:
        report_path = Path(report_path)
        logger.info(f"Reading {report_path.name}...")

        data = load_barcode_report(report_path)

        # File naming format: itsx_{barcode}_result.json
        barcode = report_path.stem.replace("itsx_", "").replace("_result", "")

        barcodes_data[barcode] = data

    # Compute totals
    totals = calculate_totals(barcodes_data)

    # Build final JSON payload
    final_report = {
        "barcodes": barcodes_data,
        "totals": totals
    }

    # Persist report
    logger.info(f"Writing consolidated report to {output_file}...")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(final_report, f, indent=2)

    # Summary banner
    logger.info("=" * 70)
    logger.info("ITSx REPORT AGGREGATION SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Barcodes processed: {len(barcodes_data)}")
    logger.info(f"Sequences processed: {totals['total_processed']}")
    logger.info(f"Subregions found: {totals['total_subregions_found']}")
    logger.info(f"Subregions inserted: {totals['total_subregions_inserted']}")
    logger.info(f"ITS reconstructed: {totals['total_reconstructed_its']}")
    logger.info(f"Non-detections: {totals['total_no_detection']}")
    logger.info("=" * 70)


def main():
    """Snakemake entry point."""
    try:
        # Gather Snakemake inputs/outputs
        input_files = snakemake.input
        output_file = Path(snakemake.output.report)

        # Configure log file handler if provided
        if hasattr(snakemake, 'log') and snakemake.log:
            log_file = snakemake.log[0]
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(
                logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            )
            logger.addHandler(file_handler)

        # Aggregate reports
        aggregate_reports(input_files, output_file)

        logger.info("ITSx aggregation complete!")

    except Exception as e:
        logger.error(f"Error while aggregating ITSx reports: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
