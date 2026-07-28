"""Aggregate per-clade SATIVA integration reports into a single report.

The integration is scattered per clade (rule sativa_integrate_task); this barrier
rule sums the per-clade statistics into sativa_integration_report.json, keeping the
same aggregate shape the downstream steps expect.
"""
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def _merge_counts(dst, src):
    """Sum int/float counters (recursively for the col_api_stats sub-dict)."""
    for key, value in src.items():
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)):
            dst[key] = dst.get(key, 0) + value
        elif isinstance(value, dict) and key == "col_api_stats":
            sub = dst.setdefault(key, {})
            _merge_counts(sub, value)


def main(snakemake):
    reports = snakemake.input.clade_reports
    if isinstance(reports, str):
        reports = [reports]

    global_stats = {}
    per_clade = []
    for path in reports:
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception as exc:  # noqa: BLE001 - a missing/corrupt clade report must not abort the aggregate
            logger.warning("Skipping unreadable clade report %s: %s", path, exc)
            continue
        _merge_counts(global_stats, data)
        meta = data.get("metadata", {})
        per_clade.append({
            "taxon": meta.get("taxon"),
            "barcode": meta.get("barcode"),
            "tax_level": meta.get("tax_level"),
            "runtime_seconds": meta.get("runtime_seconds"),
            "total_modifications": data.get("total_modifications", 0),
        })

    global_stats["metadata"] = {
        "end_time": datetime.now().isoformat(),
        "total_clades": len(reports),
        "clades_with_modifications": sum(1 for c in per_clade if c["total_modifications"]),
    }
    global_stats["per_clade"] = per_clade

    out_report = Path(snakemake.output.report)
    out_report.parent.mkdir(parents=True, exist_ok=True)
    with open(out_report, "w") as fh:
        json.dump(global_stats, fh, indent=2, default=str)

    logger.info(
        "Aggregated %d clade reports -> %s (seq_tax_created=%s, warnings_created=%s)",
        len(reports), out_report,
        global_stats.get("seq_tax_created", 0), global_stats.get("warnings_created", 0),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake)
    else:
        logger.error("This script must be run through Snakemake")
        sys.exit(1)
