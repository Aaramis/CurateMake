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

    # Cross-barcode reconciliation. Every clade was integrated independently, so ITS / ITS1 / ITS2 of
    # the SAME read can now hold contradictory SATIVA verdicts. This is the first point where the whole
    # database is visible, so it is where the most confident verdict is propagated to its siblings.
    cfg = getattr(snakemake, "config", {}) or {}
    recon_cfg = cfg.get("sativa_integration", {}).get("reconcile_barcodes", {})
    if recon_cfg.get("enabled", True):
        try:
            sys.path.insert(0, str(Path(__file__).parent))                        # reconcile_barcodes
            sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))       # same as integrate_sativa
            sys.path.insert(0, str(Path(__file__).parent.parent))
            from reconcile_barcodes import reconcile
            from db_manager import get_database_manager

            db_manager = get_database_manager(cfg)
            global_stats["barcode_reconciliation"] = reconcile(
                db_manager, apply=recon_cfg.get("apply", True),
                rank_by=recon_cfg.get("rank_by", "votes"))
        except Exception as exc:  # noqa: BLE001 - reconciliation must not lose the aggregated report
            logger.error("Barcode reconciliation failed (integration report still written): %s", exc)
            global_stats["barcode_reconciliation"] = {"error": str(exc)}
    else:
        logger.info("Barcode reconciliation disabled by config.")

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
