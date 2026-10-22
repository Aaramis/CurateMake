#!/usr/bin/env python3
"""
Aggregate SATIVA per-MSA summaries into the global analysis report.
"""

import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Any

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
)
logger = logging.getLogger("collect_sativa_results")


def cleanup_old_temp_files(sativa_output_dir: Path, retention_hours: int, cleanup_enabled: bool) -> None:
    """Remove stale per-task SATIVA scratch directories.

    Each task now writes its scratch to <its output dir>/sativa_tmp and removes it on
    success, so what is left here is the scratch of tasks that failed or timed out, kept
    deliberately for inspection. This sweep is the retention policy for those.
    """
    if not cleanup_enabled:
        logger.info("Temporary file cleanup disabled")
        return

    if not sativa_output_dir.exists():
        logger.info("SATIVA output directory not found: %s", sativa_output_dir)
        return

    logger.info("Cleaning leftover task scratch under %s (retention: %dh)",
                sativa_output_dir, retention_hours)
    now = time.time()
    cleaned = 0
    total_size = 0

    for folder in sativa_output_dir.rglob("sativa_tmp"):
        if not folder.is_dir():
            continue
        age_hours = (now - folder.stat().st_mtime) / 3600
        if age_hours <= retention_hours:
            continue

        folder_size = sum(f.stat().st_size for f in folder.rglob("*") if f.is_file())
        total_size += folder_size
        try:
            for attempt in range(3):
                try:
                    import shutil

                    shutil.rmtree(folder, ignore_errors=False)
                    break
                except Exception as exc:
                    if attempt == 2:
                        raise
                    logger.warning("Retrying cleanup of %s after error: %s", folder, exc)
                    time.sleep(1)
            cleaned += 1
        except Exception as exc:
            logger.warning("Failed to delete %s: %s", folder, exc)

    if cleaned:
        logger.info("Removed %d temp directories (%.1f MB freed)", cleaned, total_size / 1024 / 1024)
    else:
        logger.info("No temporary directories required cleaning")


def load_summary(path: Path) -> Dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            result = json.load(handle)
            if "status" not in result:
                result["status"] = "unknown"
            return result
    except FileNotFoundError:
        logger.warning("Summary missing: %s", path)
    except json.JSONDecodeError as exc:
        logger.warning("Invalid JSON in %s: %s", path, exc)
    except Exception as exc:
        logger.warning("Failed to read %s: %s", path, exc)
    return {
        "barcode": "unknown",
        "tax_level": "unknown",
        "taxon": path.parent.name,
        "n_seq": 0,
        "aln_len": 0,
        "status": "error",
        "reason": f"Missing or unreadable summary ({path})",
    }


def aggregate_results(summary_files: List[Path]) -> Dict[str, Any]:
    results_by_status: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    stats_by_barcode: Dict[str, Dict[str, int]] = defaultdict(lambda: {"total": 0, "success": 0, "skipped": 0, "error": 0, "timeout": 0})
    stats_by_tax_level: Dict[str, Dict[str, int]] = defaultdict(lambda: {"total": 0, "success": 0, "skipped": 0, "error": 0, "timeout": 0})

    for summary_path in summary_files:
        result = load_summary(summary_path)
        status = result.get("status", "unknown")
        barcode = result.get("barcode", "unknown")
        tax_level = result.get("tax_level", "unknown")

        results_by_status[status].append(result)
        stats_by_barcode[barcode]["total"] += 1
        stats_by_barcode[barcode][status] = stats_by_barcode[barcode].get(status, 0) + 1
        stats_by_tax_level[tax_level]["total"] += 1
        stats_by_tax_level[tax_level][status] = stats_by_tax_level[tax_level].get(status, 0) + 1

    report_data = {
        "total_msa_processed": len(summary_files),
        "summary": {
            "success": len(results_by_status["success"]),
            "skipped": len(results_by_status["skipped"]),
            "error": len(results_by_status["error"]),
            "timeout": len(results_by_status["timeout"]),
        },
        "by_barcode": dict(stats_by_barcode),
        "by_tax_level": dict(stats_by_tax_level),
        "results": dict(results_by_status),
    }
    return report_data


def main():
    try:
        if "snakemake" in globals():
            summary_files = [Path(p) for p in snakemake.input.summaries]
            output_path = Path(snakemake.output.report)

            cleanup_enabled = snakemake.params.cleanup_temp_files
            retention_hours = snakemake.params.temp_retention_hours
            sativa_output_dir = Path(snakemake.params.sativa_output_dir)

            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
            logger.addHandler(file_handler)
        else:
            raise RuntimeError("This script must run within Snakemake.")

        if not summary_files:
            logger.warning("No SATIVA summary files found. Producing empty report.")

        report_data = aggregate_results(summary_files)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(report_data, handle, indent=2)
        logger.info("SATIVA analysis report written to %s", output_path)

        # A per-MSA SATIVA failure is non-fatal (one bad group must not kill a 5M-sequence run), but it
        # MUST NOT be silent: a crashed group contributes ZERO corrections while the run still reports
        # success. Historically the biggest groups were exactly the ones that died (RAxML GTRGAMMA
        # numerical crash above ~10k taxa), so the failures were concentrated where the data matters.
        n_err = report_data["summary"]["error"] + report_data["summary"]["timeout"]
        if n_err:
            failed = [r.get("taxon", "?") for st in ("error", "timeout")
                      for r in report_data["results"].get(st, [])]
            logger.warning("=" * 78)
            logger.warning("%d of %d SATIVA task(s) FAILED -> they produced NO corrections.",
                           n_err, report_data["total_msa_processed"])
            logger.warning("Affected MSAs: %s", ", ".join(failed[:20]) + (" ..." if len(failed) > 20 else ""))
            logger.warning("Check logs/sativa/tasks/<barcode>/<level>/<taxon>.log for the cause.")
            logger.warning("=" * 78)

        cleanup_old_temp_files(sativa_output_dir, retention_hours, cleanup_enabled)

    except Exception as exc:
        logger.exception("Failed to aggregate SATIVA results: %s", exc)
        raise


if __name__ == "__main__":
    main()
