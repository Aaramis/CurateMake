#!/usr/bin/env python3
"""
Prepare per-clade integration tasks so the ``mafft --add`` step can be
scattered across clades by Snakemake (it used to be a single serial loop over
all clades in one job).

Reads the hmmsearch assignments + the MSA task manifest, groups the sequences
to integrate by (barcode, clade), resolves each clade's target MSA path, and
writes one task JSON per clade plus a manifest listing them. No MAFFT is run
here — that happens in the scattered worker (integrate_clade_task.py).
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

# Reuse the exact same grouping/path-resolution helpers as the legacy monolith.
from add_assignments_to_msa import normalize_clade, read_manifest

logger = logging.getLogger("prepare_integration_tasks")


def setup_logging(log_file=None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(log_file, mode="w")
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logging.getLogger().addHandler(handler)


def ensure_clean_dir(path: Path) -> None:
    if path.exists():
        for f in path.glob("*"):
            if f.is_file():
                f.unlink()
    else:
        path.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare per-clade integration tasks.")
    p.add_argument("--assignments", required=True, help="CSV with hmmsearch assignments.")
    p.add_argument("--manifest", required=True, help="MSA task manifest JSON.")
    p.add_argument("--msa-dir", required=True, help="Base directory for MSA outputs.")
    p.add_argument("--id-col", default="seq_id")
    p.add_argument("--output-manifest", required=True, help="Integration manifest JSON to write.")
    p.add_argument("--tasks-dir", required=True, help="Directory for per-clade task JSONs.")
    p.add_argument("--log-file")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(args.log_file)

    tasks_dir = Path(args.tasks_dir)
    ensure_clean_dir(tasks_dir)

    lookup = read_manifest(Path(args.manifest), Path(args.msa_dir))

    assignments = pd.read_csv(args.assignments)
    if not assignments.empty:
        assignments = assignments.dropna(subset=["assigned_msa"])

    jobs = []
    total_assignments = 0

    if not assignments.empty:
        total_assignments = int(len(assignments))
        # Same extraction as the legacy integrate_assignments().
        assignments["target_barcode"] = assignments["assigned_msa"].str.extract(r"^(ITS\d?)_")[0]
        assignments["clade_name"] = assignments["assigned_msa"].apply(normalize_clade)
        assignments = assignments.dropna(subset=["target_barcode", "clade_name"])

        index = 0
        for (barcode, clade_name), group in assignments.groupby(["target_barcode", "clade_name"]):
            info = lookup.get((barcode, clade_name))
            if not info:
                logger.warning("No MSA path for %s/%s; skipping %d sequences.",
                               barcode, clade_name, len(group))
                continue
            target_path = info["trimmed_path"] if info["trimmed_path"].exists() else info["msa_path"]
            if not Path(target_path).exists():
                logger.warning("MSA file missing for %s/%s (%s); skipping.",
                               barcode, clade_name, target_path)
                continue

            seq_ids = sorted(group[args.id_col].unique().tolist(), key=str)
            # Per-clade HMM profile (built from the trimmed MSA) — used by the
            # hmmalign integration method (mirrors hmm_inputs() in msa.smk).
            profile_path = str(Path(args.msa_dir) / barcode / "hmm" / info["level"] / f"{clade_name}.hmm")
            index += 1
            job_id = f"job_{index:05d}"
            job = {
                "job_id": job_id,
                "barcode": barcode,
                "clade": clade_name,
                "level": info["level"],
                "target_path": str(target_path),
                "profile_path": profile_path,
                "seq_ids": seq_ids,
            }
            job_path = tasks_dir / f"{job_id}.json"
            with open(job_path, "w", encoding="utf-8") as fh:
                json.dump(job, fh, indent=2)
            jobs.append({
                "job_id": job_id,
                "json_path": str(job_path),
                "barcode": barcode,
                "clade": clade_name,
                "n_seqs": len(seq_ids),
            })

    manifest = {
        "total_assignments": total_assignments,
        "n_jobs": len(jobs),
        "jobs": jobs,
    }
    out = Path(args.output_manifest)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    logger.info("Integration manifest ready: %d clade jobs (%d assignments).",
                len(jobs), total_assignments)


if __name__ == "__main__":
    main()
