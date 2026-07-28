#!/usr/bin/env python3
"""
Aggregate the bottom-up construction into a global summary + the SATIVA manifest.

Two sources:
  - PER-FAMILY manifests: residual, per-family stats.
  - CLIMBED manifest (climb_lineage): the AUTHORITATIVE list of final MSAs
    (climbed cross-family merges + passthrough), which feeds SATIVA.

SATIVA selection: every MSA with at least `sativa_min_seqs` sequences goes to SATIVA
(including genus-only / single-species MSAs — SATIVA can still flag topological
intruders there).
"""

import argparse
import json
import logging
from collections import Counter
from datetime import datetime
from pathlib import Path

logger = logging.getLogger("collect_bottomup")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate bottom-up manifests + emit SATIVA manifest.")
    p.add_argument("--family-manifests", nargs="*", default=[])
    p.add_argument("--split-manifest", required=True)
    p.add_argument("--climbed-manifest", required=True,
                   help="Final-MSA manifest (climb_lineage): source of the SATIVA tasks.")
    p.add_argument("--output-summary", required=True)
    p.add_argument("--sativa-manifest", help="SATIVA-compatible manifest to emit (optional).")
    p.add_argument("--sativa-min-seqs", type=int, default=4,
                   help="Size floor: only MSAs with >= this many sequences go to SATIVA.")
    p.add_argument("--contamination-report", default=None,
                   help="Contamination report: emptied (re-homed) MSAs are excluded from SATIVA.")
    p.add_argument("--residual-report", default=None,
                   help="Residual report: orphans RE-INSERTED by hmmsearch (per_msa) -> "
                        "re-injected per barcode into the summary (flow report).")
    p.add_argument("--log-file")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        h = logging.FileHandler(args.log_file, mode="w")
        h.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(h)

    with open(args.split_manifest, encoding="utf-8") as fh:
        split = json.load(fh)
    with open(args.climbed_manifest, encoding="utf-8") as fh:
        climbed = json.load(fh)

    # MSAs emptied by the re-homing (defensive: climb already drops empty files)
    emptied = set()
    if args.contamination_report:
        try:
            with open(args.contamination_report, encoding="utf-8") as fh:
                emptied = set(json.load(fh).get("emptied_msa", []))
        except (FileNotFoundError, json.JSONDecodeError):
            pass

    # --- residual + per-family stats (family manifests) ---
    n_residual = 0
    residual_files = []
    per_family = []
    per_barcode = {}   # barcode -> {ref_seqs, residual, n_msa}
    for mpath in args.family_manifests:
        with open(mpath, encoding="utf-8") as fh:
            m = json.load(fh)
        n_residual += m.get("n_residual", 0)
        res = Path(mpath).parent / "residual.fasta"
        if res.exists():
            residual_files.append(str(res))
        per_family.append({"fam_id": Path(mpath).parent.name,
                           "n_msa": m.get("n_msa", 0), "n_residual": m.get("n_residual", 0),
                           "stop_levels": m.get("stop_levels", {})})
        pb = per_barcode.setdefault(m.get("barcode", "ITS"), {"ref_seqs": 0, "residual": 0, "reinserted": 0, "n_msa": 0})
        pb["residual"] += m.get("n_residual", 0)

    # --- final MSAs + SATIVA tasks (climbed manifest) ---
    final_msas = [e for e in climbed.get("msas", []) if e.get("msa_id") not in emptied]
    stop_levels = Counter(e.get("level") for e in final_msas)
    gate_pass = sum(1 for e in final_msas if e.get("gate_pass"))
    sativa_tasks = []
    sativa_barcodes = set()
    for e in final_msas:
        barcode = e.get("barcode", "ITS")
        pb = per_barcode.setdefault(barcode, {"ref_seqs": 0, "residual": 0, "reinserted": 0, "n_msa": 0})
        pb["ref_seqs"] += e.get("n", 0)
        pb["n_msa"] += 1
        sativa_barcodes.add(barcode)
        if not e.get("msa_path"):
            continue
        if e.get("n", 0) < args.sativa_min_seqs:
            continue
        # taxon = msa_id (filesystem-safe: used as {taxon} in the SATIVA paths);
        # taxon_label = real label for reporting.
        sativa_tasks.append({"barcode": barcode, "tax_level": e["level"],
                             "taxon": e["msa_id"], "taxon_label": e.get("taxon", e["msa_id"]),
                             "msa_path": e["msa_path"]})

    # --- re-insertion of residual orphans (residual report) ---
    # per_msa: msa_id -> number of re-inserted orphans; msa_id = "{barcode}__{family}__{NNNN}".
    n_reinserted_total = 0
    n_unassigned_total = 0
    n_orphans_total = n_residual                # family singletons (fallback if no report)
    n_no_family_recovered = 0
    if args.residual_report:
        try:
            with open(args.residual_report, encoding="utf-8") as fh:
                residual = json.load(fh)
            n_reinserted_total = residual.get("n_integrated", 0)
            n_unassigned_total = residual.get("n_unassigned", 0)
            n_orphans_total = residual.get("n_residual", n_residual)   # family singletons + no-family
            n_no_family_recovered = residual.get("n_no_family", 0)
            for msa_id, cnt in residual.get("per_msa", {}).items():
                bc = msa_id.split("__", 1)[0]
                pb = per_barcode.setdefault(bc, {"ref_seqs": 0, "residual": 0, "reinserted": 0, "n_msa": 0})
                pb["reinserted"] += cnt
        except (FileNotFoundError, json.JSONDecodeError):
            pass
    # unassigned per barcode: source residual - re-inserted (clamped >= 0)
    for pb in per_barcode.values():
        pb["unassigned"] = max(0, pb.get("residual", 0) - pb.get("reinserted", 0))

    summary = {
        "n_families": split.get("n_families", len(per_family)),
        "total_sequences": split.get("total_sequences", 0),
        "n_no_family": split.get("n_no_family", 0),
        "n_msa_total": len(final_msas),
        "n_residual_total": n_orphans_total,
        "n_no_family_recovered": n_no_family_recovered,
        "n_reinserted_total": n_reinserted_total,
        "n_unassigned_total": n_unassigned_total,
        "gate_pass_total": gate_pass,
        "climbed": climbed.get("climbed", {}),
        "stop_levels": dict(stop_levels),
        "sativa_selected": len(sativa_tasks),
        "per_barcode": per_barcode,
        "residual_files": residual_files,
        "per_family": per_family,
    }
    Path(args.output_summary).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_summary, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    if args.sativa_manifest:
        sativa_manifest = {
            "generated_at": datetime.now().isoformat(),
            "source": "bottomup",
            "total_tasks": len(sativa_tasks),
            "min_seqs": args.sativa_min_seqs,
            "barcodes": sorted(sativa_barcodes),
            "tax_levels": sorted({t["tax_level"] for t in sativa_tasks}),
            "tasks": sativa_tasks,
        }
        Path(args.sativa_manifest).parent.mkdir(parents=True, exist_ok=True)
        with open(args.sativa_manifest, "w", encoding="utf-8") as fh:
            json.dump(sativa_manifest, fh, indent=2)
        logger.info("SATIVA manifest: %d tasks (>= %d seqs).",
                    len(sativa_tasks), args.sativa_min_seqs)

    logger.info("Bottom-up: %d families, %d final MSAs, %d residual, climbed=%s, levels=%s",
                summary["n_families"], len(final_msas), n_residual,
                climbed.get("climbed", {}), dict(stop_levels))


if __name__ == "__main__":
    main()
