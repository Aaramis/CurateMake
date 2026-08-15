#!/usr/bin/env python3
"""
Aggregate the per-barcode clustered-MSA manifests (build_msas_clustered.py) into the artifacts the
rest of the workflow expects — so the sequence-first strategy is a drop-in for the bottom-up one:

  - SATIVA manifest      (BU_SATIVA_MANIFEST): the merged {tasks:[...]} sativa.smk consumes.
  - bottomup_summary.json: the summary generate_msa_report.py + generate_sequence_flow_report.py
                           read (per_barcode ref_seqs/residual/reinserted, n_msa_total, ...). The
                           bottom-up-only tiers (residual re-insertion, cross-family climb, contam
                           re-homing) do not exist here, so their fields are 0 / empty by design.
  - residual_report.json  + contamination_report.json: compatible STUBS (the clustered path has no
                           residual-classification or coarse-contamination tier — every sequence is
                           placed at clustering, and a contaminant clusters with its true kin).

Each per-barcode manifest carries {"tasks": [...], "stats": {barcode, n_sequences, n_clusters,
n_groups, n_msa, ref_seqs, residual}}.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcode-manifests", nargs="*", default=[])
    ap.add_argument("--output-summary", required=True)
    ap.add_argument("--sativa-manifest", required=True)
    ap.add_argument("--residual-report", required=True)
    ap.add_argument("--contamination-report", required=True)
    ap.add_argument("--log-file")
    args = ap.parse_args()

    tasks = []
    per_barcode = {}
    total_seqs = n_clusters = n_residual = 0
    for mpath in args.barcode_manifests:
        with open(mpath, encoding="utf-8") as fh:
            m = json.load(fh)
        tasks.extend(m.get("tasks", []))
        st = m.get("stats", {})
        bc = st.get("barcode", "ITS")
        total_seqs += st.get("n_sequences", 0)
        n_clusters += st.get("n_clusters", 0)
        n_residual += st.get("residual", 0)
        pb = per_barcode.setdefault(bc, {"ref_seqs": 0, "n_msa": 0, "residual": 0,
                                         "reinserted": 0, "unassigned": 0})
        pb["ref_seqs"] += st.get("ref_seqs", 0)
        pb["n_msa"] += st.get("n_msa", 0)
        pb["residual"] += st.get("residual", 0)
        pb["unassigned"] += st.get("residual", 0)   # no re-insertion tier -> residual == unassigned

    n_msa_total = len(tasks)
    barcodes = sorted(per_barcode)

    summary = {
        "strategy": "clustered",
        "n_families": n_clusters,          # sequence clusters (the decomposition unit here)
        "total_sequences": total_seqs,
        "n_no_family": 0,
        "n_msa_total": n_msa_total,
        "n_residual_total": n_residual,
        "n_no_family_recovered": 0,
        "n_reinserted_total": 0,
        "n_unassigned_total": n_residual,
        "gate_pass_total": 0,
        "climbed": {},
        "stop_levels": {"clustered": n_msa_total},
        "sativa_selected": n_msa_total,
        "per_barcode": per_barcode,
        "residual_files": [],
        "per_family": [],
    }

    residual_stub = {"n_residual": n_residual, "n_classified": 0, "n_integrated": 0,
                     "n_unassigned": n_residual, "n_profiles": 0, "n_no_family": 0, "per_msa": {}}
    contamination_stub = {"n_flagged": 0, "n_msa": n_msa_total, "n_rehomed_msa": 0,
                          "flagged": [], "emptied_msa": []}
    sativa_manifest = {
        "generated_at": datetime.now().isoformat(),
        "source": "clustered",
        "total_tasks": n_msa_total,
        "barcodes": barcodes,
        "tax_levels": sorted({t.get("tax_level", "clustered") for t in tasks}),
        "tasks": tasks,
    }

    for path, obj in [(args.output_summary, summary), (args.sativa_manifest, sativa_manifest),
                      (args.residual_report, residual_stub), (args.contamination_report, contamination_stub)]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2)

    msg = (f"[collect_clustered] {len(barcodes)} barcode(s), {n_msa_total} MSAs, "
           f"{total_seqs} seqs, residual={n_residual}")
    print(msg)
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        with open(args.log_file, "w", encoding="utf-8") as fh:
            fh.write(msg + "\n")


if __name__ == "__main__":
    main()
