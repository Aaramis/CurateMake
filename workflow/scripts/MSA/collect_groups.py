#!/usr/bin/env python3
"""
Aggregate the per-barcode MSA manifests (ehmm_msa.smk) into the artifacts the
rest of the workflow expects — so the sequence-first strategy is a drop-in for the bottom-up one:

  - SATIVA manifest      (BU_SATIVA_MANIFEST): the merged {tasks:[...]} sativa.smk consumes.
  - bottomup_summary.json: the summary generate_msa_report.py + generate_sequence_flow_report.py
                           read (per_barcode ref_seqs/residual/reinserted, n_msa_total, ...). The
                           bottom-up-only tiers (residual re-insertion, cross-family climb, contam
                           re-homing) do not exist here, so their fields are 0 / empty by design.
  - residual_report.json  + contamination_report.json: compatible STUBS (this path has no
                           residual-classification or coarse-contamination tier — every sequence is
                           placed at clustering, and a contaminant clusters with its true kin).

Each per-barcode manifest carries {"tasks": [...], "stats": {barcode, n_sequences, n_clusters,
n_groups, n_msa, ref_seqs, residual}}.
"""
import argparse
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

_KINGDOM_RE = re.compile(r"kingdom=([^|\s]+)")

# An MSA narrower than this cannot support a phylogenetic placement: SATIVA's verdicts
# on it are noise. Reported, never acted upon -- the point is to make it visible.
NARROW_COLUMNS = 100
LOW_OCCUPANCY_PCT = 20.0


def measure_msa(msa_path):
    """Exact width, occupancy and kingdom mix of one aligned FASTA.

    One pass per file. Width and occupancy are what decide whether a SATIVA placement on
    this MSA means anything; the kingdom mix says how broad a reference frame the group
    provides.
    """
    width = None
    gaps = residues = n_seqs = 0
    kingdoms = Counter()
    pending = []
    try:
        with open(msa_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if pending:
                        seq = "".join(pending)
                        if width is None:
                            width = len(seq)
                        gaps += seq.count("-")
                        residues += len(seq)
                        pending = []
                    n_seqs += 1
                    m = _KINGDOM_RE.search(line)
                    kingdoms[m.group(1) if m else "NA"] += 1
                else:
                    pending.append(line.strip())
            if pending:
                seq = "".join(pending)
                if width is None:
                    width = len(seq)
                gaps += seq.count("-")
                residues += len(seq)
    except OSError:
        return None

    return {
        "n_seqs": n_seqs,
        "width": width or 0,
        "occupancy_pct": round(100.0 * (1 - gaps / residues), 1) if residues else 0.0,
        "n_kingdoms": len(kingdoms),
        "kingdoms": dict(kingdoms.most_common()),
    }


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

    # --- real per-MSA metrics (one pass over each alignment) ----------------
    per_msa = []
    for t in tasks:
        entry = {
            "msa_id": t.get("taxon"),
            "barcode": t.get("barcode"),
            "taxon_label": t.get("taxon_label"),
            "n_group": t.get("n", 0),                       # sequences in the merged group
            "n_msa_seqs": t.get("n_msa_seqs", t.get("n", 0)),  # sequences actually aligned
        }
        entry["represented_pct"] = (round(100.0 * entry["n_msa_seqs"] / entry["n_group"], 1)
                                    if entry["n_group"] else 0.0)
        m = measure_msa(t.get("msa_path", ""))
        if m:
            entry.update(m)
        per_msa.append(entry)

    n_group_total = sum(e["n_group"] for e in per_msa)
    n_aligned_total = sum(e["n_msa_seqs"] for e in per_msa)
    widths = [e["width"] for e in per_msa if e.get("width")]
    occs = [e["occupancy_pct"] for e in per_msa if "occupancy_pct" in e]
    alignment_quality = {
        "n_measured": len(widths),
        "width_min": min(widths) if widths else 0,
        "width_median": sorted(widths)[len(widths) // 2] if widths else 0,
        "width_max": max(widths) if widths else 0,
        "occupancy_min_pct": min(occs) if occs else 0.0,
        "occupancy_median_pct": sorted(occs)[len(occs) // 2] if occs else 0.0,
        "narrow_threshold": NARROW_COLUMNS,
        "low_occupancy_threshold_pct": LOW_OCCUPANCY_PCT,
        "n_narrow": sum(1 for w in widths if w < NARROW_COLUMNS),
        "n_low_occupancy": sum(1 for o in occs if o < LOW_OCCUPANCY_PCT),
        # sequences whose SATIVA verdict rests on an alignment too narrow to support it
        "seqs_in_narrow": sum(e["n_msa_seqs"] for e in per_msa
                              if e.get("width", 0) and e["width"] < NARROW_COLUMNS),
    }
    representation = {
        "n_group_total": n_group_total,
        "n_aligned_total": n_aligned_total,
        # the rest get their verdict propagated via group_members.tsv, not by placement
        "n_propagated": n_group_total - n_aligned_total,
        "propagated_pct": (round(100.0 * (n_group_total - n_aligned_total) / n_group_total, 1)
                           if n_group_total else 0.0),
        "n_msa_at_cap": sum(1 for e in per_msa if e["n_msa_seqs"] >= 10000),
    }

    summary = {
        "strategy": "ehmm",
        "n_families": n_clusters,          # sequence clusters (the decomposition unit here)
        "total_sequences": total_seqs,
        "n_msa_total": n_msa_total,
        "n_residual_total": n_residual,
        "sativa_selected": n_msa_total,
        "per_barcode": per_barcode,
        "per_msa": per_msa,
        "alignment_quality": alignment_quality,
        "representation": representation,
    }

    # This strategy has no residual-reclassification and no coarse-contamination
    # tier. These files stay (they are declared rule inputs) but say so explicitly instead
    # of reporting zeros that read as "checked, nothing found".
    not_applicable = {
        "applicable": False,
        "reason": "not part of the backbone + ensemble-HMM strategy; "
                  "every sequence is placed at clustering time",
    }
    residual_stub = dict(not_applicable, n_residual=n_residual)
    contamination_stub = dict(not_applicable, n_msa=n_msa_total)
    sativa_manifest = {
        "generated_at": datetime.now().isoformat(),
        "source": "ehmm",
        "total_tasks": n_msa_total,
        "barcodes": barcodes,
        "tax_levels": sorted({t.get("tax_level", "groups") for t in tasks}),
        "tasks": tasks,
    }

    for path, obj in [(args.output_summary, summary), (args.sativa_manifest, sativa_manifest),
                      (args.residual_report, residual_stub), (args.contamination_report, contamination_stub)]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2)

    msg = (f"[collect_groups] {len(barcodes)} barcode(s), {n_msa_total} MSAs, "
           f"{total_seqs} seqs, residual={n_residual}")
    print(msg)
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        with open(args.log_file, "w", encoding="utf-8") as fh:
            fh.write(msg + "\n")


if __name__ == "__main__":
    main()
