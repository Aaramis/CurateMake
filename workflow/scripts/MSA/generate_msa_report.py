#!/usr/bin/env python3
"""
HTML report for the sequence-first clustered MSA construction.

Consumes clustered_summary.json (written by collect_clustered.py, which measures every
alignment on disk) and reports:
  1. the construction funnel: sequences -> clusters -> MSAs, per barcode,
  2. alignment quality: width x occupancy of every MSA, with an alert list,
  3. what SATIVA actually sees: aligned vs propagated verdicts, MSAs at the size cap,
  4. taxonomic breadth of each group (the premise of the sequence-first strategy),
  5. the largest groups.
"""

import argparse
import base64
import glob
import json
import logging
import os
import sys
from datetime import datetime
from io import BytesIO
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from report_style import REPORT_CSS, CHART_COLORS, page_header, page_footer

logger = logging.getLogger("generate_msa_report")


def load_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError, TypeError):
        return default if default is not None else {}


# Alert thresholds, mirrored from collect_clustered so the report can still draw them
# when reading a summary produced by an older version.
NARROW_COLUMNS = 100
LOW_OCCUPANCY_PCT = 20.0


def collect_msa_metrics(summary, out_dir=None):
    """Per-MSA metrics as a DataFrame.

    The clustered strategy measures every alignment in collect_clustered.py and stores
    the result in the summary; the per-family manifest layout is the fallback for
    bottom-up runs.
    """
    per_msa = summary.get("per_msa") or []
    if per_msa:
        return pd.DataFrame(per_msa)
    entries = []
    if out_dir:
        for mpath in glob.glob(os.path.join(out_dir, "*", "manifest.json")):
            entries.extend(load_json(mpath).get("msas", []))
    return pd.DataFrame(entries)


def _bar_plot(counts, title, xlabel):
    """Base64 PNG bar chart (matplotlib); returns '' if matplotlib is unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""
    if not counts:
        return ""
    fig, ax = plt.subplots(figsize=(8, 4))
    labels = list(counts.keys())
    ax.bar(labels, [counts[k] for k in labels], color=CHART_COLORS["primary"])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("# MSAs")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("utf-8")


def _stats_table(df, cols=None):
    if df.empty:
        return "<p><em>No data.</em></p>"
    d = df[cols] if cols else df
    return d.to_html(classes="stats-table", index=False, border=0)


def _scatter_quality(df):
    """Width x occupancy, one point per MSA, coloured by barcode. Base64 PNG."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""
    if df.empty or "width" not in df.columns or "occupancy_pct" not in df.columns:
        return ""
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = {"ITS": CHART_COLORS["primary"], "ITS1": CHART_COLORS["secondary"],
              "ITS2": CHART_COLORS["info"]}
    for bc, sub in df.groupby(df.get("barcode", "?")):
        ax.scatter(sub["width"].clip(lower=1), sub["occupancy_pct"], s=18, alpha=0.65,
                   label=str(bc), color=colors.get(str(bc), CHART_COLORS["muted"]))
    ax.set_xscale("log")
    ax.set_xlabel("alignment width (columns, log scale)")
    ax.set_ylabel("occupancy (% non-gap)")
    ax.set_title("Alignment quality — each point is one MSA")
    ax.axvline(NARROW_COLUMNS, color=CHART_COLORS["danger"], ls="--", lw=1)
    ax.axhline(LOW_OCCUPANCY_PCT, color=CHART_COLORS["danger"], ls="--", lw=1)
    ax.legend(title="barcode", fontsize=8)
    fig.tight_layout()
    buf = BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("utf-8")


def build_html(summary, msas, residual, contamination):
    n_msa = summary.get("n_msa_total", len(msas))
    per_barcode = summary.get("per_barcode", {})
    aq = summary.get("alignment_quality", {})
    rep = summary.get("representation", {})
    narrow_cols = aq.get("narrow_threshold", NARROW_COLUMNS)
    low_occ = aq.get("low_occupancy_threshold_pct", LOW_OCCUPANCY_PCT)

    # --- A. construction funnel -------------------------------------------
    funnel_df = pd.DataFrame([
        {"Barcode": bc, "Sequences": v.get("ref_seqs", 0), "MSAs": v.get("n_msa", 0),
         "Seqs / MSA (mean)": round(v.get("ref_seqs", 0) / max(v.get("n_msa", 0), 1))}
        for bc, v in sorted(per_barcode.items())
    ])

    # --- B. alignment quality ---------------------------------------------
    q_img = _scatter_quality(msas)
    q_plot_html = (f'<div class="plot"><img src="data:image/png;base64,{q_img}" alt="quality"></div>'
                   if q_img else "")
    flagged = pd.DataFrame()
    if not msas.empty and "width" in msas.columns:
        bad = msas[(msas["width"] < narrow_cols) | (msas["occupancy_pct"] < low_occ)]
        bad = bad.sort_values("width")
        flagged = bad.head(30)[[c for c in ["msa_id", "barcode", "n_msa_seqs", "width",
                                            "occupancy_pct", "n_kingdoms"] if c in bad.columns]]
        n_flagged_total = len(bad)
    else:
        n_flagged_total = 0

    # --- C. representation --------------------------------------------------
    rep_df = pd.DataFrame()
    if not msas.empty and "represented_pct" in msas.columns:
        worst = msas.sort_values("represented_pct").head(15)
        rep_df = worst[[c for c in ["msa_id", "barcode", "n_group", "n_msa_seqs",
                                    "represented_pct"] if c in worst.columns]]

    # --- D. taxonomic composition ------------------------------------------
    king_counts = {}
    if not msas.empty and "n_kingdoms" in msas.columns:
        king_counts = {f"{int(k)} kingdom(s)": int(v)
                       for k, v in sorted(msas["n_kingdoms"].value_counts().items())}
    king_img = _bar_plot(king_counts, "Kingdoms represented per MSA", "# kingdoms in the MSA")
    king_plot_html = (f'<div class="plot"><img src="data:image/png;base64,{king_img}" alt="kingdoms"></div>'
                      if king_img else "")

    # --- E. group sizes -----------------------------------------------------
    top_df = pd.DataFrame()
    if not msas.empty and "n_group" in msas.columns:
        top = msas.sort_values("n_group", ascending=False).head(15)
        top_df = top[[c for c in ["msa_id", "barcode", "taxon_label", "n_group", "n_msa_seqs",
                                  "width", "occupancy_pct"] if c in top.columns]]

    warn = ""
    if aq.get("n_narrow"):
        warn = (f'<div class="warn"><strong>{aq.get("n_narrow", 0)} of {n_msa} alignments are '
                f'narrower than {narrow_cols} columns</strong> (min {aq.get("width_min", 0)}), '
                f'covering {aq.get("seqs_in_narrow", 0):,} sequences. A placement on that few '
                f'columns carries no phylogenetic signal, so SATIVA verdicts for those sequences '
                f'should not be trusted.</div>')

    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>MSA Analysis Report</title>
        <style>{REPORT_CSS}</style>
    </head>
    <body>
        {page_header("MSA Analysis Report", "Sequence-first clustered construction &mdash; generated on " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))}
        <div class="container">

            <div class="section">
                <h2 class="section-title">1. Construction funnel</h2>
                <p>Sequences are clustered by <strong>sequence</strong> (mmseqs2) over the whole
                   barcode, small clusters are absorbed into their nearest neighbour, then clusters
                   are merged by shared taxonomy up to a size cap. Taxonomy partitions nothing &mdash;
                   it only merges &mdash; so a mislabelled read lands next to its true kin.</p>
                <div class="summary">
                    <p>Sequence clusters: <strong>{summary.get('n_families', 0):,}</strong>
                       &nbsp;&rarr;&nbsp; MSAs sent to SATIVA: <strong>{n_msa:,}</strong></p>
                    <p>Total sequences: <strong>{summary.get('total_sequences', 0):,}</strong></p>
                </div>
                {_stats_table(funnel_df)}
            </div>

            <div class="section">
                <h2 class="section-title">2. Alignment quality</h2>
                <p>Width and occupancy of every alignment SATIVA actually receives, measured on the
                   files themselves. Dashed lines mark the alert thresholds
                   ({narrow_cols} columns, {low_occ}% occupancy).</p>
                {warn}
                <div class="summary">
                    <p>Width: median <strong>{aq.get('width_median', 0):,}</strong> columns
                       (min {aq.get('width_min', 0):,}, max {aq.get('width_max', 0):,})</p>
                    <p>Occupancy: median <strong>{aq.get('occupancy_median_pct', 0)}%</strong>
                       non-gap (min {aq.get('occupancy_min_pct', 0)}%)</p>
                    <p>Below {narrow_cols} columns: <strong>{aq.get('n_narrow', 0)}</strong> MSAs
                       &nbsp;|&nbsp; below {low_occ}% occupancy: <strong>{aq.get('n_low_occupancy', 0)}</strong> MSAs</p>
                </div>
                {q_plot_html}
                <h3>Flagged alignments ({n_flagged_total} total, worst 30 shown)</h3>
                {_stats_table(flagged)}
            </div>

            <div class="section">
                <h2 class="section-title">3. What SATIVA actually sees</h2>
                <p>A cluster contributes at most <code>max_cluster_reps</code> representatives to its
                   alignment; the remaining members receive their verdict by propagation through
                   <code>group_members.tsv</code> rather than by direct placement.</p>
                <div class="summary">
                    <p>Sequences in groups: <strong>{rep.get('n_group_total', 0):,}</strong>
                       &nbsp;|&nbsp; actually aligned: <strong>{rep.get('n_aligned_total', 0):,}</strong></p>
                    <p>Verdict by propagation, not placement:
                       <strong>{rep.get('n_propagated', 0):,}</strong> sequences
                       ({rep.get('propagated_pct', 0)}%)</p>
                    <p>MSAs at the size cap: <strong>{rep.get('n_msa_at_cap', 0)}</strong> of {n_msa}</p>
                </div>
                <h3>Least represented groups</h3>
                {_stats_table(rep_df)}
            </div>

            <div class="section">
                <h2 class="section-title">4. Taxonomic composition of the groups</h2>
                <p>The strategy claims a contaminant becomes visible because it is aligned next to
                   its true relatives. That only holds if groups are taxonomically broad &mdash;
                   this is the check.</p>
                {king_plot_html}
                <p class="note">Counted from the <code>kingdom=</code> field of the alignment headers.</p>
            </div>

            <div class="section">
                <h2 class="section-title">5. Largest groups</h2>
                {_stats_table(top_df)}
            </div>
        </div>
        {page_footer()}
    </body>
    </html>
    """


def main():
    p = argparse.ArgumentParser(description="HTML report for clustered MSA construction.")
    p.add_argument("--summary", required=True)
    p.add_argument("--out-dir", required=True, help="Per-family outputs directory (out/).")
    p.add_argument("--residual-report", default=None)
    p.add_argument("--contamination-report", default=None)
    p.add_argument("--output", required=True)
    p.add_argument("--log-file")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        logger.addHandler(logging.FileHandler(args.log_file, mode="w"))

    summary = load_json(args.summary)
    msas = collect_msa_metrics(summary, args.out_dir)
    residual = load_json(args.residual_report) if args.residual_report else {}
    contamination = load_json(args.contamination_report) if args.contamination_report else {}
    html = build_html(summary, msas, residual, contamination)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info("MSA report written: %s (%d MSAs)", args.output, len(msas))


if __name__ == "__main__":
    main()
