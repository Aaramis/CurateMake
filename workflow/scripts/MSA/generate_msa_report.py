#!/usr/bin/env python3
"""
HTML report for the bottom-up MSA construction.

Consumes bottomup_summary.json (+ per-family manifests for the per-MSA metrics
+ residual_report.json + contamination_report.json) and produces a report that
matches the styling of the other pipeline reports (same container/section/stats-table
shell). Focused on:
  - distribution of the emergent STOP LEVELS + SATIVA coverage,
  - MSA quality (ratio/occupancy, gate pass rate),
  - residual (orphans reclassified by hmmsearch),
  - coarse contamination (cross-kingdom/phylum),
  - largest MSAs.
"""

import argparse
import base64
import glob
import json
import logging
import os
from datetime import datetime
from io import BytesIO
from pathlib import Path

import pandas as pd

logger = logging.getLogger("generate_msa_report")


def load_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError, TypeError):
        return default if default is not None else {}


def collect_msa_metrics(out_dir):
    """Walk the per-family manifests + the climbed manifest -> DataFrame of MSA entries."""
    entries = []
    for mpath in glob.glob(os.path.join(out_dir, "*", "manifest.json")):
        for e in load_json(mpath).get("msas", []):
            entries.append(e)
    climbed = load_json(os.path.join(os.path.dirname(out_dir.rstrip("/")), "climbed_manifest.json"))
    for e in climbed.get("msas", []):
        if str(e.get("msa_id", "")).startswith("_climbed__"):
            entries.append(e)
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
    ax.bar(labels, [counts[k] for k in labels], color="#3498db")
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


def build_html(summary, msas, residual, contamination):
    stop = summary.get("stop_levels", {})
    n_msa = summary.get("n_msa_total", len(msas))
    n_res = summary.get("n_residual_total", 0)
    gate = summary.get("gate_pass_total", 0)
    per_barcode = summary.get("per_barcode", {})

    # --- stop-level distribution ---
    stop_df = pd.DataFrame(
        [{"Stop level": k, "# MSA": v, "%": f"{100*v/max(n_msa,1):.1f}%"}
         for k, v in sorted(stop.items(), key=lambda x: -x[1])]
    )
    stop_img = _bar_plot(dict(sorted(stop.items(), key=lambda x: -x[1])),
                         "Emergent stop levels", "stop level")

    # --- SATIVA coverage per barcode ---
    cov_df = pd.DataFrame([
        {"Barcode": bc, "Reference seqs": v.get("ref_seqs", 0), "MSAs": v.get("n_msa", 0),
         "Residual": v.get("residual", 0), "Re-inserted": v.get("reinserted", 0),
         "Unassigned": v.get("unassigned", 0)}
        for bc, v in sorted(per_barcode.items())
    ])

    # --- quality ---
    q_img = ""
    if not msas.empty and "ratio" in msas.columns:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.scatter(msas["ratio"], msas.get("occ"), s=12, alpha=0.5, color="#3498db")
            ax.set_xlabel("cols / length ratio")
            ax.set_ylabel("% columns >=50% occupied")
            ax.set_title("Alignment quality (each point = one MSA)")
            fig.tight_layout()
            buf = BytesIO(); fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
            plt.close(fig); q_img = base64.b64encode(buf.getvalue()).decode("utf-8")
        except Exception:
            q_img = ""

    # --- largest MSAs ---
    top_df = pd.DataFrame()
    if not msas.empty and "n" in msas.columns:
        top_df = (msas.sort_values("n", ascending=False)
                  .head(15)[["msa_id", "level", "taxon", "n", "n_taxa", "ratio", "occ"]]
                  if "n_taxa" in msas.columns else
                  msas.sort_values("n", ascending=False).head(15)[["msa_id", "level", "n", "ratio", "occ"]])

    # --- contamination table ---
    contam_rows = pd.DataFrame([
        {"msa_id": f["msa_id"], "rank": f["rank"], "label": f["label"],
         "homology (true)": f["homology"], "consensus": f.get("consensus", ""),
         "score": f["score"], "re-homed to": f.get("rehomed_to", "")}
        for f in contamination.get("flagged", [])[:50]
    ])

    stop_plot_html = (f'<div class="plot"><h3>Stop-level distribution</h3>'
                      f'<img src="data:image/png;base64,{stop_img}" alt="stop levels"></div>'
                      if stop_img else "")
    q_plot_html = (f'<div class="plot"><h3>Alignment quality</h3>'
                   f'<img src="data:image/png;base64,{q_img}" alt="quality"></div>'
                   if q_img else "")

    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>MSA Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 0; padding: 20px; }}
            .container {{ max-width: 1200px; margin: 0 auto; }}
            .header {{ text-align: center; margin-bottom: 30px; }}
            .section {{ margin-bottom: 40px; }}
            .section-title {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-top: 30px; }}
            .plot {{ text-align: center; margin: 20px 0; box-shadow: 0 4px 8px rgba(0,0,0,0.1); padding: 15px; background: white; border-radius: 5px; }}
            .plot img {{ max-width: 100%; height: auto; }}
            .stats-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; box-shadow: 0 2px 3px rgba(0,0,0,0.1); }}
            .stats-table th, .stats-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            .stats-table th {{ background-color: #f2f2f2; position: sticky; top: 0; }}
            .stats-table tr:nth-child(even) {{ background-color: #f9f9f9; }}
            .stats-table tr:hover {{ background-color: #f1f1f1; }}
            .summary {{ background-color: #f8f9fa; padding: 15px; border-left: 4px solid #3498db; margin: 20px 0; }}
            .subsection {{ margin: 25px 0; padding: 15px; background-color: #f8f9fa; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
            .subsection h3 {{ color: #2c3e50; margin-top: 0; padding-bottom: 10px; border-bottom: 1px solid #e1e4e8; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>MSA Analysis Report</h1>
                <p>Bottom-up construction &mdash; generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>

            <div class="section">
                <h2 class="section-title">Analysis summary</h2>
                <div class="summary">
                    <p>Families: <strong>{summary.get('n_families', '?'):,}</strong>
                       &nbsp;|&nbsp; Total sequences: <strong>{summary.get('total_sequences', 0):,}</strong></p>
                    <p>Reference MSAs: <strong>{n_msa:,}</strong>
                       &nbsp;|&nbsp; passing the quality gate: <strong>{gate:,}</strong></p>
                    <p>Cross-family climbed merges: <strong>{summary.get('climbed', {})}</strong></p>
                    <p>Residual orphans: <strong>{n_res:,}</strong>
                       (re-inserted by hmmsearch: <strong>{summary.get('n_reinserted_total', 0):,}</strong>,
                        unassigned: <strong>{summary.get('n_unassigned_total', 0):,}</strong>)</p>
                    <p>Sequences without a family: <strong>{summary.get('n_no_family', 0):,}</strong>
                       (recovered by hmmsearch into nearby MSAs at the residual step)</p>
                    <p>MSAs sent to SATIVA: <strong>{summary.get('sativa_selected', 0):,}</strong>
                       (every MSA with &ge; sativa_min_seqs sequences)</p>
                </div>
            </div>

            <div class="section">
                <h2 class="section-title">Emergent stop levels</h2>
                <p>How far each lineage climbed before the quality gate stopped it:</p>
                {_stats_table(stop_df)}
                {stop_plot_html}
            </div>

            <div class="section">
                <h2 class="section-title">SATIVA coverage by barcode</h2>
                {_stats_table(cov_df)}
            </div>

            <div class="section">
                <h2 class="section-title">Alignment quality</h2>
                {q_plot_html if q_plot_html else '<p><em>No per-MSA metrics available.</em></p>'}
            </div>

            <div class="section">
                <h2 class="section-title">Residual (reclassified by hmmsearch)</h2>
                <div class="summary">
                    <p>{residual.get('n_residual', 0):,} residual orphans &rarr;
                       {residual.get('n_classified', 0):,} classified,
                       {residual.get('n_integrated', 0):,} integrated,
                       {residual.get('n_unassigned', 0):,} unassigned
                       (against {residual.get('n_profiles', 0):,} profiles).</p>
                </div>
            </div>

            <div class="section">
                <h2 class="section-title">Coarse contamination (cross-kingdom/phylum)</h2>
                <p>Tier complementary to SATIVA: MSAs whose <strong>homology contradicts the label</strong>
                   at a high rank (a fungus labelled as a plant, etc.) &mdash; which SATIVA cannot detect
                   because they are isolated by the clustering.</p>
                <p><strong>{contamination.get('n_flagged', 0)}</strong> MSAs flagged out of
                   {contamination.get('n_msa', 0)} checked
                   ({contamination.get('n_rehomed_msa', 0)} re-homed).</p>
                {_stats_table(contam_rows) if not contam_rows.empty else ''}
            </div>

            <div class="section">
                <h2 class="section-title">Largest MSAs</h2>
                {_stats_table(top_df)}
            </div>
        </div>
    </body>
    </html>
    """


def main():
    p = argparse.ArgumentParser(description="HTML report for bottom-up MSA construction.")
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
    msas = collect_msa_metrics(args.out_dir)
    residual = load_json(args.residual_report) if args.residual_report else {}
    contamination = load_json(args.contamination_report) if args.contamination_report else {}
    html = build_html(summary, msas, residual, contamination)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info("Bottom-up report written: %s (%d MSAs)", args.output, len(msas))


if __name__ == "__main__":
    main()
