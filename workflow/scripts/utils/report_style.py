#!/usr/bin/env python3
"""
Shared visual identity for the CurateMake HTML reports.

Every report generator injects REPORT_CSS instead of carrying its own stylesheet, so the
whole set reads as one document: same palette, same type scale, same cards, tables and
callouts. Class names are the union of those the reports already used, and synonyms
(kpi-card / stat-card / metric, warn / warning / warning-box, plot / plot-container /
chart-container, summary / summary-box) are styled identically rather than renamed, so a
generator can adopt this file without touching its markup.

Usage:
    from report_style import REPORT_CSS, page_header, page_footer
    html = f"<style>{REPORT_CSS}</style>" + page_header("Title", "subtitle") + body

Tokens live in :root; a report needing more room can add the `wide` class to .container.
"""

CURATEMAKE_PRIMARY = "#667eea"
CURATEMAKE_SECONDARY = "#764ba2"

# Chart palette, so figures sit in the same colour world as the page around them.
# Mirrors the CSS tokens below; import it instead of hard-coding hex in a plot.
CHART_COLORS = {
    "primary": "#667eea",
    "secondary": "#764ba2",
    "success": "#2e9e5b",
    "warning": "#ef8c1a",
    "danger": "#dc2626",
    "info": "#3aa0d1",
    "muted": "#6b7280",
}

# Ordered sequence for categorical series (barcodes, sources, ranks...).
CHART_SEQUENCE = ["#667eea", "#764ba2", "#3aa0d1", "#2e9e5b",
                  "#ef8c1a", "#dc2626", "#6b7280"]

REPORT_CSS = """
:root {
    --cm-primary: #667eea;
    --cm-secondary: #764ba2;
    --cm-gradient: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    --cm-ink: #2c3e50;
    --cm-muted: #6b7280;
    --cm-line: #e5e7eb;
    --cm-surface: #ffffff;
    --cm-canvas: #f5f6fa;
    --cm-info: #2196f3;
    --cm-info-bg: #e3f2fd;
    --cm-ok: #2e9e5b;
    --cm-ok-bg: #e8f5e9;
    --cm-warn: #ef8c1a;
    --cm-warn-bg: #fff3e0;
    --cm-danger: #dc2626;
    --cm-danger-bg: #fef2f2;
    --cm-radius: 10px;
    --cm-shadow: 0 2px 10px rgba(0,0,0,0.06);
    --cm-font: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu,
               Cantarell, 'Helvetica Neue', sans-serif;
    --cm-mono: 'SFMono-Regular', Menlo, Consolas, 'Courier New', monospace;
}

* { box-sizing: border-box; }

body {
    font-family: var(--cm-font);
    line-height: 1.6;
    color: var(--cm-ink);
    background: var(--cm-canvas);
    margin: 0;
    padding: 0;
}

/* ---------- page furniture ---------- */
.header {
    background: var(--cm-gradient);
    color: #fff;
    padding: 48px 20px;
    text-align: center;
    box-shadow: 0 4px 6px rgba(0,0,0,0.1);
}
.header h1 { font-size: 2.2em; margin: 0 0 8px 0; font-weight: 700; }
.header p  { font-size: 1.02em; margin: 0; color: rgba(255,255,255,0.92); }

.container { max-width: 1400px; margin: 0 auto; padding: 32px 20px; }
.container.wide { max-width: 1800px; }

.section {
    background: var(--cm-surface);
    border-radius: var(--cm-radius);
    padding: 28px;
    margin-bottom: 28px;
    box-shadow: var(--cm-shadow);
}
.section h2, .section-title {
    color: var(--cm-ink);
    font-size: 1.5em;
    margin: 0 0 18px 0;
    padding-bottom: 10px;
    border-bottom: 3px solid var(--cm-primary);
}
.section h3 { color: var(--cm-secondary); font-size: 1.15em; margin: 24px 0 12px 0; }

.footer, .timestamp {
    text-align: center;
    color: var(--cm-muted);
    font-size: 0.9em;
    padding: 24px 20px;
}

/* ---------- metric cards ----------
   One hero metric carries the gradient; the rest are quiet white cards with a coloured
   rule on top. Six identical saturated blocks give the eye nothing to hold on to. */
.kpi-grid, .stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
    gap: 16px;
    margin: 24px 0;
    align-items: stretch;
}
.kpi-card, .stat-card, .stage-card, .metric {
    background: var(--cm-surface);
    color: var(--cm-ink);
    padding: 20px 18px;
    border-radius: var(--cm-radius);
    border: 1px solid var(--cm-line);
    border-top: 3px solid var(--cm-primary);
    box-shadow: var(--cm-shadow);
    display: flex;
    flex-direction: column;
    justify-content: flex-start;
    gap: 4px;
    min-height: 118px;
}
/* the one number the page is really about */
.kpi-card.primary, .stat-card.primary {
    background: var(--cm-gradient);
    border: none;
    color: #fff;
}
.kpi-card.primary .kpi-label, .kpi-card.primary .kpi-value { color: #fff; }

.kpi-icon { font-size: 1.1em; line-height: 1; opacity: 0.75; margin-bottom: 2px; }
.kpi-value, .stat-value, .metric-value, .kpi-card .number, .stat-card .number,
.stage-card .number, .metric .number {
    font-size: 1.85em;
    font-weight: 650;
    line-height: 1.15;
    margin: 0;
    font-variant-numeric: tabular-nums;
    letter-spacing: -0.01em;
    overflow-wrap: anywhere;
}
.kpi-label, .stat-label, .metric-label, .kpi-card p, .stat-card p, .stage-card p {
    font-size: 0.78em;
    line-height: 1.35;
    text-transform: none;
    letter-spacing: 0;
    color: var(--cm-muted);
    margin: 0;
}
.kpi-card.primary .kpi-label { color: rgba(255,255,255,0.88); }
.kpi-card h5, .stat-card h5, .stage-card h5 {
    margin: 0; font-size: 0.78em; font-weight: 600; color: var(--cm-muted);
}
/* accent per position, so a row of cards reads as a sequence rather than a wall */
.kpi-grid > *:nth-child(2), .stats-grid > *:nth-child(2) { border-top-color: var(--cm-secondary); }
.kpi-grid > *:nth-child(3), .stats-grid > *:nth-child(3) { border-top-color: var(--cm-info); }
.kpi-grid > *:nth-child(4), .stats-grid > *:nth-child(4) { border-top-color: var(--cm-ok); }
.kpi-grid > *:nth-child(5), .stats-grid > *:nth-child(5) { border-top-color: var(--cm-warn); }
.kpi-grid > *:nth-child(6), .stats-grid > *:nth-child(6) { border-top-color: var(--cm-muted); }

/* ---------- callouts ---------- */
.info-box, .summary, .summary-box, .interpretation, .tool-description, .note,
.warning, .warning-box, .warn, .success, .stage-details {
    border-left: 4px solid var(--cm-info);
    background: var(--cm-info-bg);
    padding: 14px 18px;
    margin: 18px 0;
    border-radius: 6px;
}
.summary, .summary-box, .stage-details { border-left-color: var(--cm-primary); background: #f4f5fd; }
.interpretation, .tool-description { border-left-color: var(--cm-secondary); background: #f7f4fb; }
.warning, .warning-box, .info-box.warning { border-left-color: var(--cm-warn); background: var(--cm-warn-bg); }
.warn { border-left-color: var(--cm-danger); background: var(--cm-danger-bg); }
.success, .info-box.success { border-left-color: var(--cm-ok); background: var(--cm-ok-bg); }
.note { border-left-color: var(--cm-line); background: #fafafa; color: var(--cm-muted); font-size: 0.9em; }

/* ---------- tables ---------- */
.table-wrapper { overflow-x: auto; margin: 18px 0; }
table, .stats-table {
    width: 100%;
    border-collapse: collapse;
    margin: 18px 0;
    background: var(--cm-surface);
    font-size: 0.95em;
}
th, td, .stats-table th, .stats-table td {
    border: 1px solid var(--cm-line);
    padding: 9px 12px;
    text-align: left;
}
th, .stats-table th {
    background: #f4f5fd;
    color: var(--cm-ink);
    font-weight: 600;
    position: sticky;
    top: 0;
}
tr:nth-child(even), .stats-table tr:nth-child(even) { background: #fbfbfd; }
tr:hover, .stats-table tr:hover { background: #f4f5fd; }
td.number, .number { font-variant-numeric: tabular-nums; }
code { font-family: var(--cm-mono); background: #f4f5fd; padding: 1px 5px; border-radius: 4px; }

/* ---------- figures ---------- */
.plot, .plot-container, .chart-container {
    background: var(--cm-surface);
    border-radius: var(--cm-radius);
    padding: 16px;
    margin: 20px 0;
    box-shadow: var(--cm-shadow);
    text-align: center;
}
.plot img, .plot-container img, .chart-container img { max-width: 100%; height: auto; }

@media (max-width: 768px) {
    .container { padding: 20px 14px; }
    .header { padding: 32px 16px; }
    .header h1 { font-size: 1.6em; }
    .kpi-grid, .stats-grid { grid-template-columns: 1fr; }
}
"""


def human_duration(seconds: float) -> str:
    """Compact, non-wrapping duration: '15 j 04 h', '3 j 19 h', '9 h 12', '47 s'.

    A KPI card is ~200px wide; 'X days, HH:MM:SS' wraps onto two lines there and the
    seconds carry no information at that scale.
    """
    s = max(0, int(round(seconds)))
    d, rem = divmod(s, 86400)
    h, rem = divmod(rem, 3600)
    m, sec = divmod(rem, 60)
    if d:
        return f"{d} j {h:02d} h"
    if h:
        return f"{h} h {m:02d}"
    if m:
        return f"{m} min {sec:02d}"
    return f"{sec} s"


def page_header(title: str, subtitle: str = "") -> str:
    """Banner shared by every report."""
    sub = f"<p>{subtitle}</p>" if subtitle else ""
    return f'<div class="header"><h1>{title}</h1>{sub}</div>'


def page_footer(note: str = "") -> str:
    """Footer shared by every report."""
    extra = f" &mdash; {note}" if note else ""
    return (f'<div class="footer">CurateMake &mdash; ITS taxonomy curation pipeline{extra}</div>')
