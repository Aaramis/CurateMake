#!/usr/bin/env python3
"""
Generates an interactive HTML report that tracks the quantitative flow of sequences
through every step of the bioinformatics pipeline.

The script builds a Plotly Sankey diagram that highlights:
- Data sources (UNITE, BOLD, PLANiTS, CALeDNA)
- The merge into the unified database
- Pre-ITSx deduplication
- ITSx processing (ITS1, ITS2, ITS detection)
- Post-ITSx deduplication
- MSA preparation and processing
- Final per-barcode results
"""

import json
import logging
import sys
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent / "utils"))
from report_style import REPORT_CSS  # shared CurateMake report styling
from datetime import datetime
import plotly.graph_objects as go
import numpy as np

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================
# UTILITY FUNCTIONS
# ============================================

def safe_get(dictionary, keys, default=0):
    """
    Retrieve a value from a nested dict structure with error handling.

    Args:
        dictionary: input dict
        keys: list of keys or dotted string (e.g. 'a.b.c')
        default: value to return when the key is missing

    Returns:
        Extracted value or the provided default
    """
    try:
        value = dictionary
        key_list = keys.split('.') if isinstance(keys, str) else keys
        for key in key_list:
            value = value[key]
        return value if value is not None else default
    except (KeyError, TypeError, AttributeError):
        return default


def load_json_safe(filepath, default=None):
    """
    Load a JSON file with proper error handling.

    Args:
        filepath: JSON file path
        default: fallback value when loading fails

    Returns:
        Parsed JSON content or the default value
    """
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
            logger.info(f"✓ Loaded: {Path(filepath).name}")
            return data
    except FileNotFoundError:
        logger.error(f"✗ File not found: {filepath}")
        return default or {}
    except json.JSONDecodeError as e:
        logger.error(f"✗ Invalid JSON in {filepath}: {e}")
        return default or {}


def extract_database_sequences(formatting_report):
    """
    Extract sequence counts per database.
    Handles name variations (e.g., caledna_mini_FITS vs caledna_FITS).

    Returns:
        dict: {db_name: sequence_count}
    """
    sequences_by_db = {}
    caledna_total = 0

    for db_info in formatting_report.get("databases_processed", []):
        db_name = db_info.get("database", "")
        processed_seqs = safe_get(db_info, ['stats', 'processed_sequences'], 0)

        # Handle CALeDNA name variations
        if "PITS" in db_name:
            sequences_by_db["caledna_PITS"] = processed_seqs
            caledna_total += processed_seqs
        elif "FITS" in db_name:
            sequences_by_db["caledna_FITS"] = processed_seqs
            caledna_total += processed_seqs
        else:
            # Keep every other source under its own name. A whitelist would silently drop
            # any database it does not know about, so the source nodes would no longer sum
            # to the merge total.
            base_name = db_name.lower()
            for known in ['unite', 'bold', 'planits']:
                if known in base_name:
                    sequences_by_db[known] = processed_seqs
                    break
            else:
                sequences_by_db[db_name] = processed_seqs

    sequences_by_db["CALeDNA"] = caledna_total
    return sequences_by_db


# ============================================
# FIGURE GENERATION
# ============================================

def hex_to_rgba(hex_color: str, alpha: float) -> str:
    """Convert a hex color to rgba."""
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f'rgba({r},{g},{b},{alpha})'


def generate_ribbon_sankey(x0, y0, x1, y1, width, num_points=100, curvature=0.4):
    """
    Generate a ribbon-style Sankey link with constant width.

    Args:
        x0, y0: start point
        x1, y1: end point
        width: ribbon width
        num_points: number of points for the curve
        curvature: bezier curvature (0-1)

    Returns:
        tuple: (ribbon_x, ribbon_y) polygon coordinates
    """
    t = np.linspace(0, 1, num_points)

    # Bezier cubic control points
    cx0 = x0 + curvature * (x1 - x0)
    cy0 = y0
    cx1 = x1 - curvature * (x1 - x0)
    cy1 = y1

    # Base curve (cubic Bezier)
    x = (1-t)**3 * x0 + 3*(1-t)**2*t * cx0 + 3*(1-t)*t**2 * cx1 + t**3 * x1
    y = (1-t)**3 * y0 + 3*(1-t)**2*t * cy0 + 3*(1-t)*t**2 * cy1 + t**3 * y1

    # Ribbon borders (constant vertical offset)
    x_upper = x
    y_upper = y + width/2
    x_lower = x
    y_lower = y - width/2

    # Close the polygon
    ribbon_x = np.concatenate([x_upper, x_lower[::-1]])
    ribbon_y = np.concatenate([y_upper, y_lower[::-1]])

    return ribbon_x, ribbon_y


def add_tool_rectangle(fig, col_start, col_end, row_start, row_end, title,
                       color="rgba(200,200,200,0.15)", border_color="rgba(150,150,150,0.4)",
                       nb_col=9, nb_row=6):
    """
    Add an annotated rectangle to highlight tool zones in the figure.

    Args:
        fig: Plotly figure
        col_start, col_end: start/end columns
        row_start, row_end: start/end rows
        title: zone label
        color: fill color
        border_color: outline color
        nb_col, nb_row: grid dimensions
    """
    x_min = (col_start - 0.3) / nb_col
    x_max = (col_end + 0.3) / nb_col
    y_min = (nb_row - row_start + 0.4) / nb_row
    y_max = (nb_row - row_end - 0.7) / nb_row

    fig.add_shape(
        type="rect", x0=x_min, y0=y_min, x1=x_max, y1=y_max,
        line=dict(color=border_color, width=2, dash="dot"),
        fillcolor=color, layer="below"
    )
    fig.add_annotation(
        x=(x_min + x_max) / 2, y=y_min - 0.06,
        text=f"<b>{title}</b>", showarrow=False,
        font=dict(size=16, color="#2c3e50", family="Arial, sans-serif"),
        align="center", xanchor="center", yanchor="bottom"
    )


# ============================================
# MAIN
# ============================================

def main():
    """Main extraction and report generation entrypoint."""

    logger.info("=" * 60)
    logger.info("Sequence Flow Report Generation")
    logger.info("=" * 60)

    # Fetch Snakemake parameters
    try:
        report_paths = {
            "formatting_summary": snakemake.input.formatting_summary,
            "itsx_deduplicate_pre": snakemake.input.itsx_deduplicate_pre,
            "itsx_deduplicate_post": snakemake.input.itsx_deduplicate_post,
            "itsx_process": snakemake.input.itsx_process,
            "itsx_prepare": snakemake.input.itsx_prepare,
            "msa_prepare": snakemake.input.msa_prepare,
            "bottomup_summary": getattr(snakemake.input, "bottomup_summary", None),
        }
        output_html = snakemake.output.report
    except NameError:
        logger.error("This script must be run through Snakemake")
        sys.exit(1)

    # ========================================
    # DATA LOADING
    # ========================================
    logger.info("\n📥 Loading JSON reports...")
    reports = {key: (load_json_safe(path) if path else {}) for key, path in report_paths.items()}

    # ========================================
    # METRIC EXTRACTION
    # ========================================
    logger.info("\n📊 Extracting sequence flow metrics...")

    # --- Data sources ---
    sequences_by_db = extract_database_sequences(reports["formatting_summary"])
    total_sequence_ini = safe_get(reports["formatting_summary"],
                                   ['overall_stats', 'total_sequences_processed'], 0)

    logger.info(f"  Sources: UNITE={sequences_by_db.get('unite', 0)}, "
                f"BOLD={sequences_by_db.get('bold', 0)}, "
                f"PLANiTS={sequences_by_db.get('planits', 0)}, "
                f"CALeDNA={sequences_by_db.get('CALeDNA', 0)}")
    # CALeDNA is reported both split (PITS/FITS) and merged; count it once.
    n_databases_integrated = len([n for n in sequences_by_db
                                  if n not in ('caledna_PITS', 'caledna_FITS')
                                  and sequences_by_db[n]])
    logger.info(f"  Total initial: {total_sequence_ini:,} from "
                f"{n_databases_integrated} database(s): {sorted(sequences_by_db)}")

    # --- Pre-ITSx deduplication ---
    unique_sequences_pre = safe_get(reports["itsx_deduplicate_pre"],
                                     ['summary', 'unique_sequences_kept'], 0)
    duplicate_sequences_pre = safe_get(reports["itsx_deduplicate_pre"],
                                        ['summary', 'duplicates_marked'], 0)
    duplicate_sequences_pre_ITS = safe_get(reports["itsx_deduplicate_pre"],
                                            ['detailed_summary', 'duplicates_by_barcode', 'ITS'], 0)
    duplicate_sequences_pre_ITS1 = safe_get(reports["itsx_deduplicate_pre"],
                                             ['detailed_summary', 'duplicates_by_barcode', 'ITS1'], 0)
    duplicate_sequences_pre_ITS2 = safe_get(reports["itsx_deduplicate_pre"],
                                             ['detailed_summary', 'duplicates_by_barcode', 'ITS2'], 0)

    # A pre-ITSx dedup report rebuilt from the database after ITSx carries post-ITSx
    # counts: it would put more sequences out of the merge node than went in. The number
    # of sequences that actually survived this pass is what ITSx was then fed, which the
    # extraction report records first-hand, so derive the split from there instead.
    dedup_pre_is_measured = "_note" not in (reports.get("itsx_deduplicate_pre") or {})
    if not dedup_pre_is_measured:
        _fed_to_itsx = safe_get(reports["itsx_prepare"], ['stats', 'total_sequences'], 0)
        if _fed_to_itsx:
            unique_sequences_pre = _fed_to_itsx
            duplicate_sequences_pre = max(0, total_sequence_ini - _fed_to_itsx)
            # per-barcode duplicate counts are not recoverable this way
            duplicate_sequences_pre_ITS = duplicate_sequences_pre_ITS1 = duplicate_sequences_pre_ITS2 = 0
            logger.warning("  pre-ITSx dedup report was reconstructed after the fact; "
                           "deriving the split from the ITSx extraction report instead "
                           f"({unique_sequences_pre:,} kept, {duplicate_sequences_pre:,} duplicates)")

    # Fallback when no pre-ITSx deduplication is available (mini datasets)
    if unique_sequences_pre == 0:
        unique_sequences_pre = total_sequence_ini
        logger.warning("  ⚠️  No pre-ITSx deduplication data, using total sequences")

    logger.info(f"  Pre-ITSx: {unique_sequences_pre:,} unique, {duplicate_sequences_pre:,} duplicates")

    # --- ITSx preparation ---
    ITSX_ITS1 = safe_get(reports["itsx_prepare"],
                         ['stats', 'barcode_details', 'ITS1', 'sequence_count'], 0)
    ITSX_ITS2 = safe_get(reports["itsx_prepare"],
                         ['stats', 'barcode_details', 'ITS2', 'sequence_count'], 0)
    ITSX_ITS = safe_get(reports["itsx_prepare"],
                        ['stats', 'barcode_details', 'ITS', 'sequence_count'], 0)

    logger.info(f"  Pre-ITSx split: ITS1={ITSX_ITS1}, ITS2={ITSX_ITS2}, ITS={ITSX_ITS}")

    # --- ITSx results ---
    itsx_barcodes = safe_get(reports["itsx_process"], ['barcodes'], {})

    # ITS1 barcode processing
    ITSX_ITS1_ITS1 = safe_get(itsx_barcodes, ['ITS1', 'subregions', 'ITS1', 'inserted'], 0)
    ITSX_ITS1_ITS2 = safe_get(itsx_barcodes, ['ITS1', 'subregions', 'ITS2', 'inserted'], 0)
    ITSX_ITS1_ITS = safe_get(itsx_barcodes, ['ITS1', 'reconstructed_its', 'inserted'], 0)
    ITSX_ITS1_no_detection = safe_get(itsx_barcodes, ['ITS1', 'no_detection'], 0)

    # ITS2 barcode processing
    ITSX_ITS2_ITS1 = safe_get(itsx_barcodes, ['ITS2', 'subregions', 'ITS1', 'inserted'], 0)
    ITSX_ITS2_ITS2 = safe_get(itsx_barcodes, ['ITS2', 'subregions', 'ITS2', 'inserted'], 0)
    ITSX_ITS2_ITS = safe_get(itsx_barcodes, ['ITS2', 'reconstructed_its', 'inserted'], 0)
    ITSX_ITS2_no_detection = safe_get(itsx_barcodes, ['ITS2', 'no_detection'], 0)

    # ITS barcode processing
    ITSX_ITS_ITS1 = safe_get(itsx_barcodes, ['ITS', 'subregions', 'ITS1', 'inserted'], 0)
    ITSX_ITS_ITS2 = safe_get(itsx_barcodes, ['ITS', 'subregions', 'ITS2', 'inserted'], 0)
    ITSX_ITS_ITS = safe_get(itsx_barcodes, ['ITS', 'reconstructed_its', 'inserted'], 0)
    ITSX_ITS_no_detection = safe_get(itsx_barcodes, ['ITS', 'no_detection'], 0)

    # Totals after ITSx processing
    ITSX_ITS1_new = ITSX_ITS1_ITS1 + ITSX_ITS2_ITS1 + ITSX_ITS_ITS1
    ITSX_ITS2_new = ITSX_ITS1_ITS2 + ITSX_ITS2_ITS2 + ITSX_ITS_ITS2
    ITSX_ITS_new = ITSX_ITS1_ITS + ITSX_ITS2_ITS + ITSX_ITS_ITS
    ITSX_ITS_no_detection_tot = ITSX_ITS1_no_detection + ITSX_ITS2_no_detection + ITSX_ITS_no_detection

    logger.info(f"  Post-ITSx: ITS1={ITSX_ITS1_new}, ITS2={ITSX_ITS2_new}, ITS={ITSX_ITS_new}, "
                f"no_detection={ITSX_ITS_no_detection_tot}")

    # --- Post-ITSx deduplication ---
    unique_sequences_post = safe_get(reports["itsx_deduplicate_post"],
                                      ['summary', 'unique_sequences_kept'], 0)
    duplicate_sequences_post = safe_get(reports["itsx_deduplicate_post"],
                                         ['detailed_summary', 'sequences_by_status', 'Duplicate'], 0)
    duplicate_sequences_post_ITS = safe_get(reports["itsx_deduplicate_post"],
                                             ['detailed_summary', 'duplicates_by_barcode', 'ITS'], 0)
    duplicate_sequences_post_ITS1 = safe_get(reports["itsx_deduplicate_post"],
                                              ['detailed_summary', 'duplicates_by_barcode', 'ITS1'], 0)
    duplicate_sequences_post_ITS2 = safe_get(reports["itsx_deduplicate_post"],
                                              ['detailed_summary', 'duplicates_by_barcode', 'ITS2'], 0)

    # Deduplication deltas. They only mean anything if the pre-ITSx report is a genuine
    # measurement: a report rebuilt from the database after the post-ITSx pass necessarily
    # carries post-ITSx counts, which would make every delta zero by construction. Such a
    # report is flagged with a `_note`, and the deltas are then reported as unavailable
    # rather than as a measured zero.
    # How many sequences the post-ITSx pass marked: the final duplicate count minus what
    # the first pass had already marked. This total is sound either way.
    duplicate_sequences_post_tot = max(0, duplicate_sequences_post - duplicate_sequences_pre)

    if dedup_pre_is_measured:
        diff_post_pre_ITS1 = duplicate_sequences_post_ITS1 - duplicate_sequences_pre_ITS1
        diff_post_pre_ITS2 = duplicate_sequences_post_ITS2 - duplicate_sequences_pre_ITS2
        diff_post_pre_ITS = duplicate_sequences_post_ITS - duplicate_sequences_pre_ITS
    else:
        # The subtraction above is the right one, but its baseline here is a post-ITSx
        # snapshot: it would return zero by construction rather than by measurement. Draw
        # no ribbon and let the table say "unavailable"; the per-barcode split comes back
        # on its own as soon as a genuine pre-ITSx report exists.
        diff_post_pre_ITS1 = diff_post_pre_ITS2 = diff_post_pre_ITS = 0

    def _fmt_delta(value):
        return f"{value:,}" if dedup_pre_is_measured else "unavailable"

    pre_dup_ITS1 = _fmt_delta(duplicate_sequences_pre_ITS1)
    pre_dup_ITS2 = _fmt_delta(duplicate_sequences_pre_ITS2)
    pre_dup_ITS = _fmt_delta(duplicate_sequences_pre_ITS)

    new_dup_ITS1 = _fmt_delta(diff_post_pre_ITS1)
    new_dup_ITS2 = _fmt_delta(diff_post_pre_ITS2)
    new_dup_ITS = _fmt_delta(diff_post_pre_ITS)
    new_dup_total = f"{duplicate_sequences_post_tot:,}"
    if not dedup_pre_is_measured:
        logger.warning("  pre-ITSx dedup report was reconstructed: new-duplicate counts "
                       "cannot be derived and are reported as unavailable")

    logger.info(f"  Post-ITSx dedup: {unique_sequences_post:,} unique, "
                f"{duplicate_sequences_post:,} duplicates (new: {duplicate_sequences_post_tot:,})")

    # --- MSA preparation ---
    MSA_unique_sequences = safe_get(reports["msa_prepare"], ['stats', 'total_sequences'], 0)
    MSA_pre_ITS = safe_get(reports["msa_prepare"], ['stats', 'barcode_distribution', 'ITS'], 0)
    MSA_pre_ITS1 = safe_get(reports["msa_prepare"], ['stats', 'barcode_distribution', 'ITS1'], 0)
    MSA_pre_ITS2 = safe_get(reports["msa_prepare"], ['stats', 'barcode_distribution', 'ITS2'], 0)

    logger.info(f"  MSA prep: Total={MSA_unique_sequences:,}, ITS1={MSA_pre_ITS1}, "
                f"ITS2={MSA_pre_ITS2}, ITS={MSA_pre_ITS}")

    # --- MSA Processing ---
    _bu_summary = reports.get("bottomup_summary") or {}
    # Clustered MSA: one clustering per barcode. A sequence either ends up in a reference
    # MSA or it does not -- there is no cross-barcode reassignment and no re-insertion
    # tier -- so the quantity worth showing is how many were dropped along the way.
    _pb = _bu_summary.get("per_barcode", {})
    MSA_post_ITS1 = safe_get(_pb, ['ITS1', 'ref_seqs'], 0)
    MSA_post_ITS2 = safe_get(_pb, ['ITS2', 'ref_seqs'], 0)
    MSA_post_ITS = safe_get(_pb, ['ITS', 'ref_seqs'], 0)

    # Dropped between extraction and the reference MSAs: clusters below min_ref that could
    # not be absorbed, groups below sativa_min_seqs, members beyond max_cluster_reps.
    MSA_dropped_ITS1 = max(0, MSA_pre_ITS1 - MSA_post_ITS1)
    MSA_dropped_ITS2 = max(0, MSA_pre_ITS2 - MSA_post_ITS2)
    MSA_dropped_ITS = max(0, MSA_pre_ITS - MSA_post_ITS)

    msa_drop_rate_ITS1 = (100 * MSA_dropped_ITS1 / MSA_pre_ITS1) if MSA_pre_ITS1 else 0.0
    msa_drop_rate_ITS2 = (100 * MSA_dropped_ITS2 / MSA_pre_ITS2) if MSA_pre_ITS2 else 0.0
    msa_drop_rate_ITS = (100 * MSA_dropped_ITS / MSA_pre_ITS) if MSA_pre_ITS else 0.0

    logger.info(f"  MSA final: ITS1={MSA_post_ITS1}, ITS2={MSA_post_ITS2}, ITS={MSA_post_ITS} "
                f"(dropped: {MSA_dropped_ITS1 + MSA_dropped_ITS2 + MSA_dropped_ITS:,})")

    # ========================================
    # CONFIGURATION DU GRAPHIQUE
    # ========================================
    logger.info("\n🎨 Generating Sankey diagram...")

    COLORS = {
        'source': '#3498db',
        'database': '#2c3e50',
        'unique': '#27ae60',
        'duplicate': '#95a5a6',
        'its1': '#e74c3c',
        'its2': '#f39c12',
        'its': '#9b59b6',
        'no_detection': '#7f8c8d',
    }
    LINK_OPACITY = 0.40
    NB_COL = 9
    NB_ROW = 6

    # ========================================
    # NODE DEFINITION
    # ========================================
    nodes_data = [
        # Colonne 0 - Sources
        (f"CALeDNA<br>{sequences_by_db.get('CALeDNA', 0):,}", 0/NB_COL, (NB_ROW-1)/NB_ROW, sequences_by_db.get('CALeDNA', 0), COLORS['source'], 50),
        (f"UNITE<br>{sequences_by_db.get('unite', 0):,}", 0/NB_COL, (NB_ROW-2)/NB_ROW, sequences_by_db.get('unite', 0), COLORS['source'], 50),
        (f"BOLD<br>{sequences_by_db.get('bold', 0):,}", 0/NB_COL, (NB_ROW-3)/NB_ROW, sequences_by_db.get('bold', 0), COLORS['source'], 50),
        (f"PLANiTS<br>{sequences_by_db.get('planits', 0):,}", 0/NB_COL, (NB_ROW-4)/NB_ROW, sequences_by_db.get('planits', 0), COLORS['source'], 50),
        # Colonne 1 - Merge
        (f"Merge<br>Database<br>{total_sequence_ini:,}", 1/NB_COL, (NB_ROW-2.5)/NB_ROW, total_sequence_ini, COLORS['database'], 70),
        # Colonne 2 - Dedup pre
        (f"Unique<br>Sequences<br>{unique_sequences_pre:,}", 2/NB_COL, (NB_ROW-4)/NB_ROW, unique_sequences_pre, COLORS['unique'], 65),
        (f"Duplicate<br>Sequences<br>{duplicate_sequences_pre:,}", 2/NB_COL, (NB_ROW-2)/NB_ROW, duplicate_sequences_pre, COLORS['duplicate'], 55),
        # Colonne 3 - Pre-ITSx
        (f"ITS1<br>{ITSX_ITS1:,}", 3/NB_COL, (NB_ROW-1)/NB_ROW, ITSX_ITS1, COLORS['its1'], 55),
        (f"ITS2<br>{ITSX_ITS2:,}", 3/NB_COL, (NB_ROW-2)/NB_ROW, ITSX_ITS2, COLORS['its2'], 55),
        (f"ITS<br>{ITSX_ITS:,}", 3/NB_COL, (NB_ROW-3)/NB_ROW, ITSX_ITS, COLORS['its'], 55),
        # Colonne 4 - Post-ITSx
        (f"ITS1<br>{ITSX_ITS1_new:,}", 4/NB_COL, (NB_ROW-1)/NB_ROW, ITSX_ITS1_new, COLORS['its1'], 55),
        (f"ITS2<br>{ITSX_ITS2_new:,}", 4/NB_COL, (NB_ROW-2)/NB_ROW, ITSX_ITS2_new, COLORS['its2'], 55),
        (f"ITS<br>{ITSX_ITS_new:,}", 4/NB_COL, (NB_ROW-3)/NB_ROW, ITSX_ITS_new, COLORS['its'], 55),
        (f"No Detection<br>{ITSX_ITS_no_detection_tot:,}", 4/NB_COL, (NB_ROW-0)/NB_ROW, ITSX_ITS_no_detection_tot, COLORS['no_detection'], 50),
        # Colonne 5 - Dedup post
        (f"Unique<br>Sequences<br>{unique_sequences_post:,}", 5/NB_COL, (NB_ROW-4)/NB_ROW, unique_sequences_post, COLORS['unique'], 65),
        (f"Duplicate<br>Sequences<br>{duplicate_sequences_post_tot:,}", 5/NB_COL, (NB_ROW-2)/NB_ROW, duplicate_sequences_post_tot, COLORS['duplicate'], 55),
        # Colonne 6 - MSA prep
        (f"ITS1<br>{MSA_pre_ITS1:,}", 6/NB_COL, (NB_ROW-1)/NB_ROW, MSA_pre_ITS1, COLORS['its1'], 55),
        (f"ITS2<br>{MSA_pre_ITS2:,}", 6/NB_COL, (NB_ROW-2)/NB_ROW, MSA_pre_ITS2, COLORS['its2'], 55),
        (f"ITS<br>{MSA_pre_ITS:,}", 6/NB_COL, (NB_ROW-3)/NB_ROW, MSA_pre_ITS, COLORS['its'], 55),
        # Colonne 7 - sequences dropped before the reference MSAs (terminal nodes)
        (f"Not aligned ITS1<br>{MSA_dropped_ITS1:,}", 7/NB_COL, (NB_ROW-1.25)/NB_ROW, MSA_dropped_ITS1, COLORS['no_detection'], 50),
        (f"Not aligned ITS2<br>{MSA_dropped_ITS2:,}", 7/NB_COL, (NB_ROW-2.25)/NB_ROW, MSA_dropped_ITS2, COLORS['no_detection'], 50),
        (f"Not aligned ITS<br>{MSA_dropped_ITS:,}", 7/NB_COL, (NB_ROW-3.25)/NB_ROW, MSA_dropped_ITS, COLORS['no_detection'], 50),
        # Colonne 8 - MSA final
        (f"ITS1<br>{MSA_post_ITS1:,}", 8/NB_COL, (NB_ROW-0.75)/NB_ROW, MSA_post_ITS1, COLORS['its1'], 50),
        (f"ITS2<br>{MSA_post_ITS2:,}", 8/NB_COL, (NB_ROW-1.75)/NB_ROW, MSA_post_ITS2, COLORS['its2'], 50),
        (f"ITS<br>{MSA_post_ITS:,}", 8/NB_COL, (NB_ROW-2.75)/NB_ROW, MSA_post_ITS, COLORS['its'], 50),
    ]

    # Any source beyond the four laid out above (CALeDNA/UNITE/BOLD/PLANiTS) is appended
    # here rather than inserted, so the indices the links below rely on stay valid.
    _extra_sources = [(n, v) for n, v in sequences_by_db.items()
                      if n not in ('CALeDNA', 'caledna_PITS', 'caledna_FITS',
                                   'unite', 'bold', 'planits') and v]
    _extra_source_idx = {}
    for _i, (_name, _val) in enumerate(_extra_sources):
        _extra_source_idx[_name] = len(nodes_data)
        nodes_data.append((f"{_name}<br>{_val:,}", 0/NB_COL, (NB_ROW-5)/NB_ROW,
                           _val, COLORS['source'], 50))

    # ========================================
    # LINK DEFINITION
    # ========================================
    links_data = [
        # Sources -> Full Database
        (0, 4, sequences_by_db.get('CALeDNA', 0), COLORS['source']),
        (1, 4, sequences_by_db.get('unite', 0), COLORS['source']),
        (2, 4, sequences_by_db.get('bold', 0), COLORS['source']),
        (3, 4, sequences_by_db.get('planits', 0), COLORS['source']),
        # Full Database -> Unique/Duplicate Pre
        (4, 5, unique_sequences_pre, COLORS['unique']),
        (4, 6, duplicate_sequences_pre, COLORS['duplicate']),
        # Unique Pre -> ITS Types
        (5, 7, ITSX_ITS1, COLORS['its1']),
        (5, 8, ITSX_ITS2, COLORS['its2']),
        (5, 9, ITSX_ITS, COLORS['its']),
        (5, 14, unique_sequences_pre, COLORS['unique']),
        # ITS Types -> ITSx Results
        (7, 10, ITSX_ITS1_ITS1, COLORS['its1']), (7, 11, ITSX_ITS1_ITS2, COLORS['its2']), (7, 12, ITSX_ITS1_ITS, COLORS['its']),
        (8, 10, ITSX_ITS2_ITS1, COLORS['its1']), (8, 11, ITSX_ITS2_ITS2, COLORS['its2']), (8, 12, ITSX_ITS2_ITS, COLORS['its']),
        (9, 10, ITSX_ITS_ITS1, COLORS['its1']), (9, 11, ITSX_ITS_ITS2, COLORS['its2']), (9, 12, ITSX_ITS_ITS, COLORS['its']),
        # No Detection
        (7, 13, ITSX_ITS1_no_detection, COLORS['no_detection']),
        (8, 13, ITSX_ITS2_no_detection, COLORS['no_detection']),
        (9, 13, ITSX_ITS_no_detection, COLORS['no_detection']),
        # ITSx Results -> Unique/Duplicate Post
        (10, 14, ITSX_ITS1_new - max(0, diff_post_pre_ITS1), COLORS['unique']),
        (11, 14, ITSX_ITS2_new - max(0, diff_post_pre_ITS2), COLORS['unique']),
        (12, 14, ITSX_ITS_new - max(0, diff_post_pre_ITS), COLORS['unique']),
        (10, 15, max(0, diff_post_pre_ITS1), COLORS['duplicate']),
        (11, 15, max(0, diff_post_pre_ITS2), COLORS['duplicate']),
        (12, 15, max(0, diff_post_pre_ITS), COLORS['duplicate']),
        # Unique Post -> MSA
        (14, 16, MSA_pre_ITS1, COLORS['its1']),
        (14, 17, MSA_pre_ITS2, COLORS['its2']),
        (14, 18, MSA_pre_ITS, COLORS['its']),
        # MSA prep -> dropped (terminal)
        (16, 19, MSA_dropped_ITS1, COLORS['no_detection']),
        (17, 20, MSA_dropped_ITS2, COLORS['no_detection']),
        (18, 21, MSA_dropped_ITS, COLORS['no_detection']),
        # MSA prep -> reference MSAs
        (16, 22, MSA_post_ITS1, COLORS['its1']),
        (17, 23, MSA_post_ITS2, COLORS['its2']),
        (18, 24, MSA_post_ITS, COLORS['its']),
    ]
    links_data += [(_idx, 4, sequences_by_db[_name], COLORS['source'])
                   for _name, _idx in _extra_source_idx.items()]

    # ========================================
    # FIGURE RENDERING
    # ========================================
    fig = go.Figure()

    fig.update_layout(
        title={
            'text': "Genomic Sequence Analysis - Processing Flow",
            'x': 0.5, 'xanchor': 'center',
            'font': {'size': 24, 'color': '#2c3e50', 'family': 'Arial, sans-serif'}
        },
        xaxis=dict(visible=False, range=[-0.05, 1.05]),
        # Bottom padding: the lowest nodes sit at y=(NB_ROW-4)/NB_ROW and carry a
        # three-line label drawn below the marker, which needs room under them.
        yaxis=dict(visible=False, range=[-0.08, 1.1]),
        paper_bgcolor='rgba(248,249,250,1)',
        plot_bgcolor='rgba(248,249,250,1)',
        margin=dict(t=100, b=50, l=50, r=50),
        width=1400,
        height=800,
        hovermode='closest'
    )

    # Scale ribbons using normalized Y units (more robust than pixels)
    max_link_val = max(max(1, v) for _, _, v, _ in links_data)

    # Maximum width as a fraction of the Y-axis height (0-1.1)
    MAX_WIDTH_Y = 0.08  # 8% de la hauteur totale pour le ruban le plus large

    # Add the ribbons (links)
    for from_idx, to_idx, value, color in links_data:
        v = max(0, value)
        if v <= 0:
            continue

        # Width expressed with normalized Y coordinates (no pixel conversion)
        width_y = (v / max_link_val) * MAX_WIDTH_Y

        x0, y0 = nodes_data[from_idx][1], nodes_data[from_idx][2]
        x1, y1 = nodes_data[to_idx][1], nodes_data[to_idx][2]

        ribbon_x, ribbon_y = generate_ribbon_sankey(
            x0, y0, x1, y1, width_y, num_points=120, curvature=0.35
        )

        fig.add_trace(go.Scatter(
            x=ribbon_x, y=ribbon_y,
            fill="toself",
            fillcolor=hex_to_rgba(color, LINK_OPACITY),
            line_color='rgba(0,0,0,0)',
            mode='lines',
            hoverinfo='none',
            showlegend=False
        ))

    # Log to make sure traces were added as expected
    logger.info(f"Traces added to figure: {len(fig.data)} (links + nodes + rectangles)")

    # Add the nodes
    max_node_val = max(max(1, n[3]) for n in nodes_data)
    MIN_SIZE, MAX_SIZE = 9, 60

    for name, x, y, value, color, _ in nodes_data:
        s = MIN_SIZE + (np.sqrt(max(0, value)) / np.sqrt(max_node_val)) * (MAX_SIZE - MIN_SIZE)

        fig.add_trace(go.Scatter(
            x=[x], y=[y],
            mode='markers+text',
            marker=dict(size=s, color=color, line=dict(width=2, color='white')),
            text=name,
            textposition='bottom center',
            textfont=dict(color='black', size=12, family='Arial, sans-serif'),
            hovertemplate=f"<b>{name.replace('<br>', ' ')}</b><br>Sequences: {value:,}<extra></extra>",
            showlegend=False
        ))

    # Add the tool frames
    add_tool_rectangle(fig, 0.9, 1.1, 2, 3, "CoL",
                       color="rgba(52, 152, 219, 0.08)",
                       border_color="rgba(52, 152, 219, 0.6)")
    add_tool_rectangle(fig, 3, 4, 0, 3, "ITSx",
                       color="rgba(39, 174, 96, 0.08)",
                       border_color="rgba(39, 174, 96, 0.6)")
    add_tool_rectangle(fig, 6, 8, 0.5, 3.3, "MSA",
                       color="rgba(231, 76, 60, 0.08)",
                       border_color="rgba(231, 76, 60, 0.6)")

    # ========================================
    # FULL HTML REPORT GENERATION
    # ========================================
    logger.info(f"\n💾 Generating comprehensive HTML report...")

    # Compute key statistics
    total_final = MSA_post_ITS1 + MSA_post_ITS2 + MSA_post_ITS
    retention_rate = (total_final / total_sequence_ini * 100) if total_sequence_ini > 0 else 0
    itsx_detection_rate = ((ITSX_ITS1_new + ITSX_ITS2_new + ITSX_ITS_new + ITSX_ITS_no_detection_tot) /
                          (ITSX_ITS1 + ITSX_ITS2 + ITSX_ITS) * 100) if (ITSX_ITS1 + ITSX_ITS2 + ITSX_ITS) > 0 else 0
    msa_processing_rate = ((MSA_post_ITS1 + MSA_post_ITS2 + MSA_post_ITS) /
                          (MSA_pre_ITS1 + MSA_pre_ITS2 + MSA_pre_ITS) * 100) if (MSA_pre_ITS1 + MSA_pre_ITS2 + MSA_pre_ITS) > 0 else 0

    # Save the Sankey diagram as a standalone HTML file
    output_path = Path(output_html)
    sankey_standalone_path = output_path.parent / "sequence_flow_sankey.html"

    logger.info(f"Saving standalone Sankey diagram to: {sankey_standalone_path}")
    fig.write_html(
        str(sankey_standalone_path),
        include_plotlyjs='cdn',
        config={'responsive': True}
    )

    # Use an iframe for embedding
    # Embed the figure in the report itself rather than pointing an iframe at the
    # standalone file: a relative src makes the report show whichever
    # sequence_flow_sankey.html sits next to it, so a report moved or shared on its own
    # silently displays a stale diagram beside up-to-date numbers.
    sankey_iframe = fig.to_html(full_html=False, include_plotlyjs='cdn',
                                default_width='100%', default_height='800px',
                                config={'responsive': True})

    # Build the full HTML layout with detailed sections
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sequence Flow Report - CurateMake Pipeline</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>{REPORT_CSS}</style>
</head>
<body>
    <div class="header">
        <h1>🌊 Sequence Flow Analysis Report</h1>
        <p>Quantitative tracking of sequences through the CurateMake pipeline</p>
        <p style="font-size: 0.9em; margin-top: 10px;">Generated on {datetime.now().strftime('%B %d, %Y at %H:%M:%S')}</p>
    </div>

    <div class="container">
        <!-- Executive Summary -->
        <div class="section">
            <h2>📊 Executive Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-label">Initial Sequences</div>
                    <div class="stat-value">{total_sequence_ini:,}</div>
                    <div class="stat-label">From 4 databases</div>
                </div>
                <div class="stat-card">
                    <div class="stat-label">Final Sequences</div>
                    <div class="stat-value">{total_final:,}</div>
                    <div class="stat-label">After full pipeline</div>
                </div>
                <div class="stat-card">
                    <div class="stat-label">Retention Rate</div>
                    <div class="stat-value">{retention_rate:.1f}%</div>
                    <div class="stat-label">Overall efficiency</div>
                </div>
                <div class="stat-card">
                    <div class="stat-label">ITSx Detection</div>
                    <div class="stat-value">{itsx_detection_rate:.1f}%</div>
                    <div class="stat-label">Successful extractions</div>
                </div>
            </div>

            <div class="info-box success">
                <h4>✓ Pipeline Overview</h4>
                <p>This report tracks the quantitative flow of ITS sequences through all major processing stages of the CurateMake pipeline:
                database merging, taxonomic harmonization, deduplication, ITSx subregion extraction, and MSA-based quality control.</p>
            </div>
        </div>

        <!-- Interactive Sankey Diagram -->
        <div class="section">
            <h2>🔄 Interactive Sequence Flow Diagram</h2>
            <p>The Sankey diagram below visualizes the quantitative flow of sequences through each pipeline stage.
            Hover over nodes and flows to see detailed numbers. Node sizes are proportional to sequence counts.</p>
            {sankey_iframe}
        </div>

        <!-- Stage 1: Database Sources & CoL -->
        <div class="section">
            <h2>1️⃣ Database Integration & Taxonomic Harmonization</h2>

            <div class="tool-description">
                <h4>🗂️ Source Databases</h4>
                <p>Sequences are aggregated from four major ITS barcode databases:</p>
                <ul>
                    <li><strong>UNITE</strong>: {sequences_by_db.get('unite', 0):,} sequences - Curated ITS database</li>
                    <li><strong>BOLD</strong>: {sequences_by_db.get('bold', 0):,} sequences - Barcode of Life Data System</li>
                    <li><strong>PLANiTS</strong>: {sequences_by_db.get('planits', 0):,} sequences - Plant ITS reference library</li>
                    <li><strong>CALeDNA</strong>: {sequences_by_db.get('CALeDNA', 0):,} sequences - California eDNA database</li>
                </ul>
            </div>

            <div class="tool-description">
                <h4>🌐 Catalogue of Life (CoL) Integration</h4>
                <p>The <strong>Catalogue of Life API</strong> is used to harmonize taxonomic nomenclature across all databases.
                This critical step ensures:</p>
                <ul>
                    <li>Standardized taxonomic ranks (Kingdom → Species)</li>
                    <li>Resolution of synonyms and outdated names</li>
                    <li>Consistent taxonomic hierarchy for downstream analyses</li>
                </ul>
            </div>

            <table>
                <tr>
                    <th>Metric</th>
                    <th>Value</th>
                </tr>
                <tr>
                    <td>Total sequences merged</td>
                    <td class="number">{total_sequence_ini:,}</td>
                </tr>
                <tr>
                    <td>Databases integrated</td>
                    <td class="number">{n_databases_integrated}</td>
                </tr>
            </table>
        </div>

        <!-- Stage 2: Pre-ITSx Deduplication -->
        <div class="section">
            <h2>2️⃣ Pre-ITSx Deduplication</h2>

            <p>Identical sequences (100% identity) are identified and marked to reduce computational load and improve downstream analysis quality.</p>

            <div class="stage-details">
                <div class="stage-card">
                    <h5>Unique Sequences</h5>
                    <p class="number">{unique_sequences_pre:,}</p>
                    <p>Retained for ITSx processing</p>
                </div>
                <div class="stage-card">
                    <h5>Duplicates Marked</h5>
                    <p class="number">{duplicate_sequences_pre:,}</p>
                    <p>Exact matches removed</p>
                </div>
            </div>

            <div class="info-box">
                <h4>ℹ️ Note on Mini Datasets</h4>
                <p>In small test datasets, pre-ITSx deduplication may show zero duplicates due to limited sequence overlap between databases.</p>
            </div>
        </div>

        <!-- Stage 3: ITSx Subregion Extraction -->
        <div class="section">
            <h2>3️⃣ ITSx Subregion Extraction</h2>

            <div class="tool-description">
                <h4>🔬 ITSx: ITS Extractor</h4>
                <p><strong>ITSx</strong> is a specialized tool for extracting and identifying ITS subregions (ITS1, 5.8S, ITS2)
                and flanking conserved regions (SSU, LSU) from ITS sequences.</p>

                <p><strong>Methodology:</strong></p>
                <ul>
                    <li>Uses Hidden Markov Model (HMM) profiles to detect conserved regions (SSU, 5.8S, LSU)</li>
                    <li>Extracts variable regions (ITS1, ITS2) based on flanking conserved sequences</li>
                    <li>Applied to all three barcode types: ITS (full), ITS1, and ITS2</li>
                    <li>Can reconstruct full ITS from detected subregions</li>
                </ul>

                <p><strong>Key Output:</strong> Sequence augmentation - one input ITS sequence can generate ITS1, ITS2, and reconstructed ITS</p>
            </div>

            <div class="info-box warning">
                <h4>⚠️ Important Considerations</h4>
                <ul>
                    <li><strong>HMM Detection Limitation:</strong> ITSx relies on detecting conserved regions (LSU, 5.8S, SSU).
                    Sequences already trimmed to ITS1 or ITS2 only may not be detected (e.g., PLANiTS database sequences).</li>
                    <li><strong>Original Sequences Retained:</strong> All original sequences are kept to preserve pre-trimmed ITS1/ITS2 data.</li>
                    <li><strong>Cross-Validation:</strong> ITS1 and ITS2 input sequences are also checked with ITSx to identify potential misannotations.</li>
                </ul>
            </div>

            <h3>Input Distribution</h3>
            <table>
                <tr>
                    <th>Barcode</th>
                    <th>Input Sequences</th>
                    <th>ITSx Detected</th>
                    <th>No Detection</th>
                    <th>Detection Rate</th>
                </tr>
                <tr>
                    <td>ITS1</td>
                    <td class="number">{ITSX_ITS1:,}</td>
                    <td class="number">{ITSX_ITS1_ITS1 + ITSX_ITS1_ITS2 + ITSX_ITS1_ITS:,}</td>
                    <td class="number">{ITSX_ITS1_no_detection:,}</td>
                    <td class="number">{((ITSX_ITS1_ITS1 + ITSX_ITS1_ITS2 + ITSX_ITS1_ITS) / ITSX_ITS1 * 100 if ITSX_ITS1 > 0 else 0):.1f}%</td>
                </tr>
                <tr>
                    <td>ITS2</td>
                    <td class="number">{ITSX_ITS2:,}</td>
                    <td class="number">{ITSX_ITS2_ITS1 + ITSX_ITS2_ITS2 + ITSX_ITS2_ITS:,}</td>
                    <td class="number">{ITSX_ITS2_no_detection:,}</td>
                    <td class="number">{((ITSX_ITS2_ITS1 + ITSX_ITS2_ITS2 + ITSX_ITS2_ITS) / ITSX_ITS2 * 100 if ITSX_ITS2 > 0 else 0):.1f}%</td>
                </tr>
                <tr>
                    <td>ITS (full)</td>
                    <td class="number">{ITSX_ITS:,}</td>
                    <td class="number">{ITSX_ITS_ITS1 + ITSX_ITS_ITS2 + ITSX_ITS_ITS:,}</td>
                    <td class="number">{ITSX_ITS_no_detection:,}</td>
                    <td class="number">{((ITSX_ITS_ITS1 + ITSX_ITS_ITS2 + ITSX_ITS_ITS) / ITSX_ITS * 100 if ITSX_ITS > 0 else 0):.1f}%</td>
                </tr>
            </table>

            <h3>Output Sequences (Augmented)</h3>
            <div class="stage-details">
                <div class="stage-card">
                    <h5>ITS1 Extracted</h5>
                    <p class="number">{ITSX_ITS1_new:,}</p>
                    <p>From all sources</p>
                </div>
                <div class="stage-card">
                    <h5>ITS2 Extracted</h5>
                    <p class="number">{ITSX_ITS2_new:,}</p>
                    <p>From all sources</p>
                </div>
                <div class="stage-card">
                    <h5>ITS Reconstructed</h5>
                    <p class="number">{ITSX_ITS_new:,}</p>
                    <p>From common parents</p>
                </div>
                <div class="stage-card">
                    <h5>No Detection</h5>
                    <p class="number">{ITSX_ITS_no_detection_tot:,}</p>
                    <p>Pre-trimmed sequences</p>
                </div>
            </div>
        </div>

        <!-- Stage 4: Post-ITSx Deduplication -->
        <div class="section">
            <h2>4️⃣ Post-ITSx Deduplication</h2>

            <p>A second deduplication pass is performed after ITSx extraction, as sequence augmentation can introduce new duplicates
            (e.g., identical ITS1 regions extracted from different full ITS sequences).</p>

            <div class="stage-details">
                <div class="stage-card">
                    <h5>Unique Sequences</h5>
                    <p class="number">{unique_sequences_post:,}</p>
                    <p>Retained for MSA</p>
                </div>
                <div class="stage-card">
                    <h5>New Duplicates</h5>
                    <p class="number">{new_dup_total}</p>
                    <p>Marked by the post-ITSx pass</p>
                </div>
            </div>

            <table>
                <tr>
                    <th>Barcode</th>
                    <th>Pre-ITSx Duplicates</th>
                    <th>Post-ITSx Duplicates</th>
                    <th>New Duplicates</th>
                </tr>
                <tr>
                    <td>ITS1</td>
                    <td class="number">{pre_dup_ITS1}</td>
                    <td class="number">{duplicate_sequences_post_ITS1:,}</td>
                    <td class="number">{new_dup_ITS1}</td>
                </tr>
                <tr>
                    <td>ITS2</td>
                    <td class="number">{pre_dup_ITS2}</td>
                    <td class="number">{duplicate_sequences_post_ITS2:,}</td>
                    <td class="number">{new_dup_ITS2}</td>
                </tr>
                <tr>
                    <td>ITS</td>
                    <td class="number">{pre_dup_ITS}</td>
                    <td class="number">{duplicate_sequences_post_ITS:,}</td>
                    <td class="number">{new_dup_ITS}</td>
                </tr>
            </table>
        </div>

        <!-- Stage 5: MSA Processing -->
        <div class="section">
            <h2>5️⃣ Multiple Sequence Alignment (MSA) & Quality Control</h2>

            <div class="tool-description">
                <h4>🧬 MSA Pipeline: Critical Quality Control Step</h4>
                <p>The MSA stage is a <strong>crucial quality control and validation step</strong> that uses multiple bioinformatics tools
                to build taxonomically coherent alignments and identify misannotated sequences.</p>

                <p><strong>Tools Used:</strong></p>
                <ul>
                    <li><strong>mmseqs2</strong>: sequence clustering — selects the backbone, and recovers reads no profile claimed</li>
                    <li><strong>MAFFT (mafft-linsi)</strong>: high-accuracy alignment, of the backbone only</li>
                    <li><strong>HMMER (hmmbuild / nhmmer / hmmalign)</strong>: builds the profile ensemble, assigns each read to a profile, then projects it</li>
                </ul>

                <p><strong>Backbone + ensemble-HMM construction (UPP-style):</strong></p>
                <ul>
                    <li>A backbone of ~1200 representative sequences is selected per barcode by clustering (mmseqs2, 90% identity) and aligned once with mafft-linsi — the only from-scratch alignment in the stage</li>
                    <li>The backbone is decomposed into an ensemble of profile HMMs (subsets of 10, bucketed by kingdom, minimum 20 sequences per bucket)</li>
                    <li>Every remaining read is assigned to its best-matching profile with nhmmer, then projected into the alignment with hmmalign; reads no profile claims are recovered by sequence identity (mmseqs2, ≥75% identity)</li>
                    <li>Because every read lands in the <strong>same column frame</strong>, groups delivered to SATIVA are a free slice of the alignment — no per-group MAFFT. The frame is fixed per barcode: 940 columns for ITS, 461 for ITS1, 476 for ITS2.</li>
                </ul>

                <p><strong>Why the reference frame spans kingdoms:</strong></p>
                <ul>
                    <li>Groups are packed kingdom-interleaved up to a size cap, so every alignment contains a bit of everything</li>
                    <li>That outgroup diversity is what lets SATIVA move a contaminant <strong>across</strong> kingdoms — a plant read mislabelled as a fungus in the source database can only be caught if plants and fungi share one column frame</li>
                    <li>Alternatives preserving per-clade columns were benchmarked and rejected: they collapse cross-kingdom detection</li>
                </ul>
            </div>

            <div class="info-box success">
                <h4>✓ MSA Quality Control Benefits</h4>
                <ul>
                    <li>Identifies sequences with taxonomic annotation errors</li>
                    <li>Reveals chimeric or contaminated sequences</li>
                    <li>Ensures high-quality input for phylogenetic analysis (SATIVA)</li>
                    <li>Improves downstream taxonomic assignment accuracy</li>
                </ul>
            </div>

            <h3>MSA Processing Statistics</h3>
            <p>Each barcode is clustered independently; a sequence either reaches a reference
               MSA or is dropped before it. Sequences drop out when their cluster stays below
               <code>min_ref</code> and cannot be absorbed into a neighbour, when the merged
               group stays below <code>sativa_min_seqs</code>, or when their cluster has
               already contributed <code>max_cluster_reps</code> representatives.</p>
            <table>
                <tr>
                    <th>Barcode</th>
                    <th>Extracted</th>
                    <th>In reference MSAs</th>
                    <th>Not aligned</th>
                    <th>Drop rate</th>
                </tr>
                <tr>
                    <td>ITS1</td>
                    <td class="number">{MSA_pre_ITS1:,}</td>
                    <td class="number">{MSA_post_ITS1:,}</td>
                    <td class="number">{MSA_dropped_ITS1:,}</td>
                    <td class="number">{msa_drop_rate_ITS1:.1f}%</td>
                </tr>
                <tr>
                    <td>ITS2</td>
                    <td class="number">{MSA_pre_ITS2:,}</td>
                    <td class="number">{MSA_post_ITS2:,}</td>
                    <td class="number">{MSA_dropped_ITS2:,}</td>
                    <td class="number">{msa_drop_rate_ITS2:.1f}%</td>
                </tr>
                <tr>
                    <td>ITS</td>
                    <td class="number">{MSA_pre_ITS:,}</td>
                    <td class="number">{MSA_post_ITS:,}</td>
                    <td class="number">{MSA_dropped_ITS:,}</td>
                    <td class="number">{msa_drop_rate_ITS:.1f}%</td>
                </tr>
            </table>
        </div>

        <!-- Final Summary -->
        <div class="section">
            <h2>📈 Pipeline Performance Summary</h2>

            <table>
                <tr>
                    <th>Stage</th>
                    <th>Sequences In</th>
                    <th>Sequences Out</th>
                    <th>Retention</th>
                </tr>
                <tr>
                    <td>Database Merge</td>
                    <td class="number">{total_sequence_ini:,}</td>
                    <td class="number">{total_sequence_ini:,}</td>
                    <td class="number">100.0%</td>
                </tr>
                <tr>
                    <td>Pre-ITSx Deduplication</td>
                    <td class="number">{total_sequence_ini:,}</td>
                    <td class="number">{unique_sequences_pre:,}</td>
                    <td class="number">{(unique_sequences_pre / total_sequence_ini * 100 if total_sequence_ini > 0 else 0):.1f}%</td>
                </tr>
                <tr>
                    <td>ITSx Extraction</td>
                    <td class="number">{ITSX_ITS1 + ITSX_ITS2 + ITSX_ITS:,}</td>
                    <td class="number">{ITSX_ITS1_new + ITSX_ITS2_new + ITSX_ITS_new + ITSX_ITS_no_detection_tot:,}</td>
                    <td class="number">{itsx_detection_rate:.1f}%</td>
                </tr>
                <tr>
                    <td>Post-ITSx Deduplication</td>
                    <td class="number">{ITSX_ITS1_new + ITSX_ITS2_new + ITSX_ITS_new + ITSX_ITS_no_detection_tot:,}</td>
                    <td class="number">{unique_sequences_post:,}</td>
                    <td class="number">{(unique_sequences_post / (ITSX_ITS1_new + ITSX_ITS2_new + ITSX_ITS_new + ITSX_ITS_no_detection_tot) * 100 if (ITSX_ITS1_new + ITSX_ITS2_new + ITSX_ITS_new) > 0 else 0):.1f}%</td>
                </tr>
                <tr>
                    <td>MSA Quality Control</td>
                    <td class="number">{MSA_pre_ITS1 + MSA_pre_ITS2 + MSA_pre_ITS:,}</td>
                    <td class="number">{total_final:,}</td>
                    <td class="number">{msa_processing_rate:.1f}%</td>
                </tr>
                <tr style="background-color: #e8f5e9; font-weight: bold;">
                    <td>Overall Pipeline</td>
                    <td class="number">{total_sequence_ini:,}</td>
                    <td class="number">{total_final:,}</td>
                    <td class="number">{retention_rate:.1f}%</td>
                </tr>
            </table>

            <div class="info-box">
                <h4>ℹ️ Understanding Retention Rates</h4>
                <p>The apparent sequence "loss" through the pipeline actually represents:</p>
                <ul>
                    <li><strong>Deduplication:</strong> Removal of redundant identical sequences</li>
                    <li><strong>Quality filtering:</strong> Removal of low-quality alignments and outliers</li>
                    <li><strong>Augmentation effects:</strong> ITSx can both increase (extraction) and consolidate (deduplication) sequences</li>
                </ul>
                <p>The final curated dataset prioritizes <strong>quality over quantity</strong> for downstream phylogenetic analyses.</p>
            </div>
        </div>
    </div>

    <div class="footer">
        <p><strong>CurateMake Pipeline</strong> - ITS Sequence Curation and Quality Control</p>
        <p>Report generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}</p>
        <p>For more information, see pipeline documentation</p>
    </div>
</body>
</html>"""

    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)

    logger.info(f"✅ Comprehensive HTML report saved to: {output_html}")
    logger.info("=" * 60)
    logger.info("✅ Sequence flow report generated successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
