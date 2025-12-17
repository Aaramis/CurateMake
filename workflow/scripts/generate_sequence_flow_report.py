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
            # Extract the base name (unite, bold, planits)
            base_name = db_name.lower()
            for known in ['unite', 'bold', 'planits']:
                if known in base_name:
                    sequences_by_db[known] = processed_seqs
                    break

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
            "msa_run": snakemake.input.msa_run,
        }
        output_html = snakemake.output.report
    except NameError:
        logger.error("This script must be run through Snakemake")
        sys.exit(1)

    # ========================================
    # DATA LOADING
    # ========================================
    logger.info("\n📥 Loading JSON reports...")
    reports = {key: load_json_safe(path) for key, path in report_paths.items()}

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
    logger.info(f"  Total initial: {total_sequence_ini:,}")

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

    # Deduplication deltas
    diff_post_pre_ITS1 = duplicate_sequences_post_ITS1 - duplicate_sequences_pre_ITS1
    diff_post_pre_ITS2 = duplicate_sequences_post_ITS2 - duplicate_sequences_pre_ITS2
    diff_post_pre_ITS = duplicate_sequences_post_ITS - duplicate_sequences_pre_ITS
    duplicate_sequences_post_tot = duplicate_sequences_post - duplicate_sequences_pre

    logger.info(f"  Post-ITSx dedup: {unique_sequences_post:,} unique, "
                f"{duplicate_sequences_post:,} duplicates (new: {duplicate_sequences_post_tot:,})")

    # --- MSA preparation ---
    MSA_unique_sequences = safe_get(reports["msa_prepare"], ['stats', 'total_sequences'], 0)
    MSA_pre_ITS = safe_get(reports["msa_prepare"], ['stats', 'barcode_distribution', 'ITS'], 0)
    MSA_pre_ITS1 = safe_get(reports["msa_prepare"], ['stats', 'barcode_distribution', 'ITS1'], 0)
    MSA_pre_ITS2 = safe_get(reports["msa_prepare"], ['stats', 'barcode_distribution', 'ITS2'], 0)

    logger.info(f"  MSA prep: Total={MSA_unique_sequences:,}, ITS1={MSA_pre_ITS1}, "
                f"ITS2={MSA_pre_ITS2}, ITS={MSA_pre_ITS}")

    # --- Traitement MSA ---
    msa_barcode_results = safe_get(reports["msa_run"], ['results', 'barcode_results'], {})

    MSA_unprocessed_ITS1 = safe_get(msa_barcode_results, ['ITS1', 'sequences_unprocessed'], 0)
    MSA_unprocessed_ITS2 = safe_get(msa_barcode_results, ['ITS2', 'sequences_unprocessed'], 0)
    MSA_unprocessed_ITS = safe_get(msa_barcode_results, ['ITS', 'sequences_unprocessed'], 0)

    MSA_treat_ITS1 = safe_get(msa_barcode_results, ['ITS1', 'sequences_processed'], 0) - MSA_unprocessed_ITS1
    MSA_treat_ITS2 = safe_get(msa_barcode_results, ['ITS2', 'sequences_processed'], 0) - MSA_unprocessed_ITS2
    MSA_treat_ITS = safe_get(msa_barcode_results, ['ITS', 'sequences_processed'], 0) - MSA_unprocessed_ITS

    # MSA reassignments
    msa_reassign = safe_get(reports["msa_run"], ['results', 'global_assignments', 'reassignment_stats'], {})

    MSA_ITS1_ITS1 = safe_get(msa_reassign, ['ITS1_to_ITS1', 'count'], 0)
    MSA_ITS2_ITS2 = safe_get(msa_reassign, ['ITS2_to_ITS2', 'count'], 0)
    MSA_ITS_ITS = safe_get(msa_reassign, ['ITS_to_ITS', 'count'], 0)

    MSA_ITS1_ITS2 = safe_get(msa_reassign, ['ITS1_to_ITS2', 'count'], 0)
    MSA_ITS1_ITS = safe_get(msa_reassign, ['ITS1_to_ITS', 'count'], 0)
    MSA_ITS2_ITS1 = safe_get(msa_reassign, ['ITS2_to_ITS1', 'count'], 0)
    MSA_ITS2_ITS = safe_get(msa_reassign, ['ITS2_to_ITS', 'count'], 0)
    MSA_ITS_ITS1 = safe_get(msa_reassign, ['ITS_to_ITS1', 'count'], 0)
    MSA_ITS_ITS2 = safe_get(msa_reassign, ['ITS_to_ITS2', 'count'], 0)

    MSA_post_ITS1 = MSA_ITS1_ITS1 + MSA_ITS2_ITS1 + MSA_ITS_ITS1 + MSA_treat_ITS1
    MSA_post_ITS2 = MSA_ITS2_ITS2 + MSA_ITS1_ITS2 + MSA_ITS_ITS2 + MSA_treat_ITS2
    MSA_post_ITS = MSA_ITS_ITS + MSA_ITS1_ITS + MSA_ITS2_ITS + MSA_treat_ITS

    logger.info(f"  MSA final: ITS1={MSA_post_ITS1}, ITS2={MSA_post_ITS2}, ITS={MSA_post_ITS}")

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
        # Colonne 7 - MSA untreated
        (f"Untreated ITS1<br>{MSA_unprocessed_ITS1:,}", 7/NB_COL, (NB_ROW-1.25)/NB_ROW, MSA_unprocessed_ITS1, COLORS['its1'], 50),
        (f"Untreated ITS2<br>{MSA_unprocessed_ITS2:,}", 7/NB_COL, (NB_ROW-2.25)/NB_ROW, MSA_unprocessed_ITS2, COLORS['its2'], 50),
        (f"Untreated ITS<br>{MSA_unprocessed_ITS:,}", 7/NB_COL, (NB_ROW-3.25)/NB_ROW, MSA_unprocessed_ITS, COLORS['its'], 50),
        # Colonne 8 - MSA final
        (f"ITS1<br>{MSA_post_ITS1:,}", 8/NB_COL, (NB_ROW-0.75)/NB_ROW, MSA_post_ITS1, COLORS['its1'], 50),
        (f"ITS2<br>{MSA_post_ITS2:,}", 8/NB_COL, (NB_ROW-1.75)/NB_ROW, MSA_post_ITS2, COLORS['its2'], 50),
        (f"ITS<br>{MSA_post_ITS:,}", 8/NB_COL, (NB_ROW-2.75)/NB_ROW, MSA_post_ITS, COLORS['its'], 50),
    ]

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
        # MSA -> Untreated
        (16, 19, MSA_unprocessed_ITS1, COLORS['its1']),
        (17, 20, MSA_unprocessed_ITS2, COLORS['its2']),
        (18, 21, MSA_unprocessed_ITS, COLORS['its']),
        # MSA -> Treated
        (16, 22, MSA_treat_ITS1, COLORS['its1']),
        (17, 23, MSA_treat_ITS2, COLORS['its2']),
        (18, 24, MSA_treat_ITS, COLORS['its']),
        # Untreated -> MSA_post (re-classification)
        (19, 22, MSA_ITS1_ITS1, COLORS['its1']),
        (19, 23, MSA_ITS1_ITS2, COLORS['its2']),
        (19, 24, MSA_ITS1_ITS, COLORS['its']),
        (20, 22, MSA_ITS2_ITS1, COLORS['its1']),
        (20, 23, MSA_ITS2_ITS2, COLORS['its2']),
        (20, 24, MSA_ITS2_ITS, COLORS['its']),
        (21, 22, MSA_ITS_ITS1, COLORS['its1']),
        (21, 23, MSA_ITS_ITS2, COLORS['its2']),
        (21, 24, MSA_ITS_ITS, COLORS['its']),
    ]

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
        yaxis=dict(visible=False, range=[0, 1.1]),
        paper_bgcolor='rgba(248,249,250,1)',
        plot_bgcolor='rgba(248,249,250,1)',
        margin=dict(t=100, b=50, l=50, r=50),
        width=1400,
        height=800,
        hovermode='closest'
    )

    # Scale ribbons using normalized Y units (more robust than pixels)
    max_link_val = max(max(1, v) for _, _, v, _ in links_data)

    # Largeur maximale en proportion de la hauteur de l'axe Y (0-1.1)
    MAX_WIDTH_Y = 0.08  # 8% de la hauteur totale pour le ruban le plus large

    # Ajouter les rubans (liens)
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

    # Ajouter les nœuds
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

    # Ajouter les cadres d'outils
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

    # Sauvegarder le graphique Sankey comme fichier HTML standalone
    output_path = Path(output_html)
    sankey_standalone_path = output_path.parent / "sequence_flow_sankey.html"

    logger.info(f"Saving standalone Sankey diagram to: {sankey_standalone_path}")
    fig.write_html(
        str(sankey_standalone_path),
        include_plotlyjs='cdn',
        config={'responsive': True}
    )

    # Use an iframe for embedding
    sankey_iframe = f'<iframe src="sequence_flow_sankey.html" style="width: 100%; height: 850px; border: none;" title="Sequence Flow Sankey Diagram"></iframe>'

    # Build the full HTML layout with detailed sections
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sequence Flow Report - CurateMake Pipeline</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #2c3e50;
            background-color: #f8f9fa;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 60px 20px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .header p {{
            font-size: 1.2em;
            opacity: 0.95;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 40px 20px;
        }}

        .section {{
            background: white;
            border-radius: 10px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        }}

        .section h2 {{
            color: #667eea;
            font-size: 1.8em;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }}

        .section h3 {{
            color: #764ba2;
            font-size: 1.4em;
            margin: 25px 0 15px 0;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}

        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .stat-value {{
            font-size: 2.5em;
            font-weight: bold;
            margin: 10px 0;
        }}

        .stat-label {{
            font-size: 0.9em;
            opacity: 0.9;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}

        .info-box {{
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 5px;
        }}

        .info-box.warning {{
            background: #fff3e0;
            border-left-color: #ff9800;
        }}

        .info-box.success {{
            background: #e8f5e9;
            border-left-color: #4caf50;
        }}

        .info-box h4 {{
            color: #1976d2;
            margin-bottom: 10px;
            font-size: 1.1em;
        }}

        .info-box.warning h4 {{
            color: #f57c00;
        }}

        .info-box.success h4 {{
            color: #388e3c;
        }}

        .tool-description {{
            background: #f5f5f5;
            padding: 20px;
            border-radius: 8px;
            margin: 15px 0;
        }}

        .tool-description h4 {{
            color: #667eea;
            margin-bottom: 10px;
            font-size: 1.2em;
        }}

        .tool-description ul {{
            margin-left: 20px;
            margin-top: 10px;
        }}

        .tool-description li {{
            margin: 8px 0;
        }}

        .stage-details {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin: 20px 0;
        }}

        .stage-card {{
            background: #fafafa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}

        .stage-card h5 {{
            color: #667eea;
            margin-bottom: 10px;
        }}

        .number {{
            font-weight: bold;
            color: #764ba2;
            font-family: 'Courier New', monospace;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}

        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}

        th {{
            background-color: #667eea;
            color: white;
            font-weight: 600;
        }}

        tr:hover {{
            background-color: #f5f5f5;
        }}

        .footer {{
            text-align: center;
            padding: 30px;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #ddd;
            margin-top: 40px;
        }}

        #sankey-diagram {{
            margin: 30px 0;
            min-height: 700px;
        }}

        @media (max-width: 768px) {{
            .stage-details {{
                grid-template-columns: 1fr;
            }}

            .header h1 {{
                font-size: 1.8em;
            }}
        }}
    </style>
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
                    <td class="number">4</td>
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
                    <p class="number">{duplicate_sequences_post_tot:,}</p>
                    <p>Introduced by ITSx augmentation</p>
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
                    <td class="number">{duplicate_sequences_pre_ITS1:,}</td>
                    <td class="number">{duplicate_sequences_post_ITS1:,}</td>
                    <td class="number">{diff_post_pre_ITS1:,}</td>
                </tr>
                <tr>
                    <td>ITS2</td>
                    <td class="number">{duplicate_sequences_pre_ITS2:,}</td>
                    <td class="number">{duplicate_sequences_post_ITS2:,}</td>
                    <td class="number">{diff_post_pre_ITS2:,}</td>
                </tr>
                <tr>
                    <td>ITS</td>
                    <td class="number">{duplicate_sequences_pre_ITS:,}</td>
                    <td class="number">{duplicate_sequences_post_ITS:,}</td>
                    <td class="number">{diff_post_pre_ITS:,}</td>
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
                    <li><strong>MAFFT</strong>: Fast multiple sequence alignment algorithm</li>
                    <li><strong>trimAl</strong>: Alignment trimming and column/sequence filtering</li>
                    <li><strong>HMMER (hmmscan)</strong>: Profile HMM searching for sequence classification</li>
                </ul>

                <p><strong>Hierarchical MSA Construction:</strong></p>
                <ul>
                    <li>Sequences are grouped by taxonomy (species → genus → family → order)</li>
                    <li>MSA built at the finest taxonomic level possible (minimum {safe_get(reports["msa_run"], ['parameters', 'min_seq_threshold'], 20)} sequences)</li>
                    <li>Alignments are quality-filtered using trimAl (overlap thresholds: seq={safe_get(reports["msa_run"], ['parameters', 'trimal_seq_overlap'], 75)}%, res={safe_get(reports["msa_run"], ['parameters', 'trimal_res_overlap'], 0.5)})</li>
                </ul>

                <p><strong>HMM-Based Sequence Classification:</strong></p>
                <ul>
                    <li>HMM profiles are built from each taxonomic MSA</li>
                    <li>Sequences removed during alignment trimming are searched against all HMM profiles</li>
                    <li>Sequences matching better to a different taxonomic group are flagged as potentially misannotated</li>
                    <li>This step is <strong>essential for optimal SATIVA performance</strong>, which relies on clean MSAs</li>
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
            <table>
                <tr>
                    <th>Barcode</th>
                    <th>Input</th>
                    <th>Processed</th>
                    <th>Unprocessed</th>
                    <th>Reassigned</th>
                    <th>Final Output</th>
                </tr>
                <tr>
                    <td>ITS1</td>
                    <td class="number">{MSA_pre_ITS1:,}</td>
                    <td class="number">{MSA_treat_ITS1:,}</td>
                    <td class="number">{MSA_unprocessed_ITS1:,}</td>
                    <td class="number">{MSA_ITS1_ITS1 + MSA_ITS2_ITS1 + MSA_ITS_ITS1:,}</td>
                    <td class="number">{MSA_post_ITS1:,}</td>
                </tr>
                <tr>
                    <td>ITS2</td>
                    <td class="number">{MSA_pre_ITS2:,}</td>
                    <td class="number">{MSA_treat_ITS2:,}</td>
                    <td class="number">{MSA_unprocessed_ITS2:,}</td>
                    <td class="number">{MSA_ITS1_ITS2 + MSA_ITS2_ITS2 + MSA_ITS_ITS2:,}</td>
                    <td class="number">{MSA_post_ITS2:,}</td>
                </tr>
                <tr>
                    <td>ITS</td>
                    <td class="number">{MSA_pre_ITS:,}</td>
                    <td class="number">{MSA_treat_ITS:,}</td>
                    <td class="number">{MSA_unprocessed_ITS:,}</td>
                    <td class="number">{MSA_ITS1_ITS + MSA_ITS2_ITS + MSA_ITS_ITS:,}</td>
                    <td class="number">{MSA_post_ITS:,}</td>
                </tr>
            </table>

            <h3>Cross-Barcode Reassignments</h3>
            <p>HMM classification can identify sequences originally annotated as one barcode type that better match another:</p>
            <table>
                <tr>
                    <th>Original → Reassigned</th>
                    <th>Count</th>
                    <th>Interpretation</th>
                </tr>
                <tr>
                    <td>ITS1 → ITS2</td>
                    <td class="number">{MSA_ITS1_ITS2:,}</td>
                    <td>Potential misannotation or extraction error</td>
                </tr>
                <tr>
                    <td>ITS1 → ITS</td>
                    <td class="number">{MSA_ITS1_ITS:,}</td>
                    <td>May contain flanking regions</td>
                </tr>
                <tr>
                    <td>ITS2 → ITS1</td>
                    <td class="number">{MSA_ITS2_ITS1:,}</td>
                    <td>Potential misannotation or extraction error</td>
                </tr>
                <tr>
                    <td>ITS2 → ITS</td>
                    <td class="number">{MSA_ITS2_ITS:,}</td>
                    <td>May contain flanking regions</td>
                </tr>
                <tr>
                    <td>ITS → ITS1</td>
                    <td class="number">{MSA_ITS_ITS1:,}</td>
                    <td>Partial ITS sequence</td>
                </tr>
                <tr>
                    <td>ITS → ITS2</td>
                    <td class="number">{MSA_ITS_ITS2:,}</td>
                    <td>Partial ITS sequence</td>
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

    # Write the HTML file
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)

    logger.info(f"✅ Comprehensive HTML report saved to: {output_html}")
    logger.info("=" * 60)
    logger.info("✅ Sequence flow report generated successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
