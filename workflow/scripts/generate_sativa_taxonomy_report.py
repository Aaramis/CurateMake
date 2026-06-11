#!/usr/bin/env python3
"""
Generate a full interactive HTML report with taxonomic sunburst visualizations.

The report includes:
- 5 KPIs (key pipeline metrics)
- Section 2: Sunburst of original databases (raw taxonomy)
- Section 3: Before/after CoL comparison (side-by-side)
- Interactive filters: database source + barcode (ITS/ITS1/ITS2)

Plotly handles all interactivity client-side (no page reloads).
"""

import sys
import json
import logging
import yaml
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Any
from collections import defaultdict
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sqlalchemy import func, case, text
from sqlalchemy.orm import sessionmaker
from dotenv import load_dotenv

# Import database models
sys.path.insert(0, str(Path(__file__).parent / "utils"))
from db_models import get_engine, get_session, Sequence, Taxonomy, SeqTax, Source
from db_manager import get_database_manager

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================
# CONFIGURATION
# ============================================

KINGDOM_COLORS = {
    'Fungi': '#7570b3',
    'Plantae': '#1b9e77',
    'Animalia': '#d95f02',
    'Chromista': '#e7298a',
    'Bacteria': '#66a61e',
    'Protozoa': '#e6ab02',
}
UNKNOWN_COLOR = '#CCCCCC'
GRADIENT_INTENSITY = 0.1
TAXO_LEVELS = ['kingdom', 'phylum', 'class_', 'order_', 'family', 'genus']
MAX_DEPTH = 3


# ============================================
# KPI QUERY FUNCTIONS
# ============================================

def get_kpi_metrics(engine, db_manager=None) -> Dict[str, int]:
    """
    Calculate all 5 KPI metrics.

    Returns dict with:
    - total_original: Original sequences (parent_seq_id IS NULL)
    - total_after: Non-duplicate sequences
    - total_sources: Number of distinct database sources
    - unique_taxa_original: Distinct lineages from original taxonomy
    - unique_taxa_after: Distinct lineages after processing (best backbone)
    """
    SessionFactory = sessionmaker(bind=engine)

    def _run(name, func):
        if db_manager:
            return db_manager._execute_with_retry(name, func)
        return func()

    def fetch_total_original():
        session = SessionFactory()
        try:
            return session.query(func.count(Sequence.seq_id)).filter(
                Sequence.parent_seq_id.is_(None)
            ).scalar()
        finally:
            session.close()

    def fetch_total_after():
        session = SessionFactory()
        try:
            return session.query(func.count(Sequence.seq_id)).filter(
                Sequence.barcode_status != 'Duplicate'
            ).scalar()
        finally:
            session.close()

    def fetch_total_sources():
        session = SessionFactory()
        try:
            return session.query(func.count(func.distinct(Source.name))).scalar()
        finally:
            session.close()

    def fetch_unique_taxa_original():
        # Step 1: seq_ids from sequences only (single-table indexed scan on parent_seq_id)
        def _fetch_seq_ids():
            session = SessionFactory()
            try:
                return [r[0] for r in session.query(Sequence.seq_id).filter(
                    Sequence.parent_seq_id.is_(None)
                ).all()]
            finally:
                session.close()

        seq_ids = _run("kpi_unique_taxa_original_seq_ids", _fetch_seq_ids)
        if not seq_ids:
            return 0

        # Step 2: collect distinct tax_ids from seq_tax in chunks (avoids cross-table subquery)
        all_tax_ids = set()
        chunk_size = 1000
        for i in range(0, len(seq_ids), chunk_size):
            chunk = seq_ids[i:i + chunk_size]

            def _fetch_chunk(chunk=chunk):
                session = SessionFactory()
                try:
                    return [r[0] for r in session.query(SeqTax.tax_id).filter(
                        SeqTax.backbone == 'Raw',
                        SeqTax.seq_id.in_(chunk)
                    ).distinct().all()]
                finally:
                    session.close()

            all_tax_ids.update(_run(f"kpi_unique_taxa_original_chunk_{i}", _fetch_chunk))

        return len(all_tax_ids)

    def fetch_unique_taxa_after():
        # Step 1: seq_ids from sequences only (single-table indexed scan on barcode_status)
        def _fetch_seq_ids():
            session = SessionFactory()
            try:
                return [r[0] for r in session.query(Sequence.seq_id).filter(
                    Sequence.barcode_status != 'Duplicate'
                ).all()]
            finally:
                session.close()

        seq_ids = _run("kpi_unique_taxa_after_seq_ids", _fetch_seq_ids)
        if not seq_ids:
            return 0

        # Step 2: collect distinct tax_ids from seq_tax in chunks (avoids cross-table subquery)
        all_tax_ids = set()
        chunk_size = 1000
        for i in range(0, len(seq_ids), chunk_size):
            chunk = seq_ids[i:i + chunk_size]

            def _fetch_chunk(chunk=chunk):
                session = SessionFactory()
                try:
                    return [r[0] for r in session.query(SeqTax.tax_id).filter(
                        SeqTax.backbone == 'CoL',
                        SeqTax.seq_id.in_(chunk)
                    ).distinct().all()]
                finally:
                    session.close()

            all_tax_ids.update(_run(f"kpi_unique_taxa_after_chunk_{i}", _fetch_chunk))

        return len(all_tax_ids)

    total_original = _run("kpi_total_original", fetch_total_original)
    total_after = _run("kpi_total_after", fetch_total_after)
    total_sources = _run("kpi_total_sources", fetch_total_sources)
    unique_taxa_original = _run("kpi_unique_taxa_original", fetch_unique_taxa_original)
    unique_taxa_after = _run("kpi_unique_taxa_after", fetch_unique_taxa_after)

    return {
        'total_original': total_original or 0,
        'total_after': total_after or 0,
        'total_sources': total_sources or 0,
        'unique_taxa_original': unique_taxa_original or 0,
        'unique_taxa_after': unique_taxa_after or 0
    }


# ============================================
# DATA QUERY HELPERS
# ============================================

def query_and_prepare_data(engine, backbone_filter: str, parent_seq_id_null: bool = False,
                           non_duplicate: bool = False, db_manager=None) -> pd.DataFrame:
    """
    Query sequences with specific filters using batched IN queries to avoid server
    connection drops on large datasets.

    Args:
        engine: Database engine
        backbone_filter: Backbone to filter ('Raw', 'CoL', 'Sativa')
        parent_seq_id_null: If True, only get original sequences (parent_seq_id IS NULL)
        non_duplicate: If True, filter out duplicates
        db_manager: Optional DatabaseManager for retry support

    Returns:
        DataFrame with sequences and taxonomy
    """
    SessionFactory = sessionmaker(bind=engine)
    batch_size = 1000

    # Step 1: fetch seq_ids from seq_tax only (single-table indexed scan, no distributed join)
    # Sequence-level filters (parent_seq_id_null, non_duplicate) are applied in Step 2
    # so that the heavy join only happens on small batches.
    def fetch_seq_ids():
        session = SessionFactory()
        try:
            return [r[0] for r in session.query(SeqTax.seq_id).filter(
                SeqTax.backbone == backbone_filter
            ).all()]
        finally:
            session.close()

    op_name_ids = f"query_prepare_data_{backbone_filter}_ids"
    if db_manager:
        seq_ids = db_manager._execute_with_retry(op_name_ids, fetch_seq_ids)
    else:
        seq_ids = fetch_seq_ids()

    if not seq_ids:
        logger.warning(f"No sequences found for backbone={backbone_filter}")
        return pd.DataFrame(columns=['seq_id', 'source_name', 'barcode', 'kingdom',
                                    'phylum', 'class_', 'order_', 'family', 'genus', 'backbone'])

    total_batches = (len(seq_ids) + batch_size - 1) // batch_size
    logger.info(f"Fetching {len(seq_ids)} sequences for backbone={backbone_filter} in {total_batches} batch(es)")

    # Step 2: fetch full taxonomy data in batches using IN on seq_ids
    all_rows = []
    for batch_num, i in enumerate(range(0, len(seq_ids), batch_size), 1):
        batch_ids = seq_ids[i:i + batch_size]

        def fetch_batch():
            session = SessionFactory()
            try:
                q = session.query(
                    Sequence.seq_id,
                    Source.name.label('source_name'),
                    Sequence.barcode,
                    Taxonomy.kingdom,
                    Taxonomy.phylum,
                    Taxonomy.class_.label('class_'),
                    Taxonomy.order_.label('order_'),
                    Taxonomy.family,
                    Taxonomy.genus,
                    SeqTax.backbone
                ).join(
                    Source, Sequence.source_id == Source.source_id
                ).join(
                    SeqTax, Sequence.seq_id == SeqTax.seq_id
                ).join(
                    Taxonomy, SeqTax.tax_id == Taxonomy.tax_id
                ).filter(
                    SeqTax.backbone == backbone_filter,
                    Sequence.seq_id.in_(batch_ids)
                )
                if parent_seq_id_null:
                    q = q.filter(Sequence.parent_seq_id.is_(None))
                if non_duplicate:
                    q = q.filter(Sequence.barcode_status != 'Duplicate')
                return q.all()
            finally:
                session.close()

        op_name = f"query_prepare_data_{backbone_filter}_batch_{i}"
        if db_manager:
            rows = db_manager._execute_with_retry(op_name, fetch_batch)
        else:
            rows = fetch_batch()

        all_rows.extend(rows)
        if batch_num % 10 == 0 or batch_num == total_batches:
            logger.info(f"  {backbone_filter}: {batch_num}/{total_batches} batches ({len(all_rows)} rows so far)")

    if not all_rows:
        logger.warning(f"No sequences found for backbone={backbone_filter}")
        return pd.DataFrame(columns=['seq_id', 'source_name', 'barcode', 'kingdom',
                                    'phylum', 'class_', 'order_', 'family', 'genus', 'backbone'])

    df = pd.DataFrame([{
        'seq_id': r.seq_id,
        'source_name': r.source_name,
        'barcode': r.barcode,
        'kingdom': r.kingdom,
        'phylum': r.phylum,
        'class_': r.class_,
        'order_': r.order_,
        'family': r.family,
        'genus': r.genus,
        'backbone': r.backbone
    } for r in all_rows])

    # Fill NULLs with "Unknown"
    for col in TAXO_LEVELS:
        df[col] = df[col].fillna('Unknown').replace('', 'Unknown')

    return df


# ============================================
# SUNBURST DATA PREPARATION
# (Same logic as test_sunburst_db.py)
# ============================================

def prepare_sunburst_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare hierarchical data for sunburst plot.
    Logic from test_sunburst_db.py - creates paths, identifies leaves, propagates values.
    """
    if df.empty:
        return pd.DataFrame()

    # Build hierarchical paths
    df['path'] = df[TAXO_LEVELS].apply(
        lambda row: '-'.join(row.astype(str)), axis=1
    )

    # Group and count sequences
    df_sunburst = df.groupby('path').size().reset_index(name='count')
    df_sunburst['kingdom'] = df_sunburst['path'].str.split('-').str[0]
    df_sunburst['level'] = df_sunburst['path'].str.split('-').str.len()

    return df_sunburst


def transform_to_plotly(df_sunburst: pd.DataFrame) -> pd.DataFrame:
    """
    Transform sunburst data to Plotly format.
    Logic from test_sunburst_db.py - identifies leaves, creates nodes, propagates values.
    """
    if df_sunburst.empty:
        return pd.DataFrame()

    # Identify leaf nodes (terminal paths)
    all_paths = set(df_sunburst['path'].tolist())
    leaf_nodes = set()

    for path in all_paths:
        is_leaf = True
        for other_path in all_paths:
            if other_path != path and other_path.startswith(path + '-'):
                is_leaf = False
                break
        if is_leaf:
            leaf_nodes.add(path)

    # Build nodes dictionary with values ONLY on true leaves
    nodes_dict = {}

    for idx, row in df_sunburst.iterrows():
        parts = row['path'].split('-')
        count = row['count']
        kingdom = parts[0]
        path = row['path']

        # Create all nodes in this path
        for lvl in range(len(parts)):
            node_id = '-'.join(parts[:lvl + 1])
            label = parts[lvl]
            parent_id = '-'.join(parts[:lvl]) if lvl > 0 else ''

            if node_id not in nodes_dict:
                nodes_dict[node_id] = {
                    'ids': node_id,
                    'labels': label,
                    'parents': parent_id,
                    'kingdom': kingdom,
                    'level': lvl + 1,
                    'values': 0,
                    'is_leaf': node_id in leaf_nodes
                }

        # Add count ONLY to the actual terminal node (the full path)
        if path in leaf_nodes:
            nodes_dict[path]['values'] += count

    # Calculate parent values bottom-up (sum of direct children)
    sorted_nodes = sorted(nodes_dict.values(), key=lambda x: x['level'], reverse=True)

    # Build children map
    children_map = {}
    for node in nodes_dict.values():
        parent_id = node['parents']
        if parent_id:
            if parent_id not in children_map:
                children_map[parent_id] = []
            children_map[parent_id].append(node['ids'])

    # Propagate values upward: parent = sum of children
    for node in sorted_nodes:
        node_id = node['ids']
        if node_id in children_map:
            children_sum = sum(nodes_dict[child_id]['values'] for child_id in children_map[node_id])
            nodes_dict[node_id]['values'] = children_sum

    # Convert to DataFrame
    result_nodes = [{k: v for k, v in node.items() if k != 'is_leaf'}
                    for node in nodes_dict.values()]
    plotly_data = pd.DataFrame(result_nodes)
    plotly_data = plotly_data.sort_values(['level', 'kingdom', 'labels'])

    return plotly_data


def prepare_sunburst_data_by_filters(df: pd.DataFrame) -> Dict:
    """
    Prepare sunburst data for all combinations of database and barcode filters.

    Returns:
        Nested dict: data[source_name][barcode] = plotly_data_df
    """
    data_by_filters = {}

    # Handle empty dataframe
    if df.empty:
        logger.warning("Empty dataframe provided to prepare_sunburst_data_by_filters")
        return {'All': {'All': pd.DataFrame()}}

    # Get unique sources and barcodes
    sources = ['All'] + sorted(df['source_name'].unique().tolist())
    barcodes = ['All'] + sorted(df['barcode'].unique().tolist())

    for source in sources:
        data_by_filters[source] = {}

        for barcode in barcodes:
            # Filter data
            df_filtered = df.copy()

            if source != 'All':
                df_filtered = df_filtered[df_filtered['source_name'] == source]

            if barcode != 'All':
                df_filtered = df_filtered[df_filtered['barcode'] == barcode]

            if len(df_filtered) == 0:
                data_by_filters[source][barcode] = pd.DataFrame()
                continue

            # Prepare sunburst data
            df_sunburst = prepare_sunburst_data(df_filtered)
            plotly_data = transform_to_plotly(df_sunburst)

            data_by_filters[source][barcode] = plotly_data

    return data_by_filters


# ============================================
# FONCTIONS DE VISUALISATION
# ============================================

def generate_colors(plotly_data: pd.DataFrame) -> List[str]:
    """Generate colors for nodes based on kingdom and level."""
    if plotly_data.empty:
        return []

    colors = []
    for idx, row in plotly_data.iterrows():
        kingdom = row['kingdom']
        level = row['level']

        if kingdom in KINGDOM_COLORS:
            base_hex = KINGDOM_COLORS[kingdom]
            r = int(base_hex[1:3], 16)
            g = int(base_hex[3:5], 16)
            b = int(base_hex[5:7], 16)

            factor = (level - 1) * GRADIENT_INTENSITY
            new_r = min(int(r + (255 - r) * factor), 255)
            new_g = min(int(g + (255 - g) * factor), 255)
            new_b = min(int(b + (255 - b) * factor), 255)

            color = f'#{new_r:02x}{new_g:02x}{new_b:02x}'
        else:
            color = UNKNOWN_COLOR

        colors.append(color)

    return colors


def create_four_way_comparison_sunburst(data_before: Dict, data_after_col: Dict,
                                        data_after_itsx: Dict, data_after_sativa: Dict) -> Tuple[go.Figure, Dict]:
    """
    Create 4-way comparison (2x2 grid):
    - Top left: Before CoL (Raw, original sequences)
    - Top right: After CoL (CoL, original sequences only)
    - Bottom left: After ITSx (CoL, all non-duplicate sequences)
    - Bottom right: After Sativa (Sativa, all non-duplicate sequences)

    All have same interactive filters.
    Returns tuple of (figure, metadata)
    """
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'sunburst'}, {'type': 'sunburst'}],
               [{'type': 'sunburst'}, {'type': 'sunburst'}]],
        subplot_titles=('Before CoL (Raw only)',
                       'After CoL (Original sequences)',
                       'After ITSx (All sequences)',
                       'After Sativa (Sativa only)'),
        vertical_spacing=0.12,
        horizontal_spacing=0.05
    )

    sources = sorted(data_before.keys())
    if not sources:
        return go.Figure(), {}
    barcodes = sorted(list(data_before[sources[0]].keys()))

    # Build trace maps dynamically as we add traces
    trace_map_before = {}
    trace_map_after_col = {}
    trace_map_after_itsx = {}
    trace_map_after_sativa = {}
    current_trace_idx = 0
    combination_index = 0

    for source in sources:
        for barcode in barcodes:
            key = f'{source}_{barcode}'

            # TOP LEFT: Before CoL
            plotly_data_before = data_before.get(source, {}).get(barcode, pd.DataFrame())
            if len(plotly_data_before) > 0:
                colors = generate_colors(plotly_data_before)
                fig.add_trace(go.Sunburst(
                    ids=plotly_data_before['ids'],
                    labels=plotly_data_before['labels'],
                    parents=plotly_data_before['parents'],
                    values=plotly_data_before['values'],
                    branchvalues='total',
                    marker=dict(colors=colors, line=dict(color='white', width=0.5)),
                    hovertemplate='<b>%{label}</b><br>Sequences: %{value:,}<br><extra></extra>',
                    textinfo='label+value',
                    visible=(combination_index == 0),
                    maxdepth=MAX_DEPTH,
                    name=f'{source}_{barcode}_before'
                ), row=1, col=1)
                trace_map_before[key] = current_trace_idx
                current_trace_idx += 1

            # TOP RIGHT: After CoL
            plotly_data_after_col = data_after_col.get(source, {}).get(barcode, pd.DataFrame())
            if len(plotly_data_after_col) > 0:
                colors = generate_colors(plotly_data_after_col)
                fig.add_trace(go.Sunburst(
                    ids=plotly_data_after_col['ids'],
                    labels=plotly_data_after_col['labels'],
                    parents=plotly_data_after_col['parents'],
                    values=plotly_data_after_col['values'],
                    branchvalues='total',
                    marker=dict(colors=colors, line=dict(color='white', width=0.5)),
                    hovertemplate='<b>%{label}</b><br>Sequences: %{value:,}<br><extra></extra>',
                    textinfo='label+value',
                    visible=(combination_index == 0),
                    maxdepth=MAX_DEPTH,
                    name=f'{source}_{barcode}_after_col'
                ), row=1, col=2)
                trace_map_after_col[key] = current_trace_idx
                current_trace_idx += 1

            # BOTTOM LEFT: After ITSx
            plotly_data_after_itsx = data_after_itsx.get(source, {}).get(barcode, pd.DataFrame())
            if len(plotly_data_after_itsx) > 0:
                colors = generate_colors(plotly_data_after_itsx)
                fig.add_trace(go.Sunburst(
                    ids=plotly_data_after_itsx['ids'],
                    labels=plotly_data_after_itsx['labels'],
                    parents=plotly_data_after_itsx['parents'],
                    values=plotly_data_after_itsx['values'],
                    branchvalues='total',
                    marker=dict(colors=colors, line=dict(color='white', width=0.5)),
                    hovertemplate='<b>%{label}</b><br>Sequences: %{value:,}<br><extra></extra>',
                    textinfo='label+value',
                    visible=(combination_index == 0),
                    maxdepth=MAX_DEPTH,
                    name=f'{source}_{barcode}_after_itsx'
                ), row=2, col=1)
                trace_map_after_itsx[key] = current_trace_idx
                current_trace_idx += 1

            # BOTTOM RIGHT: After Sativa
            plotly_data_after_sativa = data_after_sativa.get(source, {}).get(barcode, pd.DataFrame())
            if len(plotly_data_after_sativa) > 0:
                colors = generate_colors(plotly_data_after_sativa)
                fig.add_trace(go.Sunburst(
                    ids=plotly_data_after_sativa['ids'],
                    labels=plotly_data_after_sativa['labels'],
                    parents=plotly_data_after_sativa['parents'],
                    values=plotly_data_after_sativa['values'],
                    branchvalues='total',
                    marker=dict(colors=colors, line=dict(color='white', width=0.5)),
                    hovertemplate='<b>%{label}</b><br>Sequences: %{value:,}<br><extra></extra>',
                    textinfo='label+value',
                    visible=(combination_index == 0),
                    maxdepth=MAX_DEPTH,
                    name=f'{source}_{barcode}_after_sativa'
                ), row=2, col=2)
                trace_map_after_sativa[key] = current_trace_idx
                current_trace_idx += 1

            combination_index += 1

    fig.update_layout(
        title=dict(
            text=f'4-Way Taxonomic Comparison',
            x=0.5,
            xanchor='center',
            font=dict(size=20)
        ),
        height=1200,
        margin=dict(t=80, b=20, l=20, r=20),
        showlegend=False
    )

    metadata = {
        'sources': sources,
        'barcodes': barcodes,
        'trace_map_before': trace_map_before,
        'trace_map_after_col': trace_map_after_col,
        'trace_map_after_itsx': trace_map_after_itsx,
        'trace_map_after_sativa': trace_map_after_sativa
    }

    return fig, metadata


# ============================================
# HTML REPORT GENERATION
# ============================================

def generate_html_report(kpis: Dict, fig_comparison: go.Figure, metadata_comparison: Dict,
                        output_path: Path) -> None:
    """Generate comprehensive HTML report with KPIs and 4-way comparison sunburst."""

    output_path.parent.mkdir(parents=True, exist_ok=True)

    html_content = f"""
<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Taxonomic Distribution Report - CurateMake</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
            color: #2c3e50;
        }}

        .container {{ max-width: 1800px; margin: 0 auto; }}

        .header {{
            text-align: center;
            color: white;
            margin-bottom: 30px;
            padding: 30px 20px;
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }}

        .timestamp {{
            background: rgba(255,255,255,0.2);
            padding: 8px 16px;
            border-radius: 20px;
            display: inline-block;
            margin-top: 10px;
            font-size: 0.9em;
        }}

        /* KPI Cards */
        .kpi-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 40px;
        }}

        .kpi-card {{
            background: white;
            border-radius: 15px;
            padding: 25px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            text-align: center;
            transition: transform 0.3s ease;
        }}

        .kpi-card:hover {{ transform: translateY(-5px); }}

        .kpi-icon {{ font-size: 2.5em; margin-bottom: 10px; }}

        .kpi-value {{
            font-size: 2.5em;
            font-weight: bold;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin: 10px 0;
        }}

        .kpi-label {{
            color: #7f8c8d;
            font-size: 0.95em;
            font-weight: 500;
        }}

        /* Sections */
        .section {{
            background: white;
            border-radius: 15px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
        }}

        .section h2 {{
            color: #2c3e50;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
            font-size: 1.8em;
        }}

        .info-box {{
            background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
            border-left: 4px solid #2196f3;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 8px;
        }}

        .plot-container {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
        }}

        .footer {{
            text-align: center;
            color: white;
            margin-top: 40px;
            padding: 20px;
            opacity: 0.8;
        }}

        @media (max-width: 768px) {{
            .kpi-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🌳 Taxonomic Distribution Report</h1>
            <p>Comprehensive overview of the pipeline's taxonomic distribution</p>
            <div class="timestamp">Generated on {datetime.now().strftime('%d/%m/%Y at %H:%M:%S')}</div>
        </div>

        <!-- KPI Cards -->
        <div class="kpi-grid">
            <div class="kpi-card">
                <div class="kpi-icon">📦</div>
                <div class="kpi-value">{kpis['total_original']:,}</div>
                <div class="kpi-label">Total Seq in Original DB</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">✅</div>
                <div class="kpi-value">{kpis['total_after']:,}</div>
                <div class="kpi-label">Total Seq After Processes</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">🗄️</div>
                <div class="kpi-value">{kpis['total_sources']}</div>
                <div class="kpi-label">Total Database Sources</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">🧬</div>
                <div class="kpi-value">{kpis['unique_taxa_original']:,}</div>
                <div class="kpi-label">Unique Taxa Originally</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">🎯</div>
                <div class="kpi-value">{kpis['unique_taxa_after']:,}</div>
                <div class="kpi-label">Unique Taxa After Processes</div>
            </div>
        </div>

        <!-- 4-Way Comparison -->
        <div class="section">
            <h2>🔄 Pipeline Taxonomic Comparison (2x2 Grid)</h2>
            <p>Complete comparison showing how taxonomy evolves throughout the pipeline.</p>
            <div class="info-box">
                <strong>Top-Left:</strong> Before CoL (Raw only) — Original sequences with raw taxonomy<br>
                <strong>Top-Right:</strong> After CoL (Original sequences) — Original non-duplicate sequences with harmonized CoL taxonomy<br>
                <strong>Bottom-Left:</strong> After ITSx (All sequences) — All non-duplicate sequences (including ITSx extractions) using CoL taxonomy<br>
                <strong>Bottom-Right:</strong> After Sativa (Sativa only) — All non-duplicate sequences validated by SATIVA
            </div>

            <!-- Custom Dropdowns for linked filtering -->
            <div style="margin: 20px 0; display: flex; gap: 20px; justify-content: center;">
                <div>
                    <label for="db-select-comparison" style="font-weight: bold; margin-right: 10px;">Database:</label>
                    <select id="db-select-comparison" style="padding: 8px 12px; border-radius: 5px; border: 2px solid #667eea; font-size: 14px;">
                        {chr(10).join(f'<option value="{s}">{s}</option>' for s in metadata_comparison['sources'])}
                    </select>
                </div>
                <div>
                    <label for="bc-select-comparison" style="font-weight: bold; margin-right: 10px;">Barcode:</label>
                    <select id="bc-select-comparison" style="padding: 8px 12px; border-radius: 5px; border: 2px solid #667eea; font-size: 14px;">
                        {chr(10).join(f'<option value="{b}">{b}</option>' for b in metadata_comparison['barcodes'])}
                    </select>
                </div>
            </div>

            <div class="plot-container">
                <div id="sunburst-comparison"></div>
            </div>
        </div>

        <div class="footer">
            <p><strong>CurateMake Pipeline</strong> - ITS Taxonomy Curation & Analysis</p>
            <p>© 2025 - Automatically generated</p>
        </div>
    </div>

    <script>
        // Metadata for filtering
        var metadataComparison = {json.dumps(metadata_comparison)};

        // Render 4-Way Comparison
        var figComparison = {fig_comparison.to_json()};
        Plotly.newPlot('sunburst-comparison', figComparison.data, figComparison.layout, {{responsive: true}});

        // Linked filtering for 4-Way Comparison
        function updateComparisonPlot() {{
            var selectedDb = document.getElementById('db-select-comparison').value;
            var selectedBc = document.getElementById('bc-select-comparison').value;
            var key = selectedDb + '_' + selectedBc;

            var traceIndexBefore = metadataComparison.trace_map_before[key];
            var traceIndexAfterCol = metadataComparison.trace_map_after_col[key];
            var traceIndexAfterItsx = metadataComparison.trace_map_after_itsx[key];
            var traceIndexAfterSativa = metadataComparison.trace_map_after_sativa[key];

            // Create visibility array for all 4 sunbursts
            var visibility = figComparison.data.map((_, idx) =>
                idx === traceIndexBefore || idx === traceIndexAfterCol ||
                idx === traceIndexAfterItsx || idx === traceIndexAfterSativa
            );

            Plotly.restyle('sunburst-comparison', {{visible: visibility}});
        }}

        document.getElementById('db-select-comparison').addEventListener('change', updateComparisonPlot);
        document.getElementById('bc-select-comparison').addEventListener('change', updateComparisonPlot);
    </script>
</body>
</html>
"""

    output_path.write_text(html_content)


# ============================================
# MAIN FUNCTION
# ============================================

def main():
    """Entry point."""

    # Configuration
    if 'snakemake' in globals():
        output_path = Path(snakemake.output.report)
        db_manager = get_database_manager(snakemake.config)
        engine = db_manager.engine
    else:
        output_path = Path('reports/html/taxonomy_sunburst.html')
        db_path = Path('data/db/its.db')
        if not db_path.exists():
            logger.error(f"Database not found: {db_path}")
            sys.exit(1)
        db_manager = None
        engine = get_engine("sqlite", path=str(db_path))

    logger.info("Generating taxonomic report...")

    # Step 1: Calculate KPIs
    kpis = get_kpi_metrics(engine, db_manager=db_manager)

    # Step 2: Query data
    # Before CoL: Raw backbone, original sequences only
    df_before_col = query_and_prepare_data(engine, 'Raw', parent_seq_id_null=True, non_duplicate=False, db_manager=db_manager)

    # After CoL: CoL backbone, original non-duplicate sequences
    df_after_col = query_and_prepare_data(engine, 'CoL', parent_seq_id_null=True, non_duplicate=True, db_manager=db_manager)

    # After ITSx: CoL backbone, all non-duplicate sequences
    df_after_itsx = query_and_prepare_data(engine, 'CoL', parent_seq_id_null=False, non_duplicate=True, db_manager=db_manager)

    # After Sativa: Sativa backbone, all non-duplicate sequences
    df_after_sativa = query_and_prepare_data(engine, 'Sativa', parent_seq_id_null=False, non_duplicate=True, db_manager=db_manager)

    if df_before_col.empty:
        logger.error("No original sequences found with Raw backbone!")
        sys.exit(1)

    # Step 3: Prepare data for all filter combinations
    data_before_col = prepare_sunburst_data_by_filters(df_before_col)
    data_after_col = prepare_sunburst_data_by_filters(df_after_col)
    data_after_itsx = prepare_sunburst_data_by_filters(df_after_itsx)
    data_after_sativa = prepare_sunburst_data_by_filters(df_after_sativa)

    # Step 4: Create 4-way comparison visualization
    fig_comparison, metadata_comparison = create_four_way_comparison_sunburst(
        data_before_col,
        data_after_col,
        data_after_itsx,
        data_after_sativa
    )

    # Step 5: Generate HTML report
    generate_html_report(
        kpis,
        fig_comparison, metadata_comparison,
        output_path
    )

    logger.info(f"✅ Report saved: {output_path}")


if __name__ == "__main__":
    main()
