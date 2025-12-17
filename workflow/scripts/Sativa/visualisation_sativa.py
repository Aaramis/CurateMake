#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SATIVA Results Analyzer
=======================
A comprehensive tool for analyzing and visualizing taxonomic correction results from SATIVA.

Architecture:
- SativaParser: Parses .mis and .tax files
- TaxonomyAnalyzer: Computes metrics and statistics
- VisualizationEngine: Generates interactive visualizations
- ReportGenerator: Creates structured HTML reports
- SativaPipeline: Orchestrates the analysis workflow
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from datetime import datetime
import logging
import sys
from sqlalchemy.orm import sessionmaker

# Database imports
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from db_models import (
    Sequence as DBSequence,
    SeqTax,
    Taxonomy,
    Source,
    SativaWarning,
    get_engine
)
from db_manager import get_database_manager


# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)
logger = logging.getLogger(__name__)


# ==================== Data Structures ====================

class TaxonomicLevel(Enum):
    KINGDOM = 0
    PHYLUM = 1
    CLASS = 2
    ORDER = 3
    FAMILY = 4
    GENUS = 5
    SPECIES = 6


class ModificationType(Enum):
    COMPLETION = "Completion"
    DELETION = "Deletion"
    DELETION_SIMPLE = "Deletion (Simple)"
    DELETION_MULTIPLE = "Deletion (Multiple)"
    SUBSTITUTION = "Substitution"


@dataclass
class TaxonomicModification:
    seq_id: str
    barcode: str
    msa_file: str
    original_label: str
    proposed_label: str
    confidence: float
    original_path: List[str]
    proposed_path: List[str]
    modification_level: TaxonomicLevel
    modification_type: ModificationType
    per_rank_confidence: List[float] = field(default_factory=list)


@dataclass
class SequenceResult:
    seq_id: str
    barcode: str
    msa_file: str
    original_path: List[str]
    modification: Optional[TaxonomicModification] = None
    source: Optional[str] = None

    @property
    def is_modified(self) -> bool:
        return self.modification is not None


def humanize_count(n: int) -> str:
    """Format large numbers with K/M suffixes."""
    if n >= 1_000_000:
        val = n / 1_000_000
        suffix = "M"
    elif n >= 1_000:
        val = n / 1_000
        suffix = "K"
    else:
        return f"{n:,}"

    formatted = f"{val:.1f}".rstrip('0').rstrip('.')
    return f"{formatted}{suffix}"


# ==================== Parser ====================

class SativaParser:
    """Parses SATIVA output files (.mis and .tax)."""

    TAX_LEVELS = [
        TaxonomicLevel.KINGDOM,
        TaxonomicLevel.PHYLUM,
        TaxonomicLevel.CLASS,
        TaxonomicLevel.ORDER,
        TaxonomicLevel.FAMILY,
        TaxonomicLevel.GENUS,
        TaxonomicLevel.SPECIES
    ]

    @staticmethod
    def parse_mis_file(filepath: Path, barcode: str, msa_name: str) -> List[TaxonomicModification]:
        """Parse a .mis file containing SATIVA misclassification results."""
        modifications = []
        with open(filepath, 'r') as f:
            lines = f.readlines()

        if len(lines) < 6:
            return modifications

        for line in lines[5:]:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            seq_id = parts[0]
            original_label = parts[2]
            proposed_label = parts[3]
            confidence = float(parts[4])
            original_path = [x for x in parts[5].split(';') if x]
            proposed_path = [x for x in parts[6].split(';') if x]
            per_rank_conf = [float(x) for x in parts[7].split(';') if x]

            mod_level, mod_type = SativaParser._analyze_modification(original_path, proposed_path)

            modifications.append(TaxonomicModification(
                seq_id=seq_id,
                barcode=barcode,
                msa_file=msa_name,
                original_label=original_label,
                proposed_label=proposed_label,
                confidence=confidence,
                original_path=original_path,
                proposed_path=proposed_path,
                modification_level=mod_level,
                modification_type=mod_type,
                per_rank_confidence=per_rank_conf
            ))
        return modifications

    @staticmethod
    def parse_tax_file(filepath: Path) -> Dict[str, List[str]]:
        """Parse a .tax file containing taxonomic classifications."""
        taxonomies = {}
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    seq_id = parts[0]
                    tax_path = [x for x in parts[1].split(';') if x]
                    taxonomies[seq_id] = tax_path
        return taxonomies

    @staticmethod
    def _analyze_modification(original: List[str], proposed: List[str]) -> Tuple[TaxonomicLevel, ModificationType]:
        """Determine the taxonomic level and type of modification."""
        mod_level = TaxonomicLevel.KINGDOM
        min_len = min(len(original), len(proposed))

        for i in range(min_len):
            if original[i] != proposed[i]:
                mod_level = SativaParser.TAX_LEVELS[i] if i < len(SativaParser.TAX_LEVELS) else TaxonomicLevel.SPECIES
                break
        else:
            if len(proposed) > len(original):
                mod_level = SativaParser.TAX_LEVELS[len(original)] if len(original) < len(SativaParser.TAX_LEVELS) else TaxonomicLevel.SPECIES
            elif len(original) > len(proposed):
                mod_level = SativaParser.TAX_LEVELS[len(proposed)] if len(proposed) < len(SativaParser.TAX_LEVELS) else TaxonomicLevel.SPECIES

        if len(proposed) > len(original):
            mod_type = ModificationType.COMPLETION
        elif len(proposed) < len(original):
            diff = len(original) - len(proposed)
            if diff <= 1:
                mod_type = ModificationType.DELETION_SIMPLE
            else:
                mod_type = ModificationType.DELETION_MULTIPLE
        else:
            mod_type = ModificationType.SUBSTITUTION

        return mod_level, mod_type


# ==================== Analyzer ====================

class TaxonomyAnalyzer:
    """Analyzes taxonomic modifications and computes statistics."""

    def __init__(self, all_results: List[SequenceResult]):
        self.all_results = all_results
        self.df = self._to_dataframe()

    def _to_dataframe(self) -> pd.DataFrame:
        """Convert results to pandas DataFrame for analysis."""
        data = []
        for result in self.all_results:
            mod = result.modification
            row = {
                'seq_id': result.seq_id,
                'barcode': result.barcode,
                'msa_file': result.msa_file,
                'source': result.source or "Unknown",
                'original_depth': len(result.original_path),
                'original_path': result.original_path,
                'is_modified': result.is_modified
            }
            if result.is_modified:
                row.update({
                    'original_label': mod.original_label,
                    'proposed_label': mod.proposed_label,
                    'confidence': mod.confidence,
                    'modification_level': mod.modification_level.name,
                    'modification_type': mod.modification_type.value,
                    'proposed_depth': len(mod.proposed_path),
                    'proposed_path': mod.proposed_path
                })
            data.append(row)
        df = pd.DataFrame(data)
        required_columns = [
            'seq_id', 'barcode', 'msa_file', 'source',
            'original_depth', 'original_path', 'is_modified',
            'original_label', 'proposed_label', 'confidence',
            'modification_level', 'modification_type',
            'proposed_depth', 'proposed_path'
        ]
        for column in required_columns:
            if column not in df.columns:
                if column == 'is_modified':
                    df[column] = False
                else:
                    df[column] = pd.Series(dtype=object)
        return df

    def get_modification_rate_by(self, level: str) -> pd.DataFrame:
        """Calculate modification rates grouped by specified taxonomic level or category."""
        data = self.df.copy()
        if level.lower() == 'barcode':
            data['category'] = data['barcode']
        elif level.lower() == 'source':
            data['category'] = data['source'].fillna('Unknown')
        else:
            try:
                idx = TaxonomicLevel[level.upper()].value
            except KeyError:
                raise ValueError(f"Unknown level: {level}")
            data['category'] = data['original_path'].apply(
                lambda path: path[idx] if len(path) > idx else "Unclassified"
            )

        if data.empty:
            return pd.DataFrame(columns=['category', 'total', 'modified', 'percentage'])

        stats = []
        for cat in data['category'].unique():
            subset = data[data['category'] == cat]
            total = len(subset)
            modified = subset['is_modified'].sum()
            perc = (modified / total * 100) if total > 0 else 0
            stats.append({
                'category': cat,
                'total': total,
                'modified': modified,
                'percentage': round(perc, 2)
            })
        result = pd.DataFrame(stats)
        if result.empty:
            return pd.DataFrame(columns=['category', 'total', 'modified', 'percentage'])
        return result.sort_values('percentage', ascending=False)

    def get_modifications_by_level(self) -> pd.DataFrame:
        """Get distribution of modifications across taxonomic levels."""
        modified = self.df[self.df['is_modified']]
        if len(modified) == 0:
            return pd.DataFrame(columns=['level', 'count', 'percentage'])
        counts = modified['modification_level'].value_counts()
        total = len(modified)
        return pd.DataFrame({
            'level': counts.index,
            'count': counts.values,
            'percentage': [round(c/total*100, 2) for c in counts.values]
        })

    def get_modifications_by_type(self) -> pd.DataFrame:
        """Get distribution of modification types."""
        modified = self.df[self.df['is_modified']]
        if len(modified) == 0:
            return pd.DataFrame(columns=['type', 'count', 'percentage'])
        counts = modified['modification_type'].value_counts()
        total = len(modified)
        return pd.DataFrame({
            'type': counts.index,
            'count': counts.values,
            'percentage': [round(c/total*100, 2) for c in counts.values]
        })

    def get_confidence_statistics(self) -> Dict[str, float]:
        """Calculate confidence score statistics."""
        modified = self.df[self.df['is_modified']]
        if len(modified) == 0:
            return {'mean': 0, 'median': 0, 'min': 0, 'max': 0, 'std': 0}
        return {
            'mean': round(modified['confidence'].mean(), 3),
            'median': round(modified['confidence'].median(), 3),
            'min': round(modified['confidence'].min(), 3),
            'max': round(modified['confidence'].max(), 3),
            'std': round(modified['confidence'].std(), 3)
        }


# ==================== Visualization ====================

class VisualizationEngine:
    """Generates interactive visualizations using Plotly."""

    COLORS = {
        'primary': '#2563eb',
        'secondary': '#7c3aed',
        'success': '#10b981',
        'warning': '#f59e0b',
        'danger': '#ef4444',
        'info': '#06b6d4'
    }

    def __init__(self, analyzer: TaxonomyAnalyzer):
        self.analyzer = analyzer

    def plot_modification_rate(self) -> go.Figure:
        """
        Creates an interactive bar plot with dropdown menu for selecting
        taxonomic level or technical category.
        """
        levels = ["barcode", "source", "kingdom", "phylum", "class",
                  "order", "family", "genus", "species"]

        traces = []
        visibility = []
        buttons = []

        for i, level in enumerate(levels):
            data = self.analyzer.get_modification_rate_by(level)

            trace = go.Bar(
                x=data['category'],
                y=data['percentage'],
                name=level.capitalize(),
                marker_color=self.COLORS['primary'],
                text=[f"{p}%" for p in data['percentage']],
                textposition='outside',
                hovertemplate='<b>%{x}</b><br>' +
                             'Modified: %{y:.2f}%<br>' +
                             '<extra></extra>',
                visible=(i == 0)
            )
            traces.append(trace)
            visibility.append([j == i for j in range(len(levels))])

            buttons.append(dict(
                label=level.capitalize(),
                method="update",
                args=[
                    {"visible": visibility[i]},
                    {
                        "title": f"Modification Rate by {level.capitalize()}",
                        "xaxis": {"title": level.capitalize()}
                    }
                ]
            ))

        fig = go.Figure(data=traces)

        fig.update_layout(
            title=f"Modification Rate by {levels[0].capitalize()}",
            xaxis=dict(title=levels[0].capitalize()),
            yaxis_title='Percentage of Modified Sequences (%)',
            template='plotly_white',
            height=500,
            font=dict(family="Arial, sans-serif", size=12),
            updatemenus=[{
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 1.15,
                "y": 1.15,
                "xanchor": "left",
                "yanchor": "top"
            }]
        )

        return fig

    def plot_modifications_by_level(self) -> go.Figure:
        """Visualizes the distribution of modifications across taxonomic levels."""
        data = self.analyzer.get_modifications_by_level()
        if len(data) == 0:
            fig = go.Figure()
            fig.add_annotation(
                text="No modifications detected",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig

        text = [f"{p}% ({humanize_count(c)})" for p, c in zip(data['percentage'], data['count'])]

        fig = px.bar(
            data,
            x='level',
            y='percentage',
            text=text,
            color='percentage',
            color_continuous_scale='Blues',
            title="Distribution of Modifications Across Taxonomic Levels"
        )

        fig.update_traces(
            textposition='outside',
            textfont_size=12
        )

        fig.update_layout(
            xaxis_title='Taxonomic Level',
            yaxis_title='Percentage of Modified Sequences (%)',
            template='plotly_white',
            font=dict(family="Arial, sans-serif", size=12),
            showlegend=False
        )

        return fig

    def plot_modifications_by_type(self) -> go.Figure:
        """Visualizes the distribution of modification types."""
        data = self.analyzer.get_modifications_by_type()
        if len(data) == 0:
            fig = go.Figure()
            fig.add_annotation(
                text="No modifications detected",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig

        color_map = {
            "Completion": self.COLORS['success'],
            "Substitution": self.COLORS['info'],
            "Deletion": self.COLORS['danger'],
            "Deletion (Simple)": self.COLORS['warning'],
            "Deletion (Multiple)": self.COLORS['danger'],
        }
        colors = [color_map.get(t, self.COLORS['primary']) for t in data['type']]

        fig = go.Figure(data=[go.Pie(
            labels=data['type'],
            values=data['count'],
            hole=0.4,
            marker_colors=colors[:len(data)],
            textinfo='label+percent',
            textfont_size=12,
            hovertemplate='<b>%{label}</b><br>' +
                         'Count: %{value}<br>' +
                         'Percentage: %{percent}<br>' +
                         '<extra></extra>'
        )])

        fig.update_layout(
            title="Distribution of Modification Types",
            template='plotly_white',
            font=dict(family="Arial, sans-serif", size=12)
        )

        return fig

    def plot_confidence_distribution(self) -> go.Figure:
        """Visualizes the distribution of confidence scores."""
        modified = self.analyzer.df[self.analyzer.df['is_modified']]

        if len(modified) == 0:
            fig = go.Figure()
            fig.add_annotation(
                text="No modifications detected",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig

        fig = go.Figure(data=[go.Histogram(
            x=modified['confidence'],
            nbinsx=30,
            marker_color=self.COLORS['primary'],
            opacity=0.7
        )])

        stats = self.analyzer.get_confidence_statistics()

        fig.add_vline(
            x=stats['mean'],
            line_dash="dash",
            line_color=self.COLORS['danger'],
            annotation_text=f"Mean: {stats['mean']:.3f}",
            annotation_position="top"
        )

        fig.update_layout(
            title="Distribution of Confidence Scores",
            xaxis_title='Confidence Score',
            yaxis_title='Frequency',
            template='plotly_white',
            font=dict(family="Arial, sans-serif", size=12),
            showlegend=False
        )

        return fig


    def plot_taxonomic_flow_sankey(self, level: str = "genus", top_n: int = 20) -> go.Figure:
        """
        Creates a Sankey diagram showing taxonomic reassignment flows.
        Each taxon appears twice: once as source (Original) and once as target (Corrected).

        Args:
            level: Taxonomic level to analyze ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
            top_n: Number of top flows to display (to avoid overcrowding)

        Returns:
            Plotly Figure object with Sankey diagram
        """
        modified = self.analyzer.df[self.analyzer.df['is_modified']].copy()

        if len(modified) == 0:
            fig = go.Figure()
            fig.add_annotation(
                text="No modifications detected",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig

        # Get taxonomic level index
        try:
            level_idx = TaxonomicLevel[level.upper()].value
        except KeyError:
            raise ValueError(f"Unknown taxonomic level: {level}")

        # Extract source and target taxa at specified level with role labels
        def get_taxon_at_level(path, idx):
            if len(path) > idx:
                return path[idx]
            return "Unclassified"

        modified['source_taxon'] = modified['original_path'].apply(lambda x: get_taxon_at_level(x, level_idx))
        modified['target_taxon'] = modified['proposed_path'].apply(lambda x: get_taxon_at_level(x, level_idx))

        # Filter out self-loops (no change at this level)
        flows = modified[modified['source_taxon'] != modified['target_taxon']].copy()

        if len(flows) == 0:
            # Create empty figure with message
            fig = go.Figure()
            fig.add_annotation(
                text=f"No taxonomic changes detected at {level.capitalize()} level",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            fig.update_layout(
                title=f"Taxonomic Flow at {level.capitalize()} Level",
                template='plotly_white'
            )
            return fig

        # Count flows between source and target taxa
        flow_counts = flows.groupby(['source_taxon', 'target_taxon']).size().reset_index(name='count')
        flow_counts = flow_counts.sort_values('count', ascending=False).head(top_n)

        # Build separate node lists for sources and targets
        unique_sources = sorted(set(flow_counts['source_taxon']))
        unique_targets = sorted(set(flow_counts['target_taxon']))

        # Create node labels with role indicators
        source_nodes = [f"{taxon}\n(Original)" for taxon in unique_sources]
        target_nodes = [f"{taxon}\n(Corrected)" for taxon in unique_targets]

        all_node_labels = source_nodes + target_nodes

        # Create mapping from taxon name to node index
        source_dict = {taxon: idx for idx, taxon in enumerate(unique_sources)}
        target_dict = {taxon: idx + len(unique_sources) for idx, taxon in enumerate(unique_targets)}

        # Prepare Sankey data
        source_indices = [source_dict[src] for src in flow_counts['source_taxon']]
        target_indices = [target_dict[tgt] for tgt in flow_counts['target_taxon']]
        values = flow_counts['count'].tolist()

        # Calculate flow statistics for each node
        source_outflows = {taxon: flow_counts[flow_counts['source_taxon'] == taxon]['count'].sum()
                          for taxon in unique_sources}
        target_inflows = {taxon: flow_counts[flow_counts['target_taxon'] == taxon]['count'].sum()
                         for taxon in unique_targets}

        # Color nodes: sources in red/orange spectrum, targets in green/blue spectrum
        node_colors = []

        # Color source nodes (original classifications) - red to orange
        for taxon in unique_sources:
            outflow = source_outflows[taxon]
            # More intense color for higher outflow
            intensity = min(0.3 + (outflow / max(source_outflows.values())) * 0.5, 0.8)
            node_colors.append(f'rgba(239, 68, 68, {intensity})')  # Red

        # Color target nodes (corrected classifications) - green to blue
        for taxon in unique_targets:
            inflow = target_inflows[taxon]
            # More intense color for higher inflow
            intensity = min(0.3 + (inflow / max(target_inflows.values())) * 0.5, 0.8)
            node_colors.append(f'rgba(16, 185, 129, {intensity})')  # Green

        # Create hover labels with statistics
        node_hover_labels = []
        for taxon in unique_sources:
            count = source_outflows[taxon]
            node_hover_labels.append(f"<b>{taxon}</b><br>Original classification<br>Sequences reclassified: {humanize_count(count)}")

        for taxon in unique_targets:
            count = target_inflows[taxon]
            node_hover_labels.append(f"<b>{taxon}</b><br>Corrected classification<br>Sequences received: {humanize_count(count)}")

        # Create link colors with gradient effect
        link_colors = []
        for src, tgt, val in zip(flow_counts['source_taxon'], flow_counts['target_taxon'], values):
            # More opaque for larger flows
            max_val = max(values)
            opacity = 0.2 + (val / max_val) * 0.4
            link_colors.append(f'rgba(37, 99, 235, {opacity})')

        # Create Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            arrangement='snap',
            node=dict(
                pad=20,
                thickness=25,
                line=dict(color="white", width=1),
                label=all_node_labels,
                color=node_colors,
                customdata=node_hover_labels,
                hovertemplate='%{customdata}<extra></extra>'
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=link_colors,
                customdata=[f"{src} → {tgt}" for src, tgt in zip(flow_counts['source_taxon'], flow_counts['target_taxon'])],
                hovertemplate='<b>%{customdata}</b><br>' +
                             'Sequences: %{value:,}<br>' +
                             '<extra></extra>'
            )
        )])

        fig.update_layout(
            title=dict(
                text=f"Taxonomic Reassignment Flow at {level.capitalize()} Level<br>" +
                     f"<sub>Showing top {len(flow_counts)} flows | " +
                     f"Total sequences reassigned: {humanize_count(sum(values))}</sub>",
                font=dict(size=16)
            ),
            font=dict(family="Arial, sans-serif", size=11),
            height=700,
            template='plotly_white',
            margin=dict(l=20, r=20, t=80, b=20)
        )

        return fig


    def plot_taxonomic_flow_interactive(self) -> go.Figure:
        """
        Creates an interactive Sankey diagram with dropdown to select taxonomic level.
        """
        levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        default_level = "genus" if "genus" in levels else levels[0]

        figures_data = []
        available_levels = []

        for level in levels:
            fig = self.plot_taxonomic_flow_sankey(level=level, top_n=20)

            # Extract the sankey data
            if len(fig.data) > 0 and hasattr(fig.data[0], 'link'):
                sankey_data = fig.data[0]
                figures_data.append({
                    'level': level,
                    'node': sankey_data.node,
                    'link': sankey_data.link
                })
                available_levels.append(level)

        if not figures_data:
            # No data available
            fig = go.Figure()
            fig.add_annotation(
                text="No taxonomic flows available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig

        # Create figure with first dataset
        if default_level in available_levels:
            default_idx = available_levels.index(default_level)
        else:
            default_idx = 0

        fig = go.Figure()
        for idx, dataset in enumerate(figures_data):
            fig.add_trace(go.Sankey(
                node=dataset['node'],
                link=dataset['link'],
                visible=(idx == default_idx)
            ))

        buttons = []
        for idx, dataset in enumerate(figures_data):
            visible = [False] * len(figures_data)
            visible[idx] = True
            buttons.append(dict(
                label=dataset['level'].capitalize(),
                method="update",
                args=[
                    {"visible": visible},
                    {"title": f"Taxonomic Reassignment Flow at {dataset['level'].capitalize()} Level"}
                ]
            ))

        fig.update_layout(
            title=f"Taxonomic Reassignment Flow at {figures_data[default_idx]['level'].capitalize()} Level",
            font=dict(family="Arial, sans-serif", size=12),
            height=650,
            template='plotly_white',
            updatemenus=[{
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 1.15,
                "y": 1.15,
                "xanchor": "left",
                "yanchor": "top"
            }]
        )

        return fig


# ==================== Report Generator ====================

class ReportGenerator:
    """Generates comprehensive HTML reports."""

    def __init__(self, pipeline):
        if not pipeline.analyzer or not pipeline.visualizer:
            raise ValueError("Incomplete pipeline: analyzer and visualizer required")
        self.pipeline = pipeline

    def generate(self, output_file: Path = "sativa_report.html"):
        """Generate a complete HTML report with analysis and visualizations."""
        analyzer = self.pipeline.analyzer
        viz = self.pipeline.visualizer

        total = len(analyzer.df)
        modified = analyzer.df['is_modified'].sum()
        perc = round(modified/total*100, 2) if total > 0 else 0
        conf_stats = analyzer.get_confidence_statistics()

        html_content = self._generate_html(total, modified, perc, conf_stats, viz)

        with open(output_file, "w", encoding='utf-8') as f:
            f.write(html_content)

        logger.info(f"✅ Report generated: {output_file}")

    def _generate_html(self, total, modified, perc, conf_stats, viz) -> str:
        """Generate the HTML content for the report."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Get top modification level for conclusions
        mods_by_level = self.pipeline.analyzer.get_modifications_by_level()
        top_level = mods_by_level.iloc[0]['level'] if len(mods_by_level) > 0 else "unknown"

        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SATIVA Analysis Report</title>
    <style>
        body {{
            font-family: 'Arial', 'Helvetica', sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            background-color: white;
            padding: 40px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #1e40af;
            border-bottom: 3px solid #2563eb;
            padding-bottom: 10px;
            margin-bottom: 30px;
        }}
        h2 {{
            color: #1e40af;
            margin-top: 40px;
            margin-bottom: 20px;
        }}
        h3 {{
            color: #374151;
            margin-top: 30px;
        }}
        .summary-box {{
            background-color: #eff6ff;
            border-left: 4px solid #2563eb;
            padding: 20px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .metric {{
            display: inline-block;
            margin: 10px 20px 10px 0;
        }}
        .metric-label {{
            font-weight: bold;
            color: #374151;
        }}
        .metric-value {{
            font-size: 1.2em;
            color: #2563eb;
        }}
        .interpretation {{
            background-color: #f0fdf4;
            border-left: 4px solid #10b981;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .timestamp {{
            color: #6b7280;
            font-size: 0.9em;
            text-align: right;
            margin-top: 40px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e5e7eb;
        }}
        th {{
            background-color: #f9fafb;
            font-weight: bold;
            color: #374151;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>SATIVA Taxonomic Correction Analysis Report</h1>

        <div class="summary-box">
            <h2>Executive Summary</h2>
            <div class="metric">
                <span class="metric-label">Total Sequences:</span>
                <span class="metric-value">{humanize_count(total)}</span>
            </div>
            <div class="metric">
                <span class="metric-label">Modified Sequences:</span>
                <span class="metric-value">{humanize_count(modified)} ({perc}%)</span>
            </div>
            <div class="metric">
                <span class="metric-label">Mean Confidence:</span>
                <span class="metric-value">{conf_stats['mean']:.3f}</span>
            </div>
        </div>

        <h2>1. Introduction</h2>
        <p>This report presents a comprehensive analysis of taxonomic corrections performed by SATIVA
        (Semi-Automatic Taxonomy Improvement and Validation Algorithm). The analysis examines
        {humanize_count(total)} sequences, of which {humanize_count(modified)} ({perc}%) were flagged for
        taxonomic modifications.</p>

        <h2>2. Modification Rate Analysis</h2>
        <p>The following interactive visualization displays the percentage of sequences modified by SATIVA
        across different taxonomic levels or technical categories (e.g., barcode, phylum, genus).
        Use the dropdown menu to select different grouping criteria. High modification rates in specific
        categories may indicate systematic taxonomic inconsistencies or initial annotation quality issues.</p>

        {viz.plot_modification_rate().to_html(full_html=False, include_plotlyjs='cdn')}

        <div class="interpretation">
            <h3>Interpretation:</h3>
            <p>Modification rates vary across categories, reflecting differences in taxonomic stability
            and annotation reliability. Categories with higher modification rates warrant closer examination
            for potential systematic biases or data quality concerns.</p>
        </div>

        <h2>3. Taxonomic Level Distribution</h2>
        <p>This analysis reveals the taxonomic levels at which modifications predominantly occur.
        Understanding these patterns helps identify hierarchical stages where taxonomic uncertainties
        or errors are most prevalent. Percentages are calculated relative to the total number of
        modified sequences.</p>

        {viz.plot_modifications_by_level().to_html(full_html=False, include_plotlyjs=False)}

        <h2>4. Modification Type Distribution</h2>
        <p>Modifications are categorized into four types: <strong>Completion</strong> (adding taxonomic
        ranks), <strong>Deletion (Simple)</strong>, <strong>Deletion (Multiple)</strong>, and <strong>Substitution</strong>
        (replacing existing classifications). This distribution provides insights into the nature of taxonomic corrections.</p>

        <div class="interpretation">
            <h3>Deletion categories:</h3>
            <ul>
                <li><strong>Deletion (Simple)</strong>: SATIVA removes exactly one taxonomic rank compared to the original lineage.</li>
                <li><strong>Deletion (Multiple)</strong>: SATIVA removes two or more consecutive ranks, indicating deeper pruning of the lineage.</li>
            </ul>
            <p>Simple deletions typically reflect minor adjustments (e.g., removing a species label), whereas multiple deletions
            point to broader taxonomic inconsistencies in the original annotation.</p>
        </div>

        {viz.plot_modifications_by_type().to_html(full_html=False, include_plotlyjs=False)}

        <h2>5. Taxonomic Reassignment Flows</h2>
        <p>This Sankey diagram visualizes the flow of taxonomic reassignments at different hierarchical
        levels. Each flow represents sequences that were reclassified from one taxon (left) to another
        (right). The width of each flow is proportional to the number of sequences reassigned.</p>

        <div class="interpretation">
            <h3>Color Coding:</h3>
            <ul>
                <li><strong style="color: #ef4444;">Red nodes</strong>: Original taxonomic assignments (source taxa)</li>
                <li><strong style="color: #10b981;">Green nodes</strong>: Corrected taxonomic assignments (target taxa)</li>
                <li>Color intensity indicates the number of sequences involved: darker colors represent more reclassifications</li>
            </ul>
        </div>


        {viz.plot_taxonomic_flow_interactive().to_html(full_html=False, include_plotlyjs=False)}

        <div class="interpretation">
            <h3>Insights:</h3>
            <p>Hotspots of reassignment reveal systematic taxonomic issues. Large flows between specific
            taxa may indicate historical misclassifications, nomenclatural changes, or genuinely ambiguous
            taxonomic boundaries requiring expert revision.</p>
        </div>

        <h2>6. Confidence Score Analysis</h2>
        <p>SATIVA assigns confidence scores to each proposed modification. The distribution and
        statistics of these scores indicate the algorithm's certainty in its corrections.</p>

        {viz.plot_confidence_distribution().to_html(full_html=False, include_plotlyjs=False)}

        <h3>Confidence Statistics</h3>
        <table>
            <tr>
                <th>Metric</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>Mean</td>
                <td>{conf_stats['mean']:.3f}</td>
            </tr>
            <tr>
                <td>Median</td>
                <td>{conf_stats['median']:.3f}</td>
            </tr>
            <tr>
                <td>Standard Deviation</td>
                <td>{conf_stats['std']:.3f}</td>
            </tr>
            <tr>
                <td>Minimum</td>
                <td>{conf_stats['min']:.3f}</td>
            </tr>
            <tr>
                <td>Maximum</td>
                <td>{conf_stats['max']:.3f}</td>
            </tr>
        </table>

        <h2>7. Conclusions</h2>
        <p>This analysis demonstrates the utility of SATIVA in identifying and correcting taxonomic
        misclassifications. The {perc}% modification rate suggests {'substantial' if perc > 10 else 'moderate' if perc > 5 else 'limited'}
        taxonomic refinement across the dataset. Modifications predominantly occur at the
        <strong>{top_level}</strong> level,
        indicating this is a critical point for taxonomic validation.</p>

        <div class="interpretation">
            <h3>Recommendations:</h3>
            <ul>
                <li>Review high-confidence modifications for incorporation into the reference database</li>
                <li>Investigate categories with elevated modification rates for systematic issues</li>
                <li>Consider manual curation for low-confidence corrections</li>
                <li>Implement quality control measures at taxonomic levels with frequent modifications</li>
            </ul>
        </div>

        <div class="timestamp">
            Report generated on: {timestamp}
        </div>
    </div>
</body>
</html>
"""
        return html


# ==================== Pipeline ====================

class SativaPipeline:
    """Orchestrates the complete SATIVA analysis workflow."""

    def __init__(self, root_dir: Optional[Path] = None, db_engine=None):
        self.root_dir = Path(root_dir) if root_dir else None
        self.engine = db_engine
        self.SessionFactory = sessionmaker(bind=db_engine) if db_engine else None
        self.all_results: List[SequenceResult] = []
        self.modifications: List[TaxonomicModification] = []
        self.analyzer: Optional[TaxonomyAnalyzer] = None
        self.visualizer: Optional[VisualizationEngine] = None

    def load_data(self):
        """Load data from the configured source (database preferred)."""
        if self.engine:
            self._load_data_from_db()
        elif self.root_dir:
            self._load_data_from_files()
        else:
            raise ValueError("No database connection or SATIVA directory provided.")

    def _load_data_from_files(self):
        """Load and parse all SATIVA output files."""
        logger.info("📂 Loading data...")
        self.all_results = []

        for barcode_dir in self.root_dir.iterdir():
            if not barcode_dir.is_dir():
                continue
            barcode = barcode_dir.name

            for alignment_dir in barcode_dir.iterdir():
                if not alignment_dir.is_dir():
                    continue

                for msa_dir in alignment_dir.iterdir():
                    if not msa_dir.is_dir():
                        continue
                    msa_name = msa_dir.name

                    # Load taxonomies
                    tax_files = list(msa_dir.glob("*.tax"))
                    taxonomies = {}
                    if tax_files:
                        taxonomies = SativaParser.parse_tax_file(tax_files[0])

                    # Create SequenceResult objects
                    msa_results = {
                        seq_id: SequenceResult(seq_id, barcode, msa_name, original_path)
                        for seq_id, original_path in taxonomies.items()
                    }

                    # Load modifications
                    for mis_file in msa_dir.glob("*.mis"):
                        mods = SativaParser.parse_mis_file(mis_file, barcode, msa_name)
                        for mod in mods:
                            if mod.seq_id in msa_results:
                                msa_results[mod.seq_id].modification = mod

                    self.all_results.extend(msa_results.values())

        self.modifications = [r.modification for r in self.all_results if r.is_modified]
        logger.info(f"✅ Loaded {humanize_count(len(self.all_results))} sequences, "
              f"{humanize_count(len(self.modifications))} modified")

    def _load_data_from_db(self):
        """Load SATIVA analysis data directly from the database."""
        if not self.engine:
            raise ValueError("Database engine is not configured.")

        if not self.SessionFactory:
            self.SessionFactory = sessionmaker(bind=self.engine)

        session = self.SessionFactory()
        try:
            sativa_rows = session.query(
                SeqTax.seq_id.label("seq_id"),
                DBSequence.barcode.label("barcode"),
                Source.name.label("source_name"),
                SeqTax.confidence.label("confidence"),
                Taxonomy.kingdom.label("kingdom"),
                Taxonomy.phylum.label("phylum"),
                Taxonomy.class_.label("class_"),
                Taxonomy.order_.label("order_"),
                Taxonomy.family.label("family"),
                Taxonomy.genus.label("genus"),
                Taxonomy.species.label("species")
            ).join(DBSequence, SeqTax.seq_id == DBSequence.seq_id
            ).join(Source, DBSequence.source_id == Source.source_id
            ).join(Taxonomy, SeqTax.tax_id == Taxonomy.tax_id
            ).filter(SeqTax.backbone == 'Sativa').all()

            if not sativa_rows:
                logger.warning("No Sativa assignments found in the database.")

            original_rows = session.query(
                SeqTax.seq_id.label("seq_id"),
                SeqTax.backbone.label("backbone"),
                Taxonomy.kingdom.label("kingdom"),
                Taxonomy.phylum.label("phylum"),
                Taxonomy.class_.label("class_"),
                Taxonomy.order_.label("order_"),
                Taxonomy.family.label("family"),
                Taxonomy.genus.label("genus"),
                Taxonomy.species.label("species")
            ).join(Taxonomy, SeqTax.tax_id == Taxonomy.tax_id
            ).filter(SeqTax.backbone.in_(["CoL", "Raw"])).all()

            priority = {'CoL': 0, 'Raw': 1}
            original_map: Dict[str, List[str]] = {}
            for row in sorted(original_rows, key=lambda r: priority.get(r.backbone, 99)):
                if row.seq_id not in original_map:
                    original_map[row.seq_id] = self._build_taxonomy_path(row)

            seq_results: List[SequenceResult] = []
            for row in sativa_rows:
                proposed_path = self._build_taxonomy_path(row)
                original_path = original_map.get(row.seq_id, [])
                seq_results.append(
                    self._create_sequence_result(
                        seq_id=row.seq_id,
                        barcode=row.barcode,
                        source=row.source_name,
                        original_path=original_path,
                        proposed_path=proposed_path,
                        confidence=row.confidence if row.confidence is not None else 0.0
                    )
                )

            processed_ids = {res.seq_id for res in seq_results}

            warning_rows = session.query(
                SativaWarning.seq_id,
                SativaWarning.barcode,
                SativaWarning.sativa_original_taxonomy,
                SativaWarning.sativa_proposed_taxonomy,
                SativaWarning.sativa_confidence,
                Source.name.label("source_name")
            ).join(DBSequence, SativaWarning.seq_id == DBSequence.seq_id
            ).join(Source, DBSequence.source_id == Source.source_id, isouter=True).all()
            for warn in warning_rows:
                if warn.seq_id in processed_ids:
                    continue
                seq_results.append(
                    self._create_sequence_result(
                        seq_id=warn.seq_id,
                        barcode=warn.barcode,
                        source=warn.source_name,
                        original_path=self._parse_taxonomy_string(warn.sativa_original_taxonomy),
                        proposed_path=self._parse_taxonomy_string(warn.sativa_proposed_taxonomy),
                        confidence=warn.sativa_confidence or 0.0
                    )
                )

            self.all_results = seq_results
            self.modifications = [r.modification for r in self.all_results if r.is_modified]
            logger.info(
                f"✅ Loaded {humanize_count(len(self.all_results))} sequences from database, "
                f"{humanize_count(len(self.modifications))} modified"
            )
        finally:
            session.close()

    @staticmethod
    def _build_taxonomy_path(row) -> List[str]:
        """Build a taxonomy path list from ORM row fields."""
        return [
            value for value in [
                row.kingdom,
                row.phylum,
                getattr(row, 'class_', None),
                getattr(row, 'order_', None),
                row.family,
                row.genus,
                getattr(row, 'species', None)
            ] if value
        ]

    @staticmethod
    def _parse_taxonomy_string(value: Optional[str]) -> List[str]:
        """Parse semicolon-delimited taxonomy strings."""
        if not value:
            return []
        return [part for part in value.split(';') if part]

    def _create_sequence_result(
        self,
        seq_id: str,
        barcode: Optional[str],
        source: Optional[str],
        original_path: Optional[List[str]],
        proposed_path: Optional[List[str]],
        confidence: float
    ) -> SequenceResult:
        """Create a SequenceResult object from taxonomy paths."""
        orig = [v for v in (original_path or []) if v]
        prop = [v for v in (proposed_path or []) if v]

        result = SequenceResult(
            seq_id=seq_id,
            barcode=barcode or "Unknown",
            msa_file="database",
            original_path=orig,
            source=source or "Unknown"
        )

        if orig == prop:
            return result

        mod_level, mod_type = SativaParser._analyze_modification(orig, prop)
        level_idx = mod_level.value
        orig_label = orig[level_idx] if len(orig) > level_idx and orig else (orig[-1] if orig else "Unassigned")
        prop_label = prop[level_idx] if len(prop) > level_idx and prop else (prop[-1] if prop else "Unassigned")

        result.modification = TaxonomicModification(
            seq_id=seq_id,
            barcode=barcode or "Unknown",
            msa_file="database",
            original_label=orig_label,
            proposed_label=prop_label,
            confidence=confidence or 0.0,
            original_path=orig,
            proposed_path=prop,
            modification_level=mod_level,
            modification_type=mod_type,
            per_rank_confidence=[]
        )
        return result

    def analyze(self):
        """Perform statistical analysis on loaded data."""
        logger.info("📊 Analyzing data...")
        self.analyzer = TaxonomyAnalyzer(self.all_results)
        logger.info("✅ Analysis complete")

    def visualize(self):
        """Initialize visualization engine."""
        if not self.analyzer:
            raise ValueError("Must run analyze() before visualize()")
        self.visualizer = VisualizationEngine(self.analyzer)
        logger.info("✅ Visualizations prepared")

    def report(self, output_file: Path = Path("sativa_report.html")):
        """Generate comprehensive HTML report."""
        if not self.visualizer:
            raise ValueError("Must run visualize() before report()")
        generator = ReportGenerator(self)
        generator.generate(output_file)


# ==================== Main Execution ====================

def main():
    """Main function for Snakemake or standalone execution."""
    try:
        engine = None
        root_dir = None

        # Configuration from Snakemake
        if 'snakemake' in globals():
            output_file = Path(snakemake.output.report)
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            params = getattr(snakemake, 'params', {})
            config = getattr(snakemake, 'config', {})
            root_dir_param = params.get('sativa_output_dir')
            root_dir = Path(root_dir_param) if root_dir_param else None

            if config.get('database'):
                db_manager = get_database_manager(config)
                engine = db_manager.engine

            # Configure file logging
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s | %(levelname)s | %(name)s | %(message)s'
            ))
            logger.addHandler(file_handler)
        else:
            # Standalone mode
            import argparse
            parser = argparse.ArgumentParser(description="SATIVA Results Visualization")
            parser.add_argument("--sativa-dir", type=str, help="SATIVA output directory (legacy mode)")
            parser.add_argument("--database", type=str, default="data/db/its.db", help="Path to SQLite database file")
            parser.add_argument("--output", type=str, default="sativa_visualization_report.html", help="Output HTML file")
            args = parser.parse_args()

            output_file = Path(args.output)
            if args.database:
                db_path = Path(args.database)
                if not db_path.exists():
                    logger.error(f"Database not found: {db_path}")
                    sys.exit(1)
                engine = get_engine("sqlite", path=str(db_path))
            root_dir = Path(args.sativa_dir) if args.sativa_dir else None

        if not engine and not root_dir:
            logger.error("Please provide either a database configuration or a SATIVA output directory.")
            sys.exit(1)

        # Execute pipeline
        pipeline = SativaPipeline(root_dir=root_dir, db_engine=engine)
        pipeline.load_data()
        pipeline.analyze()
        pipeline.visualize()
        pipeline.report(output_file)

        logger.info(f"✅ Visualization report generated successfully: {output_file}")

    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
