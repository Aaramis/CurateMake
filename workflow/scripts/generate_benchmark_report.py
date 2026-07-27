#!/usr/bin/env python3
"""
Generate an interactive benchmarking report for the ITS pipeline.
Reads Snakemake benchmark.txt files and builds an advanced visualization dashboard.
"""

import sys
import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from pathlib import Path
from typing import Dict, List, Any, Tuple
from datetime import datetime, timedelta
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_benchmark_file(filepath: Path) -> Dict[str, Any]:
    """
    Parse a Snakemake benchmark file.

    TSV format with columns:
    s h:m:s max_rss max_vms max_uss max_pss io_in io_out mean_load cpu_time
    """
    try:
        df = pd.read_csv(filepath, sep='\t')

        # Use the last row in case Snakemake re-ran the same rule
        row = df.iloc[-1]

        return {
            'seconds': row['s'],
            'time_formatted': row['h:m:s'],
            'max_rss_mb': row['max_rss'],
            'max_vms_mb': row['max_vms'],
            'max_uss_mb': row['max_uss'],
            'max_pss_mb': row['max_pss'],
            'io_in_mb': row['io_in'],
            'io_out_mb': row['io_out'],
            'mean_load': row['mean_load'],
            'cpu_time': row.get('cpu_time', 0)
        }
    except Exception as e:
        logger.warning(f"Error while parsing {filepath}: {e}")
        return None


def collect_benchmarks(benchmark_dir: Path) -> Dict[str, Dict]:
    """Collect every benchmark file and organize them by module."""

    benchmarks = {}

    for module_dir in benchmark_dir.iterdir():
        if not module_dir.is_dir():
            continue

        module_name = module_dir.name
        benchmarks[module_name] = {}

        for bench_file in module_dir.glob('*.benchmark.txt'):
            rule_name = bench_file.stem.replace('.benchmark', '')

            data = parse_benchmark_file(bench_file)
            if data:
                benchmarks[module_name][rule_name] = data

    return benchmarks


def create_summary_stats(benchmarks: Dict) -> pd.DataFrame:
    """Build a summary DataFrame for every benchmark entry."""

    rows = []
    for module, rules in benchmarks.items():
        for rule, data in rules.items():
            total_io = data['io_in_mb'] + data['io_out_mb']
            cpu_efficiency = (data['cpu_time'] / data['seconds'] * 100) if data['seconds'] > 0 else 0

            rows.append({
                'Module': module,
                'Rule': rule,
                'Time (s)': data['seconds'],
                'Time': data['time_formatted'],
                'Memory (MB)': data['max_rss_mb'],
                'USS (MB)': data['max_uss_mb'],
                'PSS (MB)': data['max_pss_mb'],
                'I/O In (MB)': data['io_in_mb'],
                'I/O Out (MB)': data['io_out_mb'],
                'Total I/O (MB)': total_io,
                'CPU Load': data['mean_load'],
                'CPU Time (s)': data['cpu_time'],
                'CPU Efficiency (%)': cpu_efficiency
            })

    return pd.DataFrame(rows)


def create_gantt_timeline(df: pd.DataFrame) -> go.Figure:
    """
    Create a Gantt-style chart that simulates the execution order.
    Since real timestamps are unavailable, a sequential timeline is emulated.
    """

    # Build a logical execution order based on modules
    module_order = ['database', 'format_databases', 'harmonize_taxo', 'itsx', 'msa', 'sativa', 'reports']

    # Sort the DataFrame by module order
    df['Module_Order'] = df['Module'].map({m: i for i, m in enumerate(module_order)})
    df_sorted = df.sort_values(['Module_Order', 'Time (s)'], ascending=[True, False])

    # Calculate simulated cumulative timestamps
    start_time = datetime(2025, 1, 1, 0, 0, 0)
    current_time = start_time

    gantt_data = []
    for idx, row in df_sorted.iterrows():
        end_time = current_time + timedelta(seconds=row['Time (s)'])
        gantt_data.append({
            'Task': f"{row['Module']}/{row['Rule']}",
            'Start': current_time,
            'Finish': end_time,
            'Module': row['Module'],
            'Duration': row['Time'],
            'Memory': f"{row['Memory (MB)']:.0f} MB"
        })
        current_time = end_time

    gantt_df = pd.DataFrame(gantt_data)

    # Render the Gantt chart
    fig = px.timeline(
        gantt_df,
        x_start='Start',
        x_end='Finish',
        y='Task',
        color='Module',
        hover_data=['Duration', 'Memory'],
        title='Pipeline execution timeline (sequential simulation)'
    )

    fig.update_yaxes(autorange="reversed")
    fig.update_layout(
        height=max(400, len(gantt_df) * 25),
        xaxis_title='Cumulative time',
        showlegend=True,
        hovermode='closest'
    )

    return fig


def create_scatter_matrix(df: pd.DataFrame) -> go.Figure:
    """
    Create a scatter plot showing the relationship between time, memory, and I/O.
    Bubble size encodes total I/O volume.
    """

    # Normalize bubble sizes
    max_io = df['Total I/O (MB)'].max()
    df['Bubble_Size'] = (df['Total I/O (MB)'] / max_io * 50) + 5 if max_io > 0 else 10

    fig = go.Figure()

    for module in df['Module'].unique():
        module_data = df[df['Module'] == module]

        fig.add_trace(go.Scatter(
            x=module_data['Time (s)'],
            y=module_data['Memory (MB)'],
            mode='markers',
            name=module,
            marker=dict(
                size=module_data['Bubble_Size'],
                sizemode='diameter',
                opacity=0.7,
                line=dict(width=1, color='white')
            ),
            text=module_data['Rule'],
            customdata=module_data['Total I/O (MB)'],
            hovertemplate='<b>%{text}</b><br>' +
                         'Time: %{x:.1f}s<br>' +
                         'Memory: %{y:.0f} MB<br>' +
                         'Total I/O: %{customdata:.1f} MB<br>' +
                         '<extra></extra>'
        ))

    fig.update_layout(
        title='Multi-dimensional analysis: Time × Memory × I/O',
        xaxis_title='Execution time (seconds)',
        yaxis_title='Memory usage (MB)',
        hovermode='closest',
        height=600,
        showlegend=True,
        annotations=[
            dict(
                text='Bubble size = Total I/O',
                xref='paper',
                yref='paper',
                x=0.02,
                y=0.98,
                showarrow=False,
                font=dict(size=10, color='gray')
            )
        ]
    )

    return fig


def create_correlation_heatmap(df: pd.DataFrame) -> go.Figure:
    """
    Build a heatmap showing correlations between performance metrics.
    """

    # Select relevant numeric columns
    metrics = ['Time (s)', 'Memory (MB)', 'Total I/O (MB)', 'CPU Load', 'CPU Time (s)', 'CPU Efficiency (%)']
    corr_data = df[metrics].corr()

    fig = go.Figure(data=go.Heatmap(
        z=corr_data.values,
        x=metrics,
        y=metrics,
        colorscale='RdBu',
        zmid=0,
        text=corr_data.values,
        texttemplate='%{text:.2f}',
        textfont={"size": 10},
        colorbar=dict(title="Correlation")
    ))

    fig.update_layout(
        title='Correlation matrix of performance metrics',
        height=500,
        xaxis=dict(tickangle=-45)
    )

    return fig


def create_efficiency_analysis(df: pd.DataFrame) -> go.Figure:
    """
    Create a chart analyzing CPU efficiency and throughput.
    """

    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('CPU efficiency per rule', 'CPU load vs execution time'),
        vertical_spacing=0.15,
        row_heights=[0.5, 0.5]
    )

    # Chart 1: CPU efficiency
    df_sorted = df.sort_values('CPU Efficiency (%)', ascending=True)

    fig.add_trace(
        go.Bar(
            x=df_sorted['CPU Efficiency (%)'],
            y=df_sorted['Rule'],
            orientation='h',
            marker=dict(
                color=df_sorted['CPU Efficiency (%)'],
                colorscale='RdYlGn',
                showscale=False
            ),
            text=[f"{x:.1f}%" for x in df_sorted['CPU Efficiency (%)']],
            textposition='auto',
            customdata=df_sorted['Module'],
            hovertemplate='<b>%{y}</b><br>Module: %{customdata}<br>Efficiency: %{text}<br><extra></extra>'
        ),
        row=1, col=1
    )

    # Chart 2: CPU load vs time with regression
    fig.add_trace(
        go.Scatter(
            x=df['Time (s)'],
            y=df['CPU Load'],
            mode='markers',
            marker=dict(
                size=8,
                color=df['CPU Efficiency (%)'],
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(
                    title="CPU Eff. (%)",
                    y=0.25,
                    len=0.4
                )
            ),
            text=df['Rule'],
            customdata=df['Module'],
            hovertemplate='<b>%{text}</b><br>Module: %{customdata}<br>Time: %{x:.1f}s<br>Load: %{y:.1f}<br><extra></extra>'
        ),
        row=2, col=1
    )

    fig.update_xaxes(title_text="CPU efficiency (%)", row=1, col=1)
    fig.update_yaxes(title_text="Rule", row=1, col=1)
    fig.update_xaxes(title_text="Execution time (s)", row=2, col=1)
    fig.update_yaxes(title_text="Average CPU load", row=2, col=1)

    fig.update_layout(
        height=900,
        showlegend=False,
        title_text="CPU resource efficiency analysis"
    )

    return fig


def create_distribution_plots(df: pd.DataFrame) -> go.Figure:
    """
    Create distribution plots (histograms) for key metrics.
    """

    metrics = [
        ('Time (s)', 'Execution time (s)'),
        ('Memory (MB)', 'Memory (MB)'),
        ('Total I/O (MB)', 'Total I/O (MB)'),
        ('CPU Load', 'CPU load')
    ]

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[m[1] for m in metrics],
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]],
        vertical_spacing=0.12,
        horizontal_spacing=0.1
    )

    for idx, (metric, title) in enumerate(metrics):
        row = (idx // 2) + 1
        col = (idx % 2) + 1

        # Histogramme
        fig.add_trace(
            go.Histogram(
                x=df[metric],
                nbinsx=20,
                name=title,
                marker=dict(
                    color='rgba(102, 126, 234, 0.7)',
                    line=dict(color='white', width=1)
                ),
                showlegend=False
            ),
            row=row, col=col
        )

        # Add median marker
        median_val = df[metric].median()
        fig.add_vline(
            x=median_val,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Median: {median_val:.1f}",
            annotation_position="top",
            row=row, col=col
        )

        fig.update_xaxes(title_text=title, row=row, col=col)
        fig.update_yaxes(title_text="Frequency", row=row, col=col)

    fig.update_layout(
        height=700,
        title_text="Statistical distribution of performance metrics",
        showlegend=False
    )

    return fig


def create_box_plots(df: pd.DataFrame) -> go.Figure:
    """
    Create box plots comparing distributions by module.
    """

    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('Execution time', 'Memory usage', 'Total I/O'),
        horizontal_spacing=0.1
    )

    metrics = [
        ('Time (s)', 1),
        ('Memory (MB)', 2),
        ('Total I/O (MB)', 3)
    ]

    for metric, col in metrics:
        for module in sorted(df['Module'].unique()):
            module_data = df[df['Module'] == module]

            fig.add_trace(
                go.Box(
                    y=module_data[metric],
                    name=module,
                    boxmean='sd',
                    showlegend=(col == 1)
                ),
                row=1, col=col
            )

    fig.update_yaxes(title_text="Seconds", row=1, col=1)
    fig.update_yaxes(title_text="MB", row=1, col=2)
    fig.update_yaxes(title_text="MB", row=1, col=3)

    fig.update_layout(
        height=500,
        title_text="Metric distribution by module (box plots)",
        showlegend=True
    )

    return fig


def calculate_bottleneck_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute a composite bottleneck score based on normalized metrics.
    Score = weighted mean of (time_norm * 0.4 + memory_norm * 0.3 + io_norm * 0.3)
    """

    df_score = df.copy()

    # Normalize metrics to 0-100
    df_score['Time_Norm'] = (df_score['Time (s)'] / df_score['Time (s)'].max()) * 100
    df_score['Memory_Norm'] = (df_score['Memory (MB)'] / df_score['Memory (MB)'].max()) * 100
    df_score['IO_Norm'] = (df_score['Total I/O (MB)'] / df_score['Total I/O (MB)'].max()) * 100 if df_score['Total I/O (MB)'].max() > 0 else 0

    # Weighted composite score
    df_score['Bottleneck_Score'] = (
        df_score['Time_Norm'] * 0.4 +
        df_score['Memory_Norm'] * 0.3 +
        df_score['IO_Norm'] * 0.3
    )

    return df_score.sort_values('Bottleneck_Score', ascending=False)


def create_bottleneck_analysis(df: pd.DataFrame) -> Tuple[go.Figure, pd.DataFrame]:
    """
    Build a bottleneck analysis using the composite score.
    """

    df_bottleneck = calculate_bottleneck_score(df)
    top_bottlenecks = df_bottleneck.head(10)

    # Create a radar chart for the top bottlenecks
    fig = go.Figure()

    categories = ['Time<br>(norm)', 'Memory<br>(norm)', 'I/O<br>(norm)']

    for idx, row in top_bottlenecks.head(5).iterrows():
        fig.add_trace(go.Scatterpolar(
            r=[row['Time_Norm'], row['Memory_Norm'], row['IO_Norm']],
            theta=categories,
            fill='toself',
            name=f"{row['Module']}/{row['Rule']}"[:30],
            hovertemplate='<b>%{fullData.name}</b><br>%{theta}: %{r:.1f}<br><extra></extra>'
        ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 100]
            )
        ),
        showlegend=True,
        title='Top 5 bottlenecks - multi-criteria analysis',
        height=500
    )

    return fig, top_bottlenecks


def create_module_summary(df: pd.DataFrame) -> go.Figure:
    """Create a module-level summary comparing time and memory."""

    module_stats = df.groupby('Module').agg({
        'Time (s)': 'sum',
        'Memory (MB)': 'max',
        'Total I/O (MB)': 'sum',
        'CPU Efficiency (%)': 'mean'
    }).reset_index()

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('Total time per module', 'Peak memory & CPU efficiency'),
        specs=[[{"type": "bar"}, {"type": "bar"}]],
        horizontal_spacing=0.15
    )

    # Total time
    fig.add_trace(
        go.Bar(
            x=module_stats['Module'],
            y=module_stats['Time (s)'],
            name='Time (s)',
            marker=dict(color='steelblue'),
            text=[f"{x:.0f}s" for x in module_stats['Time (s)']],
            textposition='auto'
        ),
        row=1, col=1
    )

    # Memory and efficiency
    fig.add_trace(
        go.Bar(
            x=module_stats['Module'],
            y=module_stats['Memory (MB)'],
            name='Memory (MB)',
            marker=dict(color='indianred'),
            text=[f"{x:.0f} MB" for x in module_stats['Memory (MB)']],
            textposition='auto',
            yaxis='y2'
        ),
        row=1, col=2
    )

    fig.add_trace(
        go.Scatter(
            x=module_stats['Module'],
            y=module_stats['CPU Efficiency (%)'],
            name='CPU efficiency (%)',
            mode='lines+markers',
            marker=dict(size=10, color='green'),
            line=dict(width=2),
            yaxis='y3'
        ),
        row=1, col=2
    )

    fig.update_xaxes(title_text="Module", row=1, col=1)
    fig.update_yaxes(title_text="Time (s)", row=1, col=1)
    fig.update_xaxes(title_text="Module", tickangle=-45, row=1, col=2)
    fig.update_yaxes(title_text="Memory (MB)", row=1, col=2)

    fig.update_layout(
        height=500,
        showlegend=True,
        title_text="Module performance summary"
    )

    return fig


def generate_html_report(benchmarks: Dict, output_path: Path):
    """Generate the complete HTML report with advanced visualizations."""

    # Build the summary DataFrame
    df = create_summary_stats(benchmarks)

    if df.empty:
        logger.warning("No benchmark entries found!")
        return

    # Global statistics
    total_time = df['Time (s)'].sum()
    total_time_fmt = str(timedelta(seconds=int(total_time)))
    total_memory_peak = df['Memory (MB)'].max()
    total_io = df['Total I/O (MB)'].sum()
    avg_cpu_efficiency = df['CPU Efficiency (%)'].mean()

    # Create the visualizations
    logger.info("Generating visualizations...")

    # 1. Gantt timeline
    gantt_fig = create_gantt_timeline(df)
    gantt_standalone = output_path.parent / "benchmark_gantt.html"
    gantt_fig.write_html(str(gantt_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 2. Scatter matrix
    scatter_fig = create_scatter_matrix(df)
    scatter_standalone = output_path.parent / "benchmark_scatter.html"
    scatter_fig.write_html(str(scatter_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 3. Correlation heatmap
    heatmap_fig = create_correlation_heatmap(df)
    heatmap_standalone = output_path.parent / "benchmark_heatmap.html"
    heatmap_fig.write_html(str(heatmap_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 4. Efficiency analysis
    efficiency_fig = create_efficiency_analysis(df)
    efficiency_standalone = output_path.parent / "benchmark_efficiency.html"
    efficiency_fig.write_html(str(efficiency_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 5. Distribution plots
    dist_fig = create_distribution_plots(df)
    dist_standalone = output_path.parent / "benchmark_distributions.html"
    dist_fig.write_html(str(dist_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 6. Box plots
    box_fig = create_box_plots(df)
    box_standalone = output_path.parent / "benchmark_boxplots.html"
    box_fig.write_html(str(box_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 7. Bottleneck analysis
    bottleneck_fig, bottleneck_df = create_bottleneck_analysis(df)
    bottleneck_standalone = output_path.parent / "benchmark_bottlenecks.html"
    bottleneck_fig.write_html(str(bottleneck_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # 8. Module summary
    module_fig = create_module_summary(df)
    module_standalone = output_path.parent / "benchmark_modules.html"
    module_fig.write_html(str(module_standalone), include_plotlyjs='cdn', config={'responsive': True})

    # Top bottlenecks for the table
    top_bottlenecks_html = bottleneck_df.head(10)[
        ['Module', 'Rule', 'Time', 'Memory (MB)', 'Total I/O (MB)', 'Bottleneck_Score']
    ].to_html(index=False, classes='data-table', float_format=lambda x: f'{x:.1f}')

    # Build the HTML with a modern design
    html = f"""
<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Benchmark Report - Pipeline CurateMake</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
            color: #2c3e50;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}

        .header {{
            text-align: center;
            color: white;
            margin-bottom: 40px;
            padding: 40px 20px;
        }}

        .header h1 {{
            font-size: 3em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }}

        .header p {{
            font-size: 1.2em;
            opacity: 0.9;
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

        .kpi-card:hover {{
            transform: translateY(-5px);
        }}

        .kpi-icon {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .kpi-value {{
            font-size: 2.5em;
            font-weight: bold;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin: 10px 0;
        }}

        .kpi-label {{
            color: #7f8c8d;
            font-size: 0.95em;
            font-weight: 500;
        }}

        /* Section */
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

        .section h3 {{
            color: #34495e;
            margin: 25px 0 15px 0;
            font-size: 1.3em;
        }}

        .section p {{
            color: #555;
            margin-bottom: 15px;
            line-height: 1.8;
        }}

        /* Info Box */
        .info-box {{
            background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
            border-left: 4px solid #2196f3;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 8px;
        }}

        .info-box strong {{
            color: #1976d2;
        }}

        .warning-box {{
            background: linear-gradient(135deg, #fff3e0 0%, #ffe0b2 100%);
            border-left: 4px solid #ff9800;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 8px;
        }}

        .warning-box strong {{
            color: #f57c00;
        }}

        /* Chart Container */
        .chart-container {{
            margin: 30px 0;
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
        }}

        iframe {{
            border: none;
            border-radius: 10px;
        }}

        /* Tables */
        .table-wrapper {{
            width: 100%;
            overflow-x: auto; 
        }}

        .data-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            table-layout: fixed;   /* fixed-width columns prevent layout blow-up */
            font-size: 0.9em;      /* optional: slightly more compact */
        }}

        .data-table th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}

        .data-table td {{
            padding: 12px 15px;
            border-bottom: 1px solid #ecf0f1;
            word-wrap: break-word; /* prevents very long values from stretching the layout */
        }}

        .data-table tr:hover {{
            background-color: #f8f9fa;
        }}

        .data-table tr:last-child td {{
            border-bottom: none;
        }}

        /* Footer */
        .footer {{
            text-align: center;
            color: white;
            margin-top: 50px;
            padding: 20px;
            opacity: 0.8;
        }}




        /* Responsive */
        @media (max-width: 768px) {{
            .header h1 {{
                font-size: 2em;
            }}

            .kpi-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>⚡ Benchmark Performance Report</h1>
            <p>Detailed performance analysis of the CurateMake pipeline</p>
            <div class="timestamp">
                Generated on {datetime.now().strftime('%d/%m/%Y at %H:%M:%S')}
            </div>
        </div>

        <!-- KPI Cards -->
        <div class="kpi-grid">
            <div class="kpi-card">
                <div class="kpi-icon">⏱️</div>
                <div class="kpi-value">{total_time_fmt}</div>
                <div class="kpi-label">Total execution time</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">💾</div>
                <div class="kpi-value">{int(total_memory_peak)} MB</div>
                <div class="kpi-label">Peak memory usage</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">💿</div>
                <div class="kpi-value">{int(total_io)} MB</div>
                <div class="kpi-label">Total disk I/O</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-icon">🎯</div>
                <div class="kpi-value">{avg_cpu_efficiency:.1f}%</div>
                <div class="kpi-label">Average CPU efficiency</div>
            </div>
        </div>

        <!-- Executive Summary -->
        <div class="section">
            <h2>📋 Executive Summary</h2>
            <p>
                This report reviews the performance of <strong>{len(df)} rules</strong> spread across
                <strong>{df['Module'].nunique()} modules</strong> in the CurateMake pipeline. The total execution
                time is <strong>{total_time_fmt}</strong> with a peak memory consumption of
                <strong>{int(total_memory_peak)} MB</strong>.
            </p>

            <div class="info-box">
                <strong>💡 Note:</strong> The average CPU efficiency of {avg_cpu_efficiency:.1f}% indicates
                {'excellent' if avg_cpu_efficiency > 80 else 'good' if avg_cpu_efficiency > 50 else 'moderate'}
                processor utilization. The visualizations below highlight bottlenecks and optimization opportunities.
            </div>
        </div>

        <!-- Module Summary -->
        <div class="section">
            <h2>📊 Module overview</h2>
            <p>
                This section summarizes aggregated performance by pipeline module.
                Each module groups related Snakemake rules that perform focused tasks.
            </p>
            <div class="chart-container">
                <iframe src="benchmark_modules.html" width="100%" height="550"></iframe>
            </div>
        </div>

        <!-- Gantt Timeline -->
        <div class="section">
            <h2>📅 Execution timeline</h2>
            <p>
                This visualization presents the simulated execution order for pipeline rules.
                Although Snakemake can run tasks in parallel, the timeline illustrates a logical sequential
                order based on module dependencies.
            </p>
            <div class="warning-box">
                <strong>⚠️ Note:</strong> Timestamps are simulated because the benchmark files do not include
                real execution times. The timeline reflects each task duration but not the exact parallelism.
            </div>
            <div class="chart-container">
                <iframe src="benchmark_gantt.html" width="100%" height="800"></iframe>
            </div>
        </div>

        <!-- Multi-dimensional Analysis -->
        <div class="section">
            <h2>🔬 Multi-dimensional analysis</h2>
            <p>
                This bubble chart highlights the relationship between three performance axes:
                execution time (X axis), memory usage (Y axis), and I/O volume (bubble size).
                It exposes rules that are expensive across several resource categories at once.
            </p>
            <div class="info-box">
                <strong>💡 Interpretation:</strong> Bubbles in the upper-right quadrant with large sizes
                represent the most resource-intensive rules (time + memory + I/O) and are ideal optimization targets.
            </div>
            <div class="chart-container">
                <iframe src="benchmark_scatter.html" width="100%" height="650"></iframe>
            </div>
        </div>

        <!-- Correlation Analysis -->
        <div class="section">
            <h2>🔗 Correlation analysis</h2>
            <p>
                This heatmap shows the correlations between performance metrics.
                Values near +1 (red) indicate a strong positive correlation, whereas values near -1 (blue)
                show an inverse relationship.
            </p>
            <div class="info-box">
                <strong>💡 How to read:</strong> For instance, a strong correlation between "Time" and "Memory"
                suggests that long-running tasks also consume significant memory. A weak correlation between
                "I/O" and "Time" means I/O may not be the limiting factor.
            </div>
            <div class="chart-container">
                <iframe src="benchmark_heatmap.html" width="100%" height="550"></iframe>
            </div>
        </div>

        <!-- CPU Efficiency Analysis -->
        <div class="section">
            <h2>⚙️ CPU efficiency analysis</h2>
            <p>
                CPU efficiency measures the ratio between CPU time and wall-clock time. A score near 100% means the
                process fully utilizes a single CPU core. Values above 100% indicate effective multi-core usage,
                while lower values may signify I/O waits or underutilized resources.
            </p>
            <div class="warning-box">
                <strong>⚠️ Important:</strong> Low CPU efficiency is not always a problem.
                It can be normal for I/O-bound workloads, but for compute-heavy stages it may expose opportunities
                for parallelization or code optimization.
            </div>
            <div class="chart-container">
                <iframe src="benchmark_efficiency.html" width="100%" height="950"></iframe>
            </div>
        </div>

        <!-- Statistical Distributions -->
        <div class="section">
            <h2>📈 Statistical distributions</h2>
            <p>
                These histograms display the distribution of performance metrics across every pipeline rule.
                The red line marks the median value, helping to visualize variability and outliers.
            </p>
            <div class="chart-container">
                <iframe src="benchmark_distributions.html" width="100%" height="750"></iframe>
            </div>

            <h3>📊 Box plots by module</h3>
            <p>
                The box plots compare metric distributions between modules, highlighting both intra- and inter-module variability.
            </p>
            <div class="chart-container">
                <iframe src="benchmark_boxplots.html" width="100%" height="550"></iframe>
            </div>
        </div>

        <!-- Bottleneck Analysis -->
        <div class="section">
            <h2>🚨 Bottleneck analysis</h2>
            <p>
                This section identifies pipeline bottlenecks through a <strong>composite score</strong>
                that blends three dimensions:
            </p>
            <ul style="margin: 15px 0 15px 30px; line-height: 2;">
                <li><strong>Execution time (40%)</strong>: contribution to the total pipeline duration</li>
                <li><strong>Memory usage (30%)</strong>: pressure on system resources</li>
                <li><strong>I/O operations (30%)</strong>: disk/network throughput impact</li>
            </ul>
            <p>
                Rules with higher scores should be prioritized for optimization.
            </p>

            <div class="info-box">
                <strong>💡 Methodology:</strong> Each metric is normalized (0-100) then multiplied by its weight
                to compute the final score. A score of 100 would represent an extreme bottleneck on every axis.
            </div>

            <h3>🎯 Radar view of the top 5 bottlenecks</h3>
            <div class="chart-container">
                <iframe src="benchmark_bottlenecks.html" width="100%" height="550"></iframe>
            </div>

            <h3>📋 Top 10 identified bottlenecks</h3>
            {top_bottlenecks_html}
        </div>

        <!-- Detailed Data -->
        <div class="section">
            <h2>📝 Detailed metrics</h2>
            <p>
                Full performance table for every rule executed in the pipeline.
            </p>
            <div class="table-wrapper">
                {df.to_html(index=False, classes='data-table', float_format=lambda x: f'{x:.2f}')}
            </div>
        </div>


        <!-- Footer -->
        <div class="footer">
            <p><strong>Pipeline CurateMake</strong> - ITS Taxonomy Curation & Analysis</p>
            <p>Automatically generated by Snakemake | © 2025</p>
        </div>
    </div>
</body>
</html>
"""

    # Write the main report file
    output_path.write_text(html)
    logger.info(f"✅ HTML report generated: {output_path}")


def main():
    """Main entry point."""

    # Snakemake integration
    if 'snakemake' in globals():
        benchmark_dir = Path(snakemake.params.benchmarks_dir)
        output_report = Path(snakemake.output.report)
    else:
        # Standalone mode defaults
        benchmark_dir = Path('benchmarks')
        output_report = Path('reports/html/benchmark_report.html')

    logger.info(f"📊 Collecting benchmarks from: {benchmark_dir}")

    if not benchmark_dir.exists():
        logger.error(f"❌ Benchmark directory not found: {benchmark_dir}")
        logger.error("Run the pipeline first to generate benchmark files.")
        sys.exit(1)

    # Gather benchmark data
    benchmarks = collect_benchmarks(benchmark_dir)

    if not benchmarks:
        logger.warning("⚠️ No benchmark files found!")
        sys.exit(1)

    total_rules = sum(len(rules) for rules in benchmarks.values())
    logger.info(f"✅ Benchmarks collected: {total_rules} rules across {len(benchmarks)} modules")

    # Generate the report
    generate_html_report(benchmarks, output_report)

    logger.info("✅ Benchmark report generated successfully!")


if __name__ == "__main__":
    main()
