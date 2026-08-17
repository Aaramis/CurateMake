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

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent / "utils"))
from report_style import REPORT_CSS, CHART_SEQUENCE, human_duration  # shared CurateMake report styling

# One colour per pipeline module, from the shared palette.
MODULE_COLORS = {m: CHART_SEQUENCE[i % len(CHART_SEQUENCE)] for i, m in enumerate(
    ['database', 'format_databases', 'harmonize_taxo', 'itsx', 'msa', 'sativa', 'reports'])}
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
            # Benchmark files carry no timestamp; the file's mtime is when the job
            # finished, which gives a start once the duration is subtracted.
            'end_ts': filepath.stat().st_mtime,
            'start_ts': filepath.stat().st_mtime - row['s'],
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


def format_duration(seconds: float) -> str:
    """Seconds as h:mm:ss, matching the h:m:s column of a Snakemake benchmark file."""
    total = int(round(seconds))
    return f"{total // 3600}:{(total % 3600) // 60:02d}:{total % 60:02d}"


def collect_benchmarks(benchmark_dir: Path) -> Dict[str, Dict]:
    """Collect every benchmark file and organize them by module.

    Scattered rules write one benchmark per job into a sub-directory
    (benchmarks/sativa/integrate_tasks/<taxon>.benchmark.txt). Those jobs are the bulk
    of the pipeline's compute, so the walk is recursive and each family is folded into a
    single row: time and I/O summed over the jobs, memory taken as the peak of any one of
    them, and the job count kept so the row cannot be mistaken for a single execution.
    """
    benchmarks = {}

    for module_dir in benchmark_dir.iterdir():
        if not module_dir.is_dir():
            continue

        module_name = module_dir.name
        rules: Dict[str, Dict] = {}

        # Rules with one benchmark file, directly under the module directory.
        for bench_file in module_dir.glob('*.benchmark.txt'):
            data = parse_benchmark_file(bench_file)
            if data:
                data['n_jobs'] = 1
                data['wall_seconds'] = data['seconds']
                rules[bench_file.stem.replace('.benchmark', '')] = data

        # Scattered rules: one sub-directory per rule, one benchmark file per job.
        for family_dir in sorted(d for d in module_dir.iterdir() if d.is_dir()):
            jobs = [parse_benchmark_file(f)
                    for f in family_dir.rglob('*.benchmark.txt')]
            jobs = [j for j in jobs if j]
            if not jobs:
                continue
            agg = {
                'seconds': sum(j['seconds'] for j in jobs),
                'max_rss_mb': max(j['max_rss_mb'] for j in jobs),
                'max_uss_mb': max(j['max_uss_mb'] for j in jobs),
                'max_pss_mb': max(j['max_pss_mb'] for j in jobs),
                'io_in_mb': sum(j['io_in_mb'] for j in jobs),
                'io_out_mb': sum(j['io_out_mb'] for j in jobs),
                'cpu_time': sum(j['cpu_time'] for j in jobs),
                'n_jobs': len(jobs),
            }
            agg['mean_load'] = (sum(j['mean_load'] for j in jobs) / len(jobs))
            agg['time_formatted'] = format_duration(agg['seconds'])
            # Jobs of a scattered rule run concurrently: what a user waits for is the
            # span from the first start to the last finish, not the sum of durations.
            agg['start_ts'] = min(j['start_ts'] for j in jobs)
            agg['end_ts'] = max(j['end_ts'] for j in jobs)
            agg['wall_seconds'] = agg['end_ts'] - agg['start_ts']
            rules[f"{family_dir.name} ({len(jobs)} jobs)"] = agg

        benchmarks[module_name] = rules

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
                'Jobs': data.get('n_jobs', 1),
                'Start ts': data.get('start_ts'),
                'End ts': data.get('end_ts'),
                'Time (s)': data['seconds'],
                'Time': data['time_formatted'],
                'Wall-clock (s)': data.get('wall_seconds', data['seconds']),
                'Wall-clock': format_duration(data.get('wall_seconds', data['seconds'])),
                'Parallelism': round(data['seconds'] / data['wall_seconds'], 1)
                               if data.get('wall_seconds') else 1.0,
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
    """Real calendar timeline: each rule drawn on the dates it actually ran.

    Start and end come from the benchmark files themselves (mtime is the moment the file
    was written, i.e. when the job finished; subtracting the duration gives the start).
    Cross-checked against the timestamps Snakemake writes in its own log: the two agree
    to the second. A scattered rule spans from its first job's start to its last one's end,
    so overlapping bars are genuine parallelism and gaps are genuine pauses.
    """
    if 'Start ts' not in df or df['Start ts'].isna().all():
        return go.Figure()

    d = df.dropna(subset=['Start ts', 'End ts']).copy()
    # Epoch -> wall-clock in the machine's timezone. pd.to_datetime(unit='s') would
    # yield naive UTC, drawing every bar shifted by the UTC offset relative to the
    # timestamps in the logs the reader compares the chart against.
    d['Start'] = d['Start ts'].map(datetime.fromtimestamp)
    d['Finish'] = d['End ts'].map(datetime.fromtimestamp)
    d['Task'] = d['Module'] + '/' + d['Rule']
    d = d.sort_values('Start')

    # px.timeline needs a non-zero span to draw anything; sub-minute rules would vanish.
    _floor = pd.Timedelta(minutes=2)
    d.loc[d['Finish'] - d['Start'] < _floor, 'Finish'] = d['Start'] + _floor

    fig = px.timeline(
        d, x_start='Start', x_end='Finish', y='Task', color='Module',
        color_discrete_map=MODULE_COLORS,
        hover_data={'Wall-clock': True, 'Time': True, 'Jobs': True,
                    'Start': '|%Y-%m-%d %H:%M', 'Finish': '|%Y-%m-%d %H:%M'},
    )
    span_days = (d['Finish'].max() - d['Start'].min()).total_seconds() / 86400
    fig.update_yaxes(autorange="reversed", title=None)
    fig.update_xaxes(title="Date", tickformat="%d %b", dtick=86400000.0 * 2)
    fig.update_layout(
        title=(f"Pipeline timeline &#8212; {d['Start'].min():%d %b} to {d['Finish'].max():%d %b} "
               f"({span_days:.0f} days); overlapping bars ran in parallel"),
        height=max(420, len(d) * 26),
        showlegend=True, hovermode='closest', bargap=0.25,
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
    total_time_fmt = human_duration(total_time)
    # Cumulative compute counts every job; what the pipeline actually took is the span
    # from the earliest start to the latest finish (gaps between stages included, since
    # a benchmark file cannot tell a pause from a long job).
    _starts = [d['start_ts'] for m in benchmarks.values() for d in m.values() if d.get('start_ts')]
    _ends = [d['end_ts'] for m in benchmarks.values() for d in m.values() if d.get('end_ts')]
    elapsed_span = (max(_ends) - min(_starts)) if _starts and _ends else total_time
    elapsed_span_fmt = human_duration(elapsed_span)
    longest_step = df['Wall-clock (s)'].max() if 'Wall-clock (s)' in df else 0
    longest_step_rule = (df.loc[df['Wall-clock (s)'].idxmax(), 'Rule']
                         if 'Wall-clock (s)' in df and len(df) else '')
    longest_step_fmt = human_duration(longest_step)
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

    # Embed each figure in the report itself. The standalone files stay, but a report
    # that points iframes at them by relative path shows whatever sits next to it --
    # a stale chart beside fresh numbers as soon as it is moved or shared alone.
    # plotly.js is loaded once in <head>; a fragment that carried it would only define
    # the library at its own position in the document, leaving every figure above it blank.
    def _inline(fig):
        return fig.to_html(full_html=False, default_width='100%', default_height='700px',
                           include_plotlyjs=False, config={'responsive': True})

    fig_gantt_html = _inline(gantt_fig)
    fig_scatter_html = _inline(scatter_fig)
    fig_heatmap_html = _inline(heatmap_fig)
    fig_efficiency_html = _inline(efficiency_fig)
    fig_dist_html = _inline(dist_fig)
    fig_box_html = _inline(box_fig)
    fig_bottleneck_html = _inline(bottleneck_fig)
    fig_module_html = _inline(module_fig)

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
    <style>{REPORT_CSS}</style>
    <script src="https://cdn.plot.ly/plotly-3.3.0.min.js" charset="utf-8"></script>
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
            <div class="kpi-card primary">
                <div class="kpi-value">{elapsed_span_fmt}</div>
                <div class="kpi-label">Elapsed &mdash; first job to last</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{longest_step_fmt}</div>
                <div class="kpi-label">Longest single step<br><code>{longest_step_rule}</code></div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{total_time_fmt}</div>
                <div class="kpi-label">Cumulative compute<br>all {int(df['Jobs'].sum()):,} jobs summed</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{total_memory_peak/1024:.0f} GB</div>
                <div class="kpi-label">Peak memory, highest single job</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{total_io/1024:.0f} GB</div>
                <div class="kpi-label">Disk I/O, read and written</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{avg_cpu_efficiency/100:.1f}&#215;</div>
                <div class="kpi-label">Mean parallelism<br>threads busy per job</div>
            </div>
        </div>

        <!-- Executive Summary -->
        <div class="section">
            <h2>📋 Executive Summary</h2>
            <p>
                This report reviews the performance of <strong>{len(df)} rules</strong> spread across
                <strong>{df['Module'].nunique()} modules</strong> in the CurateMake pipeline. The total execution
                span is <strong>{elapsed_span_fmt}</strong> from the first job to the last, totalling <strong>{total_time_fmt}</strong> of cumulative compute across all jobs, with a peak memory consumption of
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
                {fig_module_html}
            </div>
        </div>

        <!-- Gantt Timeline -->
        <div class="section">
            <h2>📅 Execution timeline</h2>
            <p>
                Each rule is drawn on the dates it actually ran, so bars that overlap were
                genuinely running at the same time and gaps are real pauses between stages.
                order based on module dependencies.
            </p>
            <div class="warning-box">
                <strong>Source:</strong> start and end are read from the benchmark files themselves
                (file mtime = job end, minus its duration = start), cross-checked against Snakemake's
                own log timestamps.
            </div>
            <div class="chart-container">
                {fig_gantt_html}
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
                {fig_scatter_html}
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
                {fig_heatmap_html}
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
                {fig_efficiency_html}
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
                {fig_dist_html}
            </div>

            <h3>📊 Box plots by module</h3>
            <p>
                The box plots compare metric distributions between modules, highlighting both intra- and inter-module variability.
            </p>
            <div class="chart-container">
                {fig_box_html}
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
                {fig_bottleneck_html}
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
