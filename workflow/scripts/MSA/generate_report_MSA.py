#!/usr/bin/env python3
"""
Script to analyze MSA alignment metrics and generate an HTML report.
Designed for integration inside a Snakemake workflow.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import glob
import base64
from io import BytesIO
from datetime import datetime
import logging
import sys

# Snakemake-friendly logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# --- Analysis helpers (kept identical) ---
# The functions below were originally stable and remain intact.
def get_msa_metrics(msa_file, alignment_type=None, barcode=None):
    """Extract alignment metrics from an MSA file along with taxonomy metadata."""
    try:
        alignment = AlignIO.read(msa_file, 'fasta')
        
        result = {
            'file': os.path.basename(msa_file),
            'num_sequences': len(alignment),
            'alignment_length': alignment.get_alignment_length(),
            'barcode': barcode or os.path.basename(os.path.dirname(os.path.dirname(msa_file)))
        }
        
        if alignment_type:
            result['alignment_type'] = alignment_type
            if alignment_type == 'alignments_species':
                result['taxonomic_level'] = 'species'
            elif alignment_type == 'alignments_genus':
                result['taxonomic_level'] = 'genus'
            elif alignment_type == 'alignments_family':
                result['taxonomic_level'] = 'family'
            elif alignment_type == 'alignments_order':
                result['taxonomic_level'] = 'order'
        
        return result
    except Exception as e:
        logger.error(f"Error while parsing {msa_file}: {str(e)}")
        return None

def process_msa_file_with_type(msa_file, alignment_type, barcode):
    """Wrapper that processes an MSA file and preserves the alignment type."""
    try:
        result = get_msa_metrics(msa_file, alignment_type, barcode)
        return result
    except Exception as e:
        logger.error(f"process_msa_file_with_type failed for {msa_file}: {str(e)}")
        return None

def analyze_msa(msa_dir, barcodes, max_workers=10):
    """Analyze MSAs with multi-level taxonomic support."""
    logger.info("Analyzing MSA files...")

    alignment_types = ['alignments_species', 'alignments_genus', 'alignments_family', 'alignments_order']
    all_msa_files = []
    
    for barcode in barcodes:
        for alignment_type in alignment_types:
            msa_path = os.path.join(msa_dir, barcode, alignment_type)
            if os.path.exists(msa_path):
                msa_files = glob.glob(os.path.join(msa_path, '*_msa.fasta'))
                for msa_file in msa_files:
                    all_msa_files.append((msa_file, alignment_type, barcode))
    
    logger.info(f"{len(all_msa_files)} MSA files detected in total")
    
    for alignment_type in alignment_types:
        count = sum(1 for _, atype, _ in all_msa_files if atype == alignment_type)
        logger.info(f"  - {alignment_type}: {count} files")
    
    msa_metrics = []
    max_workers = max_workers or min(32, (os.cpu_count() or 1) + 4)
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_info = {
            executor.submit(process_msa_file_with_type, msa_file, alignment_type, barcode): (msa_file, alignment_type, barcode)
            for msa_file, alignment_type, barcode in all_msa_files
        }
        
        for future in tqdm(as_completed(future_to_info), 
                         total=len(future_to_info), 
                         desc="Processing MSA files"):
            try:
                result = future.result()
                if result:
                    msa_metrics.append(result)
            except Exception as e:
                msa_file, alignment_type, barcode = future_to_info[future]
                logger.error(f"Error while processing {os.path.basename(msa_file)} ({alignment_type}): {str(e)}")
    
    df_msa = pd.DataFrame(msa_metrics)
    logger.info(f"\nMetrics extracted for {len(df_msa)} alignments")
    
    if len(df_msa) > 0:
        logger.info("\nAlignment counts per type:")
        alignment_counts = df_msa['alignment_type'].value_counts()
        for alignment_type, count in alignment_counts.items():
            logger.info(f"  - {alignment_type}: {count} alignments analyzed")
    
    return df_msa

def analyze_unassigned_sequences(msa_dir, barcodes, filename):
    """Generic function that aggregates unassigned sequence information."""
    file_path = os.path.join(msa_dir, filename)
    if not os.path.exists(file_path):
        logger.info(f"Unassigned sequence file '{filename}' not found.")
        return pd.DataFrame()
        
    try:
        df = pd.read_csv(file_path)
        logger.info(f"Found {len(df)} unassigned sequences in {filename}")
        
        data = [{'unassigned_count': len(df), 'file': file_path}]
        if 'unassigned_sequences.csv' in filename:
            # Intra-barcode case, append barcode column
            data = []
            for barcode in barcodes:
                barcode_file = os.path.join(msa_dir, barcode, filename)
                if os.path.exists(barcode_file):
                    df_barcode = pd.read_csv(barcode_file)
                    data.append({
                        'barcode': barcode,
                        'unassigned_count': len(df_barcode),
                        'file': barcode_file
                    })
        
        return pd.DataFrame(data)
    except Exception as e:
        logger.error(f"Error while processing {file_path}: {str(e)}")
        return pd.DataFrame()

def plot_jointplot(df):
    # ... (function unchanged) ...
    hue_col = 'taxonomic_level' if 'taxonomic_level' in df.columns else 'barcode'
    palette = sns.color_palette("Set1", n_colors=df[hue_col].nunique())

    plt.figure()
    jp = sns.jointplot(
        data=df,
        x='alignment_length',
        y='num_sequences',
        hue=hue_col,
        palette=palette,
        kind='kde',
        fill=True,
        height=8,
        alpha=0.5
    )
    jp.fig.suptitle("Joint distribution: alignment length vs number of sequences", y=1.02)
    
    # Force axes to start at zero for readability
    if jp.ax_joint:
        jp.ax_joint.set_xlim(left=0)
        jp.ax_joint.set_ylim(bottom=0)
    if hasattr(jp, "ax_marg_x") and jp.ax_marg_x:
        jp.ax_marg_x.set_xlim(left=0)
    if hasattr(jp, "ax_marg_y") and jp.ax_marg_y:
        jp.ax_marg_y.set_ylim(bottom=0)

    return jp.fig

def plot_violinplot(df):
    # ... (function unchanged) ...
    plt.figure(figsize=(14, 6))
    
    if 'taxonomic_level' in df.columns:
        df_plot = df.copy()
        df_plot['barcode_level'] = df_plot['barcode'] + '_' + df_plot['taxonomic_level']
        palette = sns.color_palette("Set2", n_colors=df_plot['barcode_level'].nunique())
        
        ax = sns.violinplot(data=df_plot, x='barcode_level', y='alignment_length', 
                           inner='quartile', palette=palette)
        plt.title("Alignment length distribution by barcode and taxonomic level")
    else:
        palette = sns.color_palette("Set1", n_colors=df['barcode'].nunique())
        ax = sns.violinplot(data=df, x='barcode', y='alignment_length', 
                           inner='quartile', palette=palette)
        plt.title("Alignment length distribution by barcode")
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    return plt.gcf()

def plot_taxonomic_comparison(df):
    # ... (function unchanged) ...
    if 'taxonomic_level' not in df.columns:
        return None
        
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    sns.boxplot(data=df, x='taxonomic_level', y='alignment_length', ax=axes[0,0])
    axes[0,0].set_title('Alignment length distribution per taxonomic level')
    
    sns.boxplot(data=df, x='taxonomic_level', y='num_sequences', ax=axes[0,1])
    axes[0,1].set_title('Sequence count distribution per taxonomic level')
    
    count_data = df.groupby(['barcode', 'taxonomic_level']).size().reset_index(name='count')
    sns.barplot(data=count_data, x='taxonomic_level', y='count', hue='barcode', ax=axes[1,0])
    axes[1,0].set_title('Alignments per taxonomic level and barcode')
    
    for i, level in enumerate(df['taxonomic_level'].unique()):
        level_data = df[df['taxonomic_level'] == level]
        axes[1,1].scatter(level_data['alignment_length'], level_data['num_sequences'], 
                         label=level, alpha=0.6)
    axes[1,1].set_xlabel('Alignment length')
    axes[1,1].set_ylabel('Number of sequences')
    axes[1,1].set_title('Correlation by taxonomic level')
    axes[1,1].legend()
    
    plt.tight_layout()
    return fig

def create_html_report(df_msa, df_unassigned, df_cross_unassigned, output_file, barcodes):
    """Create an HTML report from the analyzed data with multi-level support."""
    # Generate the figures
    jointplot_fig = plot_jointplot(df_msa)
    violinplot_fig = plot_violinplot(df_msa)
    
    taxonomic_fig = None
    taxonomic_img = ""
    if 'taxonomic_level' in df_msa.columns:
        taxonomic_fig = plot_taxonomic_comparison(df_msa)
    
    def fig_to_base64(fig):
        if fig is None:
            return ""
        buf = BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode('utf-8')
    
    jointplot_img = fig_to_base64(jointplot_fig)
    violinplot_img = fig_to_base64(violinplot_fig)
    if taxonomic_fig:
        taxonomic_img = fig_to_base64(taxonomic_fig)
    
    if 'taxonomic_level' in df_msa.columns:
        stats = df_msa.groupby(['barcode', 'taxonomic_level']).agg({
            'file': 'count',
            'num_sequences': ['mean', 'median', 'sum'],
            'alignment_length': ['mean', 'median']
        }).round(2)
    else:
        stats = df_msa.groupby('barcode').agg({
            'file': 'count',
            'num_sequences': ['mean', 'median', 'sum'],
            'alignment_length': ['mean', 'median']
        }).round(2)
    
    taxonomic_section = ""
    if 'taxonomic_level' in df_msa.columns:
        taxonomic_section = f"""
        <div class="section">
            <h2 class="section-title">Taxonomic level analysis</h2>
            <div class="plot">
                <h3>Multi-level taxonomic comparison</h3>
                <img src="data:image/png;base64,{taxonomic_img}" alt="Taxonomic Comparison">
            </div>
            <div class="subsection">
                <h3>Statistics by taxonomic level</h3>
                {df_msa.groupby('taxonomic_level').agg({'file': 'count', 'num_sequences': ['mean', 'sum'], 'alignment_length': 'mean'}).round(2).to_html(classes='stats-table')}
            </div>
        </div>
        """
    
    html_content = f"""
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
                <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2 class="section-title">Analysis summary</h2>
                <div class="summary">
                    <p>Total alignments analyzed: <strong>{len(df_msa):,}</strong></p>
                    <p>Barcodes analyzed: <strong>{len(barcodes)}</strong> ({', '.join(barcodes)})</p>
                    {f'<p>Taxonomic levels: <strong>{", ".join(df_msa["taxonomic_level"].unique())}</strong></p>' if 'taxonomic_level' in df_msa.columns else ''}
                    {f'<p>Unassigned sequences (intra-barcode): <strong>{df_unassigned["unassigned_count"].sum():,}</strong></p>' if not df_unassigned.empty and 'barcode' in df_unassigned.columns else ''}
                    {f'<p>Unassigned sequences (cross-barcode): <strong>{df_cross_unassigned["unassigned_count"].sum():,}</strong></p>' if not df_cross_unassigned.empty else ''}
                </div>
            </div>
            
            <div class="section">
                <h2 class="section-title">Statistics by barcode{' and taxonomic level' if 'taxonomic_level' in df_msa.columns else ''}</h2>
                {stats.to_html(classes='stats-table')}
            </div>
            
            {taxonomic_section}
            
            <div class="section">
                <h2 class="section-title">Visualizations</h2>
                <div class="plot">
                    <h3>Joint distribution: alignment length vs number of sequences</h3>
                    <img src="data:image/png;base64,{jointplot_img}" alt="Jointplot">
                </div>
                <div class="plot">
                    <h3>Alignment length distribution by barcode</h3>
                    <img src="data:image/png;base64,{violinplot_img}" alt="Violin Plot">
                </div>
            </div>
            
            {f'''
            <div class="section">
                <h2 class="section-title">Unassigned sequences after intra-barcode processing</h2>
                {df_unassigned[['barcode', 'unassigned_count']].to_html(index=False, classes='stats-table')}
            </div>
            ''' if not df_unassigned.empty and 'barcode' in df_unassigned.columns else ''}
            
            {f'''
            <div class="section">
                <h2 class="section-title">Unassigned sequences after cross-barcode processing</h2>
                <p>Total unassigned sequences: {df_cross_unassigned['unassigned_count'].sum():,}</p>
            </div>
            ''' if not df_cross_unassigned.empty else ''}
            
            <div class="section">
                <h2 class="section-title">Detailed metrics</h2>
                <div class="subsection">
                    <h3>Top 10 longest alignments</h3>
                    {df_msa.nlargest(10, 'alignment_length')[['file', 'barcode'] + (['taxonomic_level'] if 'taxonomic_level' in df_msa.columns else []) + ['alignment_length', 'num_sequences']].to_html(classes='stats-table', index=False)}
                </div>
                <div class="subsection" style="margin-top: 30px;">
                    <h3>Top 10 shortest alignments</h3>
                    {df_msa.nsmallest(10, 'alignment_length')[['file', 'barcode'] + (['taxonomic_level'] if 'taxonomic_level' in df_msa.columns else []) + ['alignment_length', 'num_sequences']].to_html(classes='stats-table', index=False)}
                </div>
                <div class="subsection" style="margin-top: 30px;">
                    <h3>Top 10 alignments by sequence count</h3>
                    {df_msa.nlargest(10, 'num_sequences')[['file', 'barcode'] + (['taxonomic_level'] if 'taxonomic_level' in df_msa.columns else []) + ['num_sequences', 'alignment_length']].to_html(classes='stats-table', index=False)}
                </div>
                <div class="subsection" style="margin-top: 30px;">
                    <h3>Lowest 10 alignments by sequence count</h3>
                    {df_msa.nsmallest(10, 'num_sequences')[['file', 'barcode'] + (['taxonomic_level'] if 'taxonomic_level' in df_msa.columns else []) + ['num_sequences', 'alignment_length']].to_html(classes='stats-table', index=False)}
                </div>
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logger.info(f"Report generated successfully: {output_file}")


def main():
    """Main Snakemake entry point."""
    
    # Configuration Snakemake
    if 'snakemake' not in globals():
        logger.critical("This script must be executed via Snakemake.")
        sys.exit(1)
    
    # Configuration du logging pour Snakemake
    if snakemake.log and isinstance(snakemake.log, list) and len(snakemake.log) > 0:
        log_file = snakemake.log[0]
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    logger.info("Starting MSA analysis.")
    
    # Retrieve Snakemake parameters
    msa_dir = snakemake.params.msa_dir
    barcodes = snakemake.params.barcodes
    output_report = snakemake.output.report
    
    # Analyze MSAs
    df_msa = analyze_msa(msa_dir, barcodes)
    
    # Analyze intra-barcode unassigned sequences
    df_unassigned = analyze_unassigned_sequences(msa_dir, barcodes, 'unassigned_sequences.csv')
    
    # Analyze cross-barcode unassigned sequences
    df_cross_unassigned = analyze_unassigned_sequences(msa_dir, barcodes, 'cross_unassigned_sequences.csv')
    
    # Generate the HTML report
    create_html_report(df_msa, df_unassigned, df_cross_unassigned, output_report, barcodes)
    
    logger.info("Analysis completed.")


if __name__ == "__main__":
    main()
