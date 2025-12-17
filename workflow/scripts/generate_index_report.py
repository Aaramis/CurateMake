#!/usr/bin/env python3
"""
Generate an HTML index for every pipeline report.
Creates a single landing page to navigate between the various analyses.
"""

import os
import logging
import sys
from datetime import datetime
from pathlib import Path

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def format_file_size(size_bytes):
    """Format a file size using human-friendly units."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def generate_index(reports_dir, output_file, available_reports):
    """
    Generate index.html with links to every available report.

    Args:
        reports_dir: Path to the reports/html directory
        output_file: Output path for index.html
        available_reports: List of available report paths
    """
    reports_dir = Path(reports_dir)
    reports_info = []

    # Report metadata
    report_metadata = {
        "sequence_flow_report.html": {
            "title": "Sequence Flow Report",
            "description": "Quantitative flow of sequences through formatting, deduplication, ITSx, and MSA stages",
            "icon": "🌊"
        },
        "msa_analysis_report.html": {
            "title": "MSA Analysis Report",
            "description": "Multi-sequence alignment analysis per barcode and taxonomic level (species, genus, family, order)",
            "icon": "📊"
        },
        "sativa_visualization_report.html": {
            "title": "SATIVA Visualization Report",
            "description": "Visualization of taxonomic classification results and misidentification detection",
            "icon": "🔬"
        },
        "benchmark_report.html": {
            "title": "Benchmark Report",
            "description": "Pipeline performance metrics (time, memory, CPU)",
            "icon": "⚡"
        },
        "taxonomy_sunburst.html": {
            "title": "Taxonomic Distribution",
            "description": "Interactive sunburst showing the taxonomic distribution of final non-duplicate sequences across kingdoms",
            "icon": "🌳"
        }
    }

    # Build catalog of existing reports
    for report_path in available_reports:
        report_file = Path(report_path)
        if report_file.exists() and report_file.name != "index.html":
            file_name = report_file.name
            metadata = report_metadata.get(
                file_name,
                {
                    "title": file_name.replace("_", " ").replace(".html", "").title(),
                    "description": "Analysis report",
                    "icon": "📄"
                }
            )

            reports_info.append({
                "filename": file_name,
                "title": metadata["title"],
                "description": metadata["description"],
                "icon": metadata["icon"],
                "size": format_file_size(report_file.stat().st_size),
                "modified": datetime.fromtimestamp(report_file.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
            })

    # Sort alphabetically
    reports_info.sort(key=lambda x: x["filename"])

    # Build HTML cards for each report
    report_cards = ""
    for report in reports_info:
        report_cards += f"""
        <div class="report-card">
            <div class="report-icon">{report['icon']}</div>
            <div class="report-content">
                <h2>{report['title']}</h2>
                <p class="report-description">{report['description']}</p>
                <div class="report-meta">
                    <span>📦 {report['size']}</span>
                    <span>🕒 {report['modified']}</span>
                </div>
                <a href="{report['filename']}" class="report-link">Open report →</a>
            </div>
        </div>
        """

    if not report_cards:
        report_cards = """
        <div class="no-reports">
            <p>⚠️ No reports are available yet.</p>
            <p>Run the pipeline to generate the analysis outputs.</p>
        </div>
        """

    # Template HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pipeline Reports - CurateMake</title>
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
        }}

        .container {{
            max-width: 1000px;
            margin: 0 auto;
        }}

        .header {{
            text-align: center;
            color: white;
            margin-bottom: 50px;
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

        .reports-grid {{
            display: grid;
            gap: 30px;
        }}

        .report-card {{
            background: white;
            border-radius: 15px;
            padding: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
            display: flex;
            gap: 20px;
            align-items: flex-start;
        }}

        .report-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 15px 40px rgba(0,0,0,0.3);
        }}

        .report-icon {{
            font-size: 3em;
            flex-shrink: 0;
        }}

        .report-content {{
            flex: 1;
        }}

        .report-card h2 {{
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 1.5em;
        }}

        .report-description {{
            color: #666;
            margin-bottom: 15px;
            line-height: 1.5;
        }}

        .report-meta {{
            display: flex;
            gap: 20px;
            font-size: 0.9em;
            color: #999;
            margin-bottom: 15px;
        }}

        .report-link {{
            display: inline-block;
            padding: 12px 30px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            text-decoration: none;
            border-radius: 25px;
            font-weight: 600;
            transition: all 0.3s ease;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
        }}

        .report-link:hover {{
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(102, 126, 234, 0.6);
        }}

        .no-reports {{
            background: white;
            border-radius: 15px;
            padding: 60px 30px;
            text-align: center;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
        }}

        .no-reports p {{
            color: #666;
            font-size: 1.1em;
            margin: 10px 0;
        }}

        .footer {{
            text-align: center;
            color: white;
            margin-top: 50px;
            opacity: 0.8;
            font-size: 0.9em;
        }}

        @media (max-width: 768px) {{
            .header h1 {{
                font-size: 2em;
            }}

            .report-card {{
                flex-direction: column;
                text-align: center;
            }}

            .report-meta {{
                flex-direction: column;
                gap: 5px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>📊 Pipeline Analysis Reports</h1>
            <p>CurateMake pipeline analysis dashboards</p>
            <div class="timestamp">
                Generated on {datetime.now().strftime('%d/%m/%Y at %H:%M:%S')}
            </div>
        </div>

        <div class="reports-grid">
            {report_cards}
        </div>

        <div class="footer">
            <p>Pipeline CurateMake - ITS Taxonomy Curation</p>
            <p>© 2025 - Automatically generated by Snakemake</p>
        </div>
    </div>
</body>
</html>
"""

    # Write the file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    logger.info(f"HTML index generated successfully: {output_file}")
    logger.info(f"Number of reports referenced: {len(reports_info)}")


def main():
    """Main Snakemake entry point."""

    if 'snakemake' not in globals():
        logger.error("This script must be executed through Snakemake.")
        sys.exit(1)

    # Optional file logging
    if hasattr(snakemake, 'log') and snakemake.log:
        log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

    # Retrieve parameters
    reports_dir = Path(snakemake.params.reports_dir)
    output_file = snakemake.output.index
    available_reports = snakemake.input

    logger.info("Generating report index...")
    logger.info(f"Reports directory: {reports_dir}")
    logger.info(f"Available reports: {len(available_reports)}")

    # Generate index
    generate_index(reports_dir, output_file, available_reports)

    logger.info("✓ Report index generation complete")


if __name__ == "__main__":
    main()
