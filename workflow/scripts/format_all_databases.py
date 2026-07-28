#!/usr/bin/env python3
"""
Summary script for every formatting operation, including taxonomy deduplication.
"""

import json
import time
import logging
from pathlib import Path
from typing import Dict, Any

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_report(report_path: str) -> Dict[str, Any]:
    """Load a JSON report."""
    try:
        with open(report_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error while loading {report_path}: {e}")
        return {}

def generate_comprehensive_summary() -> Dict[str, Any]:
    """Generate a comprehensive summary for every operation."""
    
    # Configuration Snakemake
    if 'snakemake' in globals():
        formatted_reports = snakemake.params.formatted
        deduplication_report = snakemake.input.deduplication
        output_file = snakemake.output.summary
        log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    else:
        # Test mode (example only)
        # Replace these paths with your own when running outside Snakemake
        formatted_reports = [
            "results/processed/BOLD_formatted_report.json",
            "results/processed/CALeDNA_Multi_formatted_report.json",
            "results/processed/PLANITS_formatted_report.json",
            "results/processed/UNITE_formatted_report.json"
        ]
        deduplication_report = "results/processed/taxonomy_deduplication_report.json"
        output_file = "results/processed/formatting_summary.json"
    
    summary = {
        'process': 'complete_database_formatting',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'databases_processed': [],
        'taxonomy_deduplication': {},
        'overall_stats': {},
        'warnings': [],
        'success': True
    }
    
    total_sequences = 0
    total_processing_time = 0.0
    
    # 1. Process per-database formatting reports
    logger.info("Consolidating formatting reports...")
    
    all_report_paths = []
    for item in formatted_reports:
        if isinstance(item, list):
            all_report_paths.extend(item)
        else:
            all_report_paths.append(item)


    for report_path in all_report_paths:
        report = load_report(report_path)
        if report:
            db_name = Path(report_path).stem.replace('_formatted_report', '')
            
            stats = report.get('stats') or report.get('aggregated_stats', {})
            
            input_file_info = report.get('input_file')
            if not input_file_info:
                input_files_list = report.get('input_files', [])
                input_file_info = ", ".join(input_files_list) if input_files_list else 'N/A'

            db_summary = {
                'database': db_name,
                'source': report.get('source', 'Unknown'),
                'version': report.get('version', 'Unknown'),
                'input_file': input_file_info,
                'stats': stats,
                'processing_rate': report.get('processing_rate_per_second', 0),
                'optimizations': report.get('optimizations_applied', []),
                'success': True
            }
            
            summary['databases_processed'].append(db_summary)
            
            # Aggregate statistics
            total_sequences += stats.get('processed_sequences', 0)
            total_processing_time += stats.get('processing_time', 0)
        else:
            summary['warnings'].append(f"Unable to load report: {report_path}")
            summary['success'] = False
    
    # 2. Process deduplication report
    logger.info("Integrating the deduplication report...")
    
    dedup_report = load_report(deduplication_report)
    if dedup_report:
        dedup_results = dedup_report.get('results', {})
        summary['taxonomy_deduplication'] = {
            'executed': True,
            'duplicates_found': dedup_results.get('duplicates_found', False),
            'duplicates_removed': dedup_results.get('stats', {}).get('duplicates_removed', 0),
            'sequences_updated': dedup_results.get('stats', {}).get('sequences_updated', 0),
            'processing_time': dedup_results.get('processing_time', 0),
            'integrity_verified': dedup_results.get('integrity_check', {}).get('integrity_ok', False),
            'final_taxonomy_count': dedup_results.get('integrity_check', {}).get('final_taxonomy_count', 0)
        }
        
        total_processing_time += dedup_results.get('processing_time', 0)
        
        if not dedup_results.get('integrity_check', {}).get('integrity_ok', True):
            summary['warnings'].append("Integrity issue detected after deduplication")
            summary['success'] = False
    else:
        summary['warnings'].append(f"Deduplication report not available: {deduplication_report}")
        summary['taxonomy_deduplication'] = {'executed': False}
    
    # 3. Statistiques globales
    summary['overall_stats'] = {
        'total_databases_processed': len(summary['databases_processed']),
        'total_sequences_processed': total_sequences,
        'total_processing_time_seconds': round(total_processing_time, 2),
        'average_processing_rate': round(total_sequences / total_processing_time if total_processing_time > 0 else 0, 0),
        'taxonomy_optimization_applied': summary['taxonomy_deduplication'].get('executed', False),
        'duplicates_removed': summary['taxonomy_deduplication'].get('duplicates_removed', 0),
        'overall_success': summary['success'] and len(summary['warnings']) == 0
    }
    
    # 4. Recommendations
    recommendations = []
    if summary['taxonomy_deduplication'].get('duplicates_removed', 0) > 0:
        recommendations.append("Duplicate taxonomies removed successfully")
    if not summary['taxonomy_deduplication'].get('integrity_verified', True):
        recommendations.append("WARNING: Manually verify database integrity")
    if total_processing_time > 3600:  # More than one hour
        recommendations.append("Consider tuning parameters for large volumes")
    summary['recommendations'] = recommendations
    
    # 5. Persistence
    try:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Comprehensive summary saved: {output_file}")
        
        # Display the summary
        print(f"\n{'='*70}")
        print("COMPLETE SUMMARY OF DATABASE FORMATTING")
        print(f"{'='*70}")
        print(f"Databases processed: {summary['overall_stats']['total_databases_processed']}")
        print(f"Total sequences: {summary['overall_stats']['total_sequences_processed']:,}")
        print(f"Total time: {summary['overall_stats']['total_processing_time_seconds']:.1f}s")
        print(f"Average rate: {summary['overall_stats']['average_processing_rate']:,.0f} seq/s")
        print(f"Deduplication: {'✓' if summary['taxonomy_deduplication']['executed'] else '✗'}")
        if summary['taxonomy_deduplication'].get('duplicates_removed', 0) > 0:
            print(f"Duplicates removed: {summary['taxonomy_deduplication']['duplicates_removed']}")
            print(f"Relations updated: {summary['taxonomy_deduplication']['sequences_updated']}")
        print(f"Overall success: {'✓' if summary['overall_stats']['overall_success'] else '✗'}")
        if summary['warnings']:
            print(f"Warnings: {len(summary['warnings'])}")
        print(f"{'='*70}")
        
        return summary
        
    except Exception as e:
        logger.error(f"Error while saving summary: {e}")
        raise

def main():
    """Script entry point."""
    try:
        summary = generate_comprehensive_summary()
        
        if not summary['overall_stats']['overall_success']:
            logger.warning("Processing finished with warnings")
            return 1
        
        logger.info("Full formatting workflow completed successfully")
        return 0
        
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
