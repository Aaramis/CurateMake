#!/usr/bin/env python3
"""
Script de résumé pour toutes les opérations de formatage, 
incluant la déduplication taxonomique.
"""

import json
import time
import logging
from pathlib import Path
from typing import Dict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_report(report_path: str) -> Dict:
    """Charge un rapport JSON."""
    try:
        with open(report_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Erreur lors du chargement de {report_path}: {e}")
        return {}

def generate_comprehensive_summary() -> Dict:
    """Génère un résumé complet de toutes les opérations."""
    
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
        # Mode test
        formatted_reports = []
        deduplication_report = "data/processed/taxonomy_deduplication_report.json"
        output_file = "data/processed/formatting_summary.json"
    
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
    
    # 1. Traitement des rapports de formatage par base
    logger.info("Consolidation des rapports de formatage...")
    
    for report_path in formatted_reports:
        report = load_report(report_path)
        if report:
            db_name = Path(report_path).stem.replace('_formatted_report', '')
            
            db_summary = {
                'database': db_name,
                'source': report.get('source', 'Unknown'),
                'version': report.get('version', 'Unknown'),
                'input_file': report.get('input_file', ''),
                'stats': report.get('stats', {}),
                'processing_rate': report.get('processing_rate_per_second', 0),
                'optimizations': report.get('optimizations_applied', []),
                'success': True
            }
            
            summary['databases_processed'].append(db_summary)
            
            # Agrégation des statistiques
            stats = report.get('stats', {})
            total_sequences += stats.get('processed_sequences', 0)
            total_processing_time += stats.get('processing_time', 0)
        else:
            summary['warnings'].append(f"Impossible de charger le rapport: {report_path}")
            summary['success'] = False
    
    # 2. Traitement du rapport de déduplication
    logger.info("Intégration du rapport de déduplication...")
    
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
            summary['warnings'].append("Problème d'intégrité détecté après la déduplication")
            summary['success'] = False
    else:
        summary['warnings'].append("Rapport de déduplication non disponible")
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
    
    # 4. Recommandations
    recommendations = []
    
    if summary['taxonomy_deduplication'].get('duplicates_removed', 0) > 0:
        recommendations.append("Taxonomies dupliquées supprimées avec succès")
    
    if not summary['taxonomy_deduplication'].get('integrity_verified', True):
        recommendations.append("ATTENTION: Vérifier manuellement l'intégrité de la base")
    
    if total_processing_time > 3600:  # Plus d'1 heure
        recommendations.append("Considérer l'optimisation des paramètres pour les gros volumes")
    
    summary['recommendations'] = recommendations
    
    # 5. Sauvegarde
    try:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Résumé complet sauvegardé: {output_file}")
        
        # Affichage du résumé
        print(f"\n{'='*70}")
        print("RÉSUMÉ COMPLET DU FORMATAGE DES BASES DE DONNÉES")
        print(f"{'='*70}")
        print(f"Bases traitées: {summary['overall_stats']['total_databases_processed']}")
        print(f"Séquences totales: {summary['overall_stats']['total_sequences_processed']:,}")
        print(f"Temps total: {summary['overall_stats']['total_processing_time_seconds']:.1f}s")
        print(f"Taux moyen: {summary['overall_stats']['average_processing_rate']:,.0f} seq/s")
        print(f"Déduplication: {'✓' if summary['taxonomy_deduplication']['executed'] else '✗'}")
        if summary['taxonomy_deduplication'].get('duplicates_removed', 0) > 0:
            print(f"Doublons supprimés: {summary['taxonomy_deduplication']['duplicates_removed']}")
            print(f"Relations mises à jour: {summary['taxonomy_deduplication']['sequences_updated']}")
        print(f"Succès global: {'✓' if summary['overall_stats']['overall_success'] else '✗'}")
        if summary['warnings']:
            print(f"Avertissements: {len(summary['warnings'])}")
        print(f"{'='*70}")
        
        return summary
        
    except Exception as e:
        logger.error(f"Erreur lors de la sauvegarde: {e}")
        raise

def main():
    """Point d'entrée principal."""
    try:
        summary = generate_comprehensive_summary()
        
        if not summary['overall_stats']['overall_success']:
            logger.warning("Le traitement s'est terminé avec des avertissements")
            return 1
        
        logger.info("Formatage complet terminé avec succès")
        return 0
        
    except Exception as e:
        logger.error(f"Erreur fatale: {e}")
        return 1

if __name__ == "__main__":
    exit(main())