#!/usr/bin/env python3
"""
Agrège les rapports ITSx individuels de chaque barcode en un rapport final unique.

Ce script est conçu pour être utilisé avec Snakemake pour combiner les résultats
de la parallélisation par barcode.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Any

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_barcode_report(report_path: Path) -> Dict[str, Any]:
    """
    Charge un rapport JSON individuel d'un barcode.

    Args:
        report_path: Chemin vers le fichier JSON du rapport

    Returns:
        Dictionnaire contenant les données du rapport
    """
    try:
        with open(report_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Erreur lors de la lecture de {report_path}: {e}")
        raise


def calculate_totals(barcodes_data: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Calcule les totaux globaux à partir des données de tous les barcodes.

    Args:
        barcodes_data: Dictionnaire avec les données de chaque barcode

    Returns:
        Dictionnaire avec les totaux globaux
    """
    totals = {
        "total_processed": 0,
        "total_subregions_found": 0,
        "total_subregions_inserted": 0,
        "total_reconstructed_its": 0,
        "total_no_detection": 0
    }

    for barcode, stats in barcodes_data.items():
        totals["total_processed"] += stats.get("processed", 0)
        totals["total_no_detection"] += stats.get("no_detection", 0)

        # Sous-régions
        subregions = stats.get("subregions", {})
        for region, region_stats in subregions.items():
            totals["total_subregions_found"] += region_stats.get("found", 0)
            totals["total_subregions_inserted"] += region_stats.get("inserted", 0)

        # ITS reconstruit
        reconstructed = stats.get("reconstructed_its", {})
        totals["total_reconstructed_its"] += reconstructed.get("inserted", 0)

    return totals


def aggregate_reports(input_files: list, output_file: Path) -> None:
    """
    Agrège tous les rapports individuels en un rapport final.

    Args:
        input_files: Liste des chemins vers les rapports individuels
        output_file: Chemin vers le fichier de rapport final
    """
    logger.info(f"Agrégation de {len(input_files)} rapports...")

    # Charger tous les rapports individuels
    barcodes_data = {}
    for report_path in input_files:
        report_path = Path(report_path)
        logger.info(f"Lecture de {report_path.name}...")

        data = load_barcode_report(report_path)

        # Extraire le barcode du nom du fichier
        # Format: itsx_{barcode}_result.json
        barcode = report_path.stem.replace("itsx_", "").replace("_result", "")

        barcodes_data[barcode] = data

    # Calculer les totaux
    totals = calculate_totals(barcodes_data)

    # Créer le rapport final
    final_report = {
        "barcodes": barcodes_data,
        "totals": totals
    }

    # Écrire le rapport final
    logger.info(f"Écriture du rapport final vers {output_file}...")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(final_report, f, indent=2)

    # Afficher un résumé
    logger.info("=" * 70)
    logger.info("RÉSUMÉ DE L'AGRÉGATION")
    logger.info("=" * 70)
    logger.info(f"Barcodes traités: {len(barcodes_data)}")
    logger.info(f"Séquences traitées: {totals['total_processed']}")
    logger.info(f"Sous-régions trouvées: {totals['total_subregions_found']}")
    logger.info(f"Sous-régions insérées: {totals['total_subregions_inserted']}")
    logger.info(f"ITS reconstruits: {totals['total_reconstructed_its']}")
    logger.info(f"Non-détections: {totals['total_no_detection']}")
    logger.info("=" * 70)


def main():
    """Point d'entrée principal pour Snakemake."""
    try:
        # Récupérer les paramètres depuis Snakemake
        input_files = snakemake.input
        output_file = Path(snakemake.output.report)

        # Configurer le logging vers le fichier de log
        if hasattr(snakemake, 'log') and snakemake.log:
            log_file = snakemake.log[0]
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(
                logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            )
            logger.addHandler(file_handler)

        # Agréger les rapports
        aggregate_reports(input_files, output_file)

        logger.info("Agrégation terminée avec succès!")

    except Exception as e:
        logger.error(f"Erreur lors de l'agrégation: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
