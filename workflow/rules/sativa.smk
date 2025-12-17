# workflow/rules/sativa.smk
"""
Règles Snakemake pour l'analyse SATIVA
Détection des misidentifications taxonomiques sur les MSA
"""

import os
from pathlib import Path


# Règle principale pour SATIVA (batch processing)
rule sativa_all:
    """
    Lance SATIVA sur tous les MSA disponibles en parallèle.
    Utilise ThreadPoolExecutor pour traiter plusieurs MSA simultanément.
    """
    input:
        msa_report="data/processed/msa_run_report.json"  # Dépendance explicite sur MSA
    output:
        report="data/processed/sativa_analysis_report.json"
    params:
        input_msa_dir=config["sativa"]["input_msa_dir"],
        output_dir=config["sativa"]["output_dir"],
        logs_dir=config["sativa"]["logs_dir"],
        barcodes=config["sativa"]["barcodes"],
        tax_levels=config["sativa"]["tax_levels"],
        max_seqs=config["sativa"]["max_seq_count"],
        max_len=config["sativa"]["max_alignment_length"],
        sativa_script=config["sativa"]["sativa_script_path"],
        sativa_mode=config["sativa"]["sativa_mode"],
        confidence=config["sativa"]["confidence_threshold"],
        taxonomy_code=config["sativa"]["taxonomy_code"],
        threads=config["sativa"]["threads"],
        max_parallel_jobs=config["sativa"]["max_parallel_jobs"],
        timeout_seconds=config["sativa"]["timeout_seconds"]
    log:
        "logs/sativa/batch_sativa.log"
    benchmark:
        "benchmarks/sativa/sativa_all.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/batch_sativa.py"


# Règle pour générer le rapport de visualisation interactif
rule sativa_visualization:
    """
    Génère un rapport HTML interactif avec visualisations Plotly des résultats SATIVA.
    """
    input:
        report="data/processed/sativa_analysis_report.json"
    output:
        report=config["sativa"]["visualization_report"]
    params:
        sativa_output_dir=config["sativa"]["output_dir"]
    log:
        "logs/sativa/sativa_visualization.log"
    benchmark:
        "benchmarks/sativa/sativa_visualization.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/visualisation_sativa.py"


# Règle pour intégrer les résultats SATIVA dans la base de données (v2.0)
rule sativa_integrate:
    """
    Intègre les corrections taxonomiques de SATIVA dans la base de données.
    Version 2.0 avec validation CoL et gestion intelligente des modifications.

    Fonctionnalités:
    - Validation CoL pour toutes les modifications (sauf délétions simples)
    - Cache CoL pour optimiser les performances
    - Création de la table sativa_warnings
    - Traitement parallèle des fichiers .mis
    - Bulk inserts pour meilleures performances

    Gestion des modifications:
    - VALIDATED (dans .tax, pas dans .mis) → seq_tax avec status='Correct'
    - DELETION_SIMPLE (1 niveau) → warning uniquement, pas de seq_tax
    - DELETION_MULTIPLE (≥2 niveaux) → validation CoL puis intégration
    - SUBSTITUTION → validation CoL puis intégration
    - COMPLETION → validation CoL puis intégration

    Hiérarchie des décisions avec CoL:
    - CoL exact → taxonomie SATIVA, status='Correct'
    - CoL modifié → taxonomie CoL, status='Check', warning créé
    - CoL no match → taxonomie SATIVA, status='Check', warning créé
    """
    input:
        sativa_report="data/processed/sativa_analysis_report.json",
        db=get_db_path()
    output:
        report=config["sativa_integration"]["reporting"]["output_dir"] + "/sativa_integration_report.json"
    params:
        db_path=get_db_path(),
        sativa_output_dir=config["sativa"]["output_dir"],
        barcodes=config["sativa"]["barcodes"],
        tax_levels=config["sativa"]["tax_levels"]
    log:
        "logs/sativa/integrate_sativa.log"
    benchmark:
        "benchmarks/sativa/sativa_integrate.benchmark.txt"
    conda:
        "../envs/sativa.yaml"
    script:
        "../scripts/Sativa/integrate_sativa.py"
