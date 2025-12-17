# workflow/rules/itsx.smk
"""
Règles Snakemake pour l'usage de l'outil ITSx
"""

import os
from pathlib import Path

ITSX_FASTA_DIR = config.get("itsx_output_dir", "data")
ITSX_BARCODES = config["msa"].get("barcodes", ["ITS", "ITS1", "ITS2"])
ITSX_FASTA_FILES = [
    os.path.join(ITSX_FASTA_DIR, f"{barcode}.fasta") for barcode in ITSX_BARCODES
]


def get_itsx_barcodes(wildcards):
    """Récupère dynamiquement la liste des barcodes depuis data/ITSx/"""
    itsx_dir = Path("data/ITSx")
    if itsx_dir.exists():
        barcodes = [d.name for d in itsx_dir.iterdir() if d.is_dir()]
        return barcodes
    return ITSX_BARCODES  # Fallback


rule deduplicate_sequences_pre_itsx:
    """
    Déduplication des séquences avant le traitement ITSx
    """
    input:
        col_summary="data/processed/col_matching_report.json"
    params:
        db=get_db_path(),
    output:
        report="data/processed/itsx_deduplicate_pre.json"
    log:
        "logs/ITSX/itsx_deduplicate_pre.log"
    benchmark:
        "benchmarks/itsx/deduplicate_sequences_pre_itsx.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/deduplicate_sequences.py"

rule prepare_itsx:
    input:
        deduplicate_pre_report="data/processed/itsx_deduplicate_pre.json"
    params:
        db=get_db_path(),
        output_fodler="data"
    output:
        report="data/processed/itsx_prepare_report.json",
        fastas=ITSX_FASTA_FILES
    log:
        "logs/ITSX/prepare_itsx.log"
    benchmark:
        "benchmarks/itsx/prepare_itsx.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/itsx/prepare_itsx.py"


rule run_itsx_barcode:
    """
    Exécute ITSx pour un barcode spécifique (parallélisation native)
    """
    input:
        report="data/processed/itsx_prepare_report.json",
        fasta=lambda wildcards: os.path.join(ITSX_FASTA_DIR, f"{wildcards.barcode}.fasta")
    params:
        output_dir="data/ITSx/{barcode}",
        temp_dir="data/ITSx/.temp_{barcode}",
        cpu=config.get("itsx_cpu", 30)
    output:
        done="data/processed/itsx_run_{barcode}.done"
    log:
        "logs/ITSX/run_itsx_{barcode}.log"
    benchmark:
        "benchmarks/itsx/run_itsx_{barcode}.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    threads: config.get("itsx_cpu", 30)
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.temp_dir}

        echo "Processing {wildcards.barcode}..." >> {log} 2>&1

        ITSx --graphical F --preserve T --save_regions all \
            --cpu {params.cpu} --multi_thread T --fasta T --summary T \
            --temp {params.temp_dir} \
            -i {input.fasta} \
            -o {params.output_dir}/{wildcards.barcode} \
            >> {log} 2>&1

        # Nettoyer le répertoire temporaire
        rm -rf {params.temp_dir}

        touch {output.done}
        """

rule run_itsx_all:
    """
    Attend que tous les barcodes ITSx soient traités
    """
    input:
        expand("data/processed/itsx_run_{barcode}.done",
               barcode=["ITS", "ITS1", "ITS2"])
    output:
        done="data/processed/itsx_run.done"
    shell:
        "touch {output.done}"

rule process_barcode_itsx:
    """
    Traite un barcode ITSx individuellement (parallélisation Snakemake native)
    Démarre dès que son barcode ITSx est terminé
    """
    input:
        itsx_done="data/processed/itsx_run_{barcode}.done",
    params:
        db=get_db_path(),
        itsx_dir="data/ITSx",
        data_dir="data",
    output:
        report="data/processed/itsx_{barcode}_result.json"
    log:
        "logs/ITSX/process_itsx_{barcode}.log"
    benchmark:
        "benchmarks/itsx/process_barcode_{barcode}.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/itsx/process_single_barcode_itsx.py"


rule aggregate_itsx_results:
    """
    Agrège tous les rapports individuels de barcodes en un rapport final
    """
    input:
        lambda wildcards: expand(
            "data/processed/itsx_{barcode}_result.json",
            barcode=get_itsx_barcodes(wildcards)
        )
    output:
        report="data/processed/itsx_process.json"
    log:
        "logs/ITSX/aggregate_itsx_results.log"
    benchmark:
        "benchmarks/itsx/aggregate_itsx_results.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/itsx/aggregate_itsx_reports.py"

rule deduplicate_sequences_post_itsx:
    """
    Déduplication des séquences après le traitement ITSx
    """
    input:
        itsx_report="data/processed/itsx_process.json",
    params:
        db=get_db_path(),
    output:
        report="data/processed/itsx_deduplicate_post.json"
    log:
        "logs/ITSX/itsx_deduplicate_post.log"
    benchmark:
        "benchmarks/itsx/deduplicate_sequences_post_itsx.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/deduplicate_sequences.py"
