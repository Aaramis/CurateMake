# workflow/rules/itsx.smk
"""
Rules Snakemake pour l'usage de l'outil ITSx
"""

import os
from pathlib import Path

ITSX_FASTA_DIR = config.get("itsx_output_dir", "results")
ITSX_BARCODES = config["msa"].get("barcodes", ["ITS", "ITS1", "ITS2"])
ITSX_DETERMINISTIC = bool(config.get("itsx_deterministic", False))
ITSX_FASTA_FILES = [
    os.path.join(ITSX_FASTA_DIR, f"{barcode}.fasta") for barcode in ITSX_BARCODES
]


def get_itsx_barcodes(wildcards):
    """Dynamically retrieve the list of barcodes from results/ITSx/"""
    itsx_dir = Path("results/ITSx")
    if itsx_dir.exists():
        barcodes = [d.name for d in itsx_dir.iterdir() if d.is_dir()]
        return sorted(barcodes)
    return ITSX_BARCODES  # Fallback


rule deduplicate_sequences_pre_itsx:
    """
    Deduplication of sequences before ITSx processing
    """
    input:
        col_summary="results/processed/col_matching_report.json"
    params:
        db=get_db_path(),
    output:
        report="results/processed/itsx_deduplicate_pre.json"
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
        deduplicate_pre_report="results/processed/itsx_deduplicate_pre.json"
    params:
        db=get_db_path(),
        output_fodler="results"
    output:
        report="results/processed/itsx_prepare_report.json",
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
    Run ITSx for a specific barcode (native parallelisation)
    """
    input:
        report="results/processed/itsx_prepare_report.json",
        fasta=lambda wildcards: os.path.join(ITSX_FASTA_DIR, f"{wildcards.barcode}.fasta")
    params:
        output_dir="results/ITSx/{barcode}",
        temp_dir="results/ITSx/.temp_{barcode}",
        cpu=1 if ITSX_DETERMINISTIC else config.get("itsx_cpu", 30),
        multi_thread="F" if ITSX_DETERMINISTIC else "T"
    output:
        done="results/processed/itsx_run_{barcode}.done"
    log:
        "logs/ITSX/run_itsx_{barcode}.log"
    benchmark:
        "benchmarks/itsx/run_itsx_{barcode}.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    threads: 1 if ITSX_DETERMINISTIC else config.get("itsx_cpu", 30)
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.temp_dir}

        echo "Processing {wildcards.barcode}..." >> {log} 2>&1

        ITSx --graphical F --preserve T --save_regions all \
            --cpu {params.cpu} --multi_thread {params.multi_thread} --fasta T --summary T \
            --temp {params.temp_dir} \
            -i {input.fasta} \
            -o {params.output_dir}/{wildcards.barcode} \
            >> {log} 2>&1

        # Clean up the temporary directory
        rm -rf {params.temp_dir}

        touch {output.done}
        """

rule run_itsx_all:
    """
    Wait until all ITSx barcodes are processed
    """
    input:
        lambda wildcards: expand(
            "results/processed/itsx_run_{barcode}.done",
            barcode=get_itsx_barcodes(wildcards)
        )
    output:
        done="results/processed/itsx_run.done"
    shell:
        "touch {output.done}"

rule process_barcode_itsx:
    """
    Process a single ITSx barcode (native Snakemake parallelisation)
    Starts as soon as its ITSx barcode is finished
    """
    input:
        itsx_done="results/processed/itsx_run_{barcode}.done",
    params:
        db=get_db_path(),
        itsx_dir="results/ITSx",
        data_dir="results",
    output:
        report="results/processed/itsx_{barcode}_result.json"
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
    Aggregate all per-barcode reports into a final report
    """
    input:
        lambda wildcards: expand(
            "results/processed/itsx_{barcode}_result.json",
            barcode=get_itsx_barcodes(wildcards)
        )
    output:
        report="results/processed/itsx_process.json"
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
    Deduplication of sequences after ITSx processing
    """
    input:
        itsx_report="results/processed/itsx_process.json",
    params:
        db=get_db_path(),
    output:
        report="results/processed/itsx_deduplicate_post.json"
    log:
        "logs/ITSX/itsx_deduplicate_post.log"
    benchmark:
        "benchmarks/itsx/deduplicate_sequences_post_itsx.benchmark.txt"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/deduplicate_sequences.py"
