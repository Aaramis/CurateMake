# workflow/rules/itsx.smk
"""
Règles Snakemake pour l'usage de l'outil ITSx
"""

import os

rule prepare_itsx:
    input:
        col_summary="data/processed/col_matching_report.json"
    params:
        db=os.path.join(config["database"]["path"], config["database"]["name"]),
        output_fodler="data"
    output:
        report="data/processed/itsx_prepare_report.json"
    log:
        "logs/ITSX/prepare_itsx.log"
    conda:
        "../envs/itsx.yaml"
    script:
        "../scripts/prepare_itsx.py"


rule run_itsx:
    input:
        report="data/processed/itsx_prepare_report.json"
    params:
        input_dir="data",
        output_dir="data/ITSx",
        cpu=config.get("itsx_cpu", 30)
    output:
        done="data/processed/itsx_run.done"
    log:
        "logs/ITSX/run_itsx.log"
    conda:
        "../envs/itsx.yaml"
    threads: config.get("itsx_cpu", 30)
    shell:
        """
        # Traiter tous les fichiers FASTA disponibles
        for fasta_file in {params.input_dir}/*.fasta; do

            if [ -f "$fasta_file" ]; then
                barcode=$(basename "$fasta_file" .fasta)
                output_dir="{params.output_dir}/$barcode"
                mkdir -p "$output_dir"

                echo "Processing $barcode..." >> {log} 2>&1

                ITSx --graphical F --preserve T --save_regions all \
                    --cpu {params.cpu} --multi_thread T --fasta T --summary T \
                    -i "$fasta_file" \
                    -o "$output_dir/$barcode" \
                    >> {log} 2>&1
            fi
        done
        
        touch {output.done}
        """

rule process_itsx_results:
    input:
        itsx_done="data/processed/itsx_run.done",
    params:
        db=os.path.join(config["database"]["path"], config["database"]["name"]),
        itsx_dir="data/ITSx",
        min_length=config["sequence_validation"]["min_length"],
        max_n_fraction=config["sequence_validation"]["max_n_fraction"]
    output:
        report="data/processed/itsx_process.json"
    log:
        "logs/ITSX/process_itsx_results.log"
    conda:
        "../envs/itsx.yaml"
    shell:
        """
        python workflow/scripts/get_itsx_results.py \
            --database {params.db} \
            --itsx-dir {params.itsx_dir} \
            --min-length {params.min_length} \
            --max-n-fraction {params.max_n_fraction} \
            --output-report {output.report} \
            --log-file {log}
        """