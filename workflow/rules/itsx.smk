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