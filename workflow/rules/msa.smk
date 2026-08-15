# workflow/rules/msa.smk
"""
Sequence extraction for the MSA stage + the MSA analysis report.

The MSA construction itself lives in rules/clustered_msa.smk (sequence-first strategy).
This file keeps only:
  - prepare_msa: extract the sequences (+ CoL) into config["msa"]["csv_path"],
    which is the input of the bottom-up split.
  - analyze_msa_report: the bottom-up HTML report.
"""

import os


rule prepare_msa:
    input:
        deduplicate_post_report="results/processed/itsx_deduplicate_post.json"
    params:
        db=get_db_path(),
        output_fodler="results",
        csv_path=config["msa"]['csv_path']
    output:
        report="results/processed/msa_prepare_report.json",
        sequences=config["msa"]['csv_path']
    log:
        "logs/MSA/prepare_msa.log"
    benchmark:
        "benchmarks/msa/prepare_msa.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    script:
        "../scripts/MSA/prepare_msa.py"


# MSA analysis report (bottom-up): consumes the bottom-up summary + the per-family
# manifests + the residual and contamination reports.
_MSA_CL = config["msa"].get("clustered", {})
rule analyze_msa_report:
    input:
        summary=_MSA_CL.get("summary", "results/MSA/clustered/clustered_summary.json"),
        residual=_MSA_CL.get("residual_report", "results/MSA/clustered/residual_report.json"),
        contamination=_MSA_CL.get("contamination_report", "results/MSA/clustered/contamination_report.json")
    params:
        out_dir=_MSA_CL.get("out_dir", "results/MSA/clustered")
    output:
        report=config["msa"]["output_report"]
    log:
        "logs/msa_report/analyze_msa_report.log"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/generate_msa_report.py \
            --summary {input.summary} \
            --out-dir {params.out_dir} \
            --residual-report {input.residual} \
            --contamination-report {input.contamination} \
            --output {output.report} \
            --log-file {log}
        """
