# workflow/rules/msa.smk
"""
Sequence extraction for the MSA stage + the MSA analysis report.

The MSA construction itself lives in rules/bottomup_msa.smk (bottom-up strategy).
This file keeps only:
  - prepare_msa: extract the sequences (+ CoL) into config["msa"]["csv_path"],
    which is the input of the bottom-up split.
  - analyze_msa_report: the bottom-up HTML report.
"""

import os


rule prepare_msa:
    input:
        deduplicate_post_report="data/processed/itsx_deduplicate_post.json"
    params:
        db=get_db_path(),
        output_fodler="data",
        csv_path=config["msa"]['csv_path']
    output:
        report="data/processed/msa_prepare_report.json",
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
_MSA_BU = config["msa"].get("bottomup", {})
rule analyze_msa_report:
    input:
        summary=_MSA_BU.get("summary", "data/MSA/bottomup/bottomup_summary.json"),
        residual=_MSA_BU.get("residual_report", "data/MSA/bottomup/residual_report.json"),
        contamination=_MSA_BU.get("contamination_report", "data/MSA/bottomup/contamination_report.json")
    params:
        out_dir=_MSA_BU.get("out_dir", "data/MSA/bottomup/out")
    output:
        report=config["msa"]["output_report"]
    log:
        "logs/msa_report/analyze_msa_report.log"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/generate_report_bottomup.py \
            --summary {input.summary} \
            --out-dir {params.out_dir} \
            --residual-report {input.residual} \
            --contamination-report {input.contamination} \
            --output {output.report} \
            --log-file {log}
        """
