# workflow/rules/msa.smk
"""
Règles Snakemake pour la création de MSA propre
"""

import os
import json


MSA_DIR = config["msa"]["msa_dir"]
MSA_MANIFEST = config["msa"].get("task_manifest", "data/processed/msa_task_manifest.json")
_stage_markers = config["msa"].get("stage_markers", {})
MSA_ALIGNMENTS_FLAG = _stage_markers.get("alignments", "data/processed/msa_alignments_complete.flag")
MSA_TRIMMING_FLAG = _stage_markers.get("trimming", "data/processed/msa_trimming_complete.flag")
MSA_HMM_FLAG = _stage_markers.get("hmm_profiles", "data/processed/msa_hmm_complete.flag")
MSA_COMBINED_HMM_FLAG = _stage_markers.get("combined_hmm", "data/processed/msa_combined_hmm_complete.flag")
MSA_CLASSIFICATION_FLAG = _stage_markers.get("classification", "data/processed/msa_classification_complete.flag")
MSA_ASSIGNMENT_SUMMARY = config["msa"].get("assignment_integration_summary", "data/processed/msa_assignment_integration.json")

CLASSIFICATION_MANIFEST = config["msa"]["classification_manifest"]
UNPROCESSED_CSV = config["msa"]["unprocessed_sequences_csv"]
GLOBAL_HMM = config["msa"]["global_hmm"]
GLOBAL_HMM_RAW = config["msa"].get("global_hmm_raw", GLOBAL_HMM + ".raw")
_hmm_filter_cfg = config["msa"].get("hmm_filter", {})
HMM_FILTER_REPORT = _hmm_filter_cfg.get("report_json", "data/MSA/combined/hmmstat_filter_report.json")
HMM_MIN_RELENT = float(_hmm_filter_cfg.get("min_relent", 0.2))
HMM_MIN_PRELE = float(_hmm_filter_cfg.get("min_prele", 0.2))
CLASSIFICATION_BATCHES_DIR = config["msa"]["classification_batches_dir"]
CLASSIFICATION_BATCHES_JSON_DIR = config["msa"]["classification_batches_json_dir"]
CLASSIFICATION_RESULTS = config["msa"]["classification_results_csv"]
CLASSIFICATION_STATS = config["msa"]["classification_stats_json"]


def _load_classification_manifest():
    manifest_path = checkpoints.msa_classification_manifest.get().output.manifest
    with open(manifest_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def classification_batches():
    manifest = _load_classification_manifest()
    return manifest.get("batches", [])


def classification_batch_json(wildcards):
    for batch in classification_batches():
        if batch["batch_id"] == wildcards.batch:
            return batch["json_path"]
    raise ValueError(f"Batch {wildcards.batch} not found in classification manifest.")


def classification_batch_outputs(_wildcards):
    return [
        os.path.join(CLASSIFICATION_BATCHES_DIR, f"{batch['batch_id']}.csv")
        for batch in classification_batches()
    ]


def _load_manifest_from_checkpoint():
    manifest_path = checkpoints.msa_task_manifest.get().output.manifest
    with open(manifest_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _manifest_tasks():
    manifest = _load_manifest_from_checkpoint()
    tasks = []
    for barcode, payload in manifest.get("barcodes", {}).items():
        for section in ("tasks", "large_tasks"):
            for task in payload.get(section, []):
                tasks.append({
                    "barcode": barcode,
                    "level": task["level"],
                    "clade": task["clade_name"]
                })
    return tasks


def _manifest_barcodes():
    manifest = _load_manifest_from_checkpoint()
    return sorted(manifest.get("barcodes", {}).keys())


def alignment_inputs(_wildcards):
    outputs = []
    for task in _manifest_tasks():
        barcode = task["barcode"]
        level = task["level"]
        clade = task["clade"]
        outputs.append(os.path.join(MSA_DIR, barcode, f"alignments_{level}", f"{clade}_msa.fasta"))
    return outputs


def trimmed_inputs(_wildcards):
    outputs = []
    for task in _manifest_tasks():
        barcode = task["barcode"]
        level = task["level"]
        clade = task["clade"]
        outputs.append(os.path.join(MSA_DIR, barcode, f"alignments_{level}", f"{clade}_msa.trimmed.fasta"))
    return outputs


def hmm_inputs(_wildcards):
    outputs = []
    for task in _manifest_tasks():
        barcode = task["barcode"]
        level = task["level"]
        clade = task["clade"]
        outputs.append(os.path.join(MSA_DIR, barcode, "hmm", level, f"{clade}.hmm"))
    return outputs


def combined_inputs(_wildcards):
    return [os.path.join(MSA_DIR, barcode, "all_profiles.hmm") for barcode in _manifest_barcodes()]


def barcode_hmm_inputs(wildcards):
    outputs = []
    for task in _manifest_tasks():
        if task["barcode"] != wildcards.barcode:
            continue
        level = task["level"]
        clade = task["clade"]
        outputs.append(os.path.join(MSA_DIR, wildcards.barcode, "hmm", level, f"{clade}.hmm"))
    return outputs


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


checkpoint msa_task_manifest:
    input:
        prepare_msa_report="data/processed/msa_prepare_report.json",
        csv=config["msa"]['csv_path']
    params:
        barcodes=config["msa"]['barcodes'],
        hierarchy_levels=config["msa"]['hierarchy_levels'],
        tax_columns=config["msa"]['tax_columns'],
        domain=config["msa"]['domain'],
        min_domain=config["msa"]['min_domain'],
        min_seq_threshold=config["msa"]['min_seq_threshold'],
        max_seq_threshold=config["msa"]['max_seq_threshold'],
        large_task_threshold=config["msa"]['large_task_threshold'],
        id_col=config["msa"]['id_col'],
        seq_col=config["msa"]['seq_col'],
        barcode_col=config["msa"]['barcode_col'],
        verbose=False
    output:
        manifest=MSA_MANIFEST
    log:
        "logs/MSA/task_manifest.log"
    benchmark:
        "benchmarks/msa/task_manifest.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    script:
        "../scripts/MSA/task_manifest.py"


rule msa_align_task:
    input:
        manifest=MSA_MANIFEST,
        csv=config["msa"]['csv_path']
    output:
        msa=os.path.join(MSA_DIR, "{barcode}", "alignments_{level}", "{clade}_msa.fasta")
    params:
        output_folder=MSA_DIR,
        seq_col=config["msa"]['seq_col'],
        id_col=config["msa"]['id_col'],
        tax_columns=" ".join(config["msa"]['tax_columns'])
    log:
        "logs/MSA/align/{barcode}/{level}/{clade}.log"
    benchmark:
        "benchmarks/msa/align/{barcode}.{level}.{clade}.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/run_alignment_task.py \
            --manifest {input.manifest} \
            --barcode {wildcards.barcode} \
            --clade {wildcards.clade} \
            --input-csv {input.csv} \
            --output-folder {params.output_folder} \
            --seq-col {params.seq_col} \
            --id-col {params.id_col} \
            --tax-columns {params.tax_columns} \
            --log-file {log}
        """


rule msa_alignments_ready:
    input:
        alignment_inputs
    output:
        flag=MSA_ALIGNMENTS_FLAG
    shell:
        "touch {output.flag}"


rule msa_trim_task:
    input:
        msa=os.path.join(MSA_DIR, "{barcode}", "alignments_{level}", "{clade}_msa.fasta")
    output:
        trimmed=os.path.join(MSA_DIR, "{barcode}", "alignments_{level}", "{clade}_msa.trimmed.fasta"),
        removed=os.path.join(MSA_DIR, "{barcode}", "alignments_{level}", "{clade}_removed.json")
    params:
        seq_overlap=config["msa"]['trimal_seq_overlap'],
        res_overlap=config["msa"]['trimal_res_overlap']
    log:
        "logs/MSA/trimal/{barcode}/{level}/{clade}.log"
    benchmark:
        "benchmarks/msa/trimal/{barcode}.{level}.{clade}.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/trim_alignment_task.py \
            --input-msa {input.msa} \
            --trimmed-msa {output.trimmed} \
            --barcode {wildcards.barcode} \
            --clade {wildcards.clade} \
            --min-overlap {params.seq_overlap} \
            --res-overlap {params.res_overlap} \
            --removed-json {output.removed} \
            --log-file {log}
        """


rule msa_trimming_ready:
    input:
        trimmed_inputs
    output:
        flag=MSA_TRIMMING_FLAG
    shell:
        "touch {output.flag}"


rule msa_build_hmm_task:
    input:
        trimmed=os.path.join(MSA_DIR, "{barcode}", "alignments_{level}", "{clade}_msa.trimmed.fasta")
    output:
        hmm=os.path.join(MSA_DIR, "{barcode}", "hmm", "{level}", "{clade}.hmm")
    params:
        output_dir=os.path.join(MSA_DIR, "{barcode}", "hmm", "{level}")
    log:
        "logs/MSA/hmm/{barcode}/{level}/{clade}.log"
    benchmark:
        "benchmarks/msa/hmm/{barcode}.{level}.{clade}.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/build_hmm_task.py \
            --trimmed-msa {input.trimmed} \
            --output-dir {params.output_dir} \
            --barcode {wildcards.barcode} \
            --clade {wildcards.clade} \
            --log-file {log}
        """


rule msa_hmm_ready:
    input:
        hmm_inputs
    output:
        flag=MSA_HMM_FLAG
    shell:
        "touch {output.flag}"


rule msa_combine_barcode_hmms:
    input:
        hmms=barcode_hmm_inputs
    output:
        combined=os.path.join(MSA_DIR, "{barcode}", "all_profiles.hmm")
    params:
        hmm_dir=os.path.join(MSA_DIR, "{barcode}", "hmm")
    log:
        "logs/MSA/hmm_db/{barcode}.log"
    benchmark:
        "benchmarks/msa/hmm_db/{barcode}.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/combine_hmm_database.py \
            --hmm-dir {params.hmm_dir} \
            --output-hmm {output.combined} \
            --log-file {log}
        """


rule msa_combined_hmm_ready:
    input:
        GLOBAL_HMM,
        HMM_FILTER_REPORT
    output:
        flag=MSA_COMBINED_HMM_FLAG
    shell:
        "touch {output.flag}"


rule msa_combine_global_hmm:
    input:
        combined_inputs
    output:
        GLOBAL_HMM_RAW
    log:
        "logs/MSA/global_hmm.log"
    benchmark:
        "benchmarks/msa/global_hmm.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/combine_global_hmm.py \
            --input-hmms {input} \
            --output-hmm {output} \
            --log-file {log}
        """


rule msa_filter_global_hmm:
    input:
        raw=GLOBAL_HMM_RAW
    output:
        filtered=GLOBAL_HMM,
        report=HMM_FILTER_REPORT
    params:
        min_relent=HMM_MIN_RELENT,
        min_prele=HMM_MIN_PRELE
    log:
        "logs/MSA/global_hmm_filter.log"
    benchmark:
        "benchmarks/msa/global_hmm_filter.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/filter_global_hmm.py \
            --input-hmm {input.raw} \
            --output-hmm {output.filtered} \
            --report-json {output.report} \
            --min-relent {params.min_relent} \
            --min-prele {params.min_prele} \
            --log-file {log}
        """


checkpoint msa_classification_manifest:
    input:
        task_manifest=MSA_MANIFEST,
        trimmed=MSA_TRIMMING_FLAG,
        hmm_ready=MSA_HMM_FLAG,
        combined_ready=MSA_COMBINED_HMM_FLAG,
        csv=config["msa"]['csv_path']
    params:
        msa_dir=MSA_DIR,
        batch_size=config["msa"]['classification_batch_size'],
        id_col=config["msa"]['id_col'],
        barcode_col=config["msa"]['barcode_col'],
        batches_dir=CLASSIFICATION_BATCHES_JSON_DIR
    output:
        manifest=CLASSIFICATION_MANIFEST,
        unprocessed=UNPROCESSED_CSV
    log:
        "logs/MSA/classification_manifest.log"
    benchmark:
        "benchmarks/msa/classification_manifest.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/prepare_classification_batches.py \
            --input-csv {input.csv} \
            --manifest-json {input.task_manifest} \
            --msa-dir {params.msa_dir} \
            --output-csv {output.unprocessed} \
            --output-manifest {output.manifest} \
            --batches-dir {params.batches_dir} \
            --batch-size {params.batch_size} \
            --id-col {params.id_col} \
            --barcode-col {params.barcode_col} \
            --log-file {log}
        """


rule msa_classify_batch:
    input:
        csv=config["msa"]['csv_path'],
        seq_ids=classification_batch_json,
        hmm=GLOBAL_HMM,
        manifest=CLASSIFICATION_MANIFEST
    output:
        assignments=os.path.join(CLASSIFICATION_BATCHES_DIR, "{batch}.csv")
    params:
        seq_col=config["msa"]['seq_col'],
        id_col=config["msa"]['id_col'],
        tax_columns=" ".join(config["msa"]['tax_columns']),
        threads=config["msa"]['classification_threads_per_worker'],
        max_evalue=config["msa"]['classification_evalue_threshold']
    threads: config["msa"]['classification_threads_per_worker']
    log:
        "logs/MSA/classify/{batch}.log"
    benchmark:
        "benchmarks/msa/classify/{batch}.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/classify_sequences_task.py \
            --input-csv {input.csv} \
            --seq-ids-json {input.seq_ids} \
            --combined-hmm {input.hmm} \
            --output-csv {output.assignments} \
            --seq-col {params.seq_col} \
            --id-col {params.id_col} \
            --tax-columns {params.tax_columns} \
            --threads {params.threads} \
            --max-evalue {params.max_evalue} \
            --log-file {log}
        """


rule msa_collect_assignments:
    input:
        batches=classification_batch_outputs,
        unprocessed=UNPROCESSED_CSV
    output:
        assignments=CLASSIFICATION_RESULTS,
        stats=CLASSIFICATION_STATS
    params:
        id_col=config["msa"]['id_col'],
        barcode_col=config["msa"]['barcode_col']
    log:
        "logs/MSA/collect_assignments.log"
    benchmark:
        "benchmarks/msa/collect_assignments.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/collect_classification_results.py \
            --batch-csvs {input.batches} \
            --unprocessed-csv {input.unprocessed} \
            --output-assignments {output.assignments} \
            --output-json {output.stats} \
            --id-col {params.id_col} \
            --barcode-col {params.barcode_col} \
            --log-file {log}
        """


rule msa_integrate_assignments:
    input:
        assignments=CLASSIFICATION_RESULTS,
        unprocessed=UNPROCESSED_CSV,
        manifest=MSA_MANIFEST
    output:
        summary=MSA_ASSIGNMENT_SUMMARY
    params:
        msa_dir=MSA_DIR,
        id_col=config["msa"]['id_col'],
        seq_col=config["msa"]['seq_col'],
        tax_columns=config["msa"]['tax_columns'],
        max_thread_workers=config["msa"]['max_thread_workers']
    log:
        "logs/MSA/integrate_assignments.log"
    benchmark:
        "benchmarks/msa/integrate_assignments.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    script:
        "../scripts/MSA/add_assignments_to_msa.py"


rule msa_classification_ready:
    input:
        CLASSIFICATION_RESULTS,
        CLASSIFICATION_STATS,
        MSA_ASSIGNMENT_SUMMARY
    output:
        flag=MSA_CLASSIFICATION_FLAG
    shell:
        "touch {output.flag}"


rule run_msa:
    input:
        prepare_msa_report="data/processed/msa_prepare_report.json",
        task_manifest=MSA_MANIFEST,
        alignments_ready=MSA_ALIGNMENTS_FLAG,
        trimming_ready=MSA_TRIMMING_FLAG,
        hmm_ready=MSA_HMM_FLAG,
        combined_ready=MSA_COMBINED_HMM_FLAG,
        classification_manifest=CLASSIFICATION_MANIFEST,
        unprocessed_csv=UNPROCESSED_CSV,
        classification_results=CLASSIFICATION_RESULTS,
        classification_stats=CLASSIFICATION_STATS,
        classification_ready=MSA_CLASSIFICATION_FLAG
    params:
        csv_path=config["msa"]['csv_path'],
        db=get_db_path(),
        output_folder="data/MSA",
        min_seq_threshold=config["msa"]['min_seq_threshold'],
        max_seq_threshold=config["msa"]['max_seq_threshold'],
        large_task_threshold=config["msa"]['large_task_threshold'],
        min_domain=config["msa"]['min_domain'],
        trimal_seq_overlap=config["msa"]['trimal_seq_overlap'],
        trimal_res_overlap=config["msa"]['trimal_res_overlap'],
        domain=config["msa"]['domain'],
        barcodes=config["msa"]['barcodes'],
        hierarchy_levels=config["msa"]['hierarchy_levels'],
        tax_columns=config["msa"]['tax_columns'],
        id_col=config["msa"]['id_col'],
        seq_col=config["msa"]['seq_col'],
        barcode_col=config["msa"]['barcode_col'],
        max_thread_workers=config["msa"]['max_thread_workers'],
        classification_max_workers=config["msa"]['classification_max_workers'],
        classification_threads_per_worker=config["msa"]['classification_threads_per_worker'],
        classification_batch_size=config["msa"]['classification_batch_size'],
        logger_folder= "logs/MSA",
        task_manifest_path=MSA_MANIFEST,
        classification_manifest_path=CLASSIFICATION_MANIFEST

    output:
        report="data/processed/msa_run_report.json"
    log:
        "logs/MSA/run_msa.log"
    benchmark:
        "benchmarks/msa/run_msa.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    script:
        "../scripts/MSA/run_msa.py"


# Règle pour l'analyse et le reporting des MSA
rule analyze_msa_report:
    input:
        run_msa_report="data/processed/msa_run_report.json"
    params:
        msa_dir=config["msa"]["msa_dir"],
        barcodes=config["msa"]["barcodes"],
    output:
        report=config["msa"]["output_report"]
    log:
        "logs/msa_report/analyze_msa_report.log"
    benchmark:
        "benchmarks/msa/analyze_msa_report.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    script:
        "../scripts/MSA/generate_report_MSA.py"
