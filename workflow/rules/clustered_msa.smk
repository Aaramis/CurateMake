# workflow/rules/clustered_msa.smk
"""
Sequence-first MSA construction — the ONLY MSA strategy.

Per barcode: mmseqs2-cluster the WHOLE barcode -> absorb tiny clusters -> merge clusters of similar
taxonomy into kingdom-interleaved groups up to a size cap -> MAFFT + trimAl -> SATIVA manifest. No
(barcode, family) scatter, so a contaminant clusters with its true relatives (by sequence) and
cross-clade / cross-kingdom contamination becomes visible to SATIVA. Validated at scale in
10.DOI_CurateMake/simulation/SIMULATION_REPORT.md (kingdom-interleaved packing corrects 54% of
cross-kingdom contamination at 41k reads, vs ~0% for a single global naive MSA that also collapses).

Emits the artifacts sativa.smk / reports.smk / msa.smk consume:
  - CL_SATIVA_MANIFEST : {tasks:[{barcode, tax_level, taxon(=msa_id), taxon_label, msa_path, n}]}
  - CL_SUMMARY / CL_RESIDUAL_REPORT / CL_CONTAM_REPORT : the summary the reports read + stubs (the
    residual re-insertion / cross-family climb / coarse-contamination tiers do not exist here — every
    sequence is placed at clustering, and a contaminant clusters with its true kin).
"""
import os

FULL_CSV = config["msa"]["csv_path"]
BU_BARCODES = config["msa"].get("barcodes", ["ITS", "ITS1", "ITS2"])

_cl = config["msa"].get("clustered", {})
CL_OUT = _cl.get("out_dir", "results/MSA/clustered")
CL_SATIVA_MANIFEST = _cl.get("sativa_manifest", "results/MSA/clustered/sativa_manifest.json")
CL_SUMMARY = _cl.get("summary", "results/MSA/clustered/clustered_summary.json")
CL_RESIDUAL_REPORT = _cl.get("residual_report", "results/MSA/clustered/residual_report.json")
CL_CONTAM_REPORT = _cl.get("contamination_report", "results/MSA/clustered/contamination_report.json")
CL_CLUSTER_ID = float(_cl.get("cluster_id", 0.97))
CL_CLUSTER_COV = float(_cl.get("cluster_cov", 0.80))
CL_MIN_REF = int(_cl.get("min_ref", 4))
CL_MAX_GROUP = int(_cl.get("max_group", 10000))
CL_MAX_CLUSTER_REPS = int(_cl.get("max_cluster_reps", 0))
CL_SATIVA_MIN_SEQS = int(_cl.get("sativa_min_seqs", 4))
CL_MMSEQS_ENV = _cl.get("mmseqs_env", "mmseqs")
CL_TOOL_ENV = _cl.get("tool_env", "its-pipeline-dev")
CL_THREADS = int(_cl.get("threads", 32))

_SCRIPT = "workflow/scripts/MSA/build_msas_clustered.py"


rule clustered_build_barcode:
    """Cluster + merge (kingdom-interleaved) + align one barcode -> its SATIVA manifest + MSAs."""
    input:
        csv=FULL_CSV,
    output:
        manifest=os.path.join(CL_OUT, "{barcode}", "sativa_manifest.json"),
    params:
        outdir=lambda w: os.path.join(CL_OUT, w.barcode),
        script=_SCRIPT,
    threads: CL_THREADS
    log:
        "logs/clustered/{barcode}.log",
    shell:
        r"""
        python {params.script} \
            --full-csv {input.csv} --barcode {wildcards.barcode} --out-dir {params.outdir} \
            --cluster-id {CL_CLUSTER_ID} --cluster-cov {CL_CLUSTER_COV} \
            --min-ref {CL_MIN_REF} --max-group {CL_MAX_GROUP} --max-cluster-reps {CL_MAX_CLUSTER_REPS} \
            --sativa-min-seqs {CL_SATIVA_MIN_SEQS} --threads {threads} \
            --mmseqs-env {CL_MMSEQS_ENV} --tool-env {CL_TOOL_ENV} > {log} 2>&1
        """


rule clustered_collect:
    """Aggregate per-barcode manifests -> the single SATIVA manifest sativa.smk consumes + the
    summary / residual / contamination reports the downstream reports read."""
    input:
        manifests=expand(os.path.join(CL_OUT, "{barcode}", "sativa_manifest.json"),
                         barcode=BU_BARCODES),
    output:
        manifest=CL_SATIVA_MANIFEST,
        summary=CL_SUMMARY,
        residual=CL_RESIDUAL_REPORT,
        contamination=CL_CONTAM_REPORT,
    log:
        "logs/clustered/collect.log",
    shell:
        """
        python workflow/scripts/MSA/collect_clustered.py \
            --barcode-manifests {input.manifests} \
            --sativa-manifest {output.manifest} \
            --output-summary {output.summary} \
            --residual-report {output.residual} \
            --contamination-report {output.contamination} \
            --log-file {log}
        """
