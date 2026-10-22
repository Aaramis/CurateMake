# workflow/rules/msa_groups.smk
"""
Aggregation of the per-barcode MSA manifests into the single manifest SATIVA consumes.

The MSAs themselves are built by ehmm_msa.smk (backbone + ensemble-HMM projection), the only
strategy the pipeline keeps. A from-scratch per-group MAFFT variant ("clustered") lived here
until 2026-08-25; it was removed after eHMM won on both accuracy and cost. Everything below is
strategy-agnostic: it only merges manifests.

Historical note kept on purpose — the old docstring described why groups span kingdoms:

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

_cl = config["msa"].get("groups", {})
CL_OUT = _cl.get("out_dir", "results/MSA/groups")
CL_SATIVA_MANIFEST = _cl.get("sativa_manifest", "results/MSA/groups/sativa_manifest.json")
CL_SUMMARY = _cl.get("summary", "results/MSA/groups/groups_summary.json")
CL_RESIDUAL_REPORT = _cl.get("residual_report", "results/MSA/groups/residual_report.json")
CL_CONTAM_REPORT = _cl.get("contamination_report", "results/MSA/groups/contamination_report.json")
CL_CLUSTER_ID = float(_cl.get("cluster_id", 0.97))
CL_CLUSTER_COV = float(_cl.get("cluster_cov", 0.80))
CL_MIN_REF = int(_cl.get("min_ref", 4))
CL_MAX_GROUP = int(_cl.get("max_group", 10000))
CL_MAX_CLUSTER_REPS = int(_cl.get("max_cluster_reps", 0))
CL_SATIVA_MIN_SEQS = int(_cl.get("sativa_min_seqs", 4))
CL_MMSEQS_ENV = _cl.get("mmseqs_env", "mmseqs")
CL_TOOL_ENV = _cl.get("tool_env", "its-pipeline-dev")
CL_THREADS = int(_cl.get("threads", 32))
CL_KMER_PER_SEQ = int(_cl.get("cluster_kmer_per_seq", 300))  # linclust sensitivity knob



rule collect_groups:
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
        "logs/msa_groups/collect.log",
    benchmark:
        "benchmarks/msa_groups/collect.benchmark.txt"
    shell:
        """
        python workflow/scripts/MSA/collect_groups.py \
            --barcode-manifests {input.manifests} \
            --sativa-manifest {output.manifest} \
            --output-summary {output.summary} \
            --residual-report {output.residual} \
            --contamination-report {output.contamination} \
            --log-file {log}
        """
