# workflow/rules/bottomup_msa.smk
"""
Bottom-up construction of the reference MSAs.
Scatter per (barcode, family); bottom-up internal to the job (cf. INTEGRATION_DESIGN.md).
SATIVA consumes the manifest produced here (see sativa.smk).
"""

import os
import json

_bu = config["msa"].get("bottomup", {})
FULL_CSV = config["msa"]["csv_path"]
BU_FAMILIES_DIR = _bu.get("families_dir", "results/MSA/bottomup/families")
BU_MANIFEST = _bu.get("families_manifest", "results/MSA/bottomup/families_manifest.json")
BU_OUT = _bu.get("out_dir", "results/MSA/bottomup/out")
BU_SUMMARY = _bu.get("summary", "results/MSA/bottomup/bottomup_summary.json")
BU_RESIDUAL_REPORT = _bu.get("residual_report", "results/MSA/bottomup/residual_report.json")
BU_RESIDUAL_EVALUE = float(_bu.get("residual_evalue", 1e-3))
BU_SATIVA_MANIFEST = _bu.get("sativa_manifest", "results/MSA/bottomup/sativa_manifest.json")
BU_SATIVA_MIN_SEQS = int(_bu.get("sativa_min_seqs", 4))
BU_SATIVA_MIN_TAXA = int(_bu.get("sativa_min_taxa", 2))
BU_CLIMBED_MANIFEST = _bu.get("climbed_manifest", "results/MSA/bottomup/climbed_manifest.json")
BU_CLIMB_RANKS = _bu.get("climb_ranks", ["order", "class", "phylum"])
BU_CONTAM_REPORT = _bu.get("contamination_report", "results/MSA/bottomup/contamination_report.json")
BU_CONTAM_RANKS = _bu.get("contamination_flag_ranks", ["kingdom", "phylum"])
BU_CONTAM_MIN_SCORE = float(_bu.get("contamination_min_score", 80.0))
BU_CONTAM_TOPK = int(_bu.get("contamination_top_k", 5))
BU_CONTAM_DIR = _bu.get("contamination_dir", "results/MSA/bottomup/contam")
BU_CONTAM_BATCH = int(_bu.get("contamination_batch_size", 200))
BU_CONTAM_MAX_EVALUE = float(_bu.get("contamination_max_evalue", 1e-5))
BU_CLUSTER_ID = float(_bu.get("cluster_id", 0.97))
BU_GATE_RATIO = float(_bu.get("gate_ratio", 2.5))
BU_GATE_OCC = float(_bu.get("gate_occ", 40.0))
BU_MAX_GROUP = int(_bu.get("max_group", 300))
BU_MIN_REF = int(_bu.get("min_ref", 2))
BU_BARCODES = config["msa"].get("barcodes", ["ITS", "ITS1", "ITS2"])


def _bu_families():
    manifest = checkpoints.bottomup_split.get().output.manifest
    with open(manifest, "r", encoding="utf-8") as fh:
        return json.load(fh)["families"]


def _bu_lookup(field):
    def fn(wildcards):
        for fam in _bu_families():
            if fam["fam_id"] == wildcards.fam_id:
                return fam[field]
        raise ValueError(f"fam_id {wildcards.fam_id} absent from the bottom-up manifest.")
    return fn


def bu_family_outputs(_wildcards):
    return [os.path.join(BU_OUT, fam["fam_id"], "manifest.json") for fam in _bu_families()]


checkpoint bottomup_split:
    input:
        csv=FULL_CSV
    output:
        manifest=BU_MANIFEST
    params:
        families_dir=BU_FAMILIES_DIR,
        barcodes=" ".join(BU_BARCODES)
    log:
        "logs/MSA/bottomup/split.log"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/split_families.py \
            --input-csv {input.csv} \
            --families-dir {params.families_dir} \
            --output-manifest {output.manifest} \
            --barcodes {params.barcodes} \
            --log-file {log}
        """


rule build_family:
    input:
        csv=_bu_lookup("csv_path"),
        manifest=BU_MANIFEST
    output:
        manifest=os.path.join(BU_OUT, "{fam_id}", "manifest.json")
    params:
        barcode=_bu_lookup("barcode"),
        family=_bu_lookup("family"),
        out_dir=os.path.join(BU_OUT, "{fam_id}"),
        cluster_id=BU_CLUSTER_ID,
        gate_ratio=BU_GATE_RATIO,
        gate_occ=BU_GATE_OCC,
        max_group=BU_MAX_GROUP,
        min_ref=BU_MIN_REF
    threads: int(_bu.get("threads_per_family", 2))
    log:
        "logs/MSA/bottomup/build/{fam_id}.log"
    benchmark:
        "benchmarks/msa/bottomup/{fam_id}.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/build_family_msa.py \
            --input-csv {input.csv} \
            --barcode {params.barcode} \
            --family '{params.family}' \
            --fam-id {wildcards.fam_id} \
            --out-dir {params.out_dir} \
            --cluster-id {params.cluster_id} \
            --gate-ratio {params.gate_ratio} \
            --gate-occ {params.gate_occ} \
            --max-group {params.max_group} \
            --min-ref {params.min_ref} \
            --log-file {log}
        """


rule climb_lineage:
    """UNBOUNDED climb (gather): merges family MSAs -> order -> class -> phylum as long
    as the gate passes (never cross-barcode). Runs AFTER the residual step, so on
    already-enriched family MSAs (residual + re-homing integrated)."""
    input:
        families=bu_family_outputs,
        residual=BU_RESIDUAL_REPORT
    output:
        manifest=BU_CLIMBED_MANIFEST
    params:
        out_dir=BU_OUT,
        climb_ranks=" ".join(BU_CLIMB_RANKS),
        gate_ratio=BU_GATE_RATIO,
        gate_occ=BU_GATE_OCC
    log:
        "logs/MSA/bottomup/climb.log"
    benchmark:
        "benchmarks/msa/bottomup/climb.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/climb_lineage.py \
            --family-manifests {input.families} \
            --out-dir {params.out_dir} \
            --output-manifest {output.manifest} \
            --climb-ranks {params.climb_ranks} \
            --gate-ratio {params.gate_ratio} \
            --gate-occ {params.gate_occ} \
            --log-file {log}
        """


rule bottomup_collect:
    input:
        families=bu_family_outputs,
        manifest=BU_MANIFEST,
        climbed=BU_CLIMBED_MANIFEST,       # final MSAs (unbounded climb)
        contamination=BU_CONTAM_REPORT,    # after re-homing (emptied MSAs to exclude)
        residual=BU_RESIDUAL_REPORT        # re-inserted orphans (per_msa) -> summary
    output:
        summary=BU_SUMMARY,
        sativa_manifest=BU_SATIVA_MANIFEST
    params:
        min_seqs=BU_SATIVA_MIN_SEQS,
        min_taxa=BU_SATIVA_MIN_TAXA
    log:
        "logs/MSA/bottomup/collect.log"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/collect_bottomup.py \
            --family-manifests {input.families} \
            --split-manifest {input.manifest} \
            --climbed-manifest {input.climbed} \
            --output-summary {output.summary} \
            --sativa-manifest {output.sativa_manifest} \
            --sativa-min-seqs {params.min_seqs} \
            --sativa-min-taxa {params.min_taxa} \
            --contamination-report {input.contamination} \
            --residual-report {input.residual} \
            --log-file {log}
        """


# Coarse contamination (cross-kingdom/phylum) is O(M^2) — one rep per MSA searched
# against every profile. It scatters: prepare (reps.fa + profile-library shards) ->
# search (one hmmsearch per shard against all reps, parallel) -> report (flag + re-home).
# Sharding the PROFILES (hmmsearch iterates queries) processes each profile exactly once
# across shards; results are identical to the single-job version.
checkpoint contam_prepare:
    """Build reps.fa (shared target DB), profile-library shards, and an index."""
    input:
        families=bu_family_outputs
    output:
        reps=os.path.join(BU_CONTAM_DIR, "reps.fa"),
        index=os.path.join(BU_CONTAM_DIR, "index.json"),
        shards=directory(os.path.join(BU_CONTAM_DIR, "shards"))
    params:
        batch_size=BU_CONTAM_BATCH
    log:
        "logs/MSA/bottomup/contam_prepare.log"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/coarse_contamination.py prepare \
            --family-manifests {input.families} \
            --out-reps {output.reps} \
            --out-index {output.index} \
            --shards-dir {output.shards} \
            --batch-size {params.batch_size} \
            --log-file {log}
        """


def contam_hit_batches(_wildcards):
    shards_dir = checkpoints.contam_prepare.get().output.shards
    shards = glob_wildcards(os.path.join(shards_dir, "shard_{b}.hmm")).b
    return [os.path.join(BU_CONTAM_DIR, "hits", f"{b}.tbl") for b in shards]


rule contam_search:
    """One hmmsearch per profile-library shard against all reps (scatter)."""
    input:
        shard=os.path.join(BU_CONTAM_DIR, "shards", "shard_{b}.hmm"),
        reps=os.path.join(BU_CONTAM_DIR, "reps.fa")
    output:
        # transient: consumed only by the bottomup_contamination gather
        tbl=temp(os.path.join(BU_CONTAM_DIR, "hits", "{b}.tbl"))
    params:
        max_evalue=BU_CONTAM_MAX_EVALUE
    threads: int(_bu.get("contamination_search_threads", 4))
    log:
        "logs/MSA/bottomup/contam_search/{b}.log"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        hmmsearch --tblout {output.tbl} -E {params.max_evalue} --cpu {threads} --noali \
            {input.shard} {input.reps} > {log} 2>&1
        """


rule bottomup_contamination:
    """Gather: parse the batch hits, flag label-vs-homology disagreements, re-home."""
    input:
        hits=contam_hit_batches,
        index=os.path.join(BU_CONTAM_DIR, "index.json")
    output:
        report=BU_CONTAM_REPORT
    params:
        ranks=" ".join(BU_CONTAM_RANKS),
        min_score=BU_CONTAM_MIN_SCORE,
        top_k=BU_CONTAM_TOPK
    threads: int(_bu.get("contamination_threads", 8))
    log:
        "logs/MSA/bottomup/contamination.log"
    benchmark:
        "benchmarks/msa/bottomup/contamination.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/coarse_contamination.py report \
            --hits {input.hits} \
            --index {input.index} \
            --output-report {output.report} \
            --flag-ranks {params.ranks} \
            --min-score {params.min_score} \
            --top-k {params.top_k} \
            --threads {threads} \
            --log-file {log}
        """


rule bottomup_residual:
    input:
        families=bu_family_outputs,
        contamination=BU_CONTAM_REPORT   # the re-homing (MSA mutation) must run first
    output:
        report=BU_RESIDUAL_REPORT
    params:
        max_evalue=BU_RESIDUAL_EVALUE
    threads: int(_bu.get("residual_threads", 8))
    log:
        "logs/MSA/bottomup/residual.log"
    benchmark:
        "benchmarks/msa/bottomup/residual.benchmark.txt"
    conda:
        "../envs/msa.yaml"
    shell:
        """
        python workflow/scripts/MSA/residual_classify_integrate.py \
            --family-manifests {input.families} \
            --output-report {output.report} \
            --max-evalue {params.max_evalue} \
            --threads {threads} \
            --log-file {log}
        """
