# workflow/rules/ehmm_msa.smk
"""
Backbone + ensemble-HMM MSA construction — the only MSA strategy.

Emits the artifacts sativa.smk / reports.smk / msa.smk consume, exactly the same
artifacts, so sativa.smk / reports.smk are untouched:
  <CL_OUT>/{barcode}/sativa_manifest.json   then the shared collect_groups merges them.

Why (measured, 10.DOI_CurateMake/simulation, CoL-harmonised s8/r15 point, prototype in
20.eHMM_MSA_prototype):
  * a single backbone alignment gives EVERY read the same column frame, which is what lets SATIVA
    move a contaminant across kingdoms (the property the from-scratch MSA also has, but which
    per-clade schemes destroy);
  * the per-group MAFFT disappears: groups become a free slice of an already-aligned matrix, so
    re-cutting a group costs seconds instead of a realignment;
  * arm-level scores against the same input: eHMM 57.3% cross-kingdom / 7.1% false positives vs
    47.6% / 20.1% for the from-scratch MAFFT path it replaced, and a 20k task ran in ~33 min instead of 1-2 h.

Steps: build_backbone (stratified sample, aligned once) -> ehmm_project (per-clade HMMs over the
backbone frame, nhmmer assigns bucket + strand, hmmalign, keep match columns) -> recover_unassigned
(mmseqs rescue for short reads nhmmer misses) -> slice_and_pack (kingdom round-robin groups).
"""
import os

FULL_CSV = config["msa"]["csv_path"]
EH_BARCODES = config["msa"].get("barcodes", ["ITS", "ITS1", "ITS2"])

_eh = config["msa"].get("ehmm", {})
EH_WORK = _eh.get("work_dir", "results/MSA/ehmm")
EH_TARGET = int(_eh.get("backbone_target", 1200))
EH_BB_SOURCE = _eh.get("backbone_source", "cluster")   # cluster reps cover sequence space
EH_BB_CID = float(_eh.get("backbone_cluster_id", 0.90))
EH_MIN_STRATUM = int(_eh.get("backbone_min_per_stratum", 3))
EH_MAFFT = _eh.get("backbone_aligner", "mafft-linsi")
EH_MIN_COL_OCC = float(_eh.get("min_col_occ", 0.10))
EH_BUCKET_MODE = _eh.get("bucket_mode", "upp")   # upp = many small profiles (see make_buckets)
EH_SUBSET = int(_eh.get("subset_size", 25))
EH_BUCKET_RANK = _eh.get("bucket_rank", "kingdom")
EH_MIN_BUCKET = int(_eh.get("min_bucket", 20))
EH_RECOVER = bool(_eh.get("recover_unassigned", True))
EH_REC_ID = float(_eh.get("recover_min_id", 0.75))
EH_REC_COV = float(_eh.get("recover_min_qcov", 0.60))
EH_SCAFFOLD = bool(_eh.get("scaffold_every_group", True))
EH_HMMER_BIN = _eh.get("hmmer_bin", "")          # empty -> take hmmbuild/nhmmer/hmmalign from PATH
EH_THREADS = int(_eh.get("threads", 32))
EH_NHMMER_JOBS = int(_eh.get("nhmmer_jobs", 4))
EH_ALIGN_JOBS = int(_eh.get("hmmalign_jobs", 16))   # hmmalign is single-threaded
EH_NHMMER_THREADS = int(_eh.get("nhmmer_threads", 0)) or EH_THREADS
# What the rule really occupies at its widest: the nhmmer fan-out and the hmmalign fan-out are
# successive phases, and mafft/mmseqs use EH_THREADS. Declaring the max keeps Snakemake from
# scheduling anything else on top of it.
EH_RULE_THREADS = max(EH_THREADS, EH_NHMMER_JOBS * EH_NHMMER_THREADS, EH_ALIGN_JOBS)   # profiles scanned in parallel

# the size cap and the SATIVA floor live under msa.groups: same meaning, one place to change
_cl = config["msa"].get("groups", {})
EH_OUT = _cl.get("out_dir", "results/MSA/groups")
EH_MAX_GROUP = int(_cl.get("max_group", 20000))
EH_MIN_GROUP = int(_cl.get("sativa_min_seqs", 4))
EH_MMSEQS_ENV = _cl.get("mmseqs_env", "mmseqs")
EH_TOOL_ENV = _cl.get("tool_env", "its-pipeline-dev")

_SCRIPTS = "workflow/scripts/MSA/ehmm"


# ITS1 and ITS2 are exact substrings of the ITS read they were cut from (checked: 200/200
# both ways), and ITS projects ~100% of its reads where ITS1/ITS2 reach ~90%. Feeding ITS's
# projected reads to their rescue step is what keeps the ~267k reads the 2026-08-25 run
# dropped (145961 ITS1 + 121206 ITS2, written to unprojected_seq_ids.txt and never seen by
# SATIVA). ITS itself gets no cross-reference: nothing above it, and nothing to gain.
# ITS1 <-> ITS2 is deliberately absent -- they sit either side of 5.8S and share no sequence.
_EH_CROSS_SOURCE = "ITS"


def _ehmm_cross_ref(wildcards):
    if wildcards.barcode == _EH_CROSS_SOURCE:
        return []
    return [os.path.join(EH_OUT, _EH_CROSS_SOURCE, "projected_refs.fasta")]


# --- consensus-vote replicates ---------------------------------------------------------
# `sativa.replicates: N` asks for N independent splits of the SAME projection. Seed 0 is the
# native split produced above; seeds 1..N-1 are re-cuts. Only the CHEAP half of eHMM is
# repeated: the backbone, nhmmer and hmmalign work depends on the barcode alone and survives
# as the per-bucket a_*.afa files, so a replicate just re-reads them and re-slices.
# Measured 2026-08-26 on ITS2: 11 min for the re-cut against ~15 h for a full projection.
# One job per barcode emits ALL the extra seeds, because reading the projection back
# (7.8 GB for ITS2, 96 GB for ITS) is the entire cost and must not be paid N times.
EH_REPLICATES = int(config.get("sativa", {}).get("replicates", 1))
EH_REPL_ROOT = _eh.get("replicate_dir", "results/MSA/replicates")
EH_REPL_SEEDS = list(range(1, EH_REPLICATES))


if EH_REPL_SEEDS:

    rule ehmm_reslice_replicates:
        """Re-cut an already projected barcode into the N-1 extra splits, in one pass."""
        input:
            manifest=os.path.join(EH_OUT, "{barcode}", "sativa_manifest.json"),
            csv=FULL_CSV,
        output:
            manifests=expand(
                os.path.join(EH_REPL_ROOT, "{{barcode}}", "seed{seed}", "sativa_manifest.json"),
                seed=EH_REPL_SEEDS,
            ),
        params:
            work=lambda w: os.path.join(EH_WORK, w.barcode),
            root=lambda w: os.path.join(EH_REPL_ROOT, w.barcode),
            seeds=" ".join(str(s) for s in EH_REPL_SEEDS),
            scaffold=("--scaffold" if EH_SCAFFOLD else ""),
        threads: 1
        log:
            "logs/ehmm/reslice_{barcode}.log",
        benchmark:
            "benchmarks/ehmm/reslice_{barcode}.benchmark.txt"
        conda:
            "../envs/msa.yaml"
        shell:
            r"""
            python {_SCRIPTS}/reslice_ehmm.py \
                --work-dir {params.work} --full-csv {input.csv} \
                --barcode {wildcards.barcode} --out-dir {params.root} \
                --split-seed {params.seeds} \
                --max-group {EH_MAX_GROUP} --min-group {EH_MIN_GROUP} \
                {params.scaffold} > {log} 2>&1
            """


rule ehmm_build_barcode:
    """One barcode: backbone -> eHMM projection -> (rescue) -> kingdom-interleaved groups."""
    input:
        csv=FULL_CSV,
        cross=_ehmm_cross_ref,
    output:
        manifest=os.path.join(EH_OUT, "{barcode}", "sativa_manifest.json"),
        refs=os.path.join(EH_OUT, "{barcode}", "projected_refs.fasta"),
    params:
        work=lambda w: os.path.join(EH_WORK, w.barcode),
        outdir=lambda w: os.path.join(EH_OUT, w.barcode),
        recover=("--recover" if EH_RECOVER else ""),
        scaffold=("--scaffold" if EH_SCAFFOLD else ""),
        hmmer_bin=(f"--hmmer-bin {EH_HMMER_BIN}" if EH_HMMER_BIN else ""),
        cross=lambda w, input: " ".join(f"--cross-ref {p}" for p in input.cross),
    threads: EH_RULE_THREADS
    log:
        "logs/ehmm/{barcode}.log",
    benchmark:
        "benchmarks/ehmm/build_{barcode}.benchmark.txt"
    conda:
        # hmmbuild / nhmmer / hmmalign live here (envs/msa.yaml declares hmmer=3.4);
        # they are on no default PATH, so without this the first hmmbuild call fails.
        "../envs/msa.yaml"
    shell:
        r"""
        python {_SCRIPTS}/build_ehmm_msa.py \
            --full-csv {input.csv} --barcode {wildcards.barcode} \
            --work-dir {params.work} --out-dir {params.outdir} \
            --backbone-target {EH_TARGET} --backbone-min-per-stratum {EH_MIN_STRATUM} \
            --backbone-source {EH_BB_SOURCE} --backbone-cluster-id {EH_BB_CID} \
            --backbone-aligner {EH_MAFFT} --min-col-occ {EH_MIN_COL_OCC} \
            --bucket-mode {EH_BUCKET_MODE} --subset-size {EH_SUBSET} \
            --bucket-rank {EH_BUCKET_RANK} --min-bucket {EH_MIN_BUCKET} \
            --recover-min-id {EH_REC_ID} --recover-min-qcov {EH_REC_COV} \
            --max-group {EH_MAX_GROUP} --min-group {EH_MIN_GROUP} \
            --threads {threads} --nhmmer-jobs {EH_NHMMER_JOBS} \
            --hmmalign-jobs {EH_ALIGN_JOBS} --nhmmer-threads {EH_NHMMER_THREADS} \
            --mmseqs-env {EH_MMSEQS_ENV} --tool-env {EH_TOOL_ENV} \
            {params.cross} {params.recover} {params.scaffold} {params.hmmer_bin} > {log} 2>&1
        """
