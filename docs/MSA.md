# MSA Stage

The MSA rules (`workflow/rules/msa.smk`) orchestrate everything between post-ITSx deduplication and the SATIVA manifest. Jobs are derived from manifests so only the necessary alignments run.

## Manifest generation

1. `prepare_msa` exports a CSV (`config["msa"]["csv_path"]`) containing curated sequences and taxonomy metadata directly from the DB.
2. `checkpoint msa_task_manifest` reads that CSV and produces `data/processed/msa_task_manifest.json`, grouping sequences by barcode and taxonomic level according to:
   - `msa.barcodes`
   - `msa.hierarchy_levels`
   - thresholds (`min_seq_threshold`, `max_seq_threshold`, `large_task_threshold`)

Snakemake checkpoints ensure downstream rules dynamically enumerate tasks from this manifest.

## Alignment + trimming

For each `(barcode, level, clade)` triplet:

1. `msa_align_task` runs `workflow/scripts/MSA/run_alignment_task.py`, which:
   - filters the CSV for the clade
   - calls MAFFT with the configured tax columns
   - writes `<clade>_msa.fasta`
2. `msa_trim_task` trims each alignment with trimAl:
   - `msa.trimal_seq_overlap`
   - `msa.trimal_res_overlap`
   Outputs: `.trimmed.fasta` + JSON describing removed sequences.

Sentinel rules (`msa_alignments_ready`, `msa_trimming_ready`) touch the flags defined in `msa.stage_markers.*` so later stages can wait on entire phases without re-evaluating every file.

## HMM assembly + classification

1. `msa_build_hmm_task` transforms each trimmed MSA into a profile HMM (via `workflow/scripts/MSA/build_hmm_task.py`).
2. `msa_combine_barcode_hmms` aggregates clade HMMs into `<barcode>/all_profiles.hmm`.
3. `msa_combine_global_hmm` concatenates every barcode-level HMM file into `msa.global_hmm_raw`, followed by `msa_filter_global_hmm` (HMMER hmmstat filter).
4. `checkpoint msa_classification_manifest` splits sequences needing classification into batches defined by:
   - `classification_batch_size`
   - `classification_threads_per_worker`
   - `classification_batches_dir`
5. `msa_classify_batch` invokes `workflow/scripts/MSA/classify_sequences_task.py` (hmmscan) per batch.
6. `msa_collect_assignments` merges batch outputs and tracks sequences still unresolved.
7. `msa_integrate_assignments` updates the DB and summarises stats in `config["msa"]["assignment_integration_summary"]`.

## Reporting

`run_msa` is the umbrella `script:` job (`workflow/scripts/MSA/run_msa.py`) that compiles end-to-end metrics into `data/processed/msa_run_report.json`. `analyze_msa_report` renders `reports/html/msa_analysis_report.html`.

## Key configuration fields

```yaml
msa:
  csv_path: "data/full.csv"
  msa_dir: "data/MSA"
  barcodes: ["ITS","ITS1","ITS2"]
  hierarchy_levels: ["genus","family","order","class","phylum"]
  tax_columns: ["kingdom","phylum","class","order","family","genus","species"]
  min_seq_threshold: 20
  max_seq_threshold: 1000
  large_task_threshold: 2000
  trimal_seq_overlap: 75
  trimal_res_overlap: 0.5
  classification_threads_per_worker: 4
  classification_max_workers: 20
```

Adjusting these values controls manifest size, trimming aggressiveness, and hmmscan throughput.

## Tips

- Rebuild just one alignment:
  ```bash
  snakemake --use-conda data/MSA/ITS/alignments_family/Asteraceae_msa.fasta
  ```
- Inspect manifest content:
  ```bash
  jq '.' data/processed/msa_task_manifest.json | less
  ```
- Use `data/processed/msa_assignment_integration.json` to audit which sequences were classified vs deferred.
