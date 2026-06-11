# Usage Guide

This guide summarises the most common operations once installation is complete.
Run every command from the project root with the `curatemake` environment
activated.

## 1. Command templates

```bash
# Inspect the DAG without executing jobs
snakemake -np

# Run with N local cores (per-rule Conda envs are built on first use)
snakemake --use-conda --cores 32

# Build a specific target (e.g. the SATIVA taxonomy report)
snakemake --use-conda reports/html/sativa_taxonomy_report.html

# Force a rule (and its downstream) to rerun
snakemake --use-conda -R sativa_task
```

Useful global flags:

| Flag | Description |
| ---- | ----------- |
| `-n` | dry-run (plan only) |
| `-r` | show the reason for each job |
| `-p` | print shell commands |
| `--summary` | status of every output file |
| `--keep-going` | continue with independent jobs after a failure |
| `--resources sativa_jobs=10` | cap the custom `sativa_jobs` resource |

## 2. Primary targets

| Target | What it produces |
| ------ | ---------------- |
| `data/processed/formatting_summary.json` | Ingestion/formatting summary |
| `data/processed/col_matching_report.json` | Catalogue of Life harmonisation |
| `data/processed/itsx_process.json` | Aggregated ITSx summary |
| `data/processed/msa_run_report.json` | Run-level MSA metrics |
| `data/processed/sativa_analysis_report.json` | Aggregated SATIVA results |
| `reports/html/msa_analysis_report.html` | Interactive MSA analytics |
| `reports/html/sequence_flow_report.html` | Pipeline-level sequence flow |
| `reports/html/sativa_taxonomy_report.html` | 4-way taxonomy comparison |
| `reports/html/index.html` | Landing page linking all reports |

Requesting any target materialises only the upstream steps it needs.

## 3. Partial workflows

### Database ingestion only
```bash
snakemake --use-conda format_all_databases
```

### Taxonomy harmonisation (Catalogue of Life)
```bash
snakemake --use-conda data/processed/col_matching_report.json
```

### ITSx extraction for one barcode
```bash
snakemake --use-conda data/processed/itsx_run_ITS1.done
```

### MSA generation only
```bash
snakemake --use-conda data/processed/msa_run_report.json
```

### SATIVA validation without the final reports
```bash
snakemake --use-conda data/processed/sativa_analysis_report.json
```

## 4. Monitoring & troubleshooting

- Logs are stored hierarchically under `logs/<module>/...`
  (e.g. `logs/ITSX/run_itsx_ITS2.log`, `logs/MSA/align/...`).
- Benchmarks (runtime/memory) land in `benchmarks/<rule>/...` when
  `benchmarking.enabled` is true in `config/config.yaml`.
- `data/processed/*_report.json` and `data/processed/*_manifest.json` expose the
  structured metrics that feed the HTML reports.
- To debug failing tasks: `snakemake --use-conda --cores <N> --keep-going -p`
  prints shell commands; cross-reference the matching log file.

## 5. Housekeeping

- Remove processed artefacts only:
  ```bash
  snakemake clean --cores 1
  ```
- Remove everything (processed data, ITSx outputs, results, logs):
  ```bash
  snakemake clean_all --cores 1
  ```
  Use with caution; these rules delete the directories listed in the `Snakefile`.

## 6. Cluster execution

Configure a Snakemake profile (e.g. SLURM) and run:

```bash
snakemake --use-conda --profile profiles/slurm
```

Map cluster resources to the rule-level `threads:` and custom resources (e.g.
`sativa_jobs`, `db_lock`) declared in the rule files.
