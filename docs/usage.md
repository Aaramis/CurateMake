# Usage Guide

This guide summarizes the most common workflows once the installation steps are complete.

## 1. Command templates

```bash
# inspect the DAG without executing jobs
snakemake --use-conda -np

# run with N local cores
snakemake --use-conda --cores 32

# run a specific target (e.g. SATIVA taxonomy report)
snakemake --use-conda reports/html/sativa_taxonomy_report.html

# rerun a rule ignoring caches
snakemake --use-conda -R sativa_task
```

Useful global flags:

| Flag | Description |
| ---- | ----------- |
| `-n` | dry-run (plan only) |
| `-r` | show reason for each job |
| `-p` | print shell commands |
| `--summary` | status of every output file |
| `--resources sativa_jobs=10` | cap custom resources if needed |

## 2. Primary targets

| Target | What it produces |
| ------ | ---------------- |
| `config["sativa"]["visualization_report"]` | SATIVA per-MSA dashboard |
| `data/processed/itsx_process.json` | Aggregated ITSx summary |
| `data/processed/msa_run_report.json` | Run-level metrics for MSA phase |
| `reports/html/msa_analysis_report.html` | Interactive MSA analytics |
| `reports/html/sequence_flow_report.html` | Pipeline-level Sankey diagram |
| `reports/html/sativa_taxonomy_report.html` | 4-way taxonomy comparison |
| `reports/html/index.html` | Landing page linking all reports |

Invoke any of these as a Snakemake target to materialize only the necessary upstream steps.

## 3. Partial workflows

### Database formatting only
```bash
snakemake --use-conda format_all_databases
```

### ITSx extraction per barcode
```bash
snakemake --use-conda data/processed/itsx_run_ITS1.done
```

### MSA generation only
```bash
snakemake --use-conda data/processed/msa_run_report.json
```

### SATIVA integration without reports
```bash
snakemake --use-conda data/processed/sativa_analysis_report.json
```

## 4. Monitoring & troubleshooting

- Logs are stored hierarchically: `logs/<module>/...` (e.g., `logs/ITSX/run_itsx_ITS2.log`).
- Benchmarks (runtime/memory) land in `benchmarks/<rule>/` if enabled via `config["benchmarking"]["enabled"]`.
- `data/processed/*run_report.json` and `data/processed/*_manifest.json` expose structured metrics that feed the HTML reports.
- To inspect failed tasks: `snakemake --use-conda --cores <N> --keep-going -p` captures shell commands; cross-reference the matching log file.

## 5. Housekeeping

- Clean processed artifacts only:
  ```bash
  snakemake clean --cores 1
  ```
- Clean everything (processed data, logs, results):
  ```bash
  snakemake clean_all --cores 1
  ```
  Use with caution; these rules delete directories defined inside the Snakefile.

## 6. Cluster execution

Prepare a Snakemake profile (e.g., SLURM) and run:

```bash
snakemake --use-conda --profile profiles/slurm
```

Adjust resource mappings in the profile to respect rule-level `threads:` and custom resources (e.g., `sativa_jobs`).
