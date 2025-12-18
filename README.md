<div align="center">

# ITS Bioinformatics Pipeline

Snakemake workflow for curating ITS/ITS1/ITS2 barcodes, aligning them, validating taxonomy with SATIVA, and emitting reproducible reports backed by SQLite or PostgreSQL.

[Snakemake](https://snakemake.readthedocs.io) • [MkDocs Docs](docs/index.md) • [GPLv3](LICENSE)

</div>

---

## What the Pipeline Does

1. **Database bootstrap & ingestion** – loads UNITE, BOLD, PLANiTS, CALeDNA, etc. into a normalized schema with sequence/taxonomy QC (`workflow/rules/database.smk`, `format_databases.smk`).
2. **ITSx extraction** – prepares filtered FASTA files, runs ITSx per barcode, and updates barcode status tracking (`workflow/rules/itsx.smk`).
3. **MSA production** – generates per-clade MSAs, trims with Trimal, builds HMMs, and classifies remaining sequences (`workflow/rules/msa.smk`).
4. **SATIVA validation** – creates SATIVA tasks from the MSA manifest, runs the classifier, and optionally writes the results back to the DB (`workflow/rules/sativa.smk`).
5. **Reporting** – compiles JSON run summaries plus HTML dashboards (sequence flow, MSA analytics, SATIVA taxonomy) (`workflow/rules/reports.smk`, `rules/benchmarking.smk`).

All steps are orchestrated by Snakemake with per-rule Conda environments under `workflow/envs/`.

---

## Requirements

- Linux/macOS shell with Bash, GNU coreutils, and either Conda or Mamba
- Snakemake ≥7 (installed automatically by `setup_env.sh`)
- Database credentials (SQLite path or PostgreSQL DSN defined in `config/.env`)
- Optional: MkDocs + `mkdocs-material` to browse `docs/`

---

## Quick Start

1. **Clone & bootstrap**
   ```bash
   git clone https://github.com/<org>/ITS-pipeline.git
   cd ITS-pipeline/6.CurateMake
   ./setup_env.sh            # installs snakemake + rule envs
   ```
2. **Configure credentials & paths**
   ```bash
   cp config/.env.example config/.env   # edit DB + API secrets
   cp config/config.yaml config/local.yaml  # optional custom config
   ```
3. **Dry-run then execute**
   ```bash
   snakemake -np                                  # plan
   snakemake --use-conda --cores <N>              # run locally
   # or snakemake --profile <cluster_profile>
   ```

Snakemake creates rule-specific Conda envs on demand; no manual activation is necessary.

---

## Configuration Cheatsheet

| Section (config/config.yaml) | Key knobs |
| --- | --- |
| `database` | Backend selection, schema, backup policy |
| `databases.*` | Paths to UNITE/BOLD/PLANiTS/CALeDNA dumps |
| `itsx_*` / `itsx` | Thread counts, temp directories, CPU limits |
| `msa` | Thresholds for manifest generation, MAFFT/Trimal params, classification batching |
| `sativa` | Input MSA dir, timeout, max sequences, SATIVA script path, visualization outputs |
| `sativa_integration` | CoL API tuning, cache, reporting, DB write-back |

Override any parameter via `snakemake --config key=value` or by pointing to an alternate `--configfile`.

---

## Running Targeted Stages

```bash
# Only run ITSx for ITS2
snakemake --use-conda data/processed/itsx_run_ITS2.done

# Rebuild the global MSA report
snakemake --use-conda reports/html/msa_analysis_report.html

# Recompute SATIVA taxonomy report
snakemake --use-conda reports/html/sativa_taxonomy_report.html
```

Use `-n` (dry-run), `-r` (reasons), and `--summary` to inspect dependency graphs before launching long jobs.

---

## Outputs & Observability

- **Processed artifacts** live under `data/processed/` (JSON summaries, manifests, CSV exports).
- **MSA and SATIVA deliverables** are emitted to `data/MSA` and `data/SATIVA`.
- **Reports** render to `reports/html/`:
  - `msa_analysis_report.html`
  - `sequence_flow_report.html`
  - `sativa_visualization_report.html`
  - `sativa_taxonomy_report.html`
  - `benchmark_report.html` and `index.html`
- **Logs** sit in `logs/<module>/`, with English log text suitable for monitoring systems.
- **Benchmarks** are written to `benchmarks/<rule>/` when `config["benchmarking"]["enabled"]` is true.

---

## Repository Map (abridged)

```
config/          # YAML + .env configuration
data/            # Raw inputs, processed tables, ITSx/MSA/SATIVA outputs
docs/            # MkDocs pages (see docs/index.md)
workflow/
  Snakefile      # Root workflow
  rules/         # Modular Snakemake rule files
  scripts/       # Python helpers (formatting, reports, DB ops)
  envs/          # Conda env specifications
logs/, benchmarks/, reports/  # Runtime artifacts
setup_env.sh     # Helper to install Snakemake + env deps
```

---

## Documentation & Support

- Read the MkDocs portal starting at [`docs/index.md`](docs/index.md) for deeper architecture notes, database specs, and troubleshooting.
- File issues/feature requests with sample configs + logs.
- Contributions follow the typical fork/PR workflow; keep README/docs aligned with code changes.

**License:** GPLv3 (external tools such as ITSx, MAFFT, HMMER, SATIVA retain their original licenses—respect redistribution terms).
