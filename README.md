<div align="center">

# ITS Bioinformatics Pipeline

Automated Snakemake workflow for ingesting, curating, and validating ITS (Internal Transcribed Spacer) barcode sequences from heterogeneous sources.  

[Snakemake](https://snakemake.readthedocs.io) • [MkDocs Documentation](docs/) • [GPLv3](LICENSE)

</div>

---

## Overview

This repository hosts a production-ready pipeline that consolidates ITS barcode data (ITS, ITS1, ITS2) into a unified relational database, applies domain-specific cleaning, and produces harmonized analytical outputs. The workflow is orchestrated with Snakemake and leverages Conda environments for full reproducibility.

Core capabilities include:

- **Ingestion & normalization** of public datasets (BOLD, PLANITS, CALeDNA, UNITE, etc.) using high-throughput parsing scripts.
- **Quality control & deduplication** at both sequence and taxonomy levels before database insertion.
- **ITSx processing** to extract marker-specific subsequences and update barcode status tracking.
- **Multiple sequence alignment (MSA)** per taxonomic slice with MAFFT/Trimal and downstream summarization.
- **SATIVA-based taxonomic validation** plus automated reporting and optional database integration of reassigned taxa.
- **Rich observability**, including JSON run reports, HTML dashboards, and extensive logging.

The pipeline runs end-to-end against either SQLite (local prototyping) or PostgreSQL (large-scale production) backends with transparent configuration switches.

---

## Repository Layout

```
.
├── config/                 # Global configuration and environment templates
├── data/                   # Raw and processed datasets (git-ignored if large)
├── docs/                   # MkDocs documentation portal
├── resources/              # Reference files (e.g., taxonomy definitions)
├── workflow/
│   ├── Snakefile           # Primary Snakemake entry point
│   ├── envs/               # Conda environments per tool (ITSx, MAFFT, SATIVA, etc.)
│   ├── report/             # HTML templates and assets used in reports
│   ├── rules/              # Snakemake rule modules (download, itsx, msa, sativa, reports)
│   └── scripts/            # Python utilities (formatters, DB ops, report generators)
├── results/                # Generated outputs (alignments, SATIVA runs, augmented datasets)
├── tests/                  # Unit/integration tests (optional)
├── mkdocs.yml              # Documentation site configuration
└── setup_env.sh            # Helper script for local environment bootstrap
```

---

## Requirements

- Linux or macOS with Bash, GNU coreutils, and Conda (Miniconda or Mamba recommended)
- Snakemake ≥ 7.0 (automatically installed via conda/mamba)
- Access credentials for the selected database backend (SQLite path or PostgreSQL DSN)
- Optional: MkDocs + `mkdocs-material` for browsing the docs locally

---

## Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/<org>/ITS-pipeline.git
   cd ITS-pipeline/6.CurateMake
   ```

2. **Create/update Conda environments**
   ```bash
   chmod +x setup_env.sh
   ./setup_env.sh          # installs Snakemake and required tool envs
   ```

3. **Configure credentials**
   ```bash
   cp config/.env.example config/.env
   # edit config/.env to define DB connection, API tokens, etc.
   ```

4. **Review primary configuration**
   - `config/config.yaml` centralizes dataset paths, resource usage, and tool parameters.
   - `workflow/envs/*.yaml` pin the versions of ITSx, MAFFT, SATIVA, hmmer, etc.

5. **Dry-run the workflow**
   ```bash
   snakemake -np            # prints planned jobs without executing them
   ```

6. **Execute the full pipeline**
   ```bash
   snakemake --use-conda --cores <N>
   # or submit to a cluster with --profile
   ```

Snakemake automatically creates the per-tool Conda environments declared in each rule, so no manual activation is required during workflow execution.

---

## Configuration Highlights

| Section | Purpose |
|--------|---------|
| `database` | Toggle between SQLite (default) and PostgreSQL, set schema, connection string, and optional backups. |
| `data_sources` | Define file paths for BOLD/PLANITS/CALeDNA/UNITE dumps and ingestion thresholds. |
| `itsx` | CPU/thread allocation, runtime limits, and post-processing flags for ITSx batches. |
| `msa` | MAFFT/Trimal parameters, batch sizing, alignment thresholds, and reporting options. |
| `sativa` | SATIVA deployment details (tree limits, API integrations, visualization paths). |
| `reports` | Output locations for HTML dashboards and JSON metrics. |

Every parameter is documented inline within `config/config.yaml`. Override values per environment by exporting variables in `config/.env` or via `--configfile`.

---

## Running Specific Stages

```bash
# Format and ingest CALeDNA only
snakemake --use-conda results/caledna_ingestion.done

# Launch ITSx processing for ITS2
snakemake --use-conda itsx/ITS2/itsx_run.done

# Generate SATIVA visualization report
snakemake --use-conda reports/html/sativa_visualization_report.html
```

Use `-n` (dry-run) and `-r` (print reason) to inspect dependencies before executing resource-intensive jobs.

---

## Reporting & Monitoring

- **JSON run reports** (`data/processed/*_run_report.json`) capture per-stage metrics, throughput, and warnings.
- **HTML dashboards** in `reports/html/` summarize ITSx/MSA/SATIVA runs for stakeholder review.
- **Logs** are stored in `logs/<module>/` with human-readable English messaging for quicker troubleshooting.
- **Database utilities** under `workflow/scripts/utils/` (e.g., `restore_db_backup.py`) validate and restore curated DB snapshots.

---

## Development Workflow

1. Create a feature branch from `CurateMake` (or `main`, depending on your release strategy).
2. Apply linting/formatting to scripts as needed (black/isort/mypy optional).
3. Run targeted Snakemake jobs to validate changes before committing.
4. Submit pull requests with links to relevant run reports or logs.

The repository includes `.gitignore` entries for large datasets, secrets, editor files, and Conda artifacts to keep commits clean.

---

## Documentation

The MkDocs site compiles the architecture overview, database schema notes, and troubleshooting guides located in `docs/`. Serve it locally with:

```bash
pip install mkdocs mkdocs-material
mkdocs serve
# open http://127.0.0.1:8000
```

---

## Support & Contributions

- File issues or feature requests via GitHub Issues (include Snakemake logs and relevant config snippets).
- Contributions are welcome through standard fork/PR workflows. Please keep README and docs synchronized with code changes.
- For security-sensitive topics (credentials, datasets under embargo), coordinate via the appropriate private channels instead of public trackers.

---

**License:** GPLv3. This pipeline bundles third-party tools (ITSx, MAFFT, SATIVA, hmmer, etc.); ensure compliance with their respective licenses when redistributing artifacts.
