<div align="center">

# CurateMake

A reproducible **Snakemake** workflow that ingests heterogeneous ITS/ITS1/ITS2
barcode datasets, harmonises their taxonomy, extracts barcode regions with ITSx,
builds curated multiple-sequence alignments (MSAs), validates taxonomic labels
with SATIVA, and emits fully traceable HTML/JSON reports backed by SQLite or
PostgreSQL.

[Snakemake](https://snakemake.readthedocs.io) · [Documentation](docs/index.md) · [GPLv3](LICENSE)

</div>

---

## Overview

The pipeline turns raw public reference databases (UNITE, BOLD, PLANiTS,
CALeDNA) and optional user-supplied sequences into a curated, taxonomically
validated ITS reference set. Every record is stored in a normalised relational
schema so that each sequence can be traced from ingestion through to its final,
validated taxonomic assignment.

Taxonomy is refined at **two independent stages**, both recorded in the database
via the `backbone` field of the `seq_tax` table:

1. **Catalogue of Life (CoL) harmonisation** — original labels are matched and
   reconciled against the CoL taxonomic backbone (`backbone = 'CoL'`).
2. **SATIVA phylogenetic validation** — alignment-based detection of mislabelled
   sequences proposes corrections (`backbone = 'Sativa'`).

All remaining steps (ITSx, MSA, HMM classification) are curation/infrastructure
operations that assign barcodes and produce alignments; they do not alter the
taxonomy themselves.

---

## Pipeline Stages

| # | Stage | Rule file | Key outputs |
|---|-------|-----------|-------------|
| 1 | **Database bootstrap** — create the normalised schema (SQLite or PostgreSQL) | `rules/database.smk` | `data/db/*.db` or initialised schema |
| 2 | **Ingestion & formatting** — load UNITE, BOLD, PLANiTS, CALeDNA, PERSONAL; sequence QC; taxonomy deduplication | `rules/format_databases.smk` | `data/processed/*_formatted_report.json`, `formatting_summary.json` |
| 3 | **Taxonomy harmonisation** — reconcile labels against Catalogue of Life | `rules/harmonize_taxo.smk` | `data/processed/col_matching_report.json` |
| 4 | **ITSx extraction** — extract ITS/ITS1/ITS2 regions per barcode, with deduplication before and after | `rules/itsx.smk` | `data/ITSx/<barcode>/`, `itsx_process.json` |
| 5 | **MSA production** — per-clade MAFFT alignments, trimAl trimming, HMM profiles, global HMM classification of unaligned sequences | `rules/msa.smk` | `data/MSA/`, `msa_run_report.json` |
| 6 | **SATIVA validation** — per-MSA mislabel detection and write-back of corrections | `rules/sativa.smk` | `data/SATIVA/`, `sativa_integration_report.json` |
| 7 | **Reporting** — interactive HTML dashboards and benchmark summaries | `rules/reports.smk`, `rules/benchmarking.smk` | `reports/html/*.html` |

The full directed acyclic graph (DAG) is orchestrated by Snakemake, and each
rule pins its own Conda environment under `workflow/envs/` for reproducibility.

---

## Requirements

- A Linux/macOS shell with Bash and GNU coreutils
- [Conda](https://docs.conda.io) or [Mamba](https://mamba.readthedocs.io)
  (Miniforge/Mambaforge recommended)
- Several GB of free disk for intermediate FASTA/MSA/SATIVA artefacts
- PostgreSQL client libraries **only** if you use the PostgreSQL backend
- Optional: a Snakemake cluster profile (SLURM, SGE, …) for HPC execution

All bioinformatics tools (ITSx, MAFFT, trimAl, HMMER, RAxML/SATIVA) are pulled in
automatically by the per-rule Conda environments — you do **not** install them by
hand.

---

## Environment Management

The pipeline uses a two-tier environment model:

- **One orchestration environment** (`curatemake`) that you activate. It
  contains only Snakemake and its launch-time dependencies.
- **Per-rule environments** (`workflow/envs/*.yaml`) that Snakemake creates
  on demand when `--use-conda` is enabled. These hold the heavy scientific
  tooling and keep each stage isolated and reproducible.

```bash
# 1. Create and activate the orchestration environment
conda env create -f environment.yaml      # or: mamba env create -f environment.yaml
conda activate curatemake

# 2. (Optional) install the development toolchain (linters, MkDocs, pytest)
conda env create -f dev-environment.yaml
conda activate curatemake-dev
```

`setup_env.sh` automates the same steps interactively (orchestration env,
optional dev env, optional pre-commit hooks, and `.env` bootstrap).

> **Note** — never run the workflow from inside an unrelated virtualenv. The
> `Snakefile` deliberately strips a leaking `PYTHONPATH`/`VIRTUAL_ENV` so that
> each rule resolves `python` from its own Conda environment.

---

## Quick Start

```bash
# 1. Clone and enter the project
git clone https://github.com/Aaramis/CurateMake.git
cd CurateMake

# 2. Create the orchestration environment
conda env create -f environment.yaml
conda activate curatemake

# 3. Configure credentials (only needed for the PostgreSQL backend)
cp config/.env.example .env        # then edit .env
#   The Snakefile loads .env from the project root and expands {{ env.VAR }}
#   placeholders found in config/config.yaml.

# 4. Point config/config.yaml at your input data and choose a backend
#    (default backend is SQLite, requiring no credentials).

# 5. Plan, then run
snakemake -np                          # dry-run: inspect the DAG
snakemake --use-conda --cores <N>      # execute locally
# or
snakemake --use-conda --profile <profile>   # execute on a cluster
```

The first `--use-conda` run builds the per-rule environments (a few minutes);
subsequent runs reuse them.

---

## Configuration

All tunable parameters live in [`config/config.yaml`](config/config.yaml).
Database credentials are kept out of the repository in a root `.env` file
(`config/.env.example` is the template).

| Section | Purpose |
|---------|---------|
| `database` | Backend selection (`sqlite`/`postgresql`), paths/credentials, retry & backup policy |
| `database_formatting` | Bulk-insert batch sizes and query chunking |
| `databases.*` | Input paths for UNITE, BOLD, PLANiTS, CALeDNA, and the optional PERSONAL dataset |
| `sequence_validation` | Minimum length, maximum `N` fraction, allowed characters |
| `deduplication` | Sequence/taxonomy deduplication behaviour |
| `CoL_api` / `sativa_integration.col_api` | Catalogue of Life API throughput and caching |
| `itsx_cpu` / `itsx_deterministic` | ITSx thread count and reproducible mode |
| `msa` | Alignment thresholds, MAFFT/trimAl parameters, HMM filtering, classification batching |
| `sativa` | SATIVA mode, thresholds, parallelism, visualisation output |
| `sativa_integration` | Database write-back, caching, and reporting |
| `benchmarking` | Toggle and output directory for per-rule benchmarks |

Override any value at the command line:

```bash
snakemake --use-conda --config itsx_cpu=8
# or point to an alternative file
snakemake --use-conda --configfile config/local.yaml
```

To add your own sequences, see
[`docs/personal-database.md`](docs/personal-database.md).

---

## Running Targeted Stages

Because the workflow is a DAG, requesting any output materialises only the
upstream steps it requires:

```bash
# Database ingestion and harmonisation only
snakemake --use-conda data/processed/col_matching_report.json

# ITSx extraction for a single barcode
snakemake --use-conda data/processed/itsx_run_ITS2.done

# MSA phase only
snakemake --use-conda data/processed/msa_run_report.json

# Rebuild a single report
snakemake --use-conda reports/html/sativa_taxonomy_report.html
```

Useful flags: `-n` (dry-run), `-r` (reasons), `-p` (print shell commands),
`--summary` (output status). See [`docs/usage.md`](docs/usage.md) for more.

---

## Outputs

- **Processed artefacts** — `data/processed/` (JSON run reports, manifests, CSVs)
- **Alignments & validation** — `data/MSA/`, `data/SATIVA/`
- **HTML reports** — `reports/html/`:
  - `index.html` — landing page linking every report
  - `sequence_flow_report.html` — end-to-end sequence flow
  - `msa_analysis_report.html` — MSA analytics
  - `sativa_visualization_report.html` — per-MSA SATIVA dashboard
  - `sativa_taxonomy_report.html` — 4-way taxonomy comparison
  - `benchmark_report.html` — runtime/memory benchmarks
- **Logs** — `logs/<module>/...`
- **Benchmarks** — `benchmarks/<rule>/...` (when `benchmarking.enabled` is true)

---

## Repository Layout

```
CurateMake/
├── config/
│   ├── config.yaml         # master configuration
│   └── .env.example        # credentials template (copy to ../.env)
├── workflow/
│   ├── Snakefile           # root workflow + default targets
│   ├── rules/*.smk         # modular rule sets (database, ITSx, MSA, SATIVA, reports)
│   ├── scripts/            # Python helpers (formatting, QC, MSA, SATIVA, reporting, DB utils)
│   └── envs/*.yaml         # per-rule Conda environments
├── data/                   # raw inputs, processed tables, ITSx/MSA/SATIVA outputs
├── docs/                   # MkDocs documentation set (see docs/index.md)
├── reports/, logs/, benchmarks/   # runtime artefacts
├── environment.yaml        # orchestration environment (Snakemake)
├── dev-environment.yaml    # optional development toolchain
├── setup_env.sh            # interactive environment bootstrap
└── mkdocs.yml              # documentation site configuration
```

---

## Database

Records are stored in a normalised schema shared by both backends:

| Table | Role |
|-------|------|
| `sources` | One row per input dataset (UNITE, BOLD, PERSONAL, …) |
| `sequences` | Sequences with length/`N` metrics, barcode, and processing status |
| `taxonomy` | Unique taxonomic lineages (kingdom → species) |
| `seq_tax` | Sequence↔taxonomy links, keyed by `backbone` (`Raw`, `CoL`, `Sativa`) |
| `sativa_warnings` | Flagged modifications requiring manual attention |

The `seq_tax.backbone` column is what makes the two-stage curation auditable:
the original (`Raw`) assignment, the CoL-reconciled assignment, and the
SATIVA-validated assignment coexist for every sequence. See
[`docs/database.md`](docs/database.md) for the full schema and backend setup.

---

## Documentation

A complete MkDocs site lives under [`docs/`](docs/index.md):

| Page | Contents |
|------|----------|
| [Installation](docs/installation.md) | Environment setup, backend configuration |
| [Usage](docs/usage.md) | Command templates, targets, troubleshooting |
| [Database](docs/database.md) | Schema, SQLite/PostgreSQL configuration, backups |
| [ITSx](docs/ITSx.md) / [MSA](docs/MSA.md) | Stage-specific details |
| [Personal database](docs/personal-database.md) | Adding your own sequences |
| [Development](docs/dev.md) / [Contributing](docs/contributing.md) | Extending the pipeline |

Build the site locally with `mkdocs serve` (requires the dev environment).

---

## Authors

Developed and maintained by **Auguste Gardette** (UMMISCO, Institut de Recherche
pour le Développement — IRD, France). See [`AUTHORS.md`](AUTHORS.md) for the full
list of contributors.

## Citing CurateMake

If you use CurateMake in your research, please cite it. Citation metadata is
provided in [`CITATION.cff`](CITATION.cff) (GitHub renders a "Cite this
repository" button from it). For a long-term, version-pinned reference, cite the
archived release DOI (to be minted on Zenodo at publication) and the associated
*Methods in Ecology and Evolution* article once available.

## Funding & Acknowledgements

This work was carried out within the **Metaplantcode** project, funded by
**Biodiversa+**, the European Biodiversity Partnership supporting excellent
research on biodiversity with an impact for society and policy.

## License

Released under the **GNU General Public License v3.0 or later** (see
[`LICENSE`](LICENSE)).

Copyright © 2026 Auguste Gardette (IRD, UMMISCO).

External tools bundled or invoked by the pipeline — ITSx, MAFFT, trimAl, HMMER,
RAxML, SATIVA — retain their own licenses; respect their redistribution terms.
