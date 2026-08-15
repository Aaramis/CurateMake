# Installation

Follow these steps to prepare a workstation or HPC login node for the workflow.

## 1. Prerequisites

- A Linux/macOS shell with Bash, GNU coreutils, and `git`
- Conda or Mamba (Miniforge/Mambaforge recommended)
- PostgreSQL client libraries **only** if you intend to use the PostgreSQL
  backend (the default SQLite backend needs nothing extra)

```bash
# Install Miniforge (provides conda + mamba)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

## 2. Clone the repository

```bash
git clone https://github.com/Aaramis/CurateMake.git
cd CurateMake
```

## 3. Create the orchestration environment

The workflow uses a two-tier environment model:

- **One orchestration environment** (`curatemake`) that you activate — it
  contains Snakemake and its launch-time dependencies.
- **Per-rule environments** under `workflow/envs/`, created automatically by
  Snakemake the first time each rule runs with `--use-conda`.

```bash
conda env create -f environment.yaml      # or: mamba env create -f environment.yaml
conda activate curatemake
```

Alternatively, run the interactive bootstrap script, which also offers to create
the development environment and a local `.env`:

```bash
chmod +x setup_env.sh
./setup_env.sh
```

No manual installation of ITSx, MAFFT, trimAl, HMMER, or SATIVA is required —
they are provisioned per rule.

## 4. Configure paths and credentials

1. Edit `config/config.yaml` (or create an override and pass it with
   `--configfile`) to point at your raw data, set resource limits, and choose
   the database backend (`sqlite` vs `postgresql`).
2. For the **PostgreSQL backend only**, copy the credentials template to the
   project root and edit it:
   ```bash
   cp config/.env.example .env
   # edit .env: POSTGRES_HOST, POSTGRES_PORT, POSTGRES_USER, POSTGRES_PASSWORD, POSTGRES_DB
   ```
   The `Snakefile` loads `.env` from the project root and resolves
   `{{ env.VAR }}` placeholders found in `config/config.yaml`. Missing variables
   raise a clear error at start-up.

## 5. Initialise the database schema (optional)

- **SQLite** — created automatically under `results/db/` on first run.
- **PostgreSQL** — make sure the database, user, and schema already exist; the
  pipeline creates the tables on first run. You can pre-create them with:
  ```bash
  snakemake --use-conda --cores 1 create_database
  ```

## 6. Validate the installation

```bash
# Confirm Snakemake is on PATH
snakemake --version

# Dry-run the default DAG (no jobs executed)
snakemake -np
```

If the dry-run completes without errors, launch the full workflow:

```bash
snakemake --use-conda --cores <N>
```

## Troubleshooting

- **`import dotenv` / parse-time errors** — make sure the `curatemake`
  environment is activated; it provides `python-dotenv`.
- **A leaking virtualenv shadows `python`** — deactivate any active `.venv`
  before launching; the `Snakefile` also strips `VIRTUAL_ENV`/`PYTHONPATH`
  defensively.
- **Database authentication errors** — double-check `.env` and the `database`
  block in `config/config.yaml`.
- **Insufficient disk space** — the pipeline writes large FASTA/MSA/SATIVA
  artefacts under `results/`. Use symlinks or update the paths in the config.
