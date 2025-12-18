# Pipeline Installation

Follow these steps to prepare a workstation or HPC login node for the Snakemake workflow.

## 1. Prerequisites

- Linux/macOS shell with Bash, GNU coreutils, and `git`
- Conda or Mamba (Miniforge/Mambaforge recommended)
- PostgreSQL client tools if you plan to use the remote DB backend

```bash
# Install Mambaforge (recommended)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

## 2. Clone the repository

```bash
git clone https://github.com/<org>/ITS-pipeline.git
cd ITS-pipeline/6.CurateMake
```

## 3. Bootstrap Snakemake + rule environments

The repo ships with a helper script that pins Snakemake and creates the shared tooling env:

```bash
chmod +x setup_env.sh
./setup_env.sh
# This installs snakemake via Conda/Mamba and leaves envs under ~/.conda/envs/
```

Rule-specific Conda environments (ITSx, MAFFT, SATIVA, etc.) are automatically created by Snakemake the first time each rule runs. No manual activation is required.

## 4. Configure secrets and paths

1. Copy the `.env` template and define database credentials, API keys, etc.
   ```bash
   cp config/.env.example config/.env
   # edit config/.env
   ```
2. Edit `config/config.yaml` (or create a `config/local.yaml` override) to point to your raw data, tweak resource limits, and choose the database backend (`sqlite` vs `postgresql`). Environment variables referenced via `{{ env.VAR }}` will be resolved from `.env`.

## 5. Optional: initialize the database schema

- **SQLite** is created automatically in `data/db`.
- **PostgreSQL**: ensure the database, user, and schema already exist. The pipeline will create tables on first run, but you can pre-create them by running:
  ```bash
  snakemake create_database --use-conda --cores 1
  ```

## 6. Validate the installation

```bash
# Confirm snakemake is on PATH
snakemake --version

# Dry-run the default DAG (no jobs executed)
snakemake --use-conda -np
```

If the dry-run completes without errors, you're ready to launch the full workflow:

```bash
snakemake --use-conda --cores <N>
```

## Troubleshooting

- **Missing system libraries**: rerun `setup_env.sh` after updating Conda/Mamba.
- **Database authentication errors**: double-check `config/.env` and the `database` block in `config/config.yaml`.
- **Insufficient disk space**: the pipeline writes large FASTA/MSA/SATIVA artifacts under `data/`. Use symlinks or update the paths in the config if needed.
