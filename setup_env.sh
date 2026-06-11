#!/bin/bash
#
# Bootstrap script for CurateMake.
#
# Creates the orchestration Conda environment (curatemake), optionally the
# development environment (curatemake-dev), installs pre-commit hooks when
# available, and prepares a local .env file from the template.
#
# Usage:
#   ./setup_env.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Prefer mamba when available (much faster solver), fall back to conda.
if command -v mamba &> /dev/null; then
    CONDA_BIN="mamba"
elif command -v conda &> /dev/null; then
    CONDA_BIN="conda"
else
    echo "Error: neither mamba nor conda found. Install Miniforge/Mambaforge first." >&2
    exit 1
fi

echo "Using '$CONDA_BIN' to manage environments."

# 1. Orchestration environment (Snakemake)
echo "Creating the orchestration environment (curatemake)..."
"$CONDA_BIN" env create -f environment.yaml || \
    echo "Environment 'curatemake' already exists; skipping creation."

# 2. Optional development environment
read -p "Create the development environment (curatemake-dev)? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    "$CONDA_BIN" env create -f dev-environment.yaml || \
        echo "Environment 'curatemake-dev' already exists; skipping creation."
    echo "Activate it with: conda activate curatemake-dev"
fi

# 3. Pre-commit hooks (only if pre-commit is installed and a config exists)
if command -v pre-commit &> /dev/null && [ -f .pre-commit-config.yaml ]; then
    echo "Installing pre-commit hooks..."
    pre-commit install
fi

# 4. Local .env (project-root, loaded by the Snakefile)
if [ ! -f .env ]; then
    echo "Creating .env from config/.env.example..."
    cp config/.env.example .env
    echo "Edit .env to set your database credentials (PostgreSQL backend only)."
fi

echo
echo "Setup complete."
echo "Activate the pipeline with: conda activate curatemake"
echo "Then run a dry-run with:    snakemake -np"
