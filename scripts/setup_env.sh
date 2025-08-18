#!/bin/bash

# Setup script for bio-pipeline environment
# Usage: source setup_env.sh

# Exit on error
set -e

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

# # Create and activate the base environment
# echo "Creating base environment..."
# conda env create -f envs/environment.yaml
# conda activate its-pipeline

# echo "Base environment created and activated."

# Create development environment if requested
read -p "Would you like to create a development environment? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Creating development environment..."
    conda env create -f envs/dev-environment.yaml --name its-pipeline-dev
    echo "Development environment created. Activate with: conda activate its-pipeline-dev"
fi

# Set up pre-commit hooks
echo "Setting up pre-commit hooks..."
pre-commit install

# Create .env file if it doesn't exist
if [ ! -f .env ]; then
    echo "Creating .env file from template..."
    cp .env.example .env
    echo "Please edit the .env file with your configuration."
fi

echo "Setup complete!"
echo "Activate the environment with: conda activate its-pipeline"
echo "or for development: conda activate its-pipeline-dev"
