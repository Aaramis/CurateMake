# ITS Pipeline Documentation

This site describes the Snakemake workflow that ingests heterogeneous ITS barcode datasets, harmonizes them inside a relational database, and delivers curated MSAs plus SATIVA validation artifacts.

## Highlights

- **Full reproducibility** – Every rule pins its own Conda environment (`workflow/envs/`), so jobs can be replayed on HPC or cloud schedulers.
- **Traceable data model** – The pipeline writes to either SQLite (local) or PostgreSQL (production). Database schema diagrams live in `docs/database.md`.
- **Modular stages** – Ingestion/formatting, ITSx processing, MSA generation, SATIVA validation, and reporting can be run independently.
- **Observability** – JSON run reports, benchmark files, and HTML dashboards cover each major stage.

## Project Layout

```
workflow/Snakefile        # root workflow
workflow/rules/*.smk      # modular rule sets (database, ITSx, MSA, SATIVA, reports)
workflow/scripts/         # Python helpers for ingestion, QC, reporting
config/config.yaml        # master configuration (paths, resources, parameters)
config/.env               # secrets / DB credentials (not tracked)
data/                     # raw inputs + processed outputs
logs/ + benchmarks/       # runtime diagnostics
docs/                     # this documentation set
```

## Getting Started

1. **Install prerequisites** – Conda/Mamba plus Snakemake (`setup_env.sh` automates this).
2. **Configure** – Edit `config/.env` (DB credentials) and `config/config.yaml` (paths, thresholds).
3. **Dry-run** – `snakemake -np --use-conda` to inspect upcoming work.
4. **Run** – `snakemake --use-conda --cores <N>` or dispatch via a Snakemake profile.

Detailed instructions:

- [Installation](installation.md)
- [Usage](usage.md)
- [Database configuration](database.md)
- [ITSx stage](ITSx.md) and [MSA stage](MSA.md)

## Additional Resources

- [Developer guide](dev.md)
- [Contribution guidelines](contributing.md)
- [Change log](changelog.md)
- [Tooling notes](tools.md)

Happy curation!
