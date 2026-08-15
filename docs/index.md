# CurateMake — Documentation

This site documents the Snakemake workflow that ingests heterogeneous ITS
barcode datasets, harmonises their taxonomy, extracts barcode regions with ITSx,
builds curated multiple-sequence alignments (MSAs), and validates taxonomic
labels with SATIVA — storing every record in a traceable relational database.

## Highlights

- **Full reproducibility** — every rule pins its own Conda environment
  (`workflow/envs/`), so jobs can be replayed locally, on HPC, or in the cloud.
- **Traceable data model** — the pipeline writes to either SQLite (local) or
  PostgreSQL (production). The schema is described in
  [Database configuration](database.md).
- **Two-stage taxonomy curation** — labels are reconciled against the Catalogue
  of Life and then validated phylogenetically by SATIVA; both revisions coexist
  in the database through the `seq_tax.backbone` field.
- **Modular stages** — ingestion, harmonisation, ITSx, MSA, SATIVA, and
  reporting can each be run independently.
- **Observability** — JSON run reports, benchmark files, and HTML dashboards
  cover every major stage.

## Pipeline at a glance

| Stage | Rule file | Description |
|-------|-----------|-------------|
| Database bootstrap | `database.smk` | Create the normalised schema |
| Ingestion & formatting | `format_databases.smk` | Load UNITE, BOLD, PLANiTS, CALeDNA, PERSONAL; QC; deduplicate |
| Taxonomy harmonisation | `harmonize_taxo.smk` | Reconcile labels against Catalogue of Life |
| ITSx extraction | `itsx.smk` | Extract ITS/ITS1/ITS2 regions per barcode |
| MSA construction | `bottomup_msa.smk` | Bottom-up per-family build (species → family), residual reclassification, coarse cross-kingdom check |
| SATIVA validation | `sativa.smk` | Per-MSA mislabel detection (leave-one-out via EPA-ng), write back corrections |
| Reporting | `reports.smk`, `benchmarking.smk` | Interactive HTML dashboards |

## Project Layout

```
workflow/Snakefile        # root workflow
workflow/rules/*.smk      # modular rule sets (database, harmonize, ITSx, MSA, SATIVA, reports)
workflow/scripts/         # Python helpers for ingestion, QC, alignment, validation, reporting
workflow/envs/*.yaml      # per-rule Conda environments
config/config.yaml        # master configuration (paths, resources, parameters)
.env                      # secrets / DB credentials (project root, not tracked)
resources/                # raw input databases
results/                   # extracted regions, MSAs, SATIVA, processed outputs
logs/ + benchmarks/       # runtime diagnostics
docs/                     # this documentation set
```

## Getting Started

1. **Install prerequisites** — Conda/Mamba, then create the orchestration
   environment with `conda env create -f environment.yaml` (`setup_env.sh`
   automates this).
2. **Configure** — edit `config/config.yaml` (paths, thresholds, backend). For
   the PostgreSQL backend, copy `config/.env.example` to `.env` and fill in
   credentials.
3. **Dry-run** — `snakemake -np` to inspect the upcoming DAG.
4. **Run** — `snakemake --use-conda --cores <N>` or dispatch via a Snakemake
   profile.

Detailed instructions:

- [Installation](installation.md)
- [Usage](usage.md)
- [Database configuration](database.md)
- [ITSx stage](ITSx.md) and [MSA stage](MSA.md)
- [Adding your own data (PERSONAL database)](personal-database.md)

## Additional Resources

- [Developer guide](dev.md)
- [Contribution guidelines](contributing.md)
- [Change log](changelog.md)
- [Tooling notes](tools.md)

## Citing & Funding

CurateMake is developed by **Auguste Gardette** (UMMISCO, IRD). If you use it,
please cite it using the metadata in `CITATION.cff`. This work was carried out
within the **Metaplantcode** project, funded by **Biodiversa+**, the European
Biodiversity Partnership supporting excellent research on biodiversity with an
impact for society and policy. Released under the GPLv3.
