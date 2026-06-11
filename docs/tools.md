# Tooling Reference

Every Snakemake rule specifies its own Conda environment under `workflow/envs/`. Snakemake installs them automatically when `--use-conda` is enabled. Below is a quick reference to the core tools.

## Snakemake

- Version pinned via `setup_env.sh`
- Entry point: `workflow/Snakefile`
- Profiles supported: local cores or any Snakemake profile you configure (SLURM, SGE, etc.)

## Catalogue of Life harmonisation

- Env: `workflow/envs/harmonize_taxo.yaml`
- Rule: `col_taxonomic_matching`
- Script: `workflow/scripts/col_taxonomic_matcher.py`
- Config knobs: `CoL_api.*`
- Output: `data/processed/col_matching_report.json` (writes `backbone = 'CoL'`
  rows into `seq_tax`)

## ITSx

- Env: `workflow/envs/itsx.yaml`
- Rule: `run_itsx_barcode`
- Config knobs: `itsx_cpu`, `itsx_deterministic`, `sequence_validation.*`
- Outputs per barcode: `data/processed/itsx_run_<barcode>.done`, `logs/ITSX/run_itsx_<barcode>.log`

## MAFFT + trimAl

- Env: `workflow/envs/msa.yaml`
- Used in `msa_align_task` and `msa_trim_task`
- Parameters derived from `config["msa"]` (tax columns, trimming overlaps)

## HMMER

- Env: `workflow/envs/msa.yaml`
- Used for:
  - `msa_build_hmm_task`
  - `msa_combine_barcode_hmms`
  - `msa_filter_global_hmm`
  - `msa_classify_batch` (hmmscan)

## SATIVA

- Env: `workflow/envs/sativa.yaml`
- Scripts:
  - `workflow/scripts/Sativa/auto_sativa.py`
  - `workflow/scripts/Sativa/collect_sativa_results.py`
  - `workflow/scripts/Sativa/integrate_sativa.py`
- Config knobs: `sativa.*`, `sativa_integration.*`
- Outputs: per-task logs under `logs/sativa/`, JSON summaries under `data/processed/`

## Reporting stack

- `workflow/scripts/generate_sequence_flow_report.py`
- `workflow/scripts/MSA/generate_report_MSA.py`
- `workflow/scripts/generate_sativa_taxonomy_report.py`
- Dependencies (plotly, pandas, jinja2, etc.) are bundled in `workflow/envs/sativa.yaml`.

## Database utilities

- ORM models (5 tables: `sources`, `sequences`, `taxonomy`, `seq_tax`,
  `sativa_warnings`): `workflow/scripts/utils/db_models.py`.
- Shared backend-agnostic helper: `workflow/scripts/utils/db_manager.py`
  (wraps SQLAlchemy for both SQLite and PostgreSQL).
- Schema bootstrap: `workflow/scripts/create_database.py` (rule `create_database`).
- See [Database configuration](database.md) for the full schema and the
  `seq_tax.backbone` provenance model.

## Benchmarks & logging

- Enable/disable benchmarking via `config["benchmarking"]["enabled"]`.
- Each rule writes logs to `logs/<module>/...`; these directories are declared in the respective Snakemake rule `log:` sections so you can tail them during runs.
