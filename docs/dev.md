# Developer Notes

## Repository structure

```
workflow/
  Snakefile            # root DAG + default targets
  rules/               # modular rule files
  envs/                # per-rule Conda specs
  scripts/             # Python entry points used by rules
    MSA/  ITSx/  Sativa/  utils/
config/
  config.yaml          # main settings (paths, thresholds)
  .env.example         # credentials template (copy to project-root .env)
data/
  raw/  processed/  ITSx/  MSA/  SATIVA/
logs/  benchmarks/  reports/
docs/                  # MkDocs source
environment.yaml       # orchestration env (Snakemake)
dev-environment.yaml   # development toolchain
```

Rules are grouped by function: `database`, `format_databases`, `harmonize_taxo`
(Catalogue of Life), `itsx`, `msa`, `sativa`, `benchmarking`, and `reports`.
Keep heavy logic inside `workflow/scripts/` and invoke it from rules via
`script:` or `shell:` blocks to preserve reproducibility.

## Adding a new rule

1. Create `workflow/rules/<feature>.smk` (and `include:` it from the `Snakefile`)
   or extend an existing file.
2. Implement the heavy lifting in `workflow/scripts/<feature>/...py`.
3. Declare `input`, `output`, `log`, `benchmark`, and `conda` in the rule.
4. Reference configuration values via `config[...]` rather than hard-coding paths.
5. Add the new output to `rule all` in the `Snakefile` if it should be part of
   the default DAG.

Example:

```python
rule example_task:
    input:
        "data/processed/itsx_process.json"
    output:
        "data/processed/example.json"
    params:
        option=config["example"]["option"]
    log:
        "logs/example/task.log"
    conda:
        "../envs/example.yaml"
    script:
        "../scripts/example/run_example.py"
```

## Checkpoints

The MSA and SATIVA stages use Snakemake **checkpoints** (`msa_task_manifest`,
`msa_classification_manifest`, `sativa_task_manifest`) because the set of
downstream jobs is only known after a manifest is generated at runtime. When you
add per-clade or per-task rules, derive their inputs from the relevant manifest
via the checkpoint helper functions in the rule file.

## Configuration best practices

- Add new parameters to `config/config.yaml` with sensible defaults.
- Reference secrets via environment variables (`{{ env.VAR }}`), kept in the
  project-root `.env`.
- Document new keys in the README or the relevant `docs/` page.

## Testing & validation

- Use `snakemake -np` to verify the DAG after modifying rules.
- For scripts, add tests under `tests/` (pytest) or embed assertions until a
  dedicated harness exists.
- During development, build a single target (e.g. one alignment) to iterate
  quickly.

## Logging

- Always declare a `log:` file so failures are easy to trace.
- Use the `logging` module inside Python scripts with clear INFO/WARNING/ERROR
  messages; log text is kept in English to suit monitoring systems.

## Development workflow

1. Branch from `CurateMake`.
2. Make focused changes (rules + scripts + docs together).
3. Run targeted dry-runs or unit tests.
4. Commit using [Conventional Commits](https://www.conventionalcommits.org/).
5. Open a PR with reproduction instructions and relevant logs/reports.

Keeping these conventions consistent makes the pipeline easier to maintain across
its stages: ingestion → harmonisation → ITSx → MSA → SATIVA → reporting.
