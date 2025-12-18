# Developer Notes

## Repository structure

```
workflow/
  Snakefile            # root DAG
  rules/               # modular rule files
  envs/                # per-rule Conda specs
  scripts/             # Python entry points used by rules
    MSА/, ITSx/, Sativa/, utils/
config/
  config.yaml          # main settings (paths, thresholds)
  .env                 # secrets/interpolation values
data/
  raw/, processed/, ITSx/, MSА/, SATIVA/
logs/, benchmarks/, reports/
docs/
  *.md (MkDocs source)
```

Rules are grouped by function (database, ITSx, MSA, SATIVA, benchmarking, reports). Keep new logic inside `workflow/scripts/` and call it via `script:` or `shell:` blocks to preserve reproducibility.

## Adding a new rule

1. Create `workflow/rules/<feature>.smk` or extend an existing file.
2. Implement the heavy lifting in `workflow/scripts/<feature>/...py`.
3. Declare inputs/outputs/logs/benchmarks/conda env in the Snakemake rule.
4. Reference configuration values via `config[...]` to avoid hard-coding paths.
5. Add targets to `rule all` if the output should be part of the default DAG.

Example snippet:

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

## Configuration best practices

- Introduce new parameters under `config/config.yaml` with sensible defaults.
- Access secrets via environment variables (referenced as `{{ env.VAR }}`).
- Document new keys in the README or relevant doc page.

## Testing & validation

- Use `snakemake -np` to verify the DAG after modifying rules.
- For scripts, add lightweight unit tests under `tests/` or embed assertions in the script until a proper test harness exists.
- During development, run specific targets (e.g., a single alignment) to iterate quickly.

## Logging

- Always specify a `log:` file so troubleshooting is straightforward.
- Prefer structured logging inside Python scripts (`logging` module with clear INFO/WARNING/ERROR messages).

## Development workflow

1. Branch from `CurateMake`.
2. Make focused changes (rules + scripts + docs).
3. Run targeted Snakemake dry-runs or unit tests.
4. Commit using Conventional Commits.
5. Open a PR with reproduction instructions and relevant logs/reports.

## Contributions

- Branch from `CurateMake`, keep commits focused, and open PRs with:
  - Summary of rule/script changes
  - Instructions to reproduce the run (targets + config overrides)
  - Relevant log excerpts or report screenshots if behavior changed

Keeping these conventions consistent makes the pipeline easier to maintain across the different stages (database ingestion → ITSx → MSA → SATIVA → reporting).
