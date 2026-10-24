# Changelog

This document lists notable changes to the project following the [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format.

## [Unreleased]
### Added
- eHMM alignment: one `mafft-linsi` backbone per barcode, UPP-style bucketing of the
  remaining sequences, `nhmmer`/`hmmalign` projection onto the backbone's match columns,
  and an mmseqs rescue pass for sequences no profile accepts.
- SATIVA replicate voting: each group is validated several times with different seeds, and
  a mislabel is kept only on a consensus vote. The per-replicate manifests are merged into
  one.
- A test suite runner (`workflow/scripts/run_tests.sh`) and a GitHub Actions job that runs
  it on every push. The suite is offline: no test reaches ChecklistBank or a database.
- `test_auto_sativa_real.py`, which drives the packaged `sativa-epang` on a synthetic
  three-clade alignment and checks it reports a deliberately mislabelled sequence, and only
  that one.
- `test_fish_testset.py`, which replays a contributed 1,439-record marine-fish test set
  against the CoL matcher from recorded CoL responses. The data itself is not redistributed
  here; see `tests/fixtures/README.md`.
- `docs/alignment-strategy.md`: the alignment question for a hypervariable barcode, the
  per-clade proposal, the measurements taken, and why the proposal was rejected.
- Per-taxon CoL provenance (`col_id`) and config-schema validation.
### Changed
- SATIVA runs the `sativa-epang` package from bioconda instead of a copy vendored in this
  repository, in a conda environment of its own because bioconda builds it against a
  different Python. `sativa.sativa_binary` selects the executable.
- The CoL matcher was rewritten. It breaks same-kingdom homonyms on the lineage the source
  declares, distinguishes a living homonym from a synonym, and flags a pick as `Check` only
  when CoL itself has no settled opinion on the name.
- Snakemake layout: the old `data/` tree is split into `results/` for outputs and
  `resources/` for inputs.
- The bottom-up MSA path is now the only one.
- `workflow/envs/msa.yaml` declares `mmseqs2` instead of expecting a hand-made environment.
- Refreshed the README and documentation set to match the current workflow:
  documented the Catalogue of Life harmonisation stage, the PERSONAL database,
  the two-stage (`CoL` + `SATIVA`) taxonomy curation model, and the five-table
  database schema.
- Renamed the project to **CurateMake** (repository, Conda environments
  `curatemake`/`curatemake-dev`, and all documentation references) and corrected
  the `.env` location (project root, loaded by the `Snakefile`).
### Fixed
- `sativa.sativa_binary` given as a path now works. The tool calls its own siblings,
  `epa-ng` and RAxML through it, by bare name, so the binary's directory is added to the
  child process's PATH.
- CoL matcher: six spurious ambiguity flags are gone. Four were subgenus homonyms competing
  with the genus that was queried, two were an accepted name beside a synonym of the same
  spelling. Conversely, a name with no accepted usage anywhere is now flagged instead of
  being presented as settled.
- CoL matcher: a rank CoL leaves unresolved is no longer back-filled from the source, a
  misspelled genus is corrected using the family, and a binomial that contradicts the
  resolved genus is never written.
- The taxonomy report reads `seq_tax` in one pass on SQLite, keeping the hash-sliced scans
  only for YugabyteDB, where the per-RPC deadline requires them.
- The eHMM rescue step returns a mapping even when there is nothing to rescue, and tolerates
  conda being absent from PATH entirely.
- SATIVA integration drops stem placements.
- Reworked `environment.yaml` into a minimal orchestration environment and made
  `dev-environment.yaml` standalone (removed the invalid `file:` include and the
  broken `pip: -e .` install).
- Updated `setup_env.sh` to create the orchestration environment and bootstrap a
  project-root `.env`.
### Removed
- The SATIVA copy vendored at `workflow/scripts/Sativa/sativa/`, and with it `ete3` and
  `epa-ng` from `workflow/envs/sativa.yaml`: both belonged to that copy.
- The top-down MSA strategy and the `sativa_min_taxa` knob.

## [0.1.0] - 2025-08-18
### Added
- Complete MkDocs site for the project
- Restructured directory layout for maintainability
- Git configuration with comprehensive `.gitignore`
- Conda environment files (`curatemake`, `curatemake-dev`)
- Bootstrap script (`setup_env.sh`)
- Git usage documentation and best practices
- Initial GitHub repository configuration

## Changelog format

Each release should specify sections when relevant:

- **Added** for new features
- **Changed** for updates to existing behavior
- **Deprecated** for soon-to-be removed features
- **Removed** for actual removals
- **Fixed** for bug fixes
- **Security** for vulnerabilities addressed

## Versioning principles

The project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html):

- **MAJOR**: incompatible API changes
- **MINOR**: backward-compatible new functionality
- **PATCH**: backward-compatible bug fixes
