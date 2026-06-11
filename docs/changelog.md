# Changelog

This document lists notable changes to the project following the [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format.

## [Unreleased]
### Changed
- Refreshed the README and documentation set to match the current workflow:
  documented the Catalogue of Life harmonisation stage, the PERSONAL database,
  the two-stage (`CoL` + `SATIVA`) taxonomy curation model, and the five-table
  database schema.
- Renamed the project to **CurateMake** (repository, Conda environments
  `curatemake`/`curatemake-dev`, and all documentation references) and corrected
  the `.env` location (project root, loaded by the `Snakefile`).
### Fixed
- Reworked `environment.yaml` into a minimal orchestration environment and made
  `dev-environment.yaml` standalone (removed the invalid `file:` include and the
  broken `pip: -e .` install).
- Updated `setup_env.sh` to create the orchestration environment and bootstrap a
  project-root `.env`.

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
