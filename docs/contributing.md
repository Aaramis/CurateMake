# Contribution Guide

Thank you for your interest in the project! We welcome bug reports, feature ideas, and code/documentation improvements.

## Table of contents

- [How to contribute](#how-to-contribute)
- [Bug reports](#bug-reports)
- [Feature proposals](#feature-proposals)
- [Submitting changes](#submitting-changes)
- [Code style](#code-style)
- [Tests](#tests)
- [Code review](#code-review)
- [Questions & support](#questions-support)

## How to contribute

1. **Report a problem**
   - Check existing issues first
   - Open a new issue with a clear, descriptive title
   - Use the provided issue template

2. **Propose an improvement**
   - Describe the feature or enhancement in detail
   - Explain why it matters and include example use cases

3. **Submit a change**
   - Follow the [development workflow](dev.md#development-workflow)
   - Ensure tests (or Snakemake dry-runs) pass
   - Update documentation when behavior changes

## Bug reports

To report a bug:

1. Use the “Bug Report” template
2. Provide a clear, descriptive title
3. Describe the problem in detail
4. List the steps required to reproduce it
5. Specify the expected behavior
6. Attach logs or screenshots when helpful
7. Include environment details (OS, Python version, Snakemake version, etc.)

## Feature proposals

To propose a new feature:

1. Check whether a similar request already exists
2. Use the “Feature Request” template
3. Describe the desired feature clearly
4. Outline the use cases and expected benefits
5. Suggest an implementation approach if you have one

## Submitting changes

### Requirements

- Python 3.10 (matches the pipeline environments)
- Git
- A GitHub account
- Basic familiarity with Snakemake

### Workflow

1. **Fork** the repository
2. Clone your fork locally:
   ```bash
   git clone https://github.com/<your-user>/CurateMake.git
   cd CurateMake
   ```
3. Create a feature branch:
   ```bash
   git checkout -b feature/<feature-name>
   ```
4. Create the development environment (optional):
   ```bash
   conda env create -f dev-environment.yaml
   conda activate curatemake-dev
   ```
5. Implement your changes
6. Check formatting/linting:
   ```bash
   black .
   flake8
   ```
7. Run tests (or targeted Snakemake runs):
   ```bash
   pytest
   ```
8. Commit with a clear message:
   ```bash
   git commit -m "feat: add new SATIVA metrics (#123)"
   ```
9. Push to your fork:
   ```bash
   git push origin feature/<feature-name>
   ```
10. Open a Pull Request targeting the `CurateMake` branch (or whichever branch you’re asked to use)

## Code style

### Style

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Format with [Black](https://black.readthedocs.io/)
- Lint with [Flake8](https://flake8.pycqa.org/)
- Document functions with docstrings (Google-style or similar)

### Commits

- Use clear, descriptive commit messages
- Follow [Conventional Commits](https://www.conventionalcommits.org/)
- Reference issue numbers when relevant

Example:
```
feat: add SATIVA manifest diff view

Adds a new dashboard comparing SATIVA runs and closes #123.
```

Suggested prefixes:
- `feat`: new feature
- `fix`: bug fix
- `docs`: documentation updates
- `style`: formatting-only changes
- `refactor`: code restructuring
- `perf`: performance improvements
- `test`: add/modify tests
- `chore`: build/CI/config updates

## Tests

### Run tests

```bash
# The whole suite
bash workflow/scripts/run_tests.sh

# With a particular interpreter
bash workflow/scripts/run_tests.sh python3.11

# A single file
python workflow/scripts/test_col_taxonomic_matcher.py
```

The suite is offline by design. The Catalogue of Life client is either faked or replayed
from recorded responses, so no test reaches ChecklistBank, none needs a database, and no
result can drift with a CoL release. CI runs the same script on a plain Python with a
handful of pip packages and no conda.

### Two tests skip instead of failing

They print `SKIPPED` and the runner treats that as neither a pass nor a failure, so a
skipped test is visible in the log rather than silently counted as success:

- `test_auto_sativa_real.py` drives the real `sativa-epang`, which comes from bioconda and
  is not a pip dependency. It runs as soon as that tool is on PATH. The arguments we build
  around it and the per-task scratch stay covered by `test_auto_sativa.py`, which fakes the
  executable.
- `test_fish_testset.py` replays a test set contributed by someone else, which is not
  redistributed in this repository. See `tests/fixtures/README.md` for the files it wants
  and where to put them.

### Writing tests

The test files are script-style, not pytest: each one prints `PASSED n/m` or `FAILED n/m`
as its last line and exits non-zero on failure. pytest cannot collect them, because they
call `sys.exit()` at module level. Keep that shape rather than introducing a framework.

- Put the file next to the code it covers, named `test_<module>.py`.
- Add it to the `TESTS` array in `workflow/scripts/run_tests.sh`, or it will never run.
- Keep third-party imports in step with the install step of `.github/workflows/tests.yml`.
  A test that cannot import its module prints no result line and fails the run, which is
  how a missing dependency gets caught.
- Cover the failure path, not only the happy one.
- Say in a comment which regression a check exists for. Most of these checks were written
  after something broke, and that context is what stops someone deleting them later.

## Code review

Every contribution goes through a code-review process:

1. Open a Pull Request
2. Ensure tests/dry-runs pass and CI is green
3. A maintainer reviews the code and may request changes
4. Once approved, the PR is merged

Acceptance checklist:
- Code is tested or accompanied by a reproducible Snakemake command
- Documentation and changelog entries are updated when needed
- Style and conventions match the rest of the repo
- The change benefits other users or is clearly scoped

## Questions & support

If you have questions or run into issues:

1. Consult the [documentation](index.md)
2. Open a GitHub issue
3. Contact the maintainers (see AUTHORS.md) for sensitive topics

## Acknowledgements

Thanks to everyone who contributes time and ideas — your name is added to
CONTRIBUTORS.md once your first PR lands.
