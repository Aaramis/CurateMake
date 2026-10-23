# Contributed test fixtures

Some regression tests replay a test set somebody contributed rather than one we wrote.
Those files are not redistributed here until their author agrees, so this directory is
mostly empty in a fresh clone and the tests that need it skip with a message. Everything
needed to run them is in the repository; only the data is missing.

## The fish test set

`workflow/scripts/test_fish_testset.py` replays 1,439 marine-fish records plus 10 named
edge cases against the Catalogue of Life matcher, contributed by Edi Prifti with the
taxonomy expected on COL26.10 (ChecklistBank dataset 316595, DOI 10.48580/dgzh4).

To run it, put these files in `tests/fixtures/fish_testset/`:

| File | What it holds |
| --- | --- |
| `raw_records.tsv` | the input records, keyed by `record_id` |
| `expected.tsv` | one `expected_<rank>` column per rank, plus `expected_status` |
| `edge_cases.tsv` | the 10 named cases, keyed by `case_id` |
| `col_cache.json.gz` | the recorded CoL responses, gzipped JSON |

Then:

    python workflow/scripts/test_fish_testset.py     # just this test
    bash workflow/scripts/run_tests.sh               # the whole suite

The path is in `.gitignore`, so adding the files cannot commit them by accident.

### Conventions in the expected values

An expected value may list alternatives separated by `|` when either answer is
defensible, as in `Moolgarda|Crenimugil`. The marker `(or genus only)` means leaving the
finer rank empty is also correct. `expected_status` is `Correct` or `Check`, where
`Check` means the matcher must flag the record rather than present the pick as settled.

### The recorded responses

`col_cache.json.gz` maps `"<match type>|<queried name>"` to the raw response CoL gave,
so the test touches no network and cannot drift with a release. The test fails if a
query it needs is missing from the cache, rather than silently treating it as no match.

Re-record it when the pinned release changes, then check the expected values still hold:

    python workflow/scripts/record_fish_cache.py [dataset_key]

### If redistribution is agreed

Drop the `tests/fixtures/fish_testset/` line from `.gitignore` and commit the four
files. The test needs no change: it runs as soon as they are present.
