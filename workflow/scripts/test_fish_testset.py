#!/usr/bin/env python3
"""Replay Edi Prifti's fish test set against the CoL matcher, entirely offline.

The test set is 1,439 marine-fish records plus 10 named edge cases, with the taxonomy
expected on COL26.10 (ChecklistBank dataset 316595, DOI 10.48580/dgzh4). Every CoL
response the matcher needs was recorded during one live run and is replayed from
tests/fixtures/fish_testset/col_cache.json.gz, so this touches no network and cannot drift
with a CoL release. Re-record it only when the pinned release changes.

Expected values may list alternatives separated by "|" when either answer is defensible,
and "(or genus only)" means leaving the finer rank empty is also correct.
"""
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "tests" / "fixtures" / "fish_testset"

spec = importlib.util.spec_from_file_location("col", HERE / "col_taxonomic_matcher.py")
col = importlib.util.module_from_spec(spec)
sys.modules["col"] = col
spec.loader.exec_module(col)

RANKS = ("class", "order", "family", "genus", "species")
passed = failed = 0
problems = []


def check(ok, label):
    global passed, failed
    if ok:
        passed += 1
    else:
        failed += 1
        problems.append(label)


class RecordedCoL:
    """Stands in for CoLAPIClient, answering from the recorded responses only."""

    def __init__(self, cache):
        self.cache = cache
        self.api_calls = 0
        self.misses = []

    def search(self, name, match_type="EXACT", limit=20):
        key = "%s|%s" % (match_type, name)
        if key not in self.cache:
            self.misses.append(key)
        return self.cache.get(key, [])


def acceptable(raw):
    """The set of answers the test set calls correct for one rank."""
    raw = (raw or "").strip()
    alts = set()
    if "(or genus only)" in raw:
        alts.add("")
        raw = raw.replace("(or genus only)", "").strip()
    for alt in raw.split("|"):
        alts.add(alt.strip())
    return alts


def observed(result):
    tax = result.col_taxonomy or {}
    got = {r: (tax.get("class_" if r == "class" else "order_" if r == "order" else r) or "")
           for r in RANKS}
    got["status"] = "Check" if result.review else "Correct"
    return got


with gzip.open(DATA / "col_cache.json.gz", "rt", encoding="utf-8") as fh:
    client = RecordedCoL(json.load(fh))
matcher = col.TaxonomyMatcher(client)

# ---------------------------------------------------------------- the 1,439 records
raw = list(csv.DictReader(open(DATA / "raw_records.tsv", encoding="utf-8"), delimiter="\t"))
expected = {r["record_id"]: r
            for r in csv.DictReader(open(DATA / "expected.tsv", encoding="utf-8"),
                                    delimiter="\t")}
check(len(raw) == 1439, "the test set should hold 1439 records, found %d" % len(raw))
check(len(expected) == len(raw), "expected.tsv covers %d of %d records"
      % (len(expected), len(raw)))

flagged = []
for row in raw:
    exp = expected.get(row["record_id"])
    if exp is None:
        continue
    result = matcher.match_taxonomy(col.TaxonomyRecord(
        tax_id=row["record_id"],
        kingdom=row.get("kingdom") or None, class_=row.get("class") or None,
        order_=row.get("order") or None, family=row.get("family") or None,
        genus=row.get("genus") or None, species=row.get("species") or None))
    got = observed(result)
    if result.review:
        flagged.append((row["record_id"], row.get("genus") or row.get("family")))
    wrong = ["%s=%r (expected %s)" % (r, got[r], exp.get("expected_" + r))
             for r in RANKS if got[r] not in acceptable(exp.get("expected_" + r))]
    want_status = (exp.get("expected_status") or "Correct").strip()
    if got["status"] != want_status:
        wrong.append("status=%s (expected %s)" % (got["status"], want_status))
    check(not wrong, "record %s (%s): %s"
          % (row["record_id"], row.get("species") or row.get("genus"), "; ".join(wrong)))

check(not client.misses,
      "the recorded cache is missing %d queries, e.g. %s"
      % (len(client.misses), client.misses[:3]))

# Only the three Moolgarda records are undecided, and they are flagged for a reason:
# COL26.10 offers a provisionally accepted genus and an ambiguous synonym of Crenimugil,
# with no accepted usage of the name anywhere. Anything else flagged is a false positive,
# which is what Edi Prifti reported on Centrophorus, Nectamia and Pentacerotidae.
check(sorted(i for i, _ in flagged) == ["941", "981", "995"],
      "expected only records 941/981/995 flagged, got %s" % (sorted(flagged)[:8],))

# ---------------------------------------------------------------- the 10 edge cases
cases = list(csv.DictReader(open(DATA / "edge_cases.tsv", encoding="utf-8"), delimiter="\t"))
check(len(cases) == 10, "expected 10 edge cases, found %d" % len(cases))
for case in cases:
    result = matcher.match_taxonomy(col.TaxonomyRecord(
        tax_id=case["case_id"], kingdom=case.get("kingdom") or None,
        family=case.get("family") or None, genus=case.get("genus") or None,
        species=case.get("species") or None))
    got = observed(result)
    wrong = ["%s=%r (expected %s)" % (r, got[r], case.get("expected_" + r))
             for r in ("family", "genus", "species")
             if got[r] not in acceptable(case.get("expected_" + r))]
    want_status = (case.get("expected_status") or "Correct").strip()
    if got["status"] != want_status:
        wrong.append("status=%s (expected %s)" % (got["status"], want_status))
    check(not wrong, "edge case %s (%s): %s"
          % (case["case_id"], case["what_it_tests"], "; ".join(wrong)))

total = passed + failed
if failed:
    print("FAILED %d/%d" % (passed, total))
    for p in problems[:25]:
        print("  - %s" % p)
    if len(problems) > 25:
        print("  ... and %d more" % (len(problems) - 25))
    sys.exit(1)
print("PASSED %d/%d" % (passed, total))
