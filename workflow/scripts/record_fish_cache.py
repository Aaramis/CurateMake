#!/usr/bin/env python3
"""Re-record the CoL responses that tests/test_fish_testset.py replays.

This is the ONLY part of the fish regression test that touches the network. Run it when
the pinned CoL release changes, check the expected values still hold, and commit the new
cache. The test itself then stays offline and cannot drift.

    python workflow/scripts/record_fish_cache.py [dataset_key]

Default dataset key 316595 is COL26.10 (DOI 10.48580/dgzh4), the release Edi Prifti's
expected values were produced against. Roughly 860 queries, a couple of minutes at 8
requests per second.
"""
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "tests" / "fixtures" / "fish_testset"
DATASET = sys.argv[1] if len(sys.argv) > 1 else "316595"

spec = importlib.util.spec_from_file_location("col", HERE / "col_taxonomic_matcher.py")
col = importlib.util.module_from_spec(spec)
sys.modules["col"] = col
spec.loader.exec_module(col)

client = col.CoLAPIClient(max_requests_per_second=8, dataset=DATASET)
print("CoL release:", json.dumps(client.provenance))

cache = {}
_live = client.search


def recording_search(name, match_type="EXACT", limit=20):
    key = "%s|%s" % (match_type, name)
    if key not in cache:
        cache[key] = _live(name, match_type, limit)
    return cache[key]


client.search = recording_search
matcher = col.TaxonomyMatcher(client)


def run(path, id_field, with_lineage):
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    for i, row in enumerate(rows, 1):
        fields = dict(tax_id=row[id_field], kingdom=row.get("kingdom") or None,
                      family=row.get("family") or None, genus=row.get("genus") or None,
                      species=row.get("species") or None)
        if with_lineage:
            fields.update(class_=row.get("class") or None, order_=row.get("order") or None)
        matcher.match_taxonomy(col.TaxonomyRecord(**fields))
        if i % 200 == 0:
            print("  %s: %d/%d, %d API calls" % (path.name, i, len(rows), client.api_calls))
    print("  %s: %d records" % (path.name, len(rows)))


run(DATA / "raw_records.tsv", "record_id", True)
run(DATA / "edge_cases.tsv", "case_id", False)

out = DATA / "col_cache.json.gz"
with gzip.open(out, "wt", encoding="utf-8") as fh:
    json.dump(cache, fh)
print("wrote %s: %d queries, %d API calls, %.1f KB"
      % (out, len(cache), client.api_calls, out.stat().st_size / 1024))
