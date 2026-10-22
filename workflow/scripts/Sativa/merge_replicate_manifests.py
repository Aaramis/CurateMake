#!/usr/bin/env python3
"""Merge the native SATIVA task manifest with the replicate manifests.

The native split (seed 0) keeps its task names untouched, so with `sativa.replicates: 1`
this is a byte-for-byte copy of the input manifest and the DAG is exactly the historical
one. Replicate r gets its own task names -- `ITS2__groups__0042` becomes
`ITS2__groups__r1__0042` -- which means every replicate lands in its own SATIVA output
directory and `rules/sativa.smk` needs no new wildcard: one more task in the manifest is
one more `sativa_task` job.

Each task carries `replicate` (0 for the native split). Downstream:
  * `sativa_vote` reads every replicate of a barcode;
  * `sativa_integrate_task` and `sativa_analysis_report` filter on replicate == 0, so the
    database and the reports see one clade per group, never N.
"""
import argparse
import json
import logging
import re
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("merge_manifests")

_SEED_DIR = re.compile(r"^seed(\d+)$")


def replicate_taxon(taxon: str, seed: int) -> str:
    """`ITS2__groups__0042` + seed 1 -> `ITS2__groups__r1__0042`."""
    parts = taxon.split("__")
    if len(parts) < 2:
        return f"{taxon}__r{seed}"
    parts.insert(len(parts) - 1, f"r{seed}")
    return "__".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", required=True)
    ap.add_argument("--replicate-manifest", nargs="*", default=[],
                    help="paths of the form <root>/<barcode>/seed<N>/sativa_manifest.json")
    ap.add_argument("--output", required=True)
    ap.add_argument("--log-file")
    a = ap.parse_args()

    if a.log_file:
        Path(a.log_file).parent.mkdir(parents=True, exist_ok=True)
        logging.getLogger().addHandler(logging.FileHandler(a.log_file))

    manifest = json.load(open(a.native))
    tasks = manifest.get("tasks", [])
    for t in tasks:
        t["replicate"] = 0
    log.info("native split: %d task(s)", len(tasks))

    per_seed = {}
    for path in a.replicate_manifest:
        p = Path(path)
        m = _SEED_DIR.match(p.parent.name)
        if not m:
            raise SystemExit(f"cannot read a seed from {path} "
                             f"(expected .../seed<N>/sativa_manifest.json)")
        seed = int(m.group(1))
        if seed == 0:
            raise SystemExit(f"seed 0 is the native split, not a replicate: {path}")
        rep = json.load(open(p))
        n = 0
        for t in rep.get("tasks", []):
            t = dict(t)
            t["replicate"] = seed
            t["taxon"] = replicate_taxon(t["taxon"], seed)
            tasks.append(t)
            n += 1
        per_seed[seed] = n
        log.info("replicate seed %d: %d task(s) from %s", seed, n, path)

    manifest["tasks"] = tasks
    stats = manifest.setdefault("stats", {})
    stats["replicates"] = 1 + len(per_seed)
    stats["tasks_per_replicate"] = {"0": sum(1 for t in tasks if t["replicate"] == 0),
                                    **{str(k): v for k, v in sorted(per_seed.items())}}
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    json.dump(manifest, open(a.output, "w"), indent=2)
    log.info("wrote %d task(s) over %d replicate(s) -> %s",
             len(tasks), stats["replicates"], a.output)


if __name__ == "__main__":
    main()
