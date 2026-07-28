#!/usr/bin/env python3
"""
UNBOUNDED climb (cross-family gather).

build_family climbs up to the FAMILY (scatter per family, internal climb
species->genus->family). This gather continues the emergent climb
family -> order -> class -> phylum AS LONG AS the gate passes
(ratio <= gate_ratio AND occ >= gate_occ), on the already-enriched family MSAs
(residual + contamination re-homing already applied). Merge via `mafft --merge`
(preserves the sub-alignments, no de novo re-alignment). Few families per order
-> low cost (reduction).

Rules:
  - Climb candidates: level "family" MSAs only (an MSA that stayed below family
    internally will not climb any higher -> passthrough as-is).
  - NEVER a cross-barcode merge (ITS / ITS1 / ITS2 stay separate): grouping is by
    (barcode, rank-value).
  - MSAs emptied by the contamination re-homing (empty file) -> ignored.
  - The stop level stays EMERGENT (decided by the gate, no ceiling).

Output: <out_dir>/_climbed/{msa,hmm}/<id>.{fasta,hmm} for the merged MSAs, and
climbed_manifest.json listing ALL final MSAs (merged + passthrough), with the same
schema as the family manifests -> consumed by collect_bottomup + SATIVA.
"""

import argparse
import json
import logging
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from build_family_msa import (read_fasta, write_fasta, mafft_merge, metrics,
                              hmmbuild, RANKS_FULL)

logger = logging.getLogger("climb_lineage")


def parse_lineage(header: str) -> dict:
    """kingdom..species from a taxonomized header (>seq_id|kingdom=..|..)."""
    d = {r: "" for r in RANKS_FULL}
    for part in header.split("|")[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            if k in d:
                d[k] = v.replace("_", " ")
    return d


def msa_rank_value(aln, rank: str) -> str:
    """Majority value of a rank among the MSA's sequences."""
    vals = [parse_lineage(h)[rank] for h, _ in aln]
    vals = [v for v in vals if v]
    return Counter(vals).most_common(1)[0][0] if vals else ""


def count_taxa(aln) -> int:
    return len({parse_lineage(h)["species"] for h, _ in aln if parse_lineage(h)["species"]})


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Unbounded cross-family lineage climb (gather).")
    p.add_argument("--family-manifests", nargs="*", default=[])
    p.add_argument("--out-dir", required=True, help="Bottom-up out/ dir (merged MSAs go to out/_climbed/).")
    p.add_argument("--output-manifest", required=True)
    p.add_argument("--climb-ranks", nargs="*", default=["order", "class", "phylum"],
                   help="Ranks to try to climb to (empty = bounded at the family).")
    p.add_argument("--gate-ratio", type=float, default=2.5)
    p.add_argument("--gate-occ", type=float, default=40.0)
    p.add_argument("--log-file")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        h = logging.FileHandler(args.log_file, mode="w")
        h.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(h)

    # 1. collect every MSA entry (barcode attached from the family manifest)
    entries = []
    for mpath in args.family_manifests:
        with open(mpath, encoding="utf-8") as fh:
            m = json.load(fh)
        barcode = m.get("barcode", "ITS")
        for e in m.get("msas", []):
            e = dict(e)
            e.setdefault("barcode", barcode)
            entries.append(e)

    def alive(e) -> bool:
        p = e.get("msa_path")
        return bool(p and os.path.exists(p) and os.path.getsize(p) > 0)

    candidates = [e for e in entries if e.get("level") == "family" and alive(e)]
    passthrough = [e for e in entries if e.get("level") != "family" and alive(e)]

    def gate(aln) -> bool:
        mm = metrics(aln)
        return bool(mm and mm["ratio"] <= args.gate_ratio and mm["occ"] >= args.gate_occ)

    climbed_counts: Counter = Counter()
    out = Path(args.out_dir) / "_climbed"
    (out / "msa").mkdir(parents=True, exist_ok=True)
    (out / "hmm").mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="climb_") as tmp:
        # items: list of (entry, aln). "family" entry = original (has an msa_path);
        # merged entry = synthetic (_climbed=True, no file yet).
        items = [(e, read_fasta(e["msa_path"])) for e in candidates]
        uid = 0
        for rank in args.climb_ranks:
            groups = defaultdict(list)
            for e, aln in items:
                groups[(e["barcode"], msa_rank_value(aln, rank))].append((e, aln))
            new_items = []
            for key in sorted(groups):        # (barcode, rank) sorted -> deterministic
                barcode, rv = key
                members = groups[key]
                if not rv or len(members) == 1:
                    new_items.extend(members)
                    continue
                uid += 1
                merged = mafft_merge([aln for _, aln in members], tmp, f"cl{uid}")
                if merged and gate(merged):
                    climbed_counts[rank] += 1
                    new_items.append(({"level": rank, "taxon": rv, "barcode": barcode,
                                       "_climbed": True}, merged))
                else:
                    new_items.extend(members)   # gate refuses -> keep the MSAs as-is
            items = new_items

        # 2. write the final MSAs
        final = list(passthrough)              # non-family: unchanged
        i = 0
        for e, aln in items:
            if not e.get("_climbed"):
                final.append(e)                # family never merged -> keep the original file
                continue
            i += 1
            msa_id = re.sub(r"[^A-Za-z0-9._-]+", "_",
                            f"_climbed__{e['barcode']}__{e['taxon']}__{i:04d}").strip("_")
            msa_path = out / "msa" / f"{msa_id}.fasta"
            hmm_path = out / "hmm" / f"{msa_id}.hmm"
            write_fasta(aln, str(msa_path))    # headers already taxonomized (inherited from source MSAs)
            profiled = hmmbuild(str(msa_path), str(hmm_path), msa_id)
            mm = metrics(aln)
            final.append({"msa_id": msa_id, "level": e["level"], "taxon": e["taxon"],
                          "barcode": e["barcode"], "n": mm["n"], "n_taxa": count_taxa(aln),
                          "ratio": mm["ratio"], "occ": mm["occ"],
                          "gate_pass": bool(mm["ratio"] <= args.gate_ratio and mm["occ"] >= args.gate_occ),
                          "msa_path": str(msa_path),
                          "profile_path": str(hmm_path) if profiled else None})

    stop_levels = Counter(e.get("level") for e in final)
    manifest = {"n_final_msa": len(final),
                "climbed": dict(climbed_counts),
                "stop_levels": dict(stop_levels),
                "climb_ranks": args.climb_ranks,
                "msas": final}
    Path(args.output_manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_manifest, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("Climb: %d final MSAs (climbed merges: %s), stop levels=%s",
                len(final), dict(climbed_counts), dict(stop_levels))


if __name__ == "__main__":
    main()
