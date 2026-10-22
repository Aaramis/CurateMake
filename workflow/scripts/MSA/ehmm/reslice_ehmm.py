#!/usr/bin/env python3
"""Re-cut an ALREADY PROJECTED barcode into a different set of SATIVA groups.

Why this exists
---------------
The group boundaries handed to SATIVA are an artefact of the traversal order, not of biology:
pack_units_interleaved walks clusters round-robin across kingdoms and cuts at max_group. Measured
on the 2026-08-25 run, the neighbourhood drives the verdict hard -- of the 1200 backbone sequences
scaffolded into ALL 122 ITS groups (same reads, same 940-column frame, only the neighbours differ),
1101 were flagged at least once, at a median of 28 groups out of 122, and 984 of them (89.4%)
proposed a DIFFERENT lineage depending on the group, median 6 distinct proposals, up to 36.

Running SATIVA over several independent splits and voting per rank is therefore the way to tell a
real correction from an artefact of who shared the task. This script makes the extra replicates
cheap: the expensive half of eHMM -- backbone, nhmmer, hmmalign, ~15 h for ITS -- depends only on
the barcode, NOT on the split. Its output survives as the per-bucket a_*.afa files, so a replicate
only has to re-read them and re-cut, which takes minutes.

Usage
-----
    # one replicate, explicit output directory
    python reslice_ehmm.py --work-dir results/MSA/ehmm/ITS2 --full-csv results/full.csv \
        --barcode ITS2 --out-dir results/MSA/replicates/ITS2_seed1 --split-seed 1 \
        --max-group 20000 --min-group 4 --scaffold

    # several replicates in ONE pass: the projection is read once, then re-cut per seed,
    # each into <out-dir>/seed<N>. Reading the a_*.afa files is the whole cost here
    # (7.8 GB for ITS2, 96 GB for ITS), so this is what keeps N replicates affordable.
    python reslice_ehmm.py --work-dir results/MSA/ehmm/ITS2 --full-csv results/full.csv \
        --barcode ITS2 --out-dir results/MSA/replicates/ITS2 --split-seed 1 2 3 \
        --max-group 20000 --min-group 4 --scaffold

Emits the same sativa_manifest.json contract as ehmm_msa.smk, so the rest of the pipeline is
unchanged. Nothing is written to the original output directory.
"""
import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_ehmm_msa import load_barcode, read_fasta, slice_groups  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("reslice")


def load_projection(work: Path, width: int = 0):
    """Rebuild {seq_id: aligned row} from the per-bucket alignments hmmalign already wrote.

    Same filter as hmmalign_bucket: keep gaps and match (upper-case) columns, drop the
    lower-case insertions, then require the row to be exactly `width` columns wide.
    """
    rows, widths, skipped = {}, {}, 0
    files = sorted(work.glob("a_*.afa"))
    if not files:
        raise SystemExit(f"no a_*.afa under {work} -- has this barcode been projected?")
    for path in files:
        for head, seq in read_fasta(path):
            row = "".join(c for c in seq if c == "-" or c.isupper())
            widths[len(row)] = widths.get(len(row), 0) + 1
            if width and len(row) != width:
                skipped += 1
                continue
            rows[head.split()[0]] = row
    if not width:
        width = max(widths, key=widths.get)          # the frame is the dominant row length
        rows = {s: r for s, r in rows.items() if len(r) == width}
    log.info("projection reloaded: %d rows of %d columns from %d bucket alignment(s), %d skipped",
             len(rows), width, len(files), skipped)
    return rows, width


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", required=True, help="eHMM work dir holding a_*.afa")
    ap.add_argument("--full-csv", required=True)
    ap.add_argument("--barcode", required=True)
    ap.add_argument("--out-dir", required=True,
                    help="MUST differ from the original output dir. With several --split-seed "
                         "values this is the ROOT and each seed lands in <out-dir>/seed<N>.")
    ap.add_argument("--split-seed", type=int, nargs="+", required=True,
                    help="one or more non-zero seeds; each yields an independent split")
    ap.add_argument("--max-group", type=int, default=20000)
    ap.add_argument("--min-group", type=int, default=4)
    ap.add_argument("--width", type=int, default=0, help="0 = infer from the alignments")
    ap.add_argument("--scaffold", action="store_true")
    a = ap.parse_args()

    seeds = list(dict.fromkeys(a.split_seed))          # de-duplicate, keep order
    if any(s == 0 for s in seeds):
        raise SystemExit("--split-seed 0 would reproduce the original split; pick another")
    work, root = Path(a.work_dir), Path(a.out_dir)
    if work.resolve() == root.resolve():
        raise SystemExit("--out-dir must differ from --work-dir")

    # The expensive half -- reading every a_*.afa back -- is done ONCE and shared by all
    # seeds. Re-cutting is then pure CPU on an in-memory dict.
    rows, width = load_projection(work, a.width)
    tax = {sid: t for sid, _seq, t in load_barcode(a.full_csv, a.barcode)}
    log.info("taxonomy loaded for %d %s sequences", len(tax), a.barcode)

    bb_ids = set()
    bb = work / "backbone_raw.fasta"
    if bb.exists():
        bb_ids = {h.split()[0].split("|")[0] for h, _ in read_fasta(bb)}
    log.info("backbone: %d ids (scaffold=%s)", len(bb_ids), bool(a.scaffold))

    for seed in seeds:
        out = root if len(seeds) == 1 else root / f"seed{seed}"
        out.mkdir(parents=True, exist_ok=True)
        slice_groups(rows, tax, bb_ids, out, a.barcode, a.max_group, a.min_group, a.scaffold,
                     len(rows), set(), None, seed)
        man = json.load(open(out / "sativa_manifest.json"))
        log.info("seed %d -> %d task(s) in %s", seed, len(man["tasks"]), out)


if __name__ == "__main__":
    main()
