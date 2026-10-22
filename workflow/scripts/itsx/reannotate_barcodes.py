#!/usr/bin/env python3
"""Re-annotate the ITS1/ITS2 barcode of sequences ITSx could not verify, by similarity.

ITSx delimits ITS1/ITS2 from their flanking regions (SSU, 5.8S, LSU). A source database that has
ALREADY been ITSx-trimmed has no flanks left, so ITSx cannot re-detect anything and its declared
barcode is taken on trust. Measured on the full run: of 289 583 sequences declared ITS1, ITSx could
check only 197. In the unverifiable remainder the declared barcode is often wrong — 85.8 % of a
10 000-sequence sample of BOLD-style `<id>.ITS1` entries match ITSx-certified ITS2 references at a
median identity of 100 %, i.e. they ARE ITS2 (~157 000 sequences of that pool). UNITE-style
`<id>_ITS1` entries, and every ITS2 pool, came back 100 % correct — the defect is source-specific.

A wrong barcode is not cosmetic: the MSA is clustered per barcode, so those sequences are aligned
against the wrong region, which corrupts the reference frame SATIVA then places on.

Method: ITSx's own extractions are the ground truth (regions cut out of full ITS sequences, where the
flanks were present). They give a large, free, verified reference set for each region; every
unverifiable sequence is searched against both and re-labelled when the answer is unambiguous.
Validation: on sequences whose region IS ITSx-certified, this classifies 100 % correctly (n=3 489).

Only `sequences.barcode` is rewritten, and only for rows whose region ITSx did not establish. ITSx's
own children keep the region ITSx assigned them.
"""
import json
import logging
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

REGIONS = ("ITS1", "ITS2")


def _read_fasta(path: Path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def build_reference(itsx_dir: Path, out_fa: Path, per_region: int) -> Dict[str, int]:
    """One FASTA holding ITSx-certified ITS1 and ITS2, tagged REF1_/REF2_ so a hit names its region."""
    counts = {r: 0 for r in REGIONS}
    if not itsx_dir.is_dir():
        logger.warning("ITSx output directory %s is missing — no reference can be built.", itsx_dir)
        return counts
    with open(out_fa, "w") as out:
        for tag, region in (("REF1", "ITS1"), ("REF2", "ITS2")):
            for barcode_dir in sorted(p for p in itsx_dir.iterdir() if p.is_dir()
                                      and not p.name.startswith(".")):
                src = barcode_dir / f"{barcode_dir.name}.{region}.fasta"
                if not src.exists():
                    continue
                for _, seq in _read_fasta(src):
                    if counts[region] >= per_region:
                        break
                    if seq:
                        out.write(f">{tag}_{counts[region]}\n{seq}\n")
                        counts[region] += 1
    logger.info("Reference built from ITSx extractions: %s", counts)
    return counts


def classify(queries_fa: Path, refs_fa: Path, workdir: Path, threads: int,
             min_id: float, min_cov: float, env: str) -> Dict[str, Tuple[str, float]]:
    """seq_id -> (region, bitscore) for the best hit. mmseqs, one pass over all queries."""
    hits = workdir / "hits.m8"
    tmp = workdir / "mmtmp"
    tmp.mkdir(exist_ok=True)
    cmd = ["mmseqs", "easy-search", str(queries_fa), str(refs_fa), str(hits), str(tmp),
           "--search-type", "3", "--min-seq-id", str(min_id), "-c", str(min_cov),
           "--max-seqs", "5", "--threads", str(threads), "-v", "1"]
    if env:
        cmd = ["conda", "run", "-n", env] + cmd
    subprocess.run(cmd, check=True)
    best: Dict[str, Tuple[str, float]] = {}
    with open(hits) as fh:
        for line in fh:
            f = line.split("\t")
            score = float(f[11]) if len(f) > 11 else float(f[2])
            region = "ITS1" if f[1].startswith("REF1") else "ITS2"
            if f[0] not in best or score > best[f[0]][1]:
                best[f[0]] = (region, score)
    return best


def main(snakemake):
    cfg = snakemake.config
    params = cfg.get("msa", {}).get("reannotate_barcodes", {})
    itsx_dir = Path(params.get("itsx_dir", "results/ITSx"))
    per_region = int(params.get("reference_per_region", 50000))
    min_id = float(params.get("min_seq_id", 0.70))
    min_cov = float(params.get("min_coverage", 0.70))
    threads = int(getattr(snakemake, "threads", 8) or 8)
    env = params.get("mmseqs_env", "mmseqs")
    report_path = Path(snakemake.output.report)

    sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from db_manager import get_database_manager
    from sqlalchemy import text

    stats = {"candidates": 0, "classified": 0, "reassigned": 0, "confirmed": 0,
             "by_change": defaultdict(int)}

    if not params.get("enabled", True):
        logger.info("Barcode re-annotation disabled by config.")
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps({"enabled": False}, indent=2))
        return

    db = get_database_manager(cfg)
    seqs_table = getattr(db, "sequences_table", "sequences")

    with tempfile.TemporaryDirectory(prefix="reannot_") as tmpdir:
        work = Path(tmpdir)
        refs = work / "refs.fasta"
        counts = build_reference(itsx_dir, refs, per_region)
        if min(counts.values()) < 100:
            logger.warning("Not enough ITSx-certified references (%s) — skipping re-annotation.", counts)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(json.dumps({"skipped": "insufficient_reference", **counts}, indent=2))
            return

        # candidates = declared ITS1/ITS2 rows that ITSx did NOT produce (ITSx children carry the
        # region ITSx itself established, and are therefore already trustworthy)
        # Duplicates are excluded: they are already out of every downstream step, and on the
        # full run they are 318k of the 1.15M candidates -- a quarter of the mmseqs work for
        # rows nothing will ever read. no_detection rows ARE kept: ITSx finding no region is
        # precisely what makes their declared barcode unverified, so they are the target.
        where = ("barcode IN ('ITS1','ITS2') AND parent_seq_id IS NULL "
                 "AND barcode_status != 'Duplicate'")
        base = f"SELECT seq_id, barcode, sequence FROM {seqs_table} WHERE {where}"

        def fetch_candidates():
            """Walk the table's own hash space on YugabyteDB.

            The candidate set is ~1.15 M rows AND carries the sequences, so a single SELECT is far
            past the 60 s per-RPC deadline, and a server-side cursor does not survive the pooler
            (same reason prepare_msa / deduplicate_sequences slice by yb_hash_code). Each slice is an
            independent, retryable query; a slice the backend refuses is halved rather than lost.
            SQLite has neither the limit nor the function, so it reads in one go.
            """
            if db.db_type != "postgresql":
                def one():
                    with db.get_session() as session:
                        return list(session.execute(text(base)))
                yield db._execute_with_retry("reannotate_fetch_candidates", one)
                return

            def read_range(lo: int, hi: int, depth: int = 0):
                sql = (f"{base} AND yb_hash_code(seq_id) >= :lo AND yb_hash_code(seq_id) < :hi")

                def one():
                    with db.get_session() as session:
                        return list(session.execute(text(sql), {"lo": lo, "hi": hi}))
                try:
                    return db._execute_with_retry(f"reannotate_fetch_{lo}_{hi}", one)
                except Exception:
                    if hi - lo <= 1 or depth > 12:
                        raise
                    mid = (lo + hi) // 2                  # a refused slice is halved, never dropped
                    return read_range(lo, mid, depth + 1) + read_range(mid, hi, depth + 1)

            span, slices = 65536, 64
            step = span // slices
            for lo in range(0, span, step):
                yield read_range(lo, min(lo + step, span))

        # stream the slices straight into the query FASTA: holding 1.15 M sequences in memory first
        # is pointless, only the seq_id -> declared-barcode map has to survive
        qfa = work / "queries.fasta"
        declared: Dict[str, Tuple[str, str]] = {}
        n = 0
        with open(qfa, "w") as out:
            for chunk in fetch_candidates():
                for seq_id, barcode, sequence in chunk:
                    if not sequence:
                        continue
                    key = f"Q{n}"
                    declared[key] = (seq_id, barcode)
                    out.write(f">{key}\n{sequence}\n")
                    n += 1
        stats["candidates"] = n
        logger.info("Candidates whose region ITSx could not establish: %d", n)
        if not n:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(json.dumps({k: v for k, v in stats.items() if k != "by_change"},
                                              indent=2))
            return

        best = classify(qfa, refs, work, threads, min_id, min_cov, env)
        stats["classified"] = len(best)

        updates: List[Dict[str, str]] = []
        for key, (region, _score) in best.items():
            seq_id, barcode = declared[key]
            if region == barcode:
                stats["confirmed"] += 1
            else:
                updates.append({"seq_id": seq_id, "barcode": region})
                stats["by_change"][f"{barcode}->{region}"] += 1

        logger.info("Classified %d/%d candidates: %d confirmed, %d to re-annotate (%s)",
                    len(best), n, stats["confirmed"], len(updates),
                    dict(stats["by_change"]))

        if updates and params.get("apply", True):
            batch_size = 1000
            for i in range(0, len(updates), batch_size):
                batch = updates[i:i + batch_size]

                if db.db_type == "postgresql":
                    # One statement per batch. Handing SQLAlchemy a list of parameter dicts
                    # makes it fall back to cursor.executemany(), i.e. one round trip per row:
                    # measured on this cluster, 303 rows/s that way against ~17000 rows/s here.
                    ph = ",".join(f"(:s{j},:b{j})" for j in range(len(batch)))
                    bind = {}
                    for j, row in enumerate(batch):
                        bind[f"s{j}"], bind[f"b{j}"] = row["seq_id"], row["barcode"]
                    stmt = text(f"UPDATE {seqs_table} u SET barcode = v.barcode "
                                f"FROM (VALUES {ph}) AS v(seq_id, barcode) "
                                f"WHERE u.seq_id = v.seq_id")
                    payload = bind
                else:   # SQLite: local, round trips are free and UPDATE ... FROM is not portable
                    stmt = text(f"UPDATE {seqs_table} SET barcode = :barcode WHERE seq_id = :seq_id")
                    payload = batch

                def run(stmt=stmt, payload=payload, n=len(batch)):
                    with db.get_session() as session:
                        session.execute(stmt, payload)
                        session.commit()
                    return n

                stats["reassigned"] += db._execute_with_retry(f"reannotate_apply_{i}", run)
            logger.info("Re-annotated %d sequence(s).", stats["reassigned"])

    stats["by_change"] = dict(stats["by_change"])
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(stats, indent=2))


if __name__ == "__main__":
    if "snakemake" in globals():
        # NEVER block the pipeline. This step only improves barcode labels; a failure here must not
        # kill a multi-day production run, so any error is logged and recorded in the report while
        # the workflow carries on with the labels as they are.
        try:
            main(snakemake)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Barcode re-annotation failed — continuing with unchanged labels: %s", exc)
            out = Path(snakemake.output.report)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps({"status": "failed", "error": str(exc)}, indent=2))
    else:
        logger.error("This script must be run through Snakemake")
        sys.exit(1)
