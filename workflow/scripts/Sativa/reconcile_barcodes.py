#!/usr/bin/env python3
"""Propagate each SATIVA verdict across the sibling barcodes of the SAME read.

ITSx extracts ITS1 and ITS2 out of ITS, then SATIVA runs on each barcode independently. Those rows
describe ONE organism, so a correction found on ITS1 must also apply to ITS and ITS2 — otherwise the
curated database contradicts itself. Measured on the s8/r15 simulation (12 358 anchors carrying a
SATIVA verdict on several barcodes): 63.5 % already agree, 2.2 % are compatible, and 34.3 % genuinely
CONFLICT. In 90.8 % of those conflicts one barcode already holds the right answer — it is simply
trapped in the wrong row.

Arbitration is chosen by `rank_by` (config sativa_integration.reconcile_barcodes.rank_by):
  * "votes" (default since the consensus vote exists): number of independent splits that agreed,
    then lineage depth, then SATIVA confidence. A sibling carrying no voted correction scores 0
    votes, so a modification the vote retained on ONE barcode is propagated to the others.
  * "confidence": highest SATIVA confidence, lineage depth as tie-break — the pre-vote ordering,
    used automatically when no vote data exists. Measured end-to-end against the actual pipeline
    output on that point (43.9 correction / 15.0 cross-kingdom / 12.7 cross-order /
    18.3 false-positive): reconciliation gives 46.7 / 18.7 / 19.7 / 13.9 — better on all four.
    Those figures come from a SINGLE split, where a lone flag was ~63% artefact; the "votes"
    ordering revisits the same trade-off at ~15% and must be re-measured against them.
It also makes the result DETERMINISTIC: today the anchor keeps whichever sibling row the SQL scan
happened to read last, and two orderings of the same database scored 15.0 and 26.6 cross-kingdom.
Ties (an identical rank key on a DIFFERENT taxonomy) are never guessed: both rows are left in place.

Scale notes (this runs once over the whole database, 5M+ sequences):
  * agreement is decided on tax_id, never on the seven rank strings — taxonomies are de-duplicated,
    so equal tax_id means equal lineage, and the scan stays narrow;
  * depth is fetched once per DISTINCT tax_id (thousands), not per row (millions);
  * the scan walks yb_hash_code slices (same reason as deduplicate_sequences / prepare_msa: a
    server-side cursor dies through the pooler and ORDER BY has no index to lean on);
  * only rows that actually disagree with their winner are written.
"""
from __future__ import annotations

import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from db_manager import pg_hint  # noqa: E402  (one constructor for every hinted query)

logger = logging.getLogger(__name__)

RANKS = ["kingdom", "phylum", "class_", "order_", "family", "genus", "species"]
YB_HASH_SPACE = 65536          # fixed by YugabyteDB
HASH_SLICES = 64               # few wide slices beat many narrow ones (cf. deduplicate_sequences)
_NULLS = {None, "", "Na", "NA", "None", "null", "unknown", "Unknown"}

# Drive seq_tax from its primary key and batch the nested loop against the VALUES list;
# without it every VALUES-joined write spills to temp files and dies on temp_file_limit
# (1 GB). See pg_hint() in db_manager for the measurements behind every hint here.
_PK_HINT = pg_hint([("u", "seq_tax_pkey")], batched=("u", "v"))


def _table(db_manager, name: str) -> str:
    getter = getattr(db_manager, f"{name}_table", None)
    if getter:
        return getter
    schema = getattr(db_manager, "schema", None)
    return f"{schema}.{name}" if schema else name


def _iter_sativa_rows(db_manager) -> Iterable[Tuple[str, str, str, Optional[float]]]:
    """Yield (anchor, seq_id, tax_id, confidence) for every SATIVA verdict.

    anchor = the read every barcode row descends from (parent_seq_id when ITSx split it).
    """
    from sqlalchemy import text

    seqs, st = _table(db_manager, "sequences"), _table(db_manager, "seq_tax")
    select = (f"SELECT COALESCE(s.parent_seq_id, s.seq_id) AS anchor, s.seq_id, st.tax_id, "
              f"st.confidence FROM {st} st JOIN {seqs} s ON st.seq_id = s.seq_id "
              f"WHERE st.backbone = 'Sativa'")

    if db_manager.db_type != "postgresql":
        def fetch_all():
            with db_manager.get_session() as session:
                return list(session.execute(text(select)))
        yield from db_manager._execute_with_retry("reconcile_scan_all", fetch_all)
        return

    def read_range(lo: int, hi: int, depth: int = 0):
        # drive the join from seq_tax's pkey, as in bulk_get_col_taxonomy / prepare_msa
        sql = (f"{pg_hint([('st', 'seq_tax_pkey'), ('s', 'sequences_pkey')], batched=('st', 's'))}"
               f"{select} AND yb_hash_code(st.seq_id) >= :lo AND yb_hash_code(st.seq_id) < :hi")

        def fetch():
            with db_manager.get_session() as session:
                return list(session.execute(text(sql), {"lo": lo, "hi": hi}))
        try:
            return db_manager._execute_with_retry(f"reconcile_scan_{lo}_{hi}", fetch)
        except Exception:
            if hi - lo <= 1 or depth > 12:
                raise
            mid = (lo + hi) // 2                       # a refused range is halved, never dropped
            return read_range(lo, mid, depth + 1) + read_range(mid, hi, depth + 1)

    step = YB_HASH_SPACE // HASH_SLICES
    for lo in range(0, YB_HASH_SPACE, step):
        yield from read_range(lo, min(lo + step, YB_HASH_SPACE))


def _load_votes(db_manager) -> Dict[str, int]:
    """seq_id -> how many independent splits proposed this correction.

    Read from sativa_warnings, ONE table, no join -- the cluster's standing rule. A sequence
    with no row here was never a voted correction (it is either untouched or below the vote
    threshold, hence VALIDATED), and gets 0, which is what makes any voted modification
    outrank a merely-confirmed sibling.

    Returns an empty dict when sativa.replicates == 1: the column is then NULL everywhere and
    the ranking falls back to confidence exactly as before.
    """
    from sqlalchemy import text

    warn = _table(db_manager, "sativa_warnings")
    sql = f"SELECT seq_id, sativa_votes FROM {warn} WHERE sativa_votes IS NOT NULL"

    def fetch():
        with db_manager.get_session() as session:
            return list(session.execute(text(sql)))

    try:
        rows = db_manager._execute_with_retry("reconcile_load_votes", fetch)
    except Exception as exc:                       # missing column on a pre-vote database
        logger.warning("no vote information available (%s); ranking by confidence", exc)
        return {}

    votes: Dict[str, int] = {}
    for seq_id, n in rows:
        key = str(seq_id)
        if n is not None and n > votes.get(key, 0):
            votes[key] = int(n)                    # one row per sequence, but never assume it
    logger.info("Barcode reconciliation: %d voted correction(s) loaded from %s",
                len(votes), warn)
    return votes


def _tax_depth(db_manager, tax_ids: Iterable[str]) -> Dict[str, int]:
    """tax_id -> number of non-null ranks. One small query per chunk of DISTINCT tax_ids."""
    from sqlalchemy import text

    tax = _table(db_manager, "taxonomy")
    ids = list(tax_ids)
    chunk = 900 if db_manager.db_type == "sqlite" else 2000
    out: Dict[str, int] = {}
    cols = ", ".join(RANKS)
    for i in range(0, len(ids), chunk):
        part = ids[i:i + chunk]
        if db_manager.db_type == "postgresql":
            clause, params = "tax_id = ANY(:ids)", {"ids": part}
        else:
            clause = "tax_id IN (" + ", ".join(f":i{j}" for j in range(len(part))) + ")"
            params = {f"i{j}": v for j, v in enumerate(part)}
        sql = f"SELECT tax_id, {cols} FROM {tax} WHERE {clause}"

        def fetch(sql=sql, params=params):
            with db_manager.get_session() as session:
                return list(session.execute(text(sql), params))

        for row in db_manager._execute_with_retry(f"reconcile_depth_{i}", fetch):
            out[row[0]] = sum(1 for v in row[1:] if v not in _NULLS)
    return out


def _apply(db_manager, updates: List[Tuple[str, str, str]],
           deletes: List[Tuple[str, str]], batch_size: int = 1000) -> int:
    """Align each sequence on the winning tax_id: retype one Sativa row, drop the extras.

    tax_id is part of seq_tax's primary key (seq_id, tax_id, backbone), so a sequence that
    ALREADY carries the winning verdict on another row cannot have its losing row updated to
    it -- that is a duplicate key, and it takes the whole batch down. Those rows are deleted
    instead; the winner is already in place. Verified against this schema: the naive UPDATE
    raises UniqueViolation on exactly that case.

    One statement per batch rather than a list of parameter dicts: SQLAlchemy turns the latter
    into cursor.executemany(), i.e. one round trip PER ROW. Measured here, 5000 rows:
    303 rows/s that way against 11700-17000 rows/s with UPDATE ... FROM (VALUES ...).
    SQLite keeps the executemany form -- it is local, so round trips cost nothing, and neither
    UPDATE ... FROM nor DELETE ... USING is portable.
    """
    from sqlalchemy import text

    st = _table(db_manager, "seq_tax")
    done = 0

    if db_manager.db_type != "postgresql":
        upd = text(f"UPDATE {st} SET tax_id = :new WHERE seq_id = :seq_id "
                   f"AND tax_id = :old AND backbone = 'Sativa'")
        dele = text(f"DELETE FROM {st} WHERE seq_id = :seq_id AND tax_id = :old "
                    f"AND backbone = 'Sativa'")

        def run_sqlite():
            with db_manager.get_session() as session:
                if deletes:
                    session.execute(dele, [{"seq_id": s, "old": t} for s, t in deletes])
                if updates:
                    session.execute(upd, [{"seq_id": s, "old": o, "new": n} for s, o, n in updates])
                session.commit()
            return len(updates) + len(deletes)

        return db_manager._execute_with_retry("reconcile_apply_sqlite", run_sqlite)

    # Deletes first: they are what frees the winning tax_id for the updates below.
    for i in range(0, len(deletes), batch_size):
        chunk = deletes[i:i + batch_size]
        ph = ",".join(f"(:s{j},:t{j})" for j in range(len(chunk)))
        params = {}
        for j, (s, t) in enumerate(chunk):
            params[f"s{j}"], params[f"t{j}"] = s, t
        # The hint is not cosmetic: without it the planner materialises the join and
        # blows the 1 GB temp_file_limit -- measured, the bare DELETE already fails at
        # 25 rows. With it, 1000 rows delete in 0.04s (22.5k rows/s). Same hint the
        # scan above uses. Keep batch_size <= 1000: 2000 breaks the limit even hinted.
        sql = text(f"{_PK_HINT}DELETE FROM {st} u USING (VALUES {ph}) AS v(seq_id, tax_id) "
                   f"WHERE u.seq_id = v.seq_id AND u.tax_id = v.tax_id AND u.backbone = 'Sativa'")

        def run_del(sql=sql, params=params, n=len(chunk)):
            with db_manager.get_session() as session:
                session.execute(sql, params)
                session.commit()
            return n

        done += db_manager._execute_with_retry(f"reconcile_delete_{i}", run_del)

    for i in range(0, len(updates), batch_size):
        chunk = updates[i:i + batch_size]
        ph = ",".join(f"(:s{j},:o{j},:n{j})" for j in range(len(chunk)))
        params = {}
        for j, (s, o, n) in enumerate(chunk):
            params[f"s{j}"], params[f"o{j}"], params[f"n{j}"] = s, o, n
        # Same story, worse: the bare UPDATE exceeds temp_file_limit at EVERY batch size
        # tested (100/500/1000) and once closed the backend outright. Hinted: 1165 rows/s.
        sql = text(f"{_PK_HINT}UPDATE {st} u SET tax_id = v.newtax "
                   f"FROM (VALUES {ph}) AS v(seq_id, oldtax, newtax) "
                   f"WHERE u.seq_id = v.seq_id AND u.tax_id = v.oldtax AND u.backbone = 'Sativa'")

        def run_upd(sql=sql, params=params, n=len(chunk)):
            with db_manager.get_session() as session:
                session.execute(sql, params)
                session.commit()
            return n

        done += db_manager._execute_with_retry(f"reconcile_update_{i}", run_upd)

    return done


def reconcile(db_manager, apply: bool = True, rank_by: str = "votes") -> Dict[str, int]:
    """Propagate the best-supported SATIVA verdict to every sibling barcode of the same read.

    rank_by:
      "votes"      (default) — number of agreeing splits, then lineage depth, then confidence.
      "confidence" — the pre-vote ordering (confidence, then depth), kept as a fallback.
    """
    stats = {"anchors": 0, "multi_barcode_anchors": 0, "conflicting_anchors": 0,
             "rows_updated": 0, "ties_flagged": 0, "rank_by": rank_by}

    per: Dict[str, List[Tuple[str, str, float]]] = defaultdict(list)   # anchor -> [(seq_id, tax, conf)]
    for anchor, seq_id, tax_id, conf in _iter_sativa_rows(db_manager):
        if tax_id:
            per[anchor].append((seq_id, tax_id, -1.0 if conf is None else float(conf)))
    stats["anchors"] = len(per)

    votes = _load_votes(db_manager) if rank_by == "votes" else {}
    if rank_by == "votes" and not votes:
        logger.info("Barcode reconciliation: no vote data, falling back to confidence ranking.")
        rank_by = "confidence"
        stats["rank_by"] = "confidence (no vote data)"

    multi = {a: rows for a, rows in per.items() if len(rows) > 1 and len({r[1] for r in rows}) > 1}
    stats["multi_barcode_anchors"] = sum(1 for rows in per.values() if len(rows) > 1)
    stats["conflicting_anchors"] = len(multi)
    if not multi:
        logger.info("Barcode reconciliation: nothing to do (no conflicting sibling verdicts).")
        return stats

    depth = _tax_depth(db_manager, {r[1] for rows in multi.values() for r in rows})

    def rank(row):
        """Sort key. Two orderings, selected by `rank_by`.

        "votes" (default since the consensus vote exists): number of agreeing splits, then
        lineage depth, then SATIVA confidence. A sequence with no voted correction scores 0
        votes, so ANY modification that survived the vote is propagated to its siblings —
        which is the point: a flag confirmed by 7 of 10 independent neighbourhoods is not the
        single lucky draw the older ordering was calibrated against.

        "confidence": the pre-vote ordering, kept because two alternatives were built and
        measured on the s8/r15 point and both were REJECTED at the time — do not reintroduce
        either as a default without new evidence:
          * "prefer the sibling that proposes a real change" (SATIVA scores a confirmation 1.000 and
            a correction ~0.82-0.94, so confirmations do win ties): buys +0.8 pt of cross-kingdom but
            takes false positives from 13.9% to 26.1%, because forcing the bolder verdict also
            propagates SATIVA's mistakes on clean sequences to every sibling row;
          * "the most complete lineage wins": picks correctly 64.9% of the time against 92.1% for
            confidence — depth is not reliability, a species-level call is often the weaker one.
        Both measurements were made on a SINGLE split, where a lone flag carried ~63% artefact.
        The "votes" ordering revisits exactly that trade-off with the artefact rate cut to ~15%,
        and it must be re-measured against those numbers on the first full voted run.
        A fixed barcode priority (ITS1 > ITS > ITS2) scored well too, but that ordering was an
        artefact of which MSAs had crashed in that run — see sativa.cfg / the GTRGAMMA fix.
        """
        seq_id, tax_id, conf = row
        if rank_by == "votes":
            return (votes.get(seq_id, 0), depth.get(tax_id, 0), conf)
        return (conf, depth.get(tax_id, 0))

    updates: List[Tuple[str, str, str]] = []   # (seq_id, old_tax, winning_tax)
    deletes: List[Tuple[str, str]] = []        # (seq_id, tax_to_drop)
    ties: List[str] = []
    for anchor, rows in multi.items():
        best = max(rows, key=rank)
        best_key = rank(best)
        # an equally-ranked but DIFFERENT verdict is a genuine tie: flag it, never guess
        if any(r[1] != best[1] and rank(r) == best_key for r in rows):
            ties.append(anchor)
            continue
        # Each sequence must end up with exactly ONE Sativa row, on the winning tax_id.
        # tax_id is part of the primary key, so a sequence already holding the winner cannot
        # have its other rows updated onto it -- they are dropped instead. Likewise a sequence
        # with several losing rows gets one retyped and the rest dropped, because updating them
        # all would collide with each other.
        by_seq: Dict[str, List[str]] = defaultdict(list)
        for seq_id, tax_id, _ in rows:
            by_seq[seq_id].append(tax_id)
        for seq_id, taxes in by_seq.items():
            losing = [t for t in taxes if t != best[1]]
            if not losing:
                continue
            if len(losing) < len(taxes):          # the winner is already on one of its rows
                deletes.extend((seq_id, t) for t in losing)
            else:
                updates.append((seq_id, losing[0], best[1]))
                deletes.extend((seq_id, t) for t in losing[1:])

    stats["ties_flagged"] = len(ties)
    stats["rows_retyped"] = len(updates)
    stats["rows_dropped"] = len(deletes)
    logger.info("Barcode reconciliation: %d anchors, %d multi-barcode, %d conflicting -> "
                "%d row(s) retyped, %d redundant row(s) dropped, %d tie(s) left untouched.",
                stats["anchors"], stats["multi_barcode_anchors"], stats["conflicting_anchors"],
                len(updates), len(deletes), len(ties))
    if apply and (updates or deletes):
        stats["rows_updated"] = _apply(db_manager, updates, deletes)
        logger.info("Barcode reconciliation: %d sibling row(s) aligned.", stats["rows_updated"])
    return stats
