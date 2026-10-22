"""Test cases for the taxonomy report's seq_tax reader, on both backends.

Run: python workflow/scripts/test_sativa_taxonomy_report.py

The reader walks seq_tax over yb_hash_code slices, because YugabyteDB has a per-RPC
deadline and a single large scan is refused. SQLite has neither the function nor the
problem, and this reader was the only one in the pipeline that did not make that
distinction: deduplicate_sequences, prepare_msa, reconcile_barcodes and
reannotate_barcodes all already did. So the report rule, which is the LAST step of a run,
failed on every SQLite run with "no such function: yb_hash_code" after all the science
had been computed.

These cases run against a real temporary SQLite database, so they exercise the actual
failure rather than a stand-in, and they pin both halves of the behaviour: the one-pass
read on SQLite, and the sliced read staying the default for PostgreSQL/YugabyteDB.
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
for extra in (_SCRIPTS, _SCRIPTS / "utils"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

_SPEC = importlib.util.spec_from_file_location(
    "gen_report", _SCRIPTS / "generate_sativa_taxonomy_report.py")
rep = importlib.util.module_from_spec(_SPEC)
sys.modules["gen_report"] = rep
_SPEC.loader.exec_module(rep)

from db_models import (Base, Sequence, SeqTax, Source, Taxonomy,  # noqa: E402
                       get_engine, get_session)

passed = failed = 0


def check(name, ok):
    global passed, failed
    passed += bool(ok)
    failed += (not ok)
    print(f"  {name:<58} {'OK' if ok else 'ECHEC'}")


def sqlite_db(tmp, n=25):
    """A minimal real database. seq_tax has foreign keys to sequences and taxonomy, so
    those rows have to exist first."""
    engine = get_engine("sqlite", path=str(Path(tmp) / "t.db"))
    Base.metadata.create_all(engine)
    session = get_session(engine)
    session.add(Source(source_id="src", name="test-source"))
    for i in range(5):
        session.add(Taxonomy(tax_id=f"t{i}", kingdom="Fungi", genus=f"G{i}"))
    ids = [f"s{i}" for i in range(n)] + ["other"]
    for sid in ids:
        session.add(Sequence(seq_id=sid, source_id="src", original_id=sid,
                             sequence="ACGT", length=4, barcode="ITS",
                             barcode_status="Processed"))
    session.flush()
    for i in range(n):
        session.add(SeqTax(seq_id=f"s{i}", tax_id=f"t{i % 5}", backbone="CoL",
                           status="Correct"))
    session.add(SeqTax(seq_id="other", tax_id="t0", backbone="Raw", status="Correct"))
    session.commit()
    session.close()
    return engine


print("\nLecture de seq_tax sur SQLite")

with tempfile.TemporaryDirectory() as tmp:
    engine = sqlite_db(tmp)
    stream = rep._make_stream(engine, 3)

    seen = []
    rep._read_seq_tax(stream, "CoL", seen.extend, engine)
    check("toutes les lignes du backbone sont lues", len(seen) == 25)
    check("les lignes d un autre backbone ne sont pas lues",
          all(sid != "other" for sid, _ in seen))
    check("les identifiants sont ceux attendus",
          {sid for sid, _ in seen} == {f"s{i}" for i in range(25)})

    seen2 = []
    rep._read_seq_tax(stream, "Raw", seen2.extend, engine)
    check("un autre backbone est lu correctement",
          [sid for sid, _ in seen2] == ["other"])

    seen3 = []
    rep._read_seq_tax(stream, "Sativa", seen3.extend, engine)
    check("un backbone absent donne zero ligne sans erreur", seen3 == [])


print("\nLe decoupage par hachage reste le defaut")

# Checked by inspecting the SQL the reader asks for, not by running the doomed path: on
# SQLite the sliced path raises a RuntimeError only after seventeen levels of halving with
# backoff, about a hundred seconds of sleeps for a test that would prove nothing more.

def recording_stream(sql_seen, total=25):
    """Stands in for _make_stream: records the compiled SQL and returns no data rows.

    The sliced path sizes its first pass from a COUNT, so that one query has to answer
    with a number; every other query returns nothing, which is enough to see which SQL
    the reader chose to ask for.
    """
    def stream(stmt, label, max_attempts=None):
        sql_seen.append(str(stmt.compile(compile_kwargs={"literal_binds": True})))
        return [(total,)] if label.startswith("count/") else []
    return stream


class FakeDialect:
    def __init__(self, name):
        self.name = name


class FakeEngine:
    def __init__(self, name):
        self.dialect = FakeDialect(name)


seen_sql = []
rep._read_seq_tax(recording_stream(seen_sql), "CoL", lambda rows: None)
check("sans engine, le decoupage par hachage est conserve",
      any("yb_hash_code" in q for q in seen_sql))

seen_sql = []
rep._read_seq_tax(recording_stream(seen_sql), "CoL", lambda rows: None,
                  FakeEngine("postgresql"))
check("un backend dit postgresql passe par le decoupage",
      any("yb_hash_code" in q for q in seen_sql))

seen_sql = []
rep._read_seq_tax(recording_stream(seen_sql), "CoL", lambda rows: None,
                  FakeEngine("sqlite"))
check("un backend dit sqlite n emet aucun yb_hash_code",
      not any("yb_hash_code" in q for q in seen_sql))
check("et il n emet qu une seule requete", len(seen_sql) == 1)


print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
