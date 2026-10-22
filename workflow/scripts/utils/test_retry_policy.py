"""Test cases for DatabaseManager's retry classification.

Run: python workflow/scripts/utils/test_retry_policy.py

YugabyteDB reports several unrelated conditions as XX000, and they do not all deserve the
same treatment -- but the deciding factor is NOT the message. It is what the CALLER does
with the failure:
  * "Duplicate request" / leader changes -> genuinely transient; retrying IS the fix
    (PLANITS lost 75 905 sequences to one of these on 2026-08-26).
  * "Rpc timeout" -> a load overrun ("passed: 15.018s, timeout: 15.000s"), NOT a property
    of the statement: the same 1 000-row INSERT succeeded millions of times in the same run.
    It is retryable by default. Treating it as fatal killed process_barcode_itsx on ITS at
    15:15 on 2026-08-28.
  * ...unless the caller resplits or falls back, in which case a fast failure is what it
    wants: retrying turned a 15 s failure into 4 min 48 inside a recursion that doubles the
    chunk count at each level. Those three sites pass retry_deadline=False.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db_manager import DatabaseManager  # noqa: E402
from sqlalchemy.exc import (OperationalError, IntegrityError, InternalError,  # noqa: E402
                            ProgrammingError, DataError)


class Orig(Exception):
    def __init__(self, code, msg=""):
        super().__init__(msg)
        self.pgcode = code
        self._msg = msg

    def __str__(self):
        return self._msg


def mk(cls, code, msg=""):
    return cls("stmt", {}, Orig(code, msg))


RPC = "Rpc timeout, passed: 15.097s, timeout: 15.000s, now: 5011441.299s, last: 0.000s"
DUP = ("[Already present (yb/consensus/retryable_requests.cc:393): "
       "Duplicate request 5696 from client f5c4b5c4-d0b9-46bc-8956-cc102ade8290 (min running 5696)]")

# (nom, exception, reessaye par defaut, budget par defaut, reessaye si retry_deadline=False)
CASES = [
    ("40001 serialisation",       mk(OperationalError, "40001"),          True,  100, True),
    ("40P01 interblocage",        mk(OperationalError, "40P01"),          True,  100, True),
    ("XX000 Duplicate request",   mk(InternalError, "XX000", DUP),        True,    8, True),
    ("XX000 leader change",       mk(InternalError, "XX000", "Leader changed"), True, 8, True),
    ("XX000 Rpc timeout",         mk(InternalError, "XX000", RPC),        True,    8, False),
    ("XX000 Timed out waiting",   mk(InternalError, "XX000", "Timed out waiting kResponseSent"),
                                                                          True,    8, False),
    ("XX000 deadline exceeded",   mk(InternalError, "XX000", "Deadline Exceeded"), True, 8, False),
    ("IntegrityError",            mk(IntegrityError, "23505"),            False,   0, False),
    ("53400 limite config",       mk(OperationalError, "53400"),          False,   0, False),
    ("54000 limite programme",    mk(OperationalError, "54000"),          False,   0, False),
    ("42703 colonne inconnue",    mk(ProgrammingError, "42703"),          False,   0, False),
    ("22001 donnee trop longue",  mk(DataError, "22001"),                 False,   0, False),
]

passed = failed = 0
print(f"  {'cas':<28} {'reessaye':>9} {'attendu':>8} {'budget':>7} {'attendu':>8} {'sans deadline':>14}")
for name, exc, exp_retry, exp_budget, exp_nodl in CASES:
    r = DatabaseManager._is_retryable(exc)
    b = 0
    if r:
        b = min(100, DatabaseManager._INTERNAL_MAX_ATTEMPTS) \
            if DatabaseManager._is_internal_transient(exc) else 100
    nodl = DatabaseManager._is_retryable(exc, retry_deadline=False)
    ok = (r == exp_retry) and (b == exp_budget) and (nodl == exp_nodl)
    passed += ok
    failed += not ok
    print(f"  {name:<28} {str(r):>9} {str(exp_retry):>8} {b:>7} {exp_budget:>8} {str(nodl):>14}"
          f"  {'OK' if ok else 'ECHEC'}")

# 1. l'insertion (aucun repli) doit survivre a un pic de charge
rpc = mk(InternalError, "XX000", RPC)
survives = DatabaseManager._is_retryable(rpc)
passed += survives; failed += not survives
print(f"\n  bulk_insert_sequences survit a un Rpc timeout      : {'OK' if survives else 'ECHEC'}")

# 2. les appelants qui redecoupent gardent leur echec rapide
fast = not DatabaseManager._is_retryable(rpc, retry_deadline=False)
passed += fast; failed += not fast
print(f"  les 3 appelants adaptatifs echouent immediatement  : {'OK' if fast else 'ECHEC'}")

# 3. une erreur deterministe n'est jamais reessayee, quel que soit le drapeau
det = mk(IntegrityError, "23505")
never = not DatabaseManager._is_retryable(det) and not DatabaseManager._is_retryable(det, retry_deadline=False)
passed += never; failed += not never
print(f"  une IntegrityError n'est jamais reessayee          : {'OK' if never else 'ECHEC'}")

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
