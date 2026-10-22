"""Test cases for the eHMM builder's rescue step.

Run: python workflow/scripts/MSA/ehmm/test_build_ehmm_msa.py

recover() returns statistics the caller reads with .get(). It also has an early exit for
the case where every read was already projected and there is nothing to rescue. That exit
used to return 0, so the caller crashed with "'int' object has no attribute 'get'" --
never on the full database, where some read always fails to project, and always on a small
clean dataset, which is how the mini test run found it.
"""
import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "build_ehmm_msa", Path(__file__).with_name("build_ehmm_msa.py"))
mod = importlib.util.module_from_spec(_SPEC)
sys.modules["build_ehmm_msa"] = mod
_SPEC.loader.exec_module(mod)

passed = failed = 0


def check(name, ok):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    print(f"  {name:<58} {'OK' if ok else 'ECHEC'}")


print("\nSortie anticipee du repechage")

seqs = {"a": "ACGT", "b": "ACGT"}
rows = {"a": "ACGT", "b": "ACGT"}          # tout est deja projete
out = mod.recover(seqs, rows, {}, {}, Path("/tmp"), 4, "", "mmseqs", 0.75, 0.6, 1)

check("rien a repecher ne leve pas d exception", True)
check("le retour expose .get comme le reste du temps", hasattr(out, "get"))
check("et .get repond sans erreur", out.get("rescued_cross_barcode") in (None, 0))
check("il reste faux, donc le drapeau aval ne change pas", not out)

print("\nResolution de mmseqs")

try:
    mod.mmseqs_cmd(None)
    ok = False
except RuntimeError as exc:
    ok = "mmseqs" in str(exc)
except Exception:
    ok = False
check("message explicite quand mmseqs est introuvable", ok)

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
