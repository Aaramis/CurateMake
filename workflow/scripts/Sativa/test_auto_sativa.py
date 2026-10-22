"""Test cases for the SATIVA task driver, after it moved off the vendored fork.

Run: python workflow/scripts/Sativa/test_auto_sativa.py

Offline: a fake `sativa-epang` on PATH stands in for the real tool, so this checks the
contract the driver owns and not the placement itself. What matters here:

  * the per-task -tmpdir is passed. Without it the tool writes scratch INSIDE its own
    installation, shared by every concurrent task and never cleaned; a single leftover
    jplace reaches 900 MB.
  * that scratch is removed on success and KEPT on failure, where it is the only thing
    left to inspect.
  * the recorded provenance is the real version and the real command. The vendored copy
    reported a constant "SATIVA" while its internal string still claimed v0.9.1, so a
    run could not be traced back to what produced it.
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("auto_sativa", Path(__file__).with_name("auto_sativa.py"))
drv = importlib.util.module_from_spec(_SPEC)
sys.modules["auto_sativa"] = drv
_SPEC.loader.exec_module(drv)

passed = failed = 0


def check(name, ok):
    global passed, failed
    passed += bool(ok)
    failed += (not ok)
    print(f"  {name:<58} {'OK' if ok else 'ECHEC'}")


FAKE = """#!/usr/bin/env bash
if [ "${1:-}" = "-version" ]; then
  echo "sativa-epang 9.9.9 (fake, for tests)"
  exit 0
fi
# record what we were called with, then behave like a successful run
printf '%s\\n' "$*" > called_with.txt
name=""
while [ $# -gt 0 ]; do
  case "$1" in -n) name="$2"; shift 2;; *) shift;; esac
done
: > "${name:-out}.mis"
exit ${FAKE_EXIT:-0}
"""


def fake_bin(tmp):
    d = Path(tmp) / "bin"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "sativa-epang"
    f.write_text(FAKE)
    f.chmod(0o755)
    os.environ["PATH"] = f"{d}{os.pathsep}" + os.environ["PATH"]
    return str(f)


def msa(tmp):
    p = Path(tmp) / "msa.fasta"
    tax = ("kingdom=Fungi|phylum=Basidiomycota|class=Agaricomycetes|order=Russulales"
           "|family=Russulaceae|genus=Russula|species=Russula_alpha")
    p.write_text("".join(f">S{i}|{tax}\nACGTACGTACGT\n" for i in range(6)))
    return p


print("\nPilote SATIVA externalise")

with tempfile.TemporaryDirectory() as tmp:
    binary = fake_bin(tmp)
    out = Path(tmp) / "out"
    out.mkdir()
    fasta, taxf = drv.parse_msa_to_sativa_format(msa(tmp), out, "Russulaceae")
    ok = drv.run_sativa(fasta, taxf, out, binary, mode="ultrafast", confidence=0.4,
                        taxonomy_code="BOT", threads=2, timeout_seconds=60,
                        output_name="ITS_family_Russulaceae", rand_seed=42)
    check("un run reussi est rapporte comme tel", ok is True)
    called = (out / "called_with.txt").read_text()
    check("le scratch par tache est bien passe en -tmpdir", "-tmpdir" in called)
    check("le scratch pointe dans le dossier de la tache", "sativa_tmp" in called)
    check("le scratch est supprime apres un succes", not (out / "sativa_tmp").exists())
    check("la version lue est celle du binaire interroge",
          drv._tool_version(binary) == "sativa-epang 9.9.9 (fake, for tests)")

with tempfile.TemporaryDirectory() as tmp:
    binary = fake_bin(tmp)
    out = Path(tmp) / "out"
    out.mkdir()
    fasta, taxf = drv.parse_msa_to_sativa_format(msa(tmp), out, "Russulaceae")
    os.environ["FAKE_EXIT"] = "1"
    try:
        drv.run_sativa(fasta, taxf, out, binary, mode="ultrafast", confidence=0.4,
                       taxonomy_code="BOT", threads=2, timeout_seconds=60,
                       output_name="ITS_family_Russulaceae")
        raised = False
    except subprocess.CalledProcessError:
        raised = True
    finally:
        os.environ.pop("FAKE_EXIT", None)
    check("un echec remonte bien une exception", raised)
    check("le scratch d un echec est conserve pour inspection", (out / "sativa_tmp").exists())

print("\nConversion du MSA vers les entrees de l outil")

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "out"
    out.mkdir()
    fasta, taxf = drv.parse_msa_to_sativa_format(msa(tmp), out, "Russulaceae")
    lines = taxf.read_text().strip().split("\n")
    check("une ligne de taxonomie par sequence", len(lines) == 6)
    ranks = lines[0].split("\t")[1].split(";")
    check("sept rangs par ligne", len(ranks) == 7)
    check("les tirets bas de l espece redeviennent des espaces", ranks[6] == "Russula alpha")
    check("l identifiant du FASTA perd la taxonomie",
          fasta.read_text().startswith(">S0\n"))

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
