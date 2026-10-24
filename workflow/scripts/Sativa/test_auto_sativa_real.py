#!/usr/bin/env python3
"""Drive the REAL sativa-epang, not a stand-in, on a synthetic alignment.

test_auto_sativa.py checks the plumbing with a fake executable on PATH: the arguments we
build, the per-task scratch, what we keep on failure. It cannot catch anything about the
tool itself, so this file runs the packaged tool for real and checks it does the job:
given three well-separated clades and one sequence deliberately labelled into the wrong
kingdom, SATIVA must report that sequence and no other.

The input is generated here from a fixed seed, so no data file is needed and the result
is reproducible. Three clades of eight members over 240 columns is enough for placement
to be unambiguous, and the whole thing takes a couple of seconds.

The tool comes from envs/sativa_epang.yaml and is not a pip dependency, so this test
SKIPS when it is not on PATH, which is the case in CI. To run it:

    conda env create -f workflow/envs/sativa_epang.yaml   # or reuse the pipeline's
    conda activate sativa_epang
    python workflow/scripts/Sativa/test_auto_sativa_real.py

Note that sativa-epang calls its siblings (epa-ng, and raxml through them) by bare name,
so its directory must be on PATH. run_sativa adds it when the binary is given as a path,
which the last check here covers.
"""
import importlib.util
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "auto_sativa", Path(__file__).with_name("auto_sativa.py"))
auto_sativa = importlib.util.module_from_spec(_SPEC)
sys.modules["auto_sativa"] = auto_sativa
_SPEC.loader.exec_module(auto_sativa)

BINARY = os.environ.get("SATIVA_EPANG_BINARY") or "sativa-epang"
resolved = shutil.which(BINARY) if os.sep not in BINARY else BINARY
if not resolved or not Path(resolved).exists():
    print("SKIPPED: %r is not available; see this file's docstring" % BINARY)
    sys.exit(0)

passed = failed = 0


def check(name, ok):
    global passed, failed
    passed += bool(ok)
    failed += (not ok)
    print(f"  {name:<58} {'OK' if ok else 'ECHEC'}")


CLADES = [
    ("Fungi", "Basidiomycota", "Agaricomycetes", "Russulales", "Russulaceae", "Russula"),
    ("Plantae", "Tracheophyta", "Magnoliopsida", "Ericales", "Primulaceae", "Primula"),
    ("Animalia", "Arthropoda", "Insecta", "Diptera", "Culicidae", "Culex"),
]
PER_CLADE = 8
COLUMNS = 240
MISLABELLED = "russula_00"


def synthetic_input(directory):
    """Three distinct clades, one member of the first labelled as the second."""
    rng = random.Random(20261024)
    records, taxa = [], {}
    for index, clade in enumerate(CLADES):
        backbone = random.Random(100 + index)
        base = [backbone.choice("ACGT") for _ in range(COLUMNS)]
        for member in range(PER_CLADE):
            seq = list(base)
            for _ in range(COLUMNS // 25):
                seq[rng.randrange(COLUMNS)] = rng.choice("ACGT")
            seq_id = "%s_%02d" % (clade[5].lower(), member)
            records.append((seq_id, "".join(seq)))
            taxa[seq_id] = list(clade) + ["%s espece%02d" % (clade[5], member)]
    taxa[MISLABELLED] = list(CLADES[1]) + ["Primula imposteur"]

    fasta, tax = directory / "t.fasta", directory / "t.tax"
    with open(fasta, "w") as fh:
        for seq_id, seq in records:
            fh.write(">%s\n%s\n" % (seq_id, seq))
    with open(tax, "w") as fh:
        for seq_id, ranks in taxa.items():
            fh.write("%s\t%s\n" % (seq_id, ";".join(ranks)))
    return fasta, tax


def mislabels(path):
    """The .mis rows, header comment lines dropped."""
    rows = []
    for line in open(path, encoding="utf-8"):
        if line.strip() and not line.startswith(";"):
            rows.append(line.rstrip("\n").split("\t"))
    return rows


print("Le vrai sativa-epang sur une entree synthetique")
print("  binaire : %s" % resolved)

version = auto_sativa._tool_version(BINARY)
check("la version est lue sur le binaire", "sativa" in version.lower())
check("et ce n est pas la valeur de repli", "unknown" not in version)
print("  version : %s" % version)

work = Path(tempfile.mkdtemp(prefix="sativa_real_"))
try:
    fasta, tax = synthetic_input(work)
    ok = auto_sativa.run_sativa(fasta, tax, work, BINARY, mode="ultrafast",
                                confidence=0.6, taxonomy_code="BOT", threads=2,
                                output_name="t", rand_seed=1)
    check("le run reussit", ok is True)

    produced = sorted(work.glob("*.mis"))
    check("un fichier .mis est produit", len(produced) == 1)
    check("le scratch de la tache est supprime apres succes",
          not (work / "sativa_tmp").exists())

    rows = mislabels(produced[0]) if produced else []
    check("la sequence mal etiquetee est signalee",
          any(r[0] == MISLABELLED for r in rows))
    check("et elle est la seule", len(rows) == 1)

    row = next((r for r in rows if r[0] == MISLABELLED), None)
    check("le desaccord est situe au rang Kingdom", row is not None and row[1] == "Kingdom")
    check("l etiquette d origine est bien celle qu on a falsifiee",
          row is not None and row[2] == "Plantae")
    check("le rang propose est le vrai clade de la sequence",
          row is not None and row[3] == "Fungi")
    if row:
        print("    signale : %s  %s  %s -> %s" % (row[0], row[1], row[2], row[3]))

    # Same run, but the binary given as an absolute path with its directory NOT on PATH.
    # Without run_sativa putting that directory on PATH, the tool dies on "epa-ng".
    if os.sep not in BINARY:
        stripped = dict(os.environ)
        kept = [d for d in stripped.get("PATH", "").split(os.pathsep)
                if d and Path(d).resolve() != Path(resolved).resolve().parent]
        second = Path(tempfile.mkdtemp(prefix="sativa_real_path_"))
        saved = os.environ.get("PATH", "")
        try:
            os.environ["PATH"] = os.pathsep.join(kept)
            fasta2, tax2 = synthetic_input(second)
            try:
                by_path = auto_sativa.run_sativa(fasta2, tax2, second, str(Path(resolved)),
                                                 mode="ultrafast", confidence=0.6,
                                                 taxonomy_code="BOT", threads=2,
                                                 output_name="t", rand_seed=1)
            except subprocess.CalledProcessError:
                # What this looks like without the PATH handling in run_sativa: the tool
                # starts, then dies on "FileNotFoundError: 'epa-ng'".
                by_path = False
            check("un binaire donne par chemin absolu marche hors PATH", by_path is True)
            check("et il produit le meme signalement",
                  by_path and [r[0] for r in mislabels(second / "t.mis")] == [MISLABELLED])
        finally:
            os.environ["PATH"] = saved
            shutil.rmtree(second, ignore_errors=True)
finally:
    shutil.rmtree(work, ignore_errors=True)

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
