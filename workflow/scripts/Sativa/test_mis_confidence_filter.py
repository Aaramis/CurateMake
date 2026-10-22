"""Test cases for the confidence cutoff applied when a .mis file is read.

Run: python workflow/scripts/Sativa/test_mis_confidence_filter.py

Why the cutoff exists. Measured on the August 2026 full run, over 905 403 flags, the
confidence histogram carries a discrete spike of 9 181 flags at EXACTLY 0.510, against
285 at 0.508, 627 at 0.509 and 66 at 0.511. That is the stem-placement signature: the
whole placement weight sits on a branch BETWEEN two ranks, so the placement cannot choose
between them and the confidence measures the tie, not the support. Cutting just above it
removes 1.19 % of flags.

The cutoff is applied when the .mis is read, not through SATIVA's own -C, so changing it
costs a re-integration rather than a full SATIVA re-run.
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parents[1]
for extra in (_SCRIPTS, _SCRIPTS / "utils", _SCRIPTS / "Sativa"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

_SPEC = importlib.util.spec_from_file_location(
    "integrate_sativa", Path(__file__).with_name("integrate_sativa.py"))
mod = importlib.util.module_from_spec(_SPEC)
sys.modules["integrate_sativa"] = mod
_SPEC.loader.exec_module(mod)

passed = failed = 0


def check(name, ok):
    global passed, failed
    passed += bool(ok)
    failed += (not ok)
    print(f"  {name:<58} {'OK' if ok else 'ECHEC'}")


TAXA = "Fungi;Basidiomycota;Agaricomycetes;Russulales;Russulaceae;Russula;Russula alpha"
TAXB = "Fungi;Basidiomycota;Agaricomycetes;Russulales;Russulaceae;Russula;Russula beta"


def mis(tmp, confidences):
    """A .mis with one flagged sequence per confidence, plus the comment header."""
    p = Path(tmp) / "t.mis"
    rows = ["; SeqID\tMislabeledLevel\tOriginalLabel\tProposedLabel\tConfidence"]
    for i, conf in enumerate(confidences):
        rows.append("\t".join([f"s{i}", "Species", "Russula alpha", "Russula beta",
                               str(conf), TAXA, TAXB, "1;1;1;1;1;1;0.5"]))
    p.write_text("\n".join(rows) + "\n")
    return p


SPREAD = [0.418, 0.508, 0.509, 0.510, 0.511, 0.75, 1.0]

print("\nSeuil de confiance a la lecture du .mis")

with tempfile.TemporaryDirectory() as tmp:
    path = mis(tmp, SPREAD)
    parser = mod.SativaMisParser()

    kept = parser.parse_mis_file(path, "ITS", "family")
    check("sans seuil, tout est conserve", len(kept) == len(SPREAD))

    kept = parser.parse_mis_file(path, "ITS", "family", min_confidence=0.0)
    check("un seuil nul se comporte comme pas de seuil", len(kept) == len(SPREAD))

    kept = parser.parse_mis_file(path, "ITS", "family", min_confidence=0.511)
    confs = sorted(round(r.confidence, 3) for r in kept)
    check("le seuil 0.511 retire le pic a 0.510", 0.510 not in confs)
    check("et garde 0.511 lui-meme", 0.511 in confs)
    check("et garde tout ce qui est au-dessus", confs == [0.511, 0.75, 1.0])
    check("et retire la queue basse", all(c >= 0.511 for c in confs))

    kept = parser.parse_mis_file(path, "ITS", "family", min_confidence=1.0)
    check("un seuil a 1.0 ne garde que la certitude parfaite",
          [round(r.confidence, 3) for r in kept] == [1.0])

    kept = parser.parse_mis_file(path, "ITS", "family", min_confidence=2.0)
    check("un seuil hors domaine ne garde rien, sans erreur", kept == [])

print("\nRobustesse")

with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp) / "bad.mis"
    p.write_text("; en-tete\n"
                 + "\t".join(["s0", "Species", "a", "b", "pas_un_nombre", TAXA, TAXB, "x"]) + "\n"
                 + "\t".join(["s1", "Species", "a", "b", "0.9", TAXA, TAXB, "x"]) + "\n")
    parser = mod.SativaMisParser()
    kept = parser.parse_mis_file(p, "ITS", "family", min_confidence=0.511)
    check("une confiance illisible ne fait pas tomber le fichier",
          [r.seq_id for r in kept] == ["s1"])

    empty = Path(tmp) / "empty.mis"
    empty.write_text("; rien que des commentaires\n")
    check("un .mis sans ligne de donnee donne zero sans erreur",
          parser.parse_mis_file(empty, "ITS", "family", min_confidence=0.511) == [])

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
