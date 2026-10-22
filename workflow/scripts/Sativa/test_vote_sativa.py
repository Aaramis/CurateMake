"""Test cases for vote_sativa.

Run: python workflow/scripts/Sativa/test_vote_sativa.py

The invariant in test 1 is the one that would have caught the 2026-08-26 regression:
when every replicate is an identical split -- which is what happens whenever a barcode
fits in a single group -- the consensus MUST reproduce the original .mis line for line.
Observing "the vote changes nothing" is correct there; observing "the vote produces
nothing" means the chain is broken, and it silently degraded the whole pipeline to
CoL-only output (39.3% correction against 79.4% for the single-split path).
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VOTE = HERE / "vote_sativa.py"

HEADER = ";SeqID\tMislabeledLevel\tOriginalLabel\tProposedLabel\tConfidence\t" \
         "OriginalTaxonomyPath\tProposedTaxonomyPath\tPerRankConfidence\n"
ROWS = [
    ("s1", "Family", "Valsaceae", "Diaporthaceae", "1.000",
     "Fungi;Ascomycota;NA;Diaporthales;Valsaceae;NA;NA",
     "Fungi;Ascomycota;NA;Diaporthales;Diaporthaceae;Diaporthe",
     "1.000;1.000;1.000;1.000;1.000;1.000"),
    ("s2", "Kingdom", "Plantae", "Fungi", "0.900",
     "Plantae;Tracheophyta;NA;Asterales;NA;NA;NA",
     "Fungi;Basidiomycota;NA", "1.000;0.900;0.900"),
]


def build(root: Path, n_rep: int, *, write_summary=True, write_mis=True):
    """A barcode of ONE group, replicated n_rep times with an IDENTICAL split."""
    sativa = root / "SATIVA"
    tasks = []
    for rep in range(n_rep):
        taxon = "ITS2__groups__0000" if rep == 0 else f"ITS2__groups__r{rep}__0000"
        d = sativa / "ITS2" / "alignments_groups" / taxon
        d.mkdir(parents=True, exist_ok=True)
        if write_summary:
            (d / "sativa.summary.json").write_text(json.dumps({"status": "success"}))
        if write_mis:
            with open(d / f"ITS2_groups_{taxon}.mis", "w") as fh:
                fh.write(HEADER)
                for r in ROWS:
                    fh.write("\t".join(r) + "\n")
        tasks.append({"barcode": "ITS2", "tax_level": "groups", "taxon": taxon,
                      "replicate": rep, "msa_path": "", "n": 2})
    (root / "manifest.json").write_text(json.dumps({"tasks": tasks}))
    (root / "members.tsv").write_text(
        "tag\tmember_seq_id\n" + "".join(f"ITS2__groups__0000\t{r[0]}\n" for r in ROWS))
    return sativa


def run_vote(root, sativa, n_rep, threshold="2/3"):
    return subprocess.run(
        [sys.executable, str(VOTE),
         "--manifest", str(root / "manifest.json"), "--sativa-dir", str(sativa),
         "--barcode", "ITS2", "--level", "groups", "--replicates", str(n_rep),
         "--native-members", str(root / "members.tsv"), "--vote-threshold", threshold,
         "--output", str(root / "consensus.mis"), "--report", str(root / "report.json")],
        capture_output=True, text=True)


def data_lines(path: Path):
    if not path.exists():
        return []
    return [l.rstrip("\n") for l in open(path) if not l.startswith(";") and l.strip()]


passed = failed = 0


def check(name, cond, detail=""):
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL  {name}  {detail}")


# --- 1. replicats identiques -> le consensus reproduit le .mis d'origine ---------------
for n_rep in (1, 3, 10):
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        sativa = build(root, n_rep)
        r = run_vote(root, sativa, n_rep)
        check(f"identical x{n_rep}: exit 0", r.returncode == 0, r.stderr[-200:])
        got = data_lines(root / "consensus.mis")
        check(f"identical x{n_rep}: {len(ROWS)} lignes", len(got) == len(ROWS),
              f"got {len(got)}")
        ids = [l.split("\t")[0] for l in got]
        check(f"identical x{n_rep}: memes sequences", ids == [r[0] for r in ROWS], str(ids))
        for line in got:
            f = line.split("\t")
            src = next(r for r in ROWS if r[0] == f[0])
            check(f"identical x{n_rep}: lignee de {f[0]} preservee", f[6] == src[6],
                  f"{f[6]!r} != {src[6]!r}")
            check(f"identical x{n_rep}: {f[0]} a {n_rep}/{n_rep} voix",
                  f[8] == str(n_rep) and f[9] == str(n_rep), f"{f[8]}/{f[9]}")

# --- 2. le vote lance AVANT sativa doit echouer, pas produire un consensus vide --------
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    sativa = build(root, 3, write_summary=False, write_mis=False)
    r = run_vote(root, sativa, 3)
    check("sans resume: sortie en erreur", r.returncode != 0, f"rc={r.returncode}")
    check("sans resume: message explicite", "ran BEFORE SATIVA" in r.stderr, r.stderr[-200:])
    check("sans resume: aucun consensus ecrit", not (root / "consensus.mis").exists())

# --- 3. resumes presents mais aucun .mis -> echec aussi --------------------------------
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    sativa = build(root, 3, write_mis=False)
    r = run_vote(root, sativa, 3)
    check("sans .mis: sortie en erreur", r.returncode != 0, f"rc={r.returncode}")
    check("sans .mis: message explicite", "not one .mis found" in r.stderr, r.stderr[-200:])

# --- 4. seuil : 1 voix sur 3 ne doit pas passer le 2/3 ---------------------------------
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    sativa = build(root, 3)
    # ne garder le signalement de s2 que dans le replicat 0
    for rep in (1, 2):
        taxon = f"ITS2__groups__r{rep}__0000"
        p = sativa / "ITS2" / "alignments_groups" / taxon / f"ITS2_groups_{taxon}.mis"
        with open(p, "w") as fh:
            fh.write(HEADER)
            fh.write("\t".join(ROWS[0]) + "\n")
    r = run_vote(root, sativa, 3)
    got = data_lines(root / "consensus.mis")
    ids = [l.split("\t")[0] for l in got]
    check("seuil 2/3: s1 (3 voix) retenu", "s1" in ids, str(ids))
    check("seuil 2/3: s2 (1 voix) ecarte", "s2" not in ids, str(ids))


# --- 5. trou en milieu de lignee : NA au class, annotations apres -----------------------
# Le pipeline produit presque toujours des lignees a trou (class = NA). Trois regimes,
# selon combien de votants renseignent le rang troue.
def build_class_case(root: Path, n_with_class: int, n_rep: int = 10, conflicting=False):
    """n_with_class votants renseignent le class ; tous s'accordent sur le reste."""
    sativa = root / "SATIVA"
    tasks = []
    deep = "Diaporthales;Valsaceae;Valsa"
    for rep in range(n_rep):
        taxon = "ITS2__groups__0000" if rep == 0 else f"ITS2__groups__r{rep}__0000"
        d = sativa / "ITS2" / "alignments_groups" / taxon
        d.mkdir(parents=True, exist_ok=True)
        (d / "sativa.summary.json").write_text(json.dumps({"status": "success"}))
        if rep < n_with_class:
            cls = "Sordariomycetes" if not conflicting or rep % 2 == 0 else "Dothideomycetes"
        else:
            cls = "NA"
        with open(d / f"ITS2_groups_{taxon}.mis", "w") as fh:
            fh.write(HEADER)
            fh.write("\t".join(("s1", "Family", "X", "Valsa", "1.000",
                                 "Fungi;Ascomycota;NA;Diaporthales;Old;NA;NA",
                                 f"Fungi;Ascomycota;{cls};{deep};NA",
                                 "1.000;1.000;1.000;1.000;1.000;1.000;1.000")) + "\n")
        tasks.append({"barcode": "ITS2", "tax_level": "groups", "taxon": taxon,
                      "replicate": rep, "msa_path": "", "n": 1})
    (root / "manifest.json").write_text(json.dumps({"tasks": tasks}))
    (root / "members.tsv").write_text("tag\tmember_seq_id\nITS2__groups__0000\ts1\n")
    return sativa


for n_with_class, label, expect_class in (
        (0,  "aucun votant ne renseigne le class", ""),
        (4,  "minorite unanime (4/10)",            ""),
        (7,  "majorite (7/10)",                    "Sordariomycetes")):
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        sativa = build_class_case(root, n_with_class)
        r = run_vote(root, sativa, 10)
        check(f"trou/{label}: exit 0", r.returncode == 0, r.stderr[-200:])
        got = data_lines(root / "consensus.mis")
        check(f"trou/{label}: 1 ligne", len(got) == 1, str(got))
        if got:
            path = got[0].split("\t")[6].split(";")
            check(f"trou/{label}: class = {expect_class or 'NA'}",
                  (path[2] if path[2] != "NA" else "") == expect_class, str(path))
            # le point crucial : les rangs PLUS PROFONDS survivent
            check(f"trou/{label}: order conserve",   len(path) > 3 and path[3] == "Diaporthales", str(path))
            check(f"trou/{label}: family conservee", len(path) > 4 and path[4] == "Valsaceae", str(path))
            check(f"trou/{label}: genus conserve",   len(path) > 5 and path[5] == "Valsa", str(path))

# conflit reel au class (5 disent Sordariomycetes, 5 Dothideomycetes) -> troncature
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    sativa = build_class_case(root, 10, conflicting=True)
    r = run_vote(root, sativa, 10)
    got = data_lines(root / "consensus.mis")
    if got:
        path = got[0].split("\t")[6].split(";")
        check("conflit reel au class: lignee tronquee avant order", len(path) <= 2, str(path))
    else:
        check("conflit reel au class: ligne conservee", False, "aucune ligne")


# --- 6. le vote doit s'importer SANS la pile lourde ------------------------------------
# workflow/envs/sativa.yaml ne declare ni requests ni les clients HTTP. Les regles
# `script:` heritent du PYTHONPATH de l'environnement Snakemake et y trouvaient requests
# par accident ; une regle `shell:` comme sativa_vote n'herite de rien. Importer
# integrate_sativa (-> col_taxonomic_matcher -> requests) tuait donc le job AVANT qu'il
# n'ouvre son log : 0 octet, pas de consensus, et toutes les corrections SATIVA perdues
# en silence. Le vote ne doit dependre que de la bibliotheque standard.
_probe = subprocess.run(
    [sys.executable, "-c",
     "import sys; sys.path.insert(0, %r); import vote_sativa; "
     "heavy = [m for m in ('requests','sqlalchemy','psycopg2','pandas','plotly') "
     "         if m in sys.modules]; print(','.join(heavy))" % str(HERE)],
    capture_output=True, text=True)
check("import du vote: reussi", _probe.returncode == 0, _probe.stderr[-200:])
check("import du vote: aucune dependance lourde tiree",
      _probe.stdout.strip() == "", f"tire: {_probe.stdout.strip()}")

print(f"\n{'PASSED' if failed == 0 else 'FAILED'} {passed}/{passed + failed}")
sys.exit(0 if failed == 0 else 1)
