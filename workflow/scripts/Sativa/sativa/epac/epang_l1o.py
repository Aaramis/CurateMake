#!/usr/bin/env python3
"""Leave-one-out via EPA-ng, en remplacement de RAxML `-f O` dans SATIVA.

Produit une liste de placements au MÊME format que EpaJsonParser.get_placement(),
avec des numéros d'arête = numérotation B= du refjson (celle qu'attend
classify_seq via bid_taxonomy_map). La classification/décision reste 100% SATIVA.

Approche : k-fold. Pour chaque pli, on retire ses feuilles de l'arbre+alignement de
référence et on place les séquences du pli avec EPA-ng ; on remappe les arêtes EPA-ng
vers B= par bipartition de feuilles. Tous les placements sont conservés
(--filter-acc-lwr 0.99999 --filter-max 100000) pour retrouver la masse de LWR complète.
"""
import os, sys, re, json, glob, subprocess, shutil
sys.setrecursionlimit(200000)
from ete3 import Tree

# epa-ng résolu depuis le PATH (fourni par l'env conda sativa.yaml) ; surchargeable
# via CURATEMAKE_EPANG_BIN.
EPANG = os.environ.get("CURATEMAKE_EPANG_BIN") or shutil.which("epa-ng") or "epa-ng"

def _bip_map(nhx_tree_str, tag):
    if tag == "EDGE":
        nhx_tree_str = re.sub(r"\{(\d+)\}", r"[&&NHX:EDGE=\1]", nhx_tree_str)
    t = Tree(nhx_tree_str, format=1)
    allL = frozenset(t.get_leaf_names())
    m = {}
    for n in t.traverse():
        e = getattr(n, tag, None)
        if e is None:
            continue
        desc = frozenset(n.get_leaf_names())
        side = desc if len(desc) <= len(allL) - len(desc) else allL - desc
        m[frozenset(side)] = str(e)
    return m, allL

def _read_fasta(path):
    seqs, cur = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                cur = line[1:].split()[0]; seqs[cur] = []
            else:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def _write_fasta(d, path):
    with open(path, "w") as f:
        for k, v in d.items():
            f.write(">%s\n%s\n" % (k, v))

def run_epang_l1o(refjson_tree_str, refaln_path, reftree_path, raxml_outdir,
                  workdir, folds=5, threads=1, log=None):
    def _log(m):
        if log: log.info("[epang-l1o] " + m)
        else: sys.stderr.write("[epang-l1o] " + m + "\n")

    # modèle EPA-ng : RAxML_info.mfresolv si présent, sinon GTR+G (EPA-ng ré-évalue)
    info = glob.glob(os.path.join(raxml_outdir, "RAxML_info.mfresolv*"))
    model = info[0] if info else "GTR+G"
    _log("modèle EPA-ng: %s" % model)

    bidmap, allL = _bip_map(refjson_tree_str, "B")        # bipartition -> B-id
    aln = _read_fasta(refaln_path)
    # aligner les clés de l'alignement sur les noms de feuilles
    leaves = list(allL)
    aln_by_leaf = {}
    for lf in leaves:
        for cand in (lf, lf[2:] if lf.startswith("r_") else "r_"+lf):
            if cand in aln:
                aln_by_leaf[lf] = aln[cand]; break
    missing = [l for l in leaves if l not in aln_by_leaf]
    if missing:
        raise RuntimeError("epang-l1o: %d feuilles sans séquence (ex %s)" % (len(missing), missing[:3]))

    full_tree = Tree(reftree_path, format=1)
    ordered = sorted(leaves)
    K = min(folds, len(ordered))
    folds_list = [ordered[i::K] for i in range(K)]

    os.makedirs(workdir, exist_ok=True)
    placements = []
    for fi, fold in enumerate(folds_list):
        fold_set = set(fold)
        ref_leaves = [l for l in leaves if l not in fold_set]
        if len(ref_leaves) < 4 or not fold:
            continue
        wd = os.path.join(workdir, "fold_%03d" % fi); os.makedirs(wd, exist_ok=True)
        tpr = full_tree.copy(method="newick")
        tpr.prune(ref_leaves, preserve_branch_length=True)
        tpr.write(outfile=os.path.join(wd, "ref.nwk"), format=5)
        _write_fasta({l: aln_by_leaf[l] for l in ref_leaves}, os.path.join(wd, "ref.fasta"))
        _write_fasta({l: aln_by_leaf[l] for l in fold},       os.path.join(wd, "query.fasta"))
        cmd = [EPANG, "-t", os.path.join(wd, "ref.nwk"), "-s", os.path.join(wd, "ref.fasta"),
               "-q", os.path.join(wd, "query.fasta"), "-m", model, "--outdir", wd, "--redo",
               "-T", str(threads), "--filter-acc-lwr", "0.99999", "--filter-max", "100000"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        jpf = os.path.join(wd, "epa_result.jplace")
        if r.returncode != 0 or not os.path.isfile(jpf):
            _log("EPA-ng FAIL fold %d: %s" % (fi, r.stderr[-300:])); continue
        d = json.load(open(jpf))
        ie = d["fields"].index("edge_num")
        e2side, foldL = _bip_map(d["tree"], "EDGE")
        # table bipartition-restreinte -> B
        restricted = {}
        for side, b in ((s, bidmap[s]) for s in bidmap):
            s2 = side & foldL
            key = s2 if len(s2) <= len(foldL)-len(s2) else (foldL - s2)
            restricted.setdefault(frozenset(key), b)
        epa2b = {}
        for side, e in e2side.items():
            key = side if len(side) <= len(foldL)-len(side) else (foldL - side)
            b = restricted.get(frozenset(key))
            if b is not None:
                epa2b[int(e)] = b
        for pl in d["placements"]:
            name = (pl.get("n") or pl.get("nm"))[0]
            if isinstance(name, list): name = name[0]
            newp = []
            for row in pl["p"]:
                b = epa2b.get(int(row[ie]))
                if b is None: continue
                rr = list(row); rr[ie] = int(b)
                newp.append(rr)
            if newp:
                placements.append({"p": newp, "n": [name]})
    _log("placements produits: %d (K=%d folds)" % (len(placements), K))
    return placements


def run_epang_final(reftree_path, refaln_path, raxml_outdir, workdir, threads=1, log=None):
    """Passe 2 (confirmation) via EPA-ng, en remplacement du RAxML `-f v` de run_epa_once.

    Place les séquences 'suspectes' (celles élaguées de l'arbre = alignement de réf MOINS
    les feuilles de l'arbre élagué) sur l'arbre de référence débarrassé des mislabels, et
    écrit un jplace directement consommable par EpaJsonParser (arbre {N} + placements
    auto-cohérents : le bid_tax_map est reconstruit à partir de cet arbre par SATIVA).
    Renvoie le chemin du jplace, ou None si aucune séquence suspecte.
    """
    def _log(m):
        (log.info if log else (lambda x: sys.stderr.write(x + "\n")))("[epang-final] " + m)

    os.makedirs(workdir, exist_ok=True)
    info = glob.glob(os.path.join(raxml_outdir, "RAxML_info.mfresolv*"))
    model = info[0] if info else "GTR+G"

    tree = Tree(reftree_path, format=1)
    ref_leaves = set(tree.get_leaf_names())
    aln = _read_fasta(refaln_path)
    ref_seqs = {n: s for n, s in aln.items() if n in ref_leaves}
    query_seqs = {n: s for n, s in aln.items() if n not in ref_leaves}
    if not query_seqs:
        _log("aucune séquence suspecte à replacer")
        return None
    _write_fasta(ref_seqs, os.path.join(workdir, "ref.fasta"))
    _write_fasta(query_seqs, os.path.join(workdir, "query.fasta"))

    cmd = [EPANG, "-t", reftree_path, "-s", os.path.join(workdir, "ref.fasta"),
           "-q", os.path.join(workdir, "query.fasta"), "-m", model,
           "--outdir", workdir, "--redo", "-T", str(threads),
           "--filter-acc-lwr", "0.99999", "--filter-max", "100000"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    jpf = os.path.join(workdir, "epa_result.jplace")
    if r.returncode != 0 or not os.path.isfile(jpf):
        raise RuntimeError("EPA-ng final placement failed: " + (r.stderr[-300:] if r.stderr else "no jplace"))
    _log("%d suspects replacés (modèle %s)" % (len(query_seqs), os.path.basename(model)))
    return jpf
