#!/usr/bin/env python3
"""
Bottom-up construction of the reference MSAs for ONE (barcode, family).

Replaces the old top-down descent (per-clade align/trim). Logic validated in
17.msa_bottomup (see INTEGRATION_DESIGN.md):

  FLOOR (per species / genus-only / family-only):
    deterministic sort -> vsearch (homology clusters) -> per cluster:
    MAFFT -> trimAl -gappyout (trims flanks / length outliers, keeps the sequences)
  GATED CLIMB:
    merge a species' clusters -> species (or separate variants if divergent)
    -> merge species -> genus -> merge genera -> family
    add-at-level: sequences annotated EXACTLY at a level are added there (mafft --add)
    GATE at the climb (on the RAW merge): occ >= gate_occ AND ratio <= gate_ratio

Outputs:
  <out>/msa/<msa_id>.fasta   one reference MSA per final unit (n >= min_ref)
  <out>/residual.fasta       singletons + unplaced sequences -> downstream hmmsearch classif
  <out>/manifest.json        stop levels, metrics, residual (for reports + downstream)

Deterministic: sorted inputs, MAFFT FFT-NS-2 single-thread, fixed vsearch order.
"""

import argparse
import csv
import json
import logging
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger("build_family_msa")

Rec = Tuple[str, str]          # (seq_id, raw ungapped sequence)
Aln = List[Tuple[str, str]]    # [(seq_id, aligned sequence)]


# --------------------------------------------------------------------------- IO
def read_fasta(path: str) -> Aln:
    recs, h, s = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(s)))
                h, s = line[1:].split()[0], []
            else:
                s.append(line)
    if h is not None:
        recs.append((h, "".join(s)))
    return recs


def write_fasta(recs: Aln, path: str) -> None:
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n{s}\n")


RANKS_FULL = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def taxed_header(seq_id: str, tax) -> str:
    """Header in the format expected by SATIVA (auto_sativa):
    >seq_id|kingdom=..|phylum=..|...|species=Genus_species (spaces -> _)."""
    parts = [seq_id]
    for rank, val in zip(RANKS_FULL, tax or []):
        v = (val or "").strip().replace(" ", "_")
        if v:
            parts.append(f"{rank}={v}")
    return "|".join(parts)


def write_fasta_taxed(recs: Aln, path: str, tax_map) -> None:
    """Write a FASTA with taxonomy embedded in the headers (for SATIVA)."""
    with open(path, "w") as fh:
        for seq_id, s in recs:
            fh.write(f">{taxed_header(seq_id, tax_map.get(seq_id))}\n{s}\n")


# --------------------------------------------------------------- external tools
def _run(cmd: List[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                          text=True, **kw)


def mafft_align(recs: List[Rec], tmp: str, tag: str) -> Aln:
    """De novo alignment, deterministic (FFT-NS-2, single-thread, input order)."""
    if len(recs) == 1:
        return list(recs)
    inp, out = f"{tmp}/{tag}.in", f"{tmp}/{tag}.out"
    write_fasta(recs, inp)
    with open(out, "w") as fo:
        r = subprocess.run(["mafft", "--retree", "2", "--maxiterate", "0", "--quiet",
                            "--thread", "1", "--inputorder", inp],
                           stdout=fo, stderr=subprocess.DEVNULL)
    return read_fasta(out) if r.returncode == 0 else list(recs)


def mafft_merge(subalns: List[Aln], tmp: str, tag: str) -> Aln:
    """Merge sub-alignments while preserving their internal alignment (mafft --merge)."""
    subalns = [a for a in subalns if a]
    if not subalns:
        return []
    if len(subalns) == 1:
        return subalns[0]
    inp, tab, out = f"{tmp}/{tag}.in", f"{tmp}/{tag}.tab", f"{tmp}/{tag}.out"
    idx, lines = 1, []
    with open(inp, "w") as fh:
        for aln in subalns:
            grp = []
            for h, s in aln:
                fh.write(f">{h}\n{s}\n")
                grp.append(str(idx)); idx += 1
            lines.append(" ".join(grp))
    with open(tab, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(out, "w") as fo:
        r = subprocess.run(["mafft", "--merge", tab, "--quiet", "--thread", "1", inp],
                           stdout=fo, stderr=subprocess.DEVNULL)
    return read_fasta(out) if r.returncode == 0 else subalns[0]


def mafft_add(base: Aln, raw: List[Rec], tmp: str, tag: str) -> Aln:
    """Add raw sequences (annotated at this level) to an existing alignment."""
    if not raw:
        return base
    if not base:
        return mafft_align(raw, tmp, tag)
    b, a, out = f"{tmp}/{tag}.b", f"{tmp}/{tag}.a", f"{tmp}/{tag}.o"
    write_fasta(base, b)
    write_fasta(raw, a)
    with open(out, "w") as fo:
        r = subprocess.run(["mafft", "--add", a, "--reorder", "--quiet", "--thread", "1", b],
                           stdout=fo, stderr=subprocess.DEVNULL)
    return read_fasta(out) if r.returncode == 0 else base


def hmmalign_add(base: Aln, raw: List[Rec], tmp: str, tag: str) -> Aln:
    """Re-add raw sequences into a reference alignment via an ephemeral profile +
    hmmalign --mapali (much cheaper than mafft --add for large near-identical batches:
    this is the mechanism used to cover the surplus sequences beyond max_group without
    re-aligning de novo). Fallback: mafft --add, then base."""
    if not raw:
        return base
    if not base:
        return mafft_align(raw, tmp, tag)
    b, hmm, add, out = f"{tmp}/{tag}.b", f"{tmp}/{tag}.hmm", f"{tmp}/{tag}.a", f"{tmp}/{tag}.o"
    write_fasta(base, b)
    if not hmmbuild(b, hmm, tag):
        return mafft_add(base, raw, tmp, tag)
    write_fasta(raw, add)
    r = subprocess.run(["hmmalign", "--dna", "--outformat", "afa", "--mapali", b,
                        "-o", out, hmm, add],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if r.returncode != 0 or not os.path.exists(out) or os.path.getsize(out) == 0:
        return mafft_add(base, raw, tmp, tag)
    merged = [(h, s.upper().replace(".", "-")) for h, s in read_fasta(out)]  # afa -> flat FASTA
    return merged if merged else base


def gappyout(aln: Aln, tmp: str, tag: str) -> Aln:
    """trimAl -gappyout: drops sparse columns (flanks/outliers), keeps the sequences."""
    if len(aln) < 4:
        return aln
    inp, out = f"{tmp}/{tag}.gi", f"{tmp}/{tag}.go"
    write_fasta(aln, inp)
    r = _run(["trimal", "-in", inp, "-out", out, "-gappyout"])
    if r.returncode != 0 or not os.path.exists(out) or os.path.getsize(out) == 0:
        return aln
    c = read_fasta(out)
    return c if c else aln


def hmmbuild(msa_path: str, hmm_path: str, name: str) -> bool:
    """Build a (nucleotide) HMM profile from a reference MSA."""
    r = _run(["hmmbuild", "--dna", "-n", name, hmm_path, msa_path])
    return r.returncode == 0 and os.path.exists(hmm_path) and os.path.getsize(hmm_path) > 0


def vsearch_cluster(recs: List[Rec], ident: float, tmp: str, tag: str) -> List[List[Rec]]:
    """Homology clusters (deterministic: input sorted by seq_id)."""
    if len(recs) < 3:
        return [recs]
    recs = sorted(recs, key=lambda r: r[0])
    fa, uc = f"{tmp}/{tag}.cf", f"{tmp}/{tag}.uc"
    write_fasta(recs, fa)
    r = _run(["vsearch", "--cluster_fast", fa, "--id", str(ident), "--uc", uc, "--threads", "1"])
    if r.returncode != 0 or not os.path.exists(uc):
        return [recs]
    smap = dict(recs); clu: Dict[int, List[Rec]] = defaultdict(list)
    with open(uc) as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in ("S", "H"):
                sid = f[8].split()[0]
                if sid in smap:
                    clu[int(f[1])].append((sid, smap[sid]))
    return [clu[k] for k in sorted(clu)] or [recs]


# ------------------------------------------------------------------- metrics
def metrics(aln: Aln) -> Optional[Dict]:
    if not aln:
        return None
    seqs = [s for _, s in aln]; n = len(seqs); ncols = len(seqs[0])
    if ncols == 0:
        return None
    ungapped = sorted(len(s.replace("-", "")) for s in seqs)
    med = ungapped[n // 2] or 1
    occ_cols = sum(1 for c in range(ncols)
                   if sum(1 for s in seqs if s[c] != "-") / n >= 0.5)
    return {"n": n, "ncols": ncols, "med_len": med,
            "ratio": round(ncols / med, 2),
            "occ": round(100 * occ_cols / ncols, 1)}


# ------------------------------------------------------------------- algorithm
class FamilyBuilder:
    def __init__(self, cfg: Dict, tmp: str):
        self.cluster_id = cfg["cluster_id"]
        self.gate_ratio = cfg["gate_ratio"]
        self.gate_occ = cfg["gate_occ"]
        self.max_group = cfg["max_group"]
        self.tmp = tmp
        self._uid = 0

    def uid(self) -> int:
        self._uid += 1
        return self._uid

    def gate(self, aln: Aln) -> bool:
        m = metrics(aln)
        return bool(m and m["ratio"] <= self.gate_ratio and m["occ"] >= self.gate_occ)

    def floor(self, recs: List[Rec]) -> List[Aln]:
        """cluster -> per cluster (align + gappyout). Returns clean MSAs.

        Compute cap: only <= max_group representatives are aligned DE NOVO (beyond that
        the alignment does not improve). The surplus sequences are NOT discarded: they
        are re-added to the same MSA via hmmalign (cheap) -> every sequence stays covered
        by the downstream SATIVA QC (cf. design §4.3, "keep the memberships")."""
        leaves = []
        for c in vsearch_cluster(recs, self.cluster_id, self.tmp, f"cl{self.uid()}"):
            ordered = sorted(c, key=lambda r: r[0])       # deterministic
            sub, held = ordered[:self.max_group], ordered[self.max_group:]
            a = gappyout(mafft_align(sub, self.tmp, f"al{self.uid()}"), self.tmp, f"go{self.uid()}")
            if held:
                a = hmmalign_add(a, held, self.tmp, f"ha{self.uid()}")
            leaves.append(a)
        return leaves

    def species(self, recs: List[Rec]) -> List[Aln]:
        """Leaves of a species: re-merge similar clusters, keep real variants separate."""
        clusters = self.floor(recs)
        if len(clusters) <= 1:
            return clusters
        merged = mafft_merge(clusters, self.tmp, f"sm{self.uid()}")
        return [merged] if self.gate(merged) else clusters

    def climb(self, children: List[Aln], own_raw: List[Rec]) -> Optional[Aln]:
        """Merge children + add sequences annotated at this level, if the GATE passes (raw)."""
        merged = mafft_merge(children, self.tmp, f"mg{self.uid()}") if children else []
        candidate = mafft_add(merged, own_raw, self.tmp, f"ad{self.uid()}") if own_raw else merged
        return candidate if candidate and self.gate(candidate) else None


def read_hierarchy(csv_path: str):
    """Robust parse (csv, handles quoting). Returns:
    (genus->species->[recs], genus->[genus-only], [family-only], tax_map seq_id->tax)."""
    gs: Dict[str, Dict[str, List[Rec]]] = defaultdict(lambda: defaultdict(list))
    go: Dict[str, List[Rec]] = defaultdict(list)
    fo: List[Rec] = []
    tax_map: Dict[str, List[str]] = {}
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("seq_id") or "").strip()
            seq = (row.get("sequence") or "").replace("-", "").replace(".", "").upper()
            if not sid or not seq:
                continue
            gen = (row.get("genus") or "").strip()
            spe = (row.get("species") or "").strip()
            tax_map[sid] = [(row.get(r) or "").strip() for r in RANKS_FULL]
            if not gen:
                fo.append((sid, seq))
            elif not spe:
                go[gen].append((sid, seq))
            else:
                gs[gen][spe].append((sid, seq))
    return gs, go, fo, tax_map


def build(csv: str, barcode: str, family: str, cfg: Dict, tmp: str) -> Tuple[List[Tuple], List[Rec]]:
    """Returns (final_msas [(level, taxon, aln)], residual [singleton recs])."""
    b = FamilyBuilder(cfg, tmp)
    gs, go, fo, tax_map = read_hierarchy(csv)

    genus_final: List[Tuple[str, str, Aln]] = []
    for gen in sorted(set(list(gs) + list(go))):
        leaves: List[Tuple[str, str, Aln]] = []
        for spe in sorted(gs.get(gen, {})):
            for a in b.species(gs[gen][spe]):
                leaves.append(("species", f"{gen}|{spe}", a))
        for a in b.floor(go.get(gen, [])):
            leaves.append(("cluster", f"{gen}|(genus-only)", a))
        # try to climb to the genus
        climbed = b.climb([a for _, _, a in leaves], [])
        genus_final.extend([("genus", gen, climbed)] if climbed else leaves)

    # try to climb to the family (+ add the family-only sequences)
    fam_climbed = b.climb([a for _, _, a in genus_final], fo)
    if fam_climbed:
        final = [("family", family, fam_climbed)]
    else:
        final = list(genus_final)
        for a in b.floor(fo):
            final.append(("cluster", f"{family}|(family-only)", a))

    # split off the residual (1-sequence MSA = not validatable alone -> downstream hmmsearch classif)
    msas, residual = [], []
    for level, taxon, aln in final:
        if len(aln) < cfg["min_ref"]:
            residual.extend((h, s.replace("-", "")) for h, s in aln)
        else:
            msas.append((level, taxon, aln))
    return msas, residual, tax_map


# ------------------------------------------------------------------------- CLI
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Bottom-up MSA construction for one (barcode, family).")
    p.add_argument("--input-csv", required=True, help="Family CSV (seq_id,sequence,...,genus,species,barcode).")
    p.add_argument("--barcode", required=True)
    p.add_argument("--family", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cluster-id", type=float, default=0.97)
    p.add_argument("--gate-ratio", type=float, default=2.5)
    p.add_argument("--gate-occ", type=float, default=40.0)
    p.add_argument("--max-group", type=int, default=300)
    p.add_argument("--min-ref", type=int, default=2, help="Min size of a reference MSA (otherwise residual).")
    p.add_argument("--fam-id", help="Filesystem-safe family identifier (prefix of the msa_id).")
    p.add_argument("--log-file")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        h = logging.FileHandler(args.log_file, mode="w")
        h.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(h)

    cfg = dict(cluster_id=args.cluster_id, gate_ratio=args.gate_ratio, gate_occ=args.gate_occ,
               max_group=args.max_group, min_ref=args.min_ref)
    fam_id = args.fam_id or re.sub(r"[^A-Za-z0-9._-]+", "_", f"{args.barcode}__{args.family}").strip("_")
    out = Path(args.out_dir)
    (out / "msa").mkdir(parents=True, exist_ok=True)
    (out / "hmm").mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="bfm_") as tmp:
        msas, residual, tax_map = build(args.input_csv, args.barcode, args.family, cfg, tmp)

    # write the MSAs + one HMM profile per MSA + the manifest
    entries, level_counts = [], defaultdict(int)
    for i, (level, taxon, aln) in enumerate(sorted(msas, key=lambda x: (x[0], x[1])), start=1):
        msa_id = f"{fam_id}__{i:04d}"
        msa_path = out / "msa" / f"{msa_id}.fasta"
        hmm_path = out / "hmm" / f"{msa_id}.hmm"
        write_fasta_taxed(aln, str(msa_path), tax_map)   # headers with taxonomy (SATIVA)
        profiled = hmmbuild(str(msa_path), str(hmm_path), msa_id)
        m = metrics(aln)
        # n_taxa = number of DISTINCT species: the real SATIVA validation power
        # (>=2 taxa are needed to detect a mislabel, not just >=N sequences).
        sp_idx = RANKS_FULL.index("species")
        n_taxa = len({(tax_map.get(sid) or [""] * 7)[sp_idx] for sid, _ in aln
                      if (tax_map.get(sid) or [""] * 7)[sp_idx]})
        entries.append({"msa_id": msa_id, "level": level, "taxon": taxon,
                        "n": m["n"], "n_taxa": n_taxa, "ratio": m["ratio"], "occ": m["occ"],
                        "gate_pass": bool(m["ratio"] <= args.gate_ratio and m["occ"] >= args.gate_occ),
                        "msa_path": str(msa_path),
                        "profile_path": str(hmm_path) if profiled else None})
        level_counts[level] += 1

    if residual:
        # embedded taxonomy: downstream integration propagates it into the MSA -> SATIVA
        # validates the orphans (and catches contaminated singletons). The raw id stays
        # extractable (split '|') for the reclassification hmmsearch.
        write_fasta_taxed(residual, str(out / "residual.fasta"), tax_map)

    manifest = {
        "barcode": args.barcode, "family": args.family,
        "params": cfg,
        "n_msa": len(entries), "n_residual": len(residual),
        "stop_levels": dict(level_counts),
        "gate_pass": sum(e["gate_pass"] for e in entries),
        "msas": entries,
    }
    with open(out / "manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("%s/%s: %d MSA (%s), residual=%d",
                args.barcode, args.family, len(entries), dict(level_counts), len(residual))


if __name__ == "__main__":
    main()
