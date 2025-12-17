#!/usr/bin/env python3
from pathlib import Path

msa_root = Path("data/MSA/ITS2")
records = []
total = 0

for fasta in msa_root.rglob("*_msa.trimmed.fasta"):
    count = sum(1 for line in fasta.open() if line.startswith(">"))
    records.append((str(fasta), count))
    total += count

for path, n in sorted(records):
    print(f"{path}\t{n}")

print(total)
