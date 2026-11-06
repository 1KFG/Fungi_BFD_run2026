#!/usr/bin/env python3
"""Build bench.fa + labels.tsv from the UniProt set downloads (<set>.tsv/.fasta).
soluble is sampled to 1000 (seed 20261106); all other sets are used whole."""
import csv, random
random.seed(20261106)
sets = ["featexp", "multipass", "singlepass", "peripheral", "signal", "secreted", "soluble"]
labels = {}
for s in sets:
    acc = [r["Entry"] for r in csv.DictReader(open(f"{s}.tsv"), delimiter="\t")]
    if s == "soluble":
        acc = sorted(random.sample(acc, 1000))
    for a in acc:
        labels.setdefault(a, []).append(s)
seqs = {}
for s in sets:
    cur, buf = None, []
    for line in list(open(f"{s}.fasta")) + [">|END|"]:
        if line.startswith(">"):
            if cur in labels and cur not in seqs:
                seqs[cur] = "".join(buf)
            cur, buf = line.split("|")[1], []
        else:
            buf.append(line.strip())
with open("bench.fa", "w") as f:
    for a in sorted(seqs):
        f.write(f">{a}\n{seqs[a]}\n")
with open("labels.tsv", "w") as f:
    f.write("accession\tsets\tlength\n")
    for a in sorted(seqs):
        f.write(f"{a}\t{','.join(labels[a])}\t{len(seqs[a])}\n")
print(len(seqs), "proteins")
