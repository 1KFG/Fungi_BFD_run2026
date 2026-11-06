#!/usr/bin/env python3
"""sample_bfd_neg.py -- the TMHMM-negative strata of the BFD benchmark sample.

tables/tmhmm.parquet holds only proteins with PredHel >= 1 (merge_tmhmm.py skips
0-helix proteins), so negatives come from the per-genome TMHMM short files
results/function/tmhmm/<bucket>/<LOCUSTAG>.tmhmm_short.tsv.gz.

Strata (2x the target; held-out filtering happens later, as for sample_bfd.sh):
  tm_neg   2000  PredHel = 0, spread evenly over phyla
  sp_only   600  PredHel = 0 and a SignalP 6 signal peptide (tables/signalp...)
Genomes: up to 8 per phylum per round, chosen by a seeded shuffle; one random
protein per genome per stratum per round until the stratum is full.

Writes bfd_sample_neg/sample.tsv and bfd_sample_neg/bench.fa.
"""
import csv, gzip, hashlib, os, random, subprocess, sys
sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin")
from pick_representative_strain import make_sample_tag

R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
A = f"{R}/analysis/deeptmhmm_vs_tmhmm"
OUT = f"{A}/bfd_sample_neg"
SEED = 20261106
TARGET = {"tm_neg": 2000, "sp_only": 600}
rng = random.Random(SEED)
os.makedirs(OUT, exist_ok=True)

sup = {l.split(",")[0].strip() for l in open(f"{R}/suppress.txt") if l.strip() and not l.startswith("#")}
rows = [r for r in csv.DictReader(open(f"{R}/samples.csv")) if r["ASMID"].strip() not in sup]
by_phylum = {}
for r in rows:
    by_phylum.setdefault(r["PHYLUM"].strip() or "NA", []).append(r)
for v in by_phylum.values():
    rng.shuffle(v)

# SignalP SP protein ids, read once
sp_ids = set(subprocess.run(
    ["duckdb", "-csv", "-noheader", "-c",
     f"SELECT DISTINCT protein_id FROM '{R}/tables/signalp.signal_peptide.parquet'"],
    capture_output=True, text=True, check=True).stdout.split())

def short_path(tag):
    b = hashlib.sha1(tag.encode()).hexdigest()[:2]
    return f"{R}/BFD_function/results/function/tmhmm/{b}/{tag}.tmhmm_short.tsv.gz"

picked = {k: [] for k in TARGET}
used = set()
rnd = 0
while any(len(picked[k]) < TARGET[k] for k in TARGET) and rnd < 50:
    for ph, lst in sorted(by_phylum.items()):
        for r in lst[rnd * 8:(rnd + 1) * 8]:
            tag = r["LOCUSTAG"].strip()
            p = short_path(tag)
            if not os.path.exists(p):
                continue
            neg = []
            for line in gzip.open(p, "rt"):
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5 and f[4] == "PredHel=0":
                    neg.append((f[0], int(f[1].split("=")[1])))
            if not neg:
                continue
            for k in TARGET:
                if len(picked[k]) >= TARGET[k]:
                    continue
                pool = [x for x in neg if x[0] not in used and (k != "sp_only" or x[0] in sp_ids)]
                if pool:
                    pid, ln = rng.choice(pool)
                    used.add(pid)
                    picked[k].append((k, pid, tag, ph, ln, make_sample_tag(r["SPECIES"], r["STRAIN"])))
    rnd += 1

with open(f"{OUT}/sample.tsv", "w") as fh:
    fh.write("stratum\tprotein_id\tLOCUSTAG\tPHYLUM\tlen\tPredHel\tsignalp_sp\n")
    for k in TARGET:
        for (s, pid, tag, ph, ln, _o) in picked[k]:
            fh.write(f"{s}\t{pid}\t{tag}\t{ph}\t{ln}\t0\t{pid in sp_ids}\n")
want = {}
for k in TARGET:
    for (_s, pid, tag, _ph, _ln, o) in picked[k]:
        want.setdefault((tag, o), set()).add(pid)
seqs = {}
for (tag, o), ids in want.items():
    cur, buf = None, []
    for line in list(gzip.open(f"{R}/genome_annotation/{o}/predict_results/{o}.proteins.fa.gz", "rt")) + [">END"]:
        if line.startswith(">"):
            if cur in ids:
                seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
with open(f"{OUT}/bench.fa", "w") as fh:
    for pid in sorted(seqs):
        fh.write(f">{pid}\n{seqs[pid].rstrip('*')}\n")
for k in TARGET:
    print(k, len(picked[k]), "proteins,", len({x[3] for x in picked[k]}), "phyla,", len({x[2] for x in picked[k]}), "genomes")
print("bench.fa:", len(seqs), "sequences")
