#!/usr/bin/python3.12
"""Per-genome gene statistics from the CDS files: unique genes, transcripts, mean CDS length, GC, GC3 mean and sd (over genes).
Usage: gene_stats.py samples.csv cds_dir out.tsv [threads]   (reads *.cds-transcripts.fa[.gz]; genome key = makeSampleTag(SPECIES, STRAIN))"""
import sys, os, re, gzip, csv
from multiprocessing import Pool
import numpy as np

def clean(st):
    st = (st or "").strip(); st = re.sub(r"['\"]", "", st).split(";")[0].strip().replace(":", " ")
    st = re.sub(r"^\s*\*+", "", st); st = re.sub(r"\*+\s*$", "", st); st = re.sub(r"\s*\*+\s*", "-", st); st = re.sub(r"[&`$|<>()]+", " ", st).strip(); return st
def tag(sp, st):
    sp = re.sub(r"['\"]", "", (sp or "").strip()); j = "_".join(x for x in [sp, clean(st)] if x); return re.sub(r"[\s/#\[\]?{}]+", "_", j)

def one(args):
    g, d = args
    p = f"{d}/{g}.cds-transcripts.fa.gz"
    op = gzip.open
    if not os.path.exists(p):
        p = f"{d}/{g}.cds-transcripts.fa"; op = open
    if not os.path.exists(p): return (g,) + (np.nan,) * 7
    genes, lens, gc3s, gcs = set(), [], [], []
    name, buf = None, []
    def flush():
        if name is None: return
        s = "".join(buf).upper(); n = len(s) - len(s) % 3
        if n < 90: return
        s = s[:n]; lens.append(n // 3)
        gcs.append((s.count("G") + s.count("C")) / n)
        t = s[2::3]; gc3s.append((t.count("G") + t.count("C")) / max(len(t), 1))
    with op(p, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                flush(); name = l[1:].split()[0]; buf = []; genes.add(re.sub(r"-T\d+$", "", name))
            else: buf.append(l.strip())
        flush()
    n_tr = len(lens)
    if n_tr == 0: return (g, len(genes), 0, np.nan, np.nan, np.nan, np.nan, np.nan)
    return (g, len(genes), n_tr, float(np.mean(lens)), float(np.mean(gcs)), float(np.mean(gc3s)), float(np.std(gc3s)), float(np.median(lens)))

if __name__ == "__main__":
    sheet, d, out = sys.argv[1:4]; thr = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    gs = [tag(r["SPECIES"], r["STRAIN"]) for r in csv.DictReader(open(sheet))]
    with Pool(thr) as pool, open(out, "w") as o:
        o.write("genome\tgenes\ttranscripts\tmean_cds_codons\tgc\tgc3_mean\tgc3_sd\tmedian_cds_codons\n")
        for r in pool.imap(one, [(g, d) for g in gs], chunksize=4):
            o.write("\t".join("" if (isinstance(x, float) and x != x) else (f"{x:.4f}" if isinstance(x, float) else str(x)) for x in r) + "\n")
