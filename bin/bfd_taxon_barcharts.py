#!/usr/bin/env python3
"""Bar charts of BFD genome counts by taxon, all genomes vs one per species."""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--samples", default="samples.csv")
p.add_argument("--out", default="docs/BFD_genome_counts.pdf")
p.add_argument("--top", type=int, default=30)
p.add_argument("--png-dir", default=None)
a = p.parse_args()

d = pd.read_csv(a.samples, dtype=str)
for c in ["PHYLUM", "CLASS", "ORDER", "FAMILY"]:
    d[c] = d[c].fillna("Unassigned")
# collapse "Genus sp. STRAIN" (also "Chaetomiaceae sp. X") to one "Genus sp." entry
d["SPKEY"] = d["SPECIES"].str.replace(r"^(\S+) sp\..*$", r"\1 sp.", regex=True)
u = d.drop_duplicates("SPKEY")
N, NU = len(d), len(u)

ALL, UNI = "#4C78A8", "#E08E2B"
plt.rcParams.update({"font.size": 14, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})

def counts(col, phylum=None):
    dd, uu = (d, u) if phylum is None else (d[d.PHYLUM == phylum], u[u.PHYLUM == phylum])
    t = pd.DataFrame({"all": dd[col].value_counts(), "uniq": uu[col].value_counts()}).fillna(0).astype(int)
    return t.sort_values("all", ascending=False)

def paired(t, title, fname, logx=False, n=None, ratio=False, phylum=None):
    nn, nu = (N, NU) if phylum is None else ((d.PHYLUM == phylum).sum(), (u.PHYLUM == phylum).sum())
    t = t.head(n) if n else t
    t = t[::-1]
    h = max(6, 0.38 * len(t) + 2)
    fig, ax = plt.subplots(figsize=(13.33, min(h, 7.5)))
    y = range(len(t)); bh = 0.4
    b1 = ax.barh([i + bh / 2 for i in y], t["all"], bh, color=ALL, label=f"All genomes (n={nn:,})")
    b2 = ax.barh([i - bh / 2 for i in y], t["uniq"], bh, color=UNI, label=f"One per species (n={nu:,})")
    ax.set_yticks(list(y)); ax.set_yticklabels(t.index)
    if logx: ax.set_xscale("log")
    ax.bar_label(b1, fmt=lambda v: f"{int(v):,}", padding=3, fontsize=10)
    ax.bar_label(b2, fmt=lambda v: f"{int(v):,}", padding=3, fontsize=10)
    ax.set_xlabel("Number of genomes" + (" (log scale)" if logx else ""))
    ax.set_title(title, fontsize=18, fontweight="bold", loc="left")
    ax.legend(loc="lower right", frameon=False)
    ax.set_ylim(-0.7, len(t) - 0.3)
    fig.tight_layout()
    return fig

def ratio_chart(t, title, n):
    t = t.head(n).copy()
    t["r"] = t["all"] / t["uniq"].clip(lower=1)
    t = t.sort_values("r")
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    bars = ax.barh(t.index, t["r"], color="#59A14F")
    ax.bar_label(bars, labels=[f"{r:.1f}  ({a:,}/{b:,})" for r, a, b in zip(t["r"], t["all"], t["uniq"])],
                 padding=3, fontsize=10)
    ax.set_xlabel("Genomes per species (all genomes / unique species)")
    ax.set_title(title, fontsize=18, fontweight="bold", loc="left")
    ax.set_xlim(0, t["r"].max() * 1.25)
    fig.tight_layout()
    return fig

ph, cl, fa = counts("PHYLUM"), counts("CLASS"), counts("FAMILY")
fa_named = fa.drop(index="Unassigned", errors="ignore")
cl_named = cl.drop(index="Unassigned", errors="ignore")

figs = []
figs.append(paired(ph, f"Genomes per phylum: {N:,} genomes, {NU:,} species names", "", logx=True))
figs.append(paired(cl_named, f"Top {a.top} classes", "", n=a.top))
figs.append(paired(fa_named, f"Top {a.top} families by genome count", "", n=a.top))
fa_u = fa_named.sort_values("uniq", ascending=False)
figs.append(paired(fa_u, f"Top {a.top} families by number of species (one genome per species)", "", n=a.top))
od = counts("ORDER").drop(index="Unassigned", errors="ignore")
figs.insert(2, paired(od, f"Top {a.top} orders", "", n=a.top))
for ph_name in ["Ascomycota", "Basidiomycota"]:
    f_ph = counts("FAMILY", ph_name).drop(index="Unassigned", errors="ignore")
    figs.append(paired(f_ph, f"{ph_name}: top {a.top} families by genome count", "", n=a.top, phylum=ph_name))
figs.append(ratio_chart(fa_named, f"Genomes per species (redundancy), top {a.top} families", a.top))

with PdfPages(a.out) as pdf:
    for f in figs:
        note = (f"Source: samples.csv. Species = SPECIES column (exact name match); "
                f"'Xxx sp. STRAIN' names collapsed to one 'Xxx sp.'. "
                f"Missing family: {int((d.FAMILY=='Unassigned').sum()):,} genomes (excluded from family/class charts).")
        f.text(0.01, 0.005, note, fontsize=8, color="#666")
        pdf.savefig(f)
        if a.png_dir: f.savefig(f"{a.png_dir}/page{figs.index(f)+1}.png", dpi=60)
        plt.close(f)
print(ph.to_string()); print(fa_named.head(a.top).to_string())
