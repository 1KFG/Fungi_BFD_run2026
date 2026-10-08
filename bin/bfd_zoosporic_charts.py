#!/usr/bin/env python3
"""Taxon-count charts for Blastocladiomycota and Chytridiomycota: genomes, species, genera per family."""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--samples", default="samples.csv")
p.add_argument("--out", default="docs/BFD_Blastocladio_Chytridio_counts.pdf")
p.add_argument("--png-dir", default=None)
a = p.parse_args()

PHYLA = ["Chytridiomycota", "Blastocladiomycota", "Neocallimastigomycota"]
d = pd.read_csv(a.samples, dtype=str)
d = d[d.PHYLUM.isin(PHYLA)].copy()
# collapse "Xxx sp. STRAIN" to one "Xxx sp." entry, as in bfd_taxon_barcharts.py
d["SPKEY"] = d["SPECIES"].str.replace(r"^(\S+) sp\..*$", r"\1 sp.", regex=True)
d["FAMILY"] = d["FAMILY"].fillna("(no family assigned)")

COL = {"Genomes": "#4C78A8", "Species": "#E08E2B", "Genera": "#59A14F"}
plt.rcParams.update({"font.size": 14, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
NOTE = ("Source: samples.csv. Species = SPECIES column; 'Xxx sp. STRAIN' collapsed to 'Xxx sp.'. "
        "Genera = distinct GENUS values (blank GENUS not counted).")

def bars(ax, labels, vals_by_series, bh=0.27):
    n = len(vals_by_series)
    for k, (name, vals) in enumerate(vals_by_series.items()):
        pos = [i + (k - (n - 1) / 2) * bh for i in range(len(labels))]
        b = ax.barh(pos, vals, bh, color=COL[name], label=name)
        ax.bar_label(b, fmt=lambda v: f"{int(v)}", padding=3, fontsize=10)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    ax.set_ylim(-0.6, len(labels) - 0.4)

figs = []

# Page 1: taxa at each rank, one panel per phylum
ranks = [("CLASS", "Classes"), ("ORDER", "Orders"), ("FAMILY", "Families"),
         ("GENUS", "Genera"), ("SPKEY", "Species"), ("ASMID", "Genomes")]
fig, axes = plt.subplots(1, 3, figsize=(13.33, 7.5), sharex=True)
for ax, ph in zip(axes, PHYLA):
    s = d[d.PHYLUM == ph]
    fam = s[s.FAMILY != "(no family assigned)"].FAMILY.nunique()
    vals = [fam if c == "FAMILY" else s[c].nunique() for c, _ in ranks]
    labels = [l for _, l in ranks][::-1]
    b = ax.barh(labels, vals[::-1], color="#4C78A8")
    ax.bar_label(b, fmt=lambda v: f"{int(v)}", padding=3, fontsize=12)
    ax.set_title(ph, fontsize=14, fontweight="bold", loc="left")
    ax.set_xlabel("Number of taxa / genomes")
    ax.set_xlim(0, 100)
fig.suptitle("Taxonomic breadth of the zoosporic phyla in BFD",
             fontsize=18, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.02, 1, 0.95))
figs.append(fig)

# Pages 2-3: per family genomes / species / genera
for ph in PHYLA:
    s = d[d.PHYLUM == ph]
    one_family = s[s.FAMILY != "(no family assigned)"].FAMILY.nunique() <= 1
    key, rank = ("GENUS", "genus") if one_family else ("FAMILY", "family")
    if one_family:
        s = s.assign(GENUS=s.GENUS.fillna("(no genus assigned)"))
    t = s.groupby(key).agg(Genomes=("ASMID", "size"), Species=("SPKEY", "nunique"),
                           Genera=("GENUS", "nunique"))
    if one_family:
        t = t.drop(columns="Genera")
    t["noname"] = t.index.str.startswith("(no ")
    t = t.sort_values(["noname", "Genomes"], ascending=[True, False]).drop(columns="noname")
    t = t[::-1]
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    bars(ax, list(t.index), {c: t[c].tolist() for c in t.columns})
    ax.set_xlabel("Count")
    ax.set_title(f"{ph}: genomes, species{'' if one_family else ' and genera'} per {rank}\n"
                 f"({len(s)} genomes, {s.SPKEY.nunique()} species, {s.GENUS.nunique()} genera"
                 f"{', single family ' + s.FAMILY.iloc[0] if one_family else ''})",
                 fontsize=14, fontweight="bold", loc="left")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    figs.append(fig)

with PdfPages(a.out) as pdf:
    for i, f in enumerate(figs, 1):
        f.text(0.01, 0.005, NOTE, fontsize=8, color="#666")
        pdf.savefig(f)
        if a.png_dir: f.savefig(f"{a.png_dir}/z{i}.png", dpi=60)
        plt.close(f)
