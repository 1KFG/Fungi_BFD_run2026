#!/usr/bin/python3.12
"""Figures for the paralog dS survey. Every figure is written as PNG (200 dpi) and PDF (vector) to figures/."""
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"; FIG = f"{D}/figures"; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 100, "savefig.bbox": "tight", "pdf.fonttype": 42})
X = pd.read_csv(f"{D}/data/genome_table_final.tsv.gz", sep="\t"); H = pd.read_csv(f"{D}/data/ds_hist.tsv.gz", sep="\t")
H = H.rename(columns={"n_zero_tree": "h_zero_tree", "n_identical": "h_identical", "n_no_ds": "h_no_ds", "n_below_floor": "h_below_floor"})
X = X.merge(H, on="genome", how="left"); BCOLS = [f"b{k:02d}" for k in range(70)]
EDGES = np.exp(np.linspace(np.log(0.005), np.log(2.5), 71)); CEN = np.sqrt(EDGES[:-1] * EDGES[1:])
PH = {"Ascomycota": "#0072B2", "Basidiomycota": "#D55E00", "Mucoromycota": "#009E73", "Glomeromycota": "#CC79A7", "Mortierellomycota": "#56B4E9", "Kickxellomycota": "#E69F00",
      "Chytridiomycota": "#999999", "Microsporidia": "#000000"}
def pcol(p): return PH.get(p, "#BBBBBB")
def lab(x): return str(x).replace("|unassigned", " (unassigned)")
def save(fig, name):
    fig.savefig(f"{FIG}/{name}.png", dpi=200); fig.savefig(f"{FIG}/{name}.pdf"); plt.close(fig); print("wrote", name, flush=True)
def hist_norm(g):
    """per genome: log-dS histogram as a share of ALL nodes (tree nodes + identical groups), so the area is the share with a usable dS"""
    tot = (g.h_zero_tree + g.h_identical + g.h_no_ds + g[BCOLS].sum(axis=1)).replace(0, np.nan)
    return g[BCOLS].div(tot, axis=0).values
N_ALL = X.h_zero_tree + X.h_identical + X.h_no_ds + X[BCOLS].sum(axis=1)
LEG = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=7, label=p) for p, c in PH.items()] + [Line2D([0], [0], marker="o", color="w", markerfacecolor="#BBBBBB", markersize=7, label="other")]

# ---- Fig 1: ridgelines by class
CC = pd.read_csv(f"{D}/tables/clade_summary_class.tsv", sep="\t"); cl = CC[(CC.n_genomes >= 20) & (CC.eligible)].copy()
ph_order = {p: i for i, p in enumerate(["Ascomycota", "Basidiomycota", "Mucoromycota", "Mortierellomycota", "Glomeromycota", "Kickxellomycota", "Chytridiomycota", "Microsporidia"])}
cl["phylum"] = [X[X["class"] == c].phylum.iloc[0] for c in cl["class"]]; cl["po"] = cl.phylum.map(ph_order).fillna(99)
cl = cl.sort_values(["po", "med_shareAll_MID"], ascending=[True, False]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8.2, 0.5 * len(cl) + 1.5)); H_all = {}
for i, r in cl.iterrows():
    g = X[X["class"] == r["class"]]; M = np.nanmean(hist_norm(g), axis=0); y0 = (len(cl) - 1 - i)
    s = M / np.nanmax(M) * 0.9 if np.nanmax(M) > 0 else M
    ax.fill_between(np.log(CEN), y0, y0 + s, color=pcol(r.phylum), alpha=0.75, lw=0.6, edgecolor="black"); ax.plot([np.log(0.005), np.log(2.5)], [y0, y0], color="black", lw=0.3)
    zs = np.nanmean((g.h_zero_tree + g.h_identical) / N_ALL[g.index]); ax.barh(y0 + 0.02, zs * 6, left=np.log(0.0035), height=0.35, color="black")
    ax.text(np.log(2.7), y0 + 0.1, f"{lab(r['class'])}  (n={int(r.n_genomes)})", va="bottom", ha="left", fontsize=7)
ax.set_yticks([]); ax.set_ylim(-0.3, len(cl) + 1.2); ax.set_xlim(np.log(0.0025), np.log(2.5)); ax.set_xticks(np.log([0.01, 0.05, 0.1, 0.3, 0.5, 1, 1.5, 2])); ax.set_xticklabels(["0.01", "0.05", "0.1", "0.3", "0.5", "1", "1.5", "2"])
for v in (0.05, 0.3, 1.5): ax.axvline(np.log(v), color="grey", lw=0.5, ls=":")
ax.set_xlabel("node dS (yn00, log scale)"); ax.set_title("Node-dS distribution by class (genomes weighted equally; height normalised per class)\nblack bar left = share of identical / dS = 0 nodes (x6); area excludes nodes without a usable dS", fontsize=9)
ax.axvspan(np.log(1.5), np.log(2.5), color="grey", alpha=0.12, lw=0); ax.text(np.log(1.9), len(cl) + 1.1, "yn00 saturation zone\n(no position claim)", ha="center", va="top", fontsize=6, color="#444")
ax.legend(handles=LEG[:8], loc="upper left", frameon=False, fontsize=7, ncol=2); save(fig, "fig01_ridgeline_by_class")

# ---- Fig 2: identical copies
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2)); d = X[X.genes_in_identical_groups > 0]
ax[0].scatter(d.ident_genes_per_1000, d.largest_identical_group, s=7, c=[pcol(p) for p in d.phylum], alpha=0.6, lw=0)
f = d[d.flag_identical_excess]; ax[0].scatter(f.ident_genes_per_1000, f.largest_identical_group, s=26, facecolors="none", edgecolors="black", lw=0.7)
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("genes in identical groups per 1,000 genes"); ax[0].set_ylabel("largest identical group (genes)"); ax[0].set_title(f"identical-copy groups ({len(d):,} of {len(X):,} genomes have any;\ncircled = excess flag, {int(X.flag_identical_excess.sum())})", fontsize=9)
t = X[X.te_per_1000_genes.notna() & (X.ident_genes_per_1000 > 0)]; ax[1].scatter(t.te_per_1000_genes + 0.01, t.ident_genes_per_1000, s=7, c=[pcol(p) for p in t.phylum], alpha=0.6, lw=0)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("TE-domain proteins per 1,000 genes (+0.01)"); ax[1].set_ylabel("genes in identical groups per 1,000"); ax[1].set_title("identical copies vs TE-domain content", fontsize=9)
ax[1].legend(handles=LEG, frameon=False, fontsize=6, ncol=2); save(fig, "fig02_identical_copies")

# ---- Fig 3: heatmap order x bin mass
CO = pd.read_csv(f"{D}/tables/clade_summary_order.tsv", sep="\t"); co = CO[CO.eligible].copy()
co["phylum"] = [X[X["order"] == o].phylum.iloc[0] for o in co["order"]]; co["class"] = [X[X["order"] == o]["class"].iloc[0] for o in co["order"]]; co["po"] = co.phylum.map(ph_order).fillna(99)
co = co.sort_values(["po", "class", "order"]).reset_index(drop=True)
cols = [("med_shareAll_Z", "Z identical"), ("med_shareAll_H", "H 0-0.05"), ("med_shareAll_Y", "Y 0.05-0.3"), ("med_shareAll_MID", "MID 0.3-1.5"), ("med_shareAll_OLD", "OLD 1.5-2.5"), ("med_shareAll_BEYOND", "beyond window")]
M = co[[c for c, _ in cols]].values
fig, ax = plt.subplots(figsize=(7.8, 0.17 * len(co) + 1.5)); im = ax.imshow(np.sqrt(M), aspect="auto", cmap="viridis", vmin=0, vmax=np.sqrt(0.9))
ax.set_xticks(range(len(cols))); ax.set_xticklabels([n for _, n in cols], rotation=40, ha="right"); ax.set_yticks(range(len(co))); ax.set_yticklabels([f"{lab(o)}  ({int(n)}/{int(g)})" for o, n, g in zip(co["order"], co.n_genomes, co.n_genera)], fontsize=6)
for i in range(len(co)):
    for j in range(len(cols)): ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5, color="white" if np.sqrt(M[i, j]) < 0.55 else "black")
prev = None
for i, p in enumerate(co.phylum):
    if prev is not None and p != prev: ax.axhline(i - 0.5, color="white", lw=2)
    prev = p
ax.set_title("Median share of nodes per bin, orders with >= 10 genomes and >= 3 genera\n(order labels show genomes/genera; colour = sqrt of share; shares are of ALL nodes)", fontsize=9); fig.colorbar(im, ax=ax, fraction=0.03, label="sqrt(share)"); save(fig, "fig03_order_bin_heatmap")

# ---- Fig 4: variance components
if os.path.exists(f"{D}/tables/varcomp.tsv"):
    V = pd.read_csv(f"{D}/tables/varcomp.tsv", sep="\t"); lev = ["phylum", "class", "order", "family", "genus", "genome (residual)"]; colr = ["#08306b", "#2171b5", "#6baed6", "#c6dbef", "#fdae6b", "#d9d9d9"]
    resp = list(V.response.unique()); fig, ax = plt.subplots(figsize=(8.5, 4)); w = 0.38
    for k, adj in enumerate(["no covariates", "with covariates"]):
        bottom = np.zeros(len(resp))
        for l, c in zip(lev, colr):
            v = np.array([V[(V.response == r) & (V.adjustment == adj) & (V.level == l)].share.sum() for r in resp])
            ax.bar(np.arange(len(resp)) + (k - 0.5) * w, v, w, bottom=bottom, color=c, edgecolor="white", lw=0.4, label=l if k == 0 else None); bottom += v
    ax.set_xticks(np.arange(len(resp))); ax.set_xticklabels(resp, rotation=25, ha="right"); ax.set_ylabel("share of variance"); ax.set_title("Variance partition on the taxonomy (left bar: no covariates, right bar: with assembly covariates)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.28)); save(fig, "fig04_variance_components")

# ---- Fig 5: controls and examples (histograms)
def draw_hist(ax, g, title, color="#0072B2"):
    h = hist_norm(X[X.genome == g.genome].reset_index(drop=True))[0]; ax.fill_between(np.log(CEN), 0, h, color=color, alpha=0.8)
    ax.set_xticks(np.log([0.05, 0.3, 1.5])); ax.set_xticklabels(["0.05", "0.3", "1.5"], fontsize=6); ax.set_yticks([]); ax.tick_params(length=2)
    ax.set_title(title, fontsize=5.5)
C = pd.read_csv(f"{D}/tables/controls.tsv", sep="\t")
sel = []
for lab in C.control_group.unique():
    sub = C[C.control_group == lab].sort_values("n_nodes_all", ascending=False)
    sel += [(lab, r.genome) for r in sub.head(8 if "named" in lab else 5).itertuples()]
n = len(sel); nc = 6; nr = int(np.ceil(n / nc)); fig, axs = plt.subplots(nr, nc, figsize=(13, 1.9 * nr)); axs = np.atleast_2d(axs)
colmap = {"named hybrid (' x ' in name)": "#D55E00", "literature hybrid species (strain may be collapsed)": "#E69F00", "post-WGD Saccharomycetaceae": "#0072B2", "pre-WGD Saccharomycetaceae": "#56B4E9", "Mucorales WGD": "#009E73"}
for k in range(nr * nc):
    a = axs[k // nc, k % nc]
    if k < n:
        lab, gname = sel[k]; r = X[X.genome == gname].iloc[0]
        draw_hist(a, r, f"{(gname if len(gname)<=34 else gname[:16]+".."+gname[-16:])}\n2-gene {r.share_2gene_fam:.2f}  dup {r.busco_duplicated_pct:.0f}%  x{r.genes_ratio_clade:.1f}", colmap.get(lab, "#444"))
    else: a.axis("off")
fig.suptitle("Control genomes: node-dS shape (x: log dS, y: share of all nodes). Panel text: share of 2-gene families, BUSCO duplicated %, gene count / clade median.\norange = named and literature hybrids; blue = Saccharomycetaceae; green = Mucorales", fontsize=8, y=1.0)
fig.tight_layout(h_pad=1.6, w_pad=1.0); save(fig, "fig05_controls")

# ---- Fig 6: peaks
PC = pd.read_csv(f"{D}/tables/peak_catalogue.tsv", sep="\t")
if len(PC):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
    ax[0].hist(np.log(PC.peak_dS), bins=40, color="#444"); ax[0].set_xticks(np.log([0.05, 0.1, 0.3, 0.5, 1, 1.5])); ax[0].set_xticklabels(["0.05", "0.1", "0.3", "0.5", "1", "1.5"]); ax[0].set_xlabel("stable peak position (dS)"); ax[0].set_ylabel("peaks"); ax[0].set_title(f"{len(PC)} stable peaks in {PC.genome.nunique()} genomes", fontsize=9)
    sc = ax[1].scatter(PC.peak_dS, PC.gmm_weight, c=PC.share_2gene_fam, s=12, cmap="viridis", vmin=0.4, vmax=0.8); ax[1].set_xscale("log"); ax[1].set_xlabel("peak dS"); ax[1].set_ylabel("mixture component weight"); plt.colorbar(sc, ax=ax[1], label="share of 2-gene families")
    ax[1].set_title("peak mass and family structure", fontsize=9)
    h = PC[PC.control_hybrid_named]; ax[2].scatter(PC.nodes_per_1000_genes, PC.genes_ratio_clade, s=8, c="#BBBBBB"); ax[2].scatter(h.nodes_per_1000_genes, h.genes_ratio_clade, s=30, c="#D55E00", label="named hybrids")
    ax[2].set_xlabel("nodes per 1,000 genes"); ax[2].set_ylabel("gene count / clade median"); ax[2].legend(frameon=False, fontsize=7); ax[2].set_title("duplication density vs gene count", fontsize=9); fig.tight_layout(); save(fig, "fig06_peaks")

# ---- Fig 7: young mass vs BUSCO duplication and N50
X["young"] = X.shareAll_H + X.shareAll_Y
fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
ax[0].scatter(X.busco_duplicated_pct + 0.2, X.shareAll_H, s=5, c=[pcol(p) for p in X.phylum], alpha=0.5, lw=0); ax[0].set_xscale("log"); ax[0].set_xlabel("BUSCO duplicated % (+0.2)"); ax[0].set_ylabel("H share (dS 0-0.05)")
ax[1].scatter(X.busco_duplicated_pct + 0.2, X.shareAll_Y, s=5, c=[pcol(p) for p in X.phylum], alpha=0.5, lw=0); ax[1].set_xscale("log"); ax[1].set_xlabel("BUSCO duplicated % (+0.2)"); ax[1].set_ylabel("Y share (dS 0.05-0.3)")
ax[2].scatter(X.N50_bp, X.young, s=5, c=[pcol(p) for p in X.phylum], alpha=0.5, lw=0); ax[2].set_xscale("log"); ax[2].set_xlabel("assembly N50 (bp)"); ax[2].set_ylabel("H + Y share")
ax[2].set_ylim(0, 1.15); ax[2].legend(handles=LEG, frameon=False, fontsize=6, ncol=4, loc="upper center"); fig.suptitle("Young-node mass against assembly duplication (BUSCO) and contiguity (N50)", fontsize=9); fig.tight_layout(); save(fig, "fig07_young_vs_assembly")

# ---- Fig 8: nodes per 1000 genes vs gene count
fig, ax = plt.subplots(figsize=(6.5, 4.6)); ax.scatter(X.genes, X.nodes_per_1000_genes, s=5, c=[pcol(p) for p in X.phylum], alpha=0.45, lw=0); ax.set_xscale("log")
for lab, col, mk in [("post-WGD Saccharomycetaceae", "#0072B2", "^"), ("Mucorales WGD", "#009E73", "s"), ("named hybrid (' x ' in name)", "#D55E00", "D")]:
    c = C[C.control_group == lab]; ax.scatter(c.genes, c.nodes_per_1000_genes, s=28, marker=mk, facecolors="none", edgecolors=col, lw=1.0, label=lab)
ax.set_xlabel("genes"); ax.set_ylabel("nodes (copy events) per 1,000 genes"); ax.set_title("Duplication density against gene count", fontsize=9)
l1 = ax.legend(handles=LEG[:6], frameon=False, fontsize=6, loc="upper left"); ax.add_artist(l1); ax.legend(frameon=False, fontsize=6, loc="upper right"); save(fig, "fig08_nodes_density")

# ---- Fig 9: sensitivity
S = X[(X.n_nodes_usable >= 100) & (X.s100_n_nodes >= 50)]
fig, ax = plt.subplots(1, 4, figsize=(12, 3.2))
for a, nm in zip(ax[:3], ["H", "Y", "MID"]):
    den = S[[f"s100_n_{k}" for k in ["Z", "H", "Y", "MID", "OLD"]]].sum(axis=1)
    a.scatter(S["n_" + nm] / S.n_nodes_usable, S["s100_n_" + nm] / den, s=4, alpha=0.4, c="#444", lw=0); a.plot([0, 1], [0, 1], color="red", lw=0.7); a.set_xlabel("share, >= 30 codons"); a.set_ylabel("share, >= 100 codons"); a.set_title(f"bin {nm}", fontsize=9)
PKall = pd.read_csv(f"{D}/data/peaks2.tsv.gz", sep="\t"); w = PKall[PKall["in_window_0.05_1.5"]]
crit = {"in window": len(w), "node bootstrap >= 0.8": int((w.support_node_boot >= 0.8).sum()), "+ stable at >= 100 codons": int(((w.support_node_boot >= 0.8) & (w.stable_codon100 == True)).sum()),
        "+ stable under pair resampling": int(((w.support_node_boot >= 0.8) & (w.stable_codon100 == True) & (w.stable_resample >= 0.8)).sum())}
ax[3].barh(list(crit.keys())[::-1], list(crit.values())[::-1], color="#555"); ax[3].set_xlabel("peaks"); ax[3].set_title("peak stability filters", fontsize=9)
fig.tight_layout(); save(fig, "fig09_sensitivity")

# ---- Fig 10: yn00 vs codeml
if os.path.exists(f"{D}/tables/yn00_vs_codeml_nodes_binned.tsv"):
    B = pd.read_csv(f"{D}/tables/yn00_vs_codeml_nodes_binned.tsv", sep="\t"); Q = pd.read_csv(f"{D}/tables/yn00_codeml_qq.tsv", sep="\t"); G = pd.read_csv(f"{D}/tables/yn00_vs_codeml_genomes.tsv", sep="\t")
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6)); x = np.arange(len(B))
    ax[0].errorbar(x, B.median_ratio_codeml_over_yn00, yerr=[B.median_ratio_codeml_over_yn00 - B.ratio_p10, B.ratio_p90 - B.median_ratio_codeml_over_yn00], fmt="o", color="#333", capsize=2, ms=3)
    ax[0].axhline(1, color="red", lw=0.6); ax[0].set_yscale("log"); ax[0].set_xticks(x); ax[0].set_xticklabels(B.yn00_dS_bin, rotation=60, ha="right", fontsize=6); ax[0].set_ylabel("codeml dS / yn00 dS (median, p10-p90)"); ax[0].set_xlabel("yn00 node dS bin"); ax[0].set_title("same nodes, two estimators", fontsize=9)
    ax[1].plot(Q.yn00, Q.codeml, color="#0072B2"); ax[1].plot([0, 2.5], [0, 2.5], color="red", lw=0.6); ax[1].set_xlabel("yn00 dS quantile value"); ax[1].set_ylabel("codeml dS at the same quantile"); ax[1].set_title("quantile map (genome-weighted)", fontsize=9)
    ax[2].scatter(G.yn00_median_pos, G.codeml_median_pos, s=8, c="#444"); ax[2].plot([0, 2], [0, 2], color="red", lw=0.6); ax[2].set_xlabel("genome median node dS, yn00"); ax[2].set_ylabel("genome median node dS, codeml"); ax[2].set_title(f"{len(G)} genomes", fontsize=9)
    fig.tight_layout(); save(fig, "fig10_yn00_vs_codeml")

# ---- Fig 11: coverage by phylum
ph = X.groupby("phylum").filter(lambda g: len(g) >= 9); order = ph.groupby("phylum").size().sort_values(ascending=False).index.tolist()
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4)); ax[0].barh(order[::-1], [int((X.phylum == p).sum()) for p in order][::-1], color=[pcol(p) for p in order][::-1]); ax[0].set_xscale("log"); ax[0].set_xlabel("genomes")
ax[1].boxplot([X[X.phylum == p].frac_nodes_with_ds.dropna() for p in order[::-1]], vert=False, tick_labels=order[::-1], showfliers=False); ax[1].set_xlabel("share of nodes with a usable dS (others lie beyond the identity window)")
fig.tight_layout(); save(fig, "fig11_coverage_by_phylum")

# ---- Fig 12: example gallery
groups = [("identical-copy excess", X[X.flag_identical_excess].sort_values("z_ident_genes_per_1000", ascending=False), "#000000"), ("H excess", X[X.flag_H_excess].sort_values("z_shareAll_H", ascending=False), "#D55E00"),
          ("MID excess", X[X.flag_MID_excess].sort_values("z_shareAll_MID", ascending=False), "#009E73")]
fig, axs = plt.subplots(len(groups), 6, figsize=(11, 1.9 * len(groups)))
for i, (lab, g, col) in enumerate(groups):
    for j in range(6):
        a = axs[i, j]
        if j < len(g): r = g.iloc[j]; draw_hist(a, r, f"{r.genome[:28]}\nZ {r.shareAll_Z:.2f} H {r.shareAll_H:.2f} Y {r.shareAll_Y:.2f}", col)
        else: a.axis("off")
    axs[i, 0].set_ylabel(lab, fontsize=7)
fig.tight_layout(); save(fig, "fig12_example_gallery")
