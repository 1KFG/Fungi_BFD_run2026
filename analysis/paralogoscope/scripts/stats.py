#!/usr/bin/python3.12
"""Statistics for the paralog dS survey (follows fable_context.md). Writes tables/*.tsv.
- clade summaries (phylum, class, order): medians per genome-level metric, genus-cluster bootstrap CIs for the key metrics, n genomes and n genera;
- class effect sizes against the rest of the phylum (ratio of medians with genus-cluster bootstrap CI, Cliff's delta);
- one-genome-per-genus resampling of class medians (100 draws);
- outlier flags from RAW values (no shrinkage) as robust z within the most specific clade with >= 10 genomes, with evidence columns;
- control genomes (named hybrids, literature hybrids, old-WGD clades);
- peak catalogue (peaks stable under node bootstrap, >= 100-codon pairs and pair resampling, inside the 0.05-1.5 window);
- sensitivity of the bin shares to the 100-codon gate."""
import re, sys, numpy as np, pandas as pd
D = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
X = pd.read_csv(f"{D}/data/genome_table.tsv.gz", sep="\t")
PK = pd.read_csv(f"{D}/data/peaks2.tsv.gz", sep="\t"); GM = pd.read_csv(f"{D}/data/gmm_components2.tsv.gz", sep="\t")
rng = np.random.default_rng(2026)
KEY = ["frac_nodes_with_ds", "shareAll_Z", "shareAll_H", "shareAll_Y", "shareAll_MID", "shareAll_OLD", "shareAll_BEYOND", "nodes_per_1000_genes", "paralog_gene_frac",
       "share_2gene_fam", "ident_genes_per_1000", "genes", "busco_duplicated_pct", "te_per_1000_genes"]
BOOTM = ["shareAll_Z", "shareAll_H", "shareAll_Y", "shareAll_MID", "nodes_per_1000_genes", "paralog_gene_frac", "share_2gene_fam", "ident_genes_per_1000"]
for c in KEY: assert c in X.columns, c

def boot_median(groups, B=500):
    """cluster bootstrap by genus: groups = list of arrays (one per genus)"""
    k = len(groups); out = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, k, k); out[b] = np.nanmedian(np.concatenate([groups[i] for i in idx]))
    return np.percentile(out, [2.5, 97.5])

def clade_table(level, minn=10, ming=3):
    rows = []
    for name, g in X.groupby(level):
        r = {level: name, "n_genomes": len(g), "n_genera": g.genus.nunique(), "eligible": len(g) >= minn and g.genus.nunique() >= ming}
        for c in KEY: r["med_" + c] = g[c].median()
        if r["eligible"]:
            gen = [grp[c].values for _, grp in g.groupby("genus") for c in [None]] if False else None
            for c in BOOTM:
                groups = [grp[c].dropna().values for _, grp in g.groupby("genus")]; groups = [a for a in groups if len(a)]
                lo, hi = boot_median(groups, 300); r[c + "_ci_lo"], r[c + "_ci_hi"] = lo, hi
        rows.append(r)
    T = pd.DataFrame(rows).sort_values("n_genomes", ascending=False); T.to_csv(f"{D}/tables/clade_summary_{level}.tsv", sep="\t", index=False); return T

def cliffs_delta(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0: return np.nan
    ranks = pd.Series(np.concatenate([a, b])).rank().values; ra = ranks[:len(a)].sum()
    u = ra - len(a) * (len(a) + 1) / 2; return 2 * u / (len(a) * len(b)) - 1

print("clade summaries ...", flush=True)
CP, CC, CO = clade_table("phylum"), clade_table("class"), clade_table("order")
print(" phyla", len(CP), "| classes", len(CC), "(eligible", int(CC.eligible.sum()), ") | orders", len(CO), "(eligible", int(CO.eligible.sum()), ")")

# ---- class effects vs the rest of the phylum
rows = []
for cl in CC[CC.eligible]["class"]:
    g = X[X["class"] == cl]; ph = g.phylum.iloc[0]; rest = X[(X.phylum == ph) & (X["class"] != cl)]
    if len(rest) < 10: continue
    for c in BOOTM:
        a = g[c].dropna().values; b = rest[c].dropna().values
        if len(a) < 5 or len(b) < 5: continue
        groups = [grp[c].dropna().values for _, grp in g.groupby("genus")]; groups = [x for x in groups if len(x)]
        lo, hi = boot_median(groups, 300); mb = np.median(b)
        rows.append({"class": cl, "phylum": ph, "metric": c, "n_class": len(a), "n_genera": g.genus.nunique(), "n_rest": len(b), "median_class": np.median(a), "median_rest": mb,
                     "diff": np.median(a) - mb, "ratio": (np.median(a) / mb) if mb > 0 else np.nan, "diff_ci_lo": lo - mb, "diff_ci_hi": hi - mb, "cliffs_delta": cliffs_delta(a, b)})
EF = pd.DataFrame(rows); EF.to_csv(f"{D}/tables/class_effects_vs_phylum.tsv", sep="\t", index=False)

# ---- one genome per genus, 100 draws, class medians
rows = []
for cl in CC[CC.eligible]["class"]:
    g = X[X["class"] == cl]
    for c in BOOTM:
        meds = []
        for _ in range(100):
            s = g.groupby("genus").sample(1, random_state=int(rng.integers(1e9)))[c]; meds.append(s.median())
        rows.append({"class": cl, "metric": c, "n_genera": g.genus.nunique(), "median_all_genomes": g[c].median(), "median_1_per_genus_mean": np.mean(meds), "median_1_per_genus_sd": np.std(meds)})
pd.DataFrame(rows).to_csv(f"{D}/tables/class_medians_one_per_genus.tsv", sep="\t", index=False)

# ---- outlier flags (raw values, robust z in the most specific clade with >= 10 genomes)
def ref_group(r):
    for lvl in ["order", "class", "phylum"]:
        if (X[lvl] == r[lvl]).sum() >= 10: return lvl, r[lvl]
    return "all", "all"
def transform(col):
    if col in ("ident_genes_per_1000", "te_per_1000_genes"): return np.log1p(X[col])
    if col == "nodes_per_1000_genes": return np.log(X[col])
    s = X[col].clip(1e-4, 1 - 1e-4); return np.log(s / (1 - s))
Z = {}
refs = X.apply(ref_group, axis=1, result_type="expand"); X["ref_level"], X["ref_name"] = refs[0], refs[1]
FLAG_COLS = {"ident_genes_per_1000": 0.10, "shareAll_H": 0.30, "shareAll_Y": 0.30, "shareAll_MID": 0.30, "nodes_per_1000_genes": 0.05, "share_2gene_fam": 0.10}
for c, floor in FLAG_COLS.items():
    t = transform(c); z = np.full(len(X), np.nan)
    for (lvl, nm), idx in X.groupby(["ref_level", "ref_name"]).groups.items():
        v = t.loc[idx].values; med = np.nanmedian(v); mad = 1.4826 * np.nanmedian(np.abs(v - med)); mad = max(mad, floor)   # floor: MAD is 0 when most genomes have none
        z[X.index.get_indexer(idx)] = (v - med) / mad
    X["z_" + c] = z
X["flag_identical_excess"] = (X.z_ident_genes_per_1000 >= 3.5) & (X.ident_genes_per_1000 >= 5)
X["flag_H_excess"] = (X.z_shareAll_H >= 3.5) & (X.n_nodes_usable >= 100) & (X.shareAll_H >= 0.05)
X["flag_Y_excess"] = (X.z_shareAll_Y >= 3.5) & (X.n_nodes_usable >= 100) & (X.shareAll_Y >= 0.15)
X["flag_MID_excess"] = (X.z_shareAll_MID >= 3.5) & (X.n_nodes_usable >= 100) & (X.shareAll_MID >= 0.25)
X["flag_dup_dense"] = (X.z_nodes_per_1000_genes >= 3.5) & (X.nodes_per_1000_genes >= 450)
X["flag_2gene_excess"] = (X.z_share_2gene_fam >= 3.5) & (X.share_2gene_fam >= 0.65)
FL = [c for c in X.columns if c.startswith("flag_")]
keep = ["genome", "species_prefix", "phylum", "class", "order", "family_tax", "genus", "ref_level", "ref_name", "hybrid", "genes", "genes_ratio_clade", "n_nodes_all", "n_nodes_usable", "frac_nodes_with_ds",
        "shareAll_Z", "shareAll_H", "shareAll_Y", "shareAll_MID", "shareAll_OLD", "shareAll_BEYOND", "ident_genes_per_1000", "genes_in_identical_groups", "largest_identical_group", "n_identical_groups_ge5",
        "nodes_per_1000_genes", "share_2gene_fam", "paralog_gene_frac", "max_family_genes", "busco_duplicated_pct", "contig_count", "N50_bp", "total_length_bp", "te_per_1000_genes",
        "ci_Z_lo", "ci_Z_hi", "ci_H_lo", "ci_H_hi", "ci_Y_lo", "ci_Y_hi", "ci_MID_lo", "ci_MID_hi"] + [c for c in X.columns if c.startswith("z_")] + FL
X.loc[X[FL].any(axis=1), keep].sort_values("z_ident_genes_per_1000", ascending=False).to_csv(f"{D}/tables/outlier_flags.tsv", sep="\t", index=False)
X[keep].to_csv(f"{D}/data/genome_table_flags.tsv.gz", sep="\t", index=False)
print("flags:", {c: int(X[c].sum()) for c in FL})

# ---- controls
LIT = r"pastorianus|parabailii|orthopsilosis|metapsilosis|Aspergillus latus|longisporum|sorbitophila|Millerozyma|Zygosaccharomyces rouxii|Cryptococcus.*(deneoformans|AD hybrid)|Candida inconspicua|Pichia sorbitophila"
OLDWGD = {"post-WGD Saccharomycetaceae": r"^(Saccharomyces|Nakaseomyces|Kazachstania|Naumovozyma|Tetrapisispora|Vanderwaltozyma|Candida glabrata|Nakaseomyces glabratus)",
          "pre-WGD Saccharomycetaceae": r"^(Kluyveromyces|Lachancea|Eremothecium|Zygosaccharomyces|Torulaspora)", "Mucorales WGD": r"^(Rhizopus|Phycomyces|Mucor )"}
ctrl = []
X["SPECIES_NAME"] = X.genome.str.replace("_", " ")
def add(mask, label):
    g = X[mask].copy(); g["control_group"] = label; ctrl.append(g)
add(X.hybrid == 1, "named hybrid (' x ' in name)")
add(X.SPECIES_NAME.str.contains(LIT, regex=True, case=False) & (X.hybrid == 0), "literature hybrid species (strain may be collapsed)")
for k, v in OLDWGD.items(): add(X.SPECIES_NAME.str.contains(v, regex=True) & (X.hybrid == 0), k)
C = pd.concat(ctrl) if ctrl else pd.DataFrame()
ccols = ["control_group", "genome", "genes", "n_nodes_all", "frac_nodes_with_ds", "shareAll_Z", "shareAll_H", "shareAll_Y", "shareAll_MID", "shareAll_OLD", "nodes_per_1000_genes", "share_2gene_fam", "busco_duplicated_pct", "genes_ratio_clade", "n_peaks", "n_supported_peaks_window", "N50_bp"]
C[ccols].to_csv(f"{D}/tables/controls.tsv", sep="\t", index=False); print("controls:", C.groupby("control_group").size().to_dict())

# ---- peak catalogue
pk = PK[PK["in_window_0.05_1.5"] & (PK.support_node_boot >= 0.8) & (PK.stable_codon100 == True) & (PK.stable_resample >= 0.8)].copy()
comp = GM[GM["in_window_0.05_1.5"]].copy()
def nearest(r):
    c = comp[comp.genome == r.genome]
    if c.empty: return pd.Series({"gmm_weight": np.nan, "gmm_sd_log": np.nan, "gmm_mean_dS": np.nan})
    j = (np.log(c.mean_dS) - np.log(r.peak_dS)).abs().idxmin(); return pd.Series({"gmm_weight": c.loc[j, "weight"], "gmm_sd_log": c.loc[j, "sd_log"], "gmm_mean_dS": c.loc[j, "mean_dS"]})
pk = pd.concat([pk, pk.apply(nearest, axis=1)], axis=1)
PC = pk.merge(X[["genome", "phylum", "class", "order", "hybrid", "genes", "genes_ratio_clade", "n_nodes_usable", "frac_nodes_with_ds", "share_2gene_fam", "nodes_per_1000_genes", "busco_duplicated_pct", "N50_bp",
                 "contig_count", "te_per_1000_genes", "shareAll_H", "shareAll_Y", "shareAll_MID", "z_share_2gene_fam", "z_nodes_per_1000_genes"]], on="genome", how="left")
PC["control_hybrid_named"] = PC.hybrid == 1
PC.sort_values(["peak_dS"]).to_csv(f"{D}/tables/peak_catalogue.tsv", sep="\t", index=False)
print("peak catalogue:", len(PC), "peaks in", PC.genome.nunique(), "genomes")

# ---- sensitivity to the 100-codon gate
S = X[(X.n_nodes_usable >= 100) & (X.s100_n_nodes >= 50)].copy()
rows = []
for nm in ["Z", "H", "Y", "MID", "OLD"]:
    p30 = S["n_" + nm] / S.n_nodes_usable; den100 = S[[f"s100_n_{k}" for k in ["Z", "H", "Y", "MID", "OLD"]]].sum(axis=1); p100 = S["s100_n_" + nm] / den100
    rows.append({"bin": nm, "genomes": len(S), "median_share_30codons": p30.median(), "median_share_100codons": p100.median(), "spearman": p30.corr(p100, method="spearman"), "median_abs_diff": (p30 - p100).abs().median()})
pd.DataFrame(rows).to_csv(f"{D}/tables/sensitivity_codon100.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).round(3).to_string(index=False))
X.drop(columns=["SPECIES_NAME"]).to_csv(f"{D}/data/genome_table_final.tsv.gz", sep="\t", index=False)
