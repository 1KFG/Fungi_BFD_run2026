#!/usr/bin/python3.12
"""Wave-1 genomes (300) run with both estimators: codeml (nf_neuro/tables, node sampling) and yn00 (tables_yn00, node sampling).
1. node-level paired comparison (the same family tree gives the same node ids): dS_codeml vs dS_yn00 per (genome, family, node);
2. genome-level comparison: median node dS and bin shares;
3. a monotone yn00 -> codeml mapping from genome-weighted pooled quantiles.
Outputs: tables/yn00_vs_codeml_nodes_binned.tsv, tables/yn00_vs_codeml_genomes.tsv, tables/yn00_to_codeml_map.tsv, data/yn00_vs_codeml_pairs.tsv.gz"""
import duckdb, numpy as np, pandas as pd, csv, re, sys
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; P = f"{R}/paralogoscope_pilot"; E = f"{P}/eval/explore"
CODEML = f"{P}/nf_neuro/tables/wgd.ks.nodes.parquet"; YN00 = f"{P}/nf_full/tables_yn00/wgd.ks.nodes.parquet"
tax = pd.read_csv(f"{E}/data/taxonomy.tsv", sep="\t")
w1 = pd.read_csv(f"{P}/nf_wave1/samples_wave1.csv", dtype=str, keep_default_na=False)
def clean(st):
    st = (st or "").strip(); st = re.sub(r"['\"]", "", st).split(";")[0].strip().replace(":", " ")
    st = re.sub(r"^\s*\*+", "", st); st = re.sub(r"\*+\s*$", "", st); st = re.sub(r"\s*\*+\s*", "-", st); st = re.sub(r"[&`$|<>()]+", " ", st).strip(); return st
def tag(sp, st):
    sp = re.sub(r"['\"]", "", (sp or "").strip()); j = "_".join(x for x in [sp, clean(st)] if x); return re.sub(r"[\s/#\[\]?{}]+", "_", j)
genomes = [tag(a, b) for a, b in zip(w1.SPECIES, w1.STRAIN)]
con = duckdb.connect(); con.execute("SET memory_limit='12GB'; SET threads=4;")
con.register("g", pd.DataFrame({"genome": genomes}))
q = lambda f: f"""SELECT genome, family, node, dS_mean FROM read_parquet('{f}') WHERE kind='tree' AND n_valid>0 AND dS_mean IS NOT NULL AND genome IN (SELECT genome FROM g)"""
C = con.execute(q(CODEML)).df().rename(columns={"dS_mean": "dS_codeml"}); Y = con.execute(q(YN00)).df().rename(columns={"dS_mean": "dS_yn00"})
M = C.merge(Y, on=["genome", "family", "node"], how="inner")
print("codeml nodes with dS", len(C), "| yn00", len(Y), "| in both", len(M), "| genomes", M.genome.nunique())
M.to_csv(f"{E}/data/yn00_vs_codeml_pairs.tsv.gz", sep="\t", index=False)
# --- node-level, binned by the yn00 dS (what we have after the first pass)
bins = [0, 1e-9, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 1.6, 2.0, 3.0]
M["ybin"] = pd.cut(M.dS_yn00, bins, right=True, include_lowest=True)
rows = []
for b, g in M.groupby("ybin", observed=True):
    r = g.dS_codeml / g.dS_yn00.replace(0, np.nan)
    rows.append((str(b), len(g), g.dS_yn00.median(), g.dS_codeml.median(), r.median(), r.quantile(.1), r.quantile(.9), float(((r > .8) & (r < 1.25)).mean()), float(np.corrcoef(g.dS_yn00, g.dS_codeml)[0, 1]) if len(g) > 10 else np.nan))
B = pd.DataFrame(rows, columns=["yn00_dS_bin", "n_nodes", "median_yn00", "median_codeml", "median_ratio_codeml_over_yn00", "ratio_p10", "ratio_p90", "share_within_20pct", "pearson"])
B.to_csv(f"{E}/tables/yn00_vs_codeml_nodes_binned.tsv", sep="\t", index=False)
print(B.round(3).to_string(index=False))
# --- genome-level
def gstats(df, col):
    d = df[col].values; pos = d[(d > 0) & (d <= 2.0)]
    g = np.full(len(d), 4); g[d <= 1.5] = 3; g[d <= 0.3] = 2; g[d <= 0.05] = 1; g[d == 0] = 0
    out = {"median_pos": np.median(pos) if len(pos) else np.nan, "n": len(d)}
    for k, nm in enumerate(["Z", "H", "Y", "MID", "OLD"]): out["share_" + nm] = float((g == k).mean())
    return out
rows = []
for gname, g in M.groupby("genome"):
    c = gstats(g, "dS_codeml"); y = gstats(g, "dS_yn00")
    rows.append({"genome": gname, "n_nodes_both": len(g), **{"codeml_" + k: v for k, v in c.items()}, **{"yn00_" + k: v for k, v in y.items()}})
G = pd.DataFrame(rows); G["median_ratio"] = G.codeml_median_pos / G.yn00_median_pos
G = G.merge(tax[["genome", "PHYLUM", "CLASS"]], on="genome", how="left"); G.to_csv(f"{E}/tables/yn00_vs_codeml_genomes.tsv", sep="\t", index=False)
print("genomes", len(G), "| median codeml/yn00 genome-median ratio %.2f (IQR %.2f-%.2f)" % (G.median_ratio.median(), G.median_ratio.quantile(.25), G.median_ratio.quantile(.75)))
# --- monotone map yn00 -> codeml (genome-weighted quantiles of the paired nodes)
w = 1.0 / M.groupby("genome").dS_yn00.transform("size")
def wq(x, w, qs):
    o = np.argsort(x); x = x[o]; c = np.cumsum(w[o]) / w.sum(); return np.interp(qs, c, x)
qs = np.linspace(0.02, 0.98, 49)
mapq = pd.DataFrame({"quantile": qs, "yn00": wq(M.dS_yn00.values, w.values, qs), "codeml": wq(M.dS_codeml.values, w.values, qs)})
pts = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.5, 2.0]
mp = pd.DataFrame({"yn00_dS": pts, "codeml_dS_at_same_quantile": np.interp(pts, mapq.yn00, mapq.codeml)})
mp.to_csv(f"{E}/tables/yn00_to_codeml_map.tsv", sep="\t", index=False); mapq.to_csv(f"{E}/tables/yn00_codeml_qq.tsv", sep="\t", index=False); print(mp.round(3).to_string(index=False))
