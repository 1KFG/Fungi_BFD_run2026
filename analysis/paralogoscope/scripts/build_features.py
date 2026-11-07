#!/usr/bin/python3.12
"""Per-genome features of the paralog dS distribution from the paralogoscope node table (one parquet row group per genome).

Unit of analysis: a node = one copy event in a family tree (tree nodes carry the node-mean yn00 dS of up to 5 sampled pairs; identical nodes are
groups of genes with a byte-identical CDS, dS = 0 by construction). dS window: 0 < dS <= 2 (yn00 saturates near 2 and the gate is 2.5).

Outputs (data/): features.tsv.gz (one row per genome), gmm_components.tsv.gz, peaks.tsv.gz, famsize.tsv.gz (family-size histogram per genome).
Usage: build_features.py nodes.parquet summary.parquet outdir [threads]
"""
import sys, json, warnings
import numpy as np, pandas as pd, pyarrow.parquet as pq
from multiprocessing import Pool
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from sklearn.mixture import GaussianMixture
warnings.filterwarnings("ignore")

DS_BINS = [("zero", 0, 0), ("lt0.01", 0, 0.01), ("0.01-0.05", 0.01, 0.05), ("0.05-0.2", 0.05, 0.2), ("0.2-0.5", 0.2, 0.5),
           ("0.5-1", 0.5, 1.0), ("1-1.5", 1.0, 1.5), ("1.5-2", 1.5, 2.0), ("gt2", 2.0, 2.5)]
FAM_BINS = [("f2", 2, 2), ("f3", 3, 3), ("f4_5", 4, 5), ("f6_10", 6, 10), ("f11_50", 11, 50), ("f51_200", 51, 200), ("fgt200", 201, 10 ** 9)]
GRID = np.linspace(np.log(0.005), np.log(2.0), 80)   # log-dS grid for the KDE-like histogram
MIN_NODES_MIX = 200
RNG = np.random.default_rng(12345)

def peaks_of(h):
    s = gaussian_filter1d(h.astype(float), 1.8)
    if s.max() <= 0: return np.array([], int), s
    pk, pr = find_peaks(s, prominence=0.12 * s.max())
    return pk, s

def genome_features(args):
    i, path = args
    t = pq.ParquetFile(path).read_row_group(i, columns=["genome", "family", "kind", "node_pairs", "n_valid", "dS_mean"]).to_pandas()
    g = t.genome.iloc[0]
    tree = t[t.kind == "tree"]; ident = t[t.kind == "identical"]
    out = {"genome": g, "n_tree_nodes": len(tree), "n_identical_groups": len(ident)}
    # --- family size from the node table: reps = tree nodes + 1 (or the identical group of a family without tree nodes)
    fam_tree = tree.groupby("family").size()
    copies = ((1 + np.sqrt(1 + 8 * ident.node_pairs.astype(float))) / 2).round().astype(int)
    fam_extra = (copies - 1).groupby(ident.family.values).sum() if len(ident) else pd.Series(dtype=int)
    fams = pd.DataFrame({"tn": fam_tree}).join(fam_extra.rename("extra"), how="outer").fillna(0)
    ident_only = set(ident.family) - set(fam_tree.index)
    fams["reps"] = np.where(fams.tn > 0, fams.tn + 1, 1)
    fams["genes"] = (fams.reps + fams.extra).astype(int)
    for nm, lo, hi in FAM_BINS: out[nm] = int(((fams.genes >= lo) & (fams.genes <= hi)).sum())
    out["n_families"] = len(fams); out["genes_in_families"] = int(fams.genes.sum()); out["max_family_genes"] = int(fams.genes.max()) if len(fams) else 0
    out["genes_in_big_families_gt50"] = int(fams.genes[fams.genes > 50].sum())
    # --- node-level dS distribution (window 0 < dS <= 2.5 for binning, 0 < dS <= 2 for the mixture)
    v = tree[(tree.n_valid > 0) & tree.dS_mean.notna()]
    ds = v.dS_mean.values
    out["n_nodes_ds"] = len(ds)
    zero_tree = int((ds == 0).sum())
    out["zero_tree_nodes"] = zero_tree
    out["zero_nodes_total"] = zero_tree + len(ident)          # identical groups + tree nodes whose sampled dS is exactly 0
    out["nodes_for_share"] = len(ds) + len(ident)              # denominator of the dS == 0 share
    for nm, lo, hi in DS_BINS:
        if nm == "zero": out["n_zero"] = zero_tree + len(ident)
        else: out["n_" + nm] = int(((ds > lo) & (ds <= hi)).sum()) if nm != "lt0.01" else int(((ds > 0) & (ds < 0.01)).sum())
    pos = ds[(ds > 0) & (ds <= 2.0)]
    out["n_pos_ds"] = len(pos)
    for q in (10, 25, 50, 75, 90): out[f"ds_q{q}"] = float(np.percentile(pos, q)) if len(pos) >= 5 else np.nan
    out["ds_mean_log"] = float(np.mean(np.log(pos))) if len(pos) else np.nan
    out["ds_sd_log"] = float(np.std(np.log(pos))) if len(pos) > 2 else np.nan
    comps, peaks = [], []
    if len(pos) >= MIN_NODES_MIX:
        x = np.log(pos).reshape(-1, 1); best = None
        for k in (1, 2, 3, 4):
            gm = GaussianMixture(k, n_init=3, random_state=1).fit(x); b = gm.bic(x)
            if best is None or b < best[0]: best = (b, k, gm)
        _, k, gm = best; o = np.argsort(gm.means_.ravel())
        out["gmm_k"] = k
        for r, j in enumerate(o):
            comps.append({"genome": g, "k": k, "comp": r + 1, "mean_dS": float(np.exp(gm.means_.ravel()[j])), "sd_log": float(np.sqrt(gm.covariances_.ravel()[j])), "weight": float(gm.weights_[j])})
        # --- peak support by bootstrap of the log-dS histogram
        h, _ = np.histogram(np.log(pos), bins=np.append(GRID, GRID[-1] + (GRID[1] - GRID[0])))
        pk, s = peaks_of(h)
        p = h / h.sum(); B = RNG.multinomial(len(pos), p, size=100)
        bpk = [set(peaks_of(b)[0]) for b in B]
        for q in pk:
            sup = float(np.mean([any(abs(q - z) <= 3 for z in bs) for bs in bpk]))
            peaks.append({"genome": g, "peak_dS": float(np.exp(GRID[q])), "height_rel": float(s[q] / s.max()), "support": sup})
        out["n_peaks"] = len(peaks); out["n_supported_peaks"] = sum(1 for z in peaks if z["support"] >= 0.8)
    return out, comps, peaks, fams.assign(genome=g)

if __name__ == "__main__":
    nodes, summ, od = sys.argv[1:4]; thr = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    n = pq.ParquetFile(nodes).metadata.num_row_groups
    rows, comps, peaks, fam = [], [], [], []
    with Pool(thr) as pool:
        for r, c, p, f in pool.imap(genome_features, [(i, nodes) for i in range(n)], chunksize=8):
            rows.append(r); comps += c; peaks += p
    F = pd.DataFrame(rows); S = pd.read_parquet(summ)
    F = F.merge(S[["genome", "species_prefix", "n_pairs_raw", "hybrid", "identical_pair_share", "n_genes_in_identical_groups", "largest_identical_group"]], on="genome", how="left")
    F.to_csv(f"{od}/features.tsv.gz", sep="\t", index=False)
    pd.DataFrame(comps).to_csv(f"{od}/gmm_components.tsv.gz", sep="\t", index=False)
    pd.DataFrame(peaks).to_csv(f"{od}/peaks.tsv.gz", sep="\t", index=False)
    print("genomes", len(F), "with mixture", int(F.gmm_k.notna().sum()), "components", len(comps), "peaks", len(peaks))
