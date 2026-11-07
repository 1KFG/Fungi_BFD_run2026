#!/usr/bin/python3.12
"""Per-genome features of the paralog dS distribution (version 2, follows the Fable review in fable_context.md).

Unit: a node = one copy event of a family tree. Tree nodes carry the mean yn00 dS of up to 5 sampled pairs; identical nodes are groups of genes with a
byte-identical CDS (dS = 0 by construction). Bins (do not split finer than 0.02 below dS 0.1):
  Z    identical groups + tree nodes with dS exactly 0
  H    0 < dS <= 0.05        (heterozygosity scale: haplotigs, dikaryons, near-identical copies)
  Y    0.05 < dS <= 0.3      (young duplicates, RIP products)
  MID  0.3 < dS <= 1.5       (candidate recent events)
  OLD  1.5 < dS <= 2.5       (saturated or unreliable for yn00; no position claim)
Two denominators: all nodes (tree nodes + identical groups) and usable nodes (nodes with a dS, plus identical groups).
Robustness: family-level cluster bootstrap CIs for the bin shares; the same bins with pairs of >= 100 codons only; a pair-resampling replicate (each pair kept
with p = 0.5); peaks of the log-dS histogram are called with a node bootstrap and re-called under both sensitivity settings.
Usage: build_features2.py nodes.parquet pairs.parquet summary.parquet outdir [threads]
"""
import sys, warnings
import numpy as np, pandas as pd, pyarrow.parquet as pq
from multiprocessing import Pool
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from sklearn.mixture import GaussianMixture
warnings.filterwarnings("ignore")

BINS = [("Z", None, None), ("H", 0, 0.05), ("Y", 0.05, 0.3), ("MID", 0.3, 1.5), ("OLD", 1.5, 2.5)]
FINE = [("lt0.02", 0, 0.02), ("0.02-0.05", 0.02, 0.05), ("0.05-0.1", 0.05, 0.1), ("0.1-0.3", 0.1, 0.3), ("0.3-0.7", 0.3, 0.7), ("0.7-1.5", 0.7, 1.5), ("1.5-2", 1.5, 2.0), ("2-2.5", 2.0, 2.5)]
FAM_BINS = [("f2", 2, 2), ("f3", 3, 3), ("f4_5", 4, 5), ("f6_10", 6, 10), ("f11_50", 11, 50), ("f51_200", 51, 200), ("fgt200", 201, 10 ** 9)]
GRID = np.linspace(np.log(0.005), np.log(2.0), 80)
EDGES = np.append(GRID, GRID[-1] + (GRID[1] - GRID[0]))
MIN_NODES_MIX = 200
BOOT = 200

def peaks_of(h):
    s = gaussian_filter1d(h.astype(float), 1.8)
    if s.max() <= 0: return np.array([], int), s
    pk, _ = find_peaks(s, prominence=0.12 * s.max())
    return pk, s

def group_of(ds):
    """bin index 1..4 for H, Y, MID, OLD (0 = zero)"""
    g = np.full(len(ds), 4)
    g[ds <= 1.5] = 3; g[ds <= 0.3] = 2; g[ds <= 0.05] = 1; g[ds == 0] = 0
    return g

def peak_set(vals):
    """peak positions (grid index) of the log-dS histogram of vals (0 < dS <= 2)"""
    v = vals[(vals > 0) & (vals <= 2.0)]
    if len(v) < 50: return set(), None, None
    h, _ = np.histogram(np.log(v), bins=EDGES); pk, s = peaks_of(h)
    return set(pk.tolist()), h, s

def one(args):
    i, npath, ppath = args
    t = pq.ParquetFile(npath).read_row_group(i, columns=["genome", "family", "kind", "node", "node_pairs", "n_valid", "dS_mean"]).to_pandas()
    p = pq.ParquetFile(ppath).read_row_group(i, columns=["genome", "family", "node", "status", "dS", "pair_codons"]).to_pandas()
    g = t.genome.iloc[0]; assert p.genome.iloc[0] == g, (i, g, p.genome.iloc[0])
    rng = np.random.default_rng(1000 + i)
    tree = t[t.kind == "tree"]; ident = t[t.kind == "identical"]
    out = {"genome": g, "n_tree_nodes": len(tree), "n_identical_groups": len(ident)}
    # ---- family structure from the node table
    fam_tree = tree.groupby("family").size()
    copies = ((1 + np.sqrt(1 + 8 * ident.node_pairs.astype(float))) / 2).round().astype(int)
    fam_extra = (copies - 1).groupby(ident.family.values).sum() if len(ident) else pd.Series(dtype=int)
    fams = pd.DataFrame({"tn": fam_tree}).join(fam_extra.rename("extra"), how="outer").fillna(0)
    fams["reps"] = np.where(fams.tn > 0, fams.tn + 1, 1); fams["genes"] = (fams.reps + fams.extra).astype(int)
    for nm, lo, hi in FAM_BINS: out[nm] = int(((fams.genes >= lo) & (fams.genes <= hi)).sum())
    out["n_families"] = len(fams); out["genes_in_families"] = int(fams.genes.sum()); out["max_family_genes"] = int(fams.genes.max()) if len(fams) else 0
    out["genes_in_big_families_gt50"] = int(fams.genes[fams.genes > 50].sum())
    out["genes_in_identical_groups"] = int(copies.sum()); out["largest_identical_group"] = int(copies.max()) if len(copies) else 0
    out["n_identical_groups_ge5"] = int((copies >= 5).sum())
    # ---- node dS (primary): tree nodes with a dS
    v = tree[(tree.n_valid > 0) & tree.dS_mean.notna()]
    ds = v.dS_mean.values; fam_v = v.family.values
    out["n_nodes_ds"] = len(ds); out["n_nodes_all"] = len(tree) + len(ident); out["n_nodes_usable"] = len(ds) + len(ident)
    # ---- bin counts by node (Z includes identical groups)
    zero_t = ds == 0
    gi = group_of(ds)
    cnt = {nm: 0 for nm, _, _ in BINS}
    cnt["Z"] = int(zero_t.sum()) + len(ident)
    for k, (nm, lo, hi) in enumerate(BINS[1:], start=1): cnt[nm] = int((gi == k).sum())
    for nm in cnt: out["n_" + nm] = cnt[nm]
    out["n_zero_tree"] = int(zero_t.sum())
    for nm, lo, hi in FINE: out["fine_" + nm] = int(((ds > lo) & (ds <= hi)).sum())
    pos = ds[(ds > 0) & (ds <= 2.0)]; out["n_pos_ds"] = len(pos)
    for q in (10, 25, 50, 75, 90): out[f"ds_q{q}"] = float(np.percentile(pos, q)) if len(pos) >= 5 else np.nan
    # ---- family-level cluster bootstrap of the bin shares (usable denominator)
    allfam = pd.Index(sorted(set(fam_v) | set(ident.family)))
    if len(allfam) >= 5 and (len(ds) + len(ident)) > 0:
        fi = {f: j for j, f in enumerate(allfam)}; C = np.zeros((len(allfam), 5))
        for f_, gg in zip(fam_v, gi): C[fi[f_], gg if gg > 0 else 0] += 1   # gg 0 -> Z
        for f_ in ident.family: C[fi[f_], 0] += 1
        # columns: Z,H,Y,MID,OLD
        idx = rng.integers(0, len(allfam), (BOOT, len(allfam))); S = C[idx].sum(1); den = S.sum(1, keepdims=True)
        sh = S / np.maximum(den, 1)
        for k, (nm, _, _) in enumerate(BINS):
            out[f"ci_{nm}_lo"] = float(np.percentile(sh[:, k], 2.5)); out[f"ci_{nm}_hi"] = float(np.percentile(sh[:, k], 97.5))
    # ---- sensitivity: pairs with >= 100 codons, and pair resampling
    pp = p[(p.status == "ok") & p.dS.notna()]
    nm100 = pp[pp.pair_codons >= 100].groupby(["family", "node"]).dS.mean()
    d100 = nm100.values; g100 = group_of(d100) if len(d100) else np.array([], int)
    for k, (nm, _, _) in enumerate(BINS[1:], start=1): out["s100_n_" + nm] = int((g100 == k).sum())
    out["s100_n_Z"] = int((d100 == 0).sum()) + len(ident); out["s100_n_nodes"] = len(d100)
    key = pp.family.astype(str) + "|" + pp.node.astype(str)
    comps, peaks = [], []
    # ---- mixture and peaks (window 0 < dS <= 2); nodes >= MIN_NODES_MIX
    if len(pos) >= MIN_NODES_MIX:
        x = np.log(pos).reshape(-1, 1); best = None
        for k in (1, 2, 3, 4):
            gm = GaussianMixture(k, n_init=3, random_state=1).fit(x); b = gm.bic(x)
            if best is None or b < best[0]: best = (b, k, gm)
        _, k, gm = best; o = np.argsort(gm.means_.ravel()); out["gmm_k"] = k
        for r, j in enumerate(o):
            mu = float(np.exp(gm.means_.ravel()[j])); comps.append({"genome": g, "k": k, "comp": r + 1, "mean_dS": mu, "sd_log": float(np.sqrt(gm.covariances_.ravel()[j])),
                                                              "weight": float(gm.weights_[j]), "in_window_0.05_1.5": bool(0.05 <= mu <= 1.5)})
        pk0, h0, s0 = peak_set(pos)
        # node bootstrap
        pr = h0 / h0.sum(); Bt = rng.multinomial(len(pos), pr, size=100); bpk = [set(peaks_of(b)[0].tolist()) for b in Bt]
        # >=100 codon peaks
        pk100 = peak_set(d100)[0] if len(d100) >= 100 else None
        # resample replicates (each pair kept with p = 0.5, at least one pair per node)
        rp = []
        if len(pp):
            grp = pp.groupby(key).indices
            vals = pp.dS.values
            for _ in range(20):
                keep = rng.random(len(pp)) < 0.5
                m = []
                for ix in grp.values():
                    kk = ix[keep[ix]] if keep[ix].any() else ix[:1]
                    m.append(vals[kk].mean())
                rp.append(peak_set(np.array(m))[0])
        for q in sorted(pk0):
            sup = float(np.mean([any(abs(q - z) <= 3 for z in bs) for bs in bpk]))
            st100 = (None if pk100 is None else any(abs(q - z) <= 3 for z in pk100))
            strs = float(np.mean([any(abs(q - z) <= 3 for z in rs) for rs in rp])) if rp else np.nan
            peaks.append({"genome": g, "peak_dS": float(np.exp(GRID[q])), "height_rel": float(s0[q] / s0.max()), "support_node_boot": sup,
                          "stable_codon100": st100, "stable_resample": strs, "in_window_0.05_1.5": bool(0.05 <= np.exp(GRID[q]) <= 1.5)})
        out["n_peaks"] = len(peaks)
        out["n_supported_peaks_window"] = sum(1 for z in peaks if z["in_window_0.05_1.5"] and z["support_node_boot"] >= 0.8 and z["stable_codon100"] is True and z["stable_resample"] >= 0.8)
    return out, comps, peaks

if __name__ == "__main__":
    nodes, pairs, summ, od = sys.argv[1:5]; thr = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    n = pq.ParquetFile(nodes).metadata.num_row_groups; assert n == pq.ParquetFile(pairs).metadata.num_row_groups
    rows, comps, peaks = [], [], []
    with Pool(thr) as pool:
        for r, c, p in pool.imap(one, [(i, nodes, pairs) for i in range(n)], chunksize=8):
            rows.append(r); comps += c; peaks += p
    F = pd.DataFrame(rows); S = pd.read_parquet(summ)[["genome", "species_prefix", "hybrid", "n_pairs_raw"]]
    F = F.merge(S, on="genome", how="left")
    F.to_csv(f"{od}/features2.tsv.gz", sep="\t", index=False)
    pd.DataFrame(comps).to_csv(f"{od}/gmm_components2.tsv.gz", sep="\t", index=False)
    pd.DataFrame(peaks).to_csv(f"{od}/peaks2.tsv.gz", sep="\t", index=False)
    print("genomes", len(F), "with mixture", int(F.gmm_k.notna().sum()), "peaks", len(peaks), "genomes with a supported in-window peak", int((F.n_supported_peaks_window.fillna(0) > 0).sum()))
