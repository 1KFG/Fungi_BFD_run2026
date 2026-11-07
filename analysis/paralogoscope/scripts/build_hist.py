#!/usr/bin/python3.12
"""Per-genome histogram of node dS on a log grid (70 bins, 0.005-2.5; tree nodes with a dS > 0), plus counts of dS == 0 tree nodes, identical groups
and nodes without a usable dS. Output data/ds_hist.tsv.gz (one row per genome). Usage: build_hist.py nodes.parquet out [threads]"""
import sys, numpy as np, pandas as pd, pyarrow.parquet as pq
from multiprocessing import Pool
EDGES = np.exp(np.linspace(np.log(0.005), np.log(2.5), 71))
def one(a):
    i, path = a
    t = pq.ParquetFile(path).read_row_group(i, columns=["genome", "kind", "n_valid", "dS_mean"]).to_pandas()
    tree = t[t.kind == "tree"]; v = tree[(tree.n_valid > 0) & tree.dS_mean.notna()].dS_mean.values
    pos = v[v > 0]; h, _ = np.histogram(np.clip(pos, EDGES[0] + 1e-12, EDGES[-1] - 1e-12), bins=EDGES)
    r = {"genome": t.genome.iloc[0], "n_zero_tree": int((v == 0).sum()), "n_identical": int((t.kind == "identical").sum()), "n_no_ds": int(len(tree) - len(v)), "n_below_floor": int((pos < EDGES[0]).sum())}
    r.update({f"b{k:02d}": int(c) for k, c in enumerate(h)}); return r
if __name__ == "__main__":
    path, out = sys.argv[1:3]; thr = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    n = pq.ParquetFile(path).metadata.num_row_groups
    with Pool(thr) as pool: rows = list(pool.imap(one, [(i, path) for i in range(n)], chunksize=16))
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False); print("genomes", len(rows))
