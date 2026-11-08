#!/usr/bin/python3.12
"""Select genomes for the codeml second pass (Option E): genomes with a possible interior dS peak (shape tiers A, B, C; peak <= 1.0), excluding giants and genomes already run with codeml.
Cost model (CLAUDE.md sizing guideline: about 1-1.5 h per task): codeml 0.2 s per run, runs about 2.2 x nodes, aligner 7 ms per pair; bundles of <= 40,000 CPU-s.
Output: option_e_codeml/genomes.tsv (genome, transl_table, est_cpu_s, bundle) and bundle_<n>.txt"""
import numpy as np, pandas as pd, duckdb
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; P=f"{R}/paralogoscope_pilot"; E=f"{P}/eval/explore"; O=f"{P}/option_e_codeml"
sp=pd.read_csv(f"{E}/tables/shape_peaks.tsv",sep="\t",keep_default_na=False,na_values=[""])
sel=sp[~sp.tier.str.startswith("D")][["genome","tier","peak1_dS","busco_dup"]].drop_duplicates("genome")
tax=pd.read_csv(f"{E}/data/taxonomy.tsv",sep="\t",keep_default_na=False)[["genome","TRANSL_TABLE"]]
S=duckdb.sql(f"select genome, n_pairs, n_nodes, n_nodes_with_ds from '{R}/tables/wgd.ks.summary.parquet'").df()
d=sel.merge(tax,on="genome",how="left").merge(S,on="genome",how="left")
done=set(pd.read_csv(f"{E}/data/yn00_vs_codeml_genomes.tsv",sep="\t").genome) if False else set(pd.read_csv(f"{E}/tables/yn00_vs_codeml_genomes.tsv",sep="\t").genome)
d["already_codeml"]=d.genome.isin(done)
d["est_cpu_s"]=0.2*2.2*d.n_nodes.fillna(0)+0.007*d.n_pairs.fillna(0)
print("selected tiers A-C:",len(d),"; already codeml:",int(d.already_codeml.sum()),"; est_cpu_s max",round(d.est_cpu_s.max()))
d=d[~d.already_codeml & d.est_cpu_s.between(1,60000)].sort_values("est_cpu_s",ascending=False).reset_index(drop=True)
# first-fit decreasing into bundles of <= 40,000 CPU-s
bins=[]; 
for r in d.itertuples():
    for b in bins:
        if b["load"]+r.est_cpu_s<=40000: b["g"].append(r.genome); b["load"]+=r.est_cpu_s; break
    else: bins.append({"g":[r.genome],"load":r.est_cpu_s})
bm={g:i+1 for i,b in enumerate(bins) for g in b["g"]}; d["bundle"]=d.genome.map(bm)
d.to_csv(f"{O}/genomes.tsv",sep="\t",index=False)
for i,b in enumerate(bins,1): open(f"{O}/bundle_{i:03d}.txt","w").write("\n".join(b["g"])+"\n")
print("genomes to run:",len(d),"; bundles:",len(bins),"; total est CPU-h:",round(d.est_cpu_s.sum()/3600,1)); print(d.tier.str[0].value_counts().to_dict()); print(d.head(5)[["genome","tier","n_nodes","est_cpu_s","TRANSL_TABLE"]].to_string())
