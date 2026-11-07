#!/usr/bin/python3.12
"""Figure 15: Pfam enrichment volcano plots (recent dS<=0.1 and Y 0.1-0.3, clean genomes) from tables/gene_enrich_*_pf_clean.tsv,
and protein length by dS bin from data/gene_features_by_genome_bin.tsv.gz. Writes figures/fig15_gene_enrichment.{png,pdf}."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
fig,ax=plt.subplots(1,3,figsize=(15,4.8))
for a,(t,ttl) in zip(ax[:2],[("recent","dS <= 0.1 (identical + H1) vs all other paralog genes"),("Y","dS 0.1-0.3 vs all other paralog genes")]):
    d=pd.read_csv(f"{E}/tables/gene_enrich_{t}_pf_clean.tsv",sep="\t"); d["lq"]=(-np.log10(d.q_BH.clip(lower=1e-300))).clip(upper=60); d["l2"]=np.log2(d.OR_MH)
    sig=(d.q_BH<0.01)&(d.OR_MH>=2)
    a.scatter(d.l2[~sig],d.lq[~sig],s=4,c="#bbb"); a.scatter(d.l2[sig],d.lq[sig],s=np.clip(np.sqrt(d.aB_tot[sig])/4,3,30),c="#D55E00",alpha=.6)
    top=d[sig&(d.aB_tot>=300)].sort_values("OR_MH",ascending=False).head(14)
    for r in top.itertuples(): a.annotate(r.name if isinstance(r.name,str) and r.name else r.dom,(r.l2,r.lq),fontsize=6,xytext=(2,2),textcoords="offset points")
    a.axvline(1,color="k",lw=.5,ls=":"); a.set_xlabel("log2 odds ratio (MH, genome x length strata)"); a.set_ylabel("-log10 q (BH), capped at 60"); a.set_title(f"Pfam, clean genomes: {ttl}",fontsize=7)
F=pd.read_csv(f"{E}/data/gene_features_by_genome_bin.tsv.gz",sep="\t",keep_default_na=False,na_values=[""])
F=F[F.clean.astype(str)=="True"]; order=["Z","H1","Y","MID","OLD"]; lab=["dS=0","0-0.1","0.1-0.3","0.3-1.5",">1.5"]
for m,c,off in [("med_len","#0072B2",0),("frac_short150","#D55E00",1)]:
    pass
vals=[F[(F.bin==b)&(F.n>=10)].med_len.dropna() for b in order]; v2=[F[(F.bin==b)&(F.n>=10)].frac_short150.dropna() for b in order]
a=ax[2]; bp=a.boxplot(vals,tick_labels=lab,showfliers=False,patch_artist=True); [p.set_facecolor("#9ec5e6") for p in bp["boxes"]]
a.set_ylabel("per-genome median protein length (aa)"); a.set_xlabel("dS bin of the gene (minimum over sampled pairs)"); a2=a.twinx(); a2.plot(range(1,6),[v.median() for v in v2],"o-",c="#D55E00"); a2.set_ylabel("median share of proteins < 150 aa",color="#D55E00")
a.set_title("Protein length by dS bin (clean genomes, >=10 genes per bin)",fontsize=7)
fig.tight_layout(); fig.savefig(f"{E}/figures/fig15_gene_enrichment.png",dpi=200); fig.savefig(f"{E}/figures/fig15_gene_enrichment.pdf")
