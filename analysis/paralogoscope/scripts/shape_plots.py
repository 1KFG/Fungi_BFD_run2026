#!/usr/bin/python3.12
"""Figures 13-14: gallery of shape-detected interior dS peaks (tiers A/B/C, named hybrids first) and peak position / share of genomes with a peak
against BUSCO duplication. Reads data/ds_hist.tsv.gz and tables/shape_peaks*.tsv. Writes figures/fig13_*, fig14_*."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
EDG=np.exp(np.linspace(np.log(0.005),np.log(2.5),71)); CEN=np.sqrt(EDG[:-1]*EDG[1:]); bc=[f"b{i:02d}" for i in range(70)]
H=pd.read_csv(f"{E}/data/ds_hist.tsv.gz",sep="\t",keep_default_na=False).set_index("genome")
P=pd.read_csv(f"{E}/tables/shape_peaks.tsv",sep="\t")
col={"A":"#D55E00","B":"#E69F00","C":"#0072B2","D":"#999999"}
def panel(ax,g,r):
    c=H.loc[g,bc].astype(float).values; p=c/c.sum(); ax.fill_between(np.log(CEN),0,p,color=col[r.tier[0]],alpha=.85)
    ax.axvspan(np.log(1.5),np.log(2.5),color="grey",alpha=.12,lw=0)
    for j in (1,2):
        v=r.get(f"peak{j}_dS")
        if pd.notna(v): ax.axvline(np.log(v),color="k",lw=.6,ls=":")
    ax.set_xticks(np.log([0.05,0.3,1.0])); ax.set_xticklabels(["0.05","0.3","1"],fontsize=6); ax.set_yticks([]); 
    nm=g if len(g)<=30 else g[:14]+".."+g[-14:]
    ax.set_title(f"{nm}\npeak {r.peak1_dS:.2f}  mass {r.peak1_mass:.2f}  BUSCO dup {r.busco_dup:.0f}%",fontsize=5.5)
sel=[]
for t,n in [("A",18),("B",6),("C",12)]:
    d=P[P.tier.str.startswith(t)].sort_values("peak1_mass",ascending=False).drop_duplicates("genus")
    sel+=list(d.head(n).itertuples(index=False)) if False else [r for _,r in d.head(n).iterrows()]
named=P[P.genome.str.contains("Zygosaccharomyces|Verticillium_longisporum|Ogataea_polymorpha_x|Cryptococcus_neoformans_x")]
sel=[r for _,r in named.iterrows()]+[r for r in sel if r.genome not in set(named.genome)]
nc=6; nr=int(np.ceil(len(sel)/nc)); fig,axs=plt.subplots(nr,nc,figsize=(13,1.9*nr)); axs=np.atleast_2d(axs)
for k in range(nr*nc):
    a=axs[k//nc,k%nc]
    if k<len(sel): panel(a,sel[k].genome,sel[k])
    else: a.axis("off")
fig.suptitle("Interior dS peaks (centre <= 1.0) detected from histogram shape. Orange: sharp peak + BUSCO dup >= 40%; amber: 10-40%; blue: sharp peak, BUSCO dup < 10%.\nFirst panels: named/literature hybrids and Zygosaccharomyces. Dotted line: detected peak. Grey band: yn00 saturation zone.",fontsize=8,y=1.0)
fig.tight_layout(h_pad=1.5); fig.savefig(f"{E}/figures/fig13_shape_peaks_gallery.png",dpi=200,bbox_inches="tight"); fig.savefig(f"{E}/figures/fig13_shape_peaks_gallery.pdf",bbox_inches="tight")
# summary scatter: peak position vs BUSCO dup
fig,ax=plt.subplots(1,2,figsize=(11,4))
for t in "DCBA":
    d=P[P.tier.str.startswith(t)]; ax[0].scatter(d.busco_dup+0.2,d.peak1_dS,s=8,c=col[t],alpha=.6,label=f"{d.tier.iloc[0]} (n={len(d)})")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("BUSCO duplicated % (+0.2)"); ax[0].set_ylabel("peak dS"); ax[0].legend(fontsize=6,frameon=False)
Z=pd.read_csv(f"{E}/tables/shape_peaks_all.tsv",sep="\t"); Z["has"]=Z.n_peaks>0
Z["bin"]=pd.cut(Z.busco_dup,[-1,1,5,10,20,40,101],labels=["<1","1-5","5-10","10-20","20-40",">=40"])
s=Z.groupby("bin",observed=True).has.agg(["mean","size"]); ax[1].bar(s.index.astype(str),s["mean"],color="#555"); 
for i,(m,n) in enumerate(zip(s["mean"],s["size"])): ax[1].text(i,m+.01,f"n={n}",ha="center",fontsize=6)
ax[1].set_xlabel("BUSCO duplicated %"); ax[1].set_ylabel("share of genomes with an interior peak <= 1.0"); fig.tight_layout()
fig.savefig(f"{E}/figures/fig14_shape_peaks_vs_busco.png",dpi=200); fig.savefig(f"{E}/figures/fig14_shape_peaks_vs_busco.pdf")
print(s)
