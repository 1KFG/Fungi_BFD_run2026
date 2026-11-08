#!/usr/bin/python3.12
"""Option E: compare codeml (second pass, option_e_codeml/ks) with yn00 (tables/wgd.ks.parquet) for genomes with a possible interior dS peak.
Node dS = mean dS of the gated sampled pairs of a node (status ok, dS <= 2.5, aligned codons >= 30), computed the same way for both estimators.
(1) Per genome: histogram of tree-node dS (70 log bins, 0.005-2.5) and the shape-peak detector of shape_peaks.py (same thresholds) -> peak position and mass with each estimator.
(2) Node level: codeml / yn00 ratio in yn00 dS bins, on nodes present with a dS in both runs (node sampling differs between runs, so single-node ratios are noisy).
Outputs: tables/optionE_peak_shift.tsv, optionE_node_ratio.tsv, figures/fig18_optionE.*"""
import glob, os, numpy as np, pandas as pd, duckdb, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; P=f"{R}/paralogoscope_pilot"; E=f"{P}/eval/explore"; O=f"{P}/option_e_codeml"
EDG=np.exp(np.linspace(np.log(0.005),np.log(2.5),71)); CEN=np.sqrt(EDG[:-1]*EDG[1:]); PROM,MINMASS,MINN,CMAX,CMIN=0.35,0.08,300,1.0,0.03
def detect(v):
    v=v[v>0]
    if len(v)<MINN: return None
    h,_=np.histogram(np.clip(v,EDG[0]+1e-12,EDG[-1]-1e-12),bins=EDG); p=h/h.sum(); s=np.convolve(np.pad(p,1,mode="edge"),np.ones(3)/3,mode="valid")
    pk,pr=find_peaks(np.pad(s,1,constant_values=0),prominence=PROM*s.max()*0.5); pk=pk-1; best=None
    for i,k in enumerate(pk):
        if k<0 or k>=70 or not (CMIN<=CEN[k]<=CMAX): continue
        if pr["prominences"][i]<PROM*s[k]: continue
        w=peak_widths(np.pad(s,1,constant_values=0),[k+1],rel_height=0.5); lo,hi=int(max(0,np.floor(w[2][0])-1)),int(min(69,np.ceil(w[3][0])-1)); mass=p[lo:hi+1].sum()
        if mass<MINMASS: continue
        inwin=v[(v>=EDG[lo])&(v<EDG[min(hi+1,70)])]; cand=(float(np.median(inwin)) if len(inwin) else CEN[k],mass,np.log(EDG[min(hi+1,70)]/EDG[lo]))   # refined position: median of node dS inside the peak window
        if best is None or cand[1]>best[1]: best=cand
    return best
def node_means(pairs):
    q=pairs[(pairs.status=="ok")&pairs.dS.notna()&(pairs.dS<=2.5)&(pairs.strippedalignmentlength>=30)]
    return q.groupby(["family","node"]).dS.mean().rename("dS")
files=sorted(glob.glob(f"{O}/ks/*.cds-transcripts.fa.tsv.ks.tsv.gz")); files=[f for f in files if not f.endswith(".nodes.tsv.gz")]
genomes=[os.path.basename(f).replace(".cds-transcripts.fa.tsv.ks.tsv.gz","") for f in files]; print("codeml genomes:",len(genomes),flush=True)
con=duckdb.connect()
Y=con.sql(f"select genome, family, node, status, dS, strippedalignmentlength from '{R}/tables/wgd.ks.parquet' where genome in ({','.join(repr(g) for g in genomes)}) and status<>'identical'").df()
rows=[]; nodes=[]
for g,f in zip(genomes,files):
    c=pd.read_csv(f,sep="\t",usecols=["family","node","status","dS","strippedalignmentlength"]); cm=node_means(c)
    ym=node_means(Y[Y.genome==g]); j=pd.concat([ym.rename("yn00"),cm.rename("codeml")],axis=1)
    both=j.dropna(); d1,d2=detect(ym.values),detect(cm.values)
    rows.append(dict(genome=g,n_nodes_yn00=len(ym),n_nodes_codeml=len(cm),n_both=len(both),yn00_peak=d1[0] if d1 else np.nan,yn00_mass=d1[1] if d1 else np.nan,yn00_logwidth=d1[2] if d1 else np.nan,
                     codeml_peak=d2[0] if d2 else np.nan,codeml_mass=d2[1] if d2 else np.nan,codeml_logwidth=d2[2] if d2 else np.nan,median_ratio_both=(both.codeml/both.yn00).replace([np.inf],np.nan).median() if len(both) else np.nan))
    both=both[both.yn00>0].copy(); both["genome"]=g; nodes.append(both.reset_index(drop=True)[["genome","yn00","codeml"]])
T=pd.DataFrame(rows); sh=pd.read_csv(f"{E}/tables/shape_peaks.tsv",sep="\t")[["genome","tier","peak1_dS"]].rename(columns={"peak1_dS":"shape_peak_tableyn00"}); T=T.merge(sh,on="genome",how="left")
T["shift_ratio"]=T.codeml_peak/T.yn00_peak; T.to_csv(f"{E}/tables/optionE_peak_shift.tsv",sep="\t",index=False)
both_peak=T.dropna(subset=["yn00_peak","codeml_peak"]); print("genomes:",len(T),"| peak detected with yn00:",int(T.yn00_peak.notna().sum()),"| with codeml:",int(T.codeml_peak.notna().sum()),"| both:",len(both_peak))
print("peak shift ratio codeml/yn00 (both): median",round(both_peak.shift_ratio.median(),3),"IQR",both_peak.shift_ratio.quantile([.25,.75]).round(3).tolist(),flush=True)
N=pd.concat(nodes); N["bin"]=pd.cut(N.yn00,[0,0.05,0.1,0.2,0.3,0.5,0.8,1.2,1.6,2.5]); N["ratio"]=N.codeml/N.yn00
NR=N.groupby("bin",observed=True).agg(n_nodes=("ratio","size"),median_yn00=("yn00","median"),median_codeml=("codeml","median"),median_ratio=("ratio","median"),p10=("ratio",lambda x:x.quantile(.1)),p90=("ratio",lambda x:x.quantile(.9)),share_within_20pct=("ratio",lambda x:((x>0.8)&(x<1.2)).mean())).reset_index()
NR.to_csv(f"{E}/tables/optionE_node_ratio.tsv",sep="\t",index=False); print(NR.round(3).to_string())
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
ax[0].scatter(both_peak.yn00_peak,both_peak.codeml_peak,s=14,c="#0072B2",alpha=.7); lim=[0.03,1.6]; ax[0].plot(lim,lim,"k:",lw=.7); ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("peak dS, yn00"); ax[0].set_ylabel("peak dS, codeml"); ax[0].set_title(f"Peak position, {len(both_peak)} genomes with a peak in both",fontsize=8)
for r in both_peak[both_peak.genome.str.contains("Verticillium_longisporum|Zygosaccharomyces_bailii|Ogataea_polymorpha_x")].itertuples(): ax[0].annotate(r.genome.split("_")[0][:5]+"_"+r.genome.split("_")[1][:6],(r.yn00_peak,r.codeml_peak),fontsize=6,xytext=(3,3),textcoords="offset points")
ax[1].plot(range(len(NR)),NR.median_ratio,"o-",c="#D55E00"); ax[1].fill_between(range(len(NR)),NR.p10,NR.p90,color="#D55E00",alpha=.2); ax[1].axhline(1,color="k",lw=.6,ls=":"); ax[1].set_xticks(range(len(NR))); ax[1].set_xticklabels([str(b) for b in NR.bin],rotation=45,fontsize=6,ha="right"); ax[1].set_ylabel("codeml / yn00 node dS (median, p10-p90)"); ax[1].set_title("Node dS ratio by yn00 bin",fontsize=8)
ax[2].hist(both_peak.shift_ratio.dropna(),bins=25,color="#999"); ax[2].axvline(1,color="k",lw=.6,ls=":"); ax[2].set_xlabel("peak position ratio codeml / yn00"); ax[2].set_ylabel("genomes")
fig.tight_layout(); fig.savefig(f"{E}/figures/fig18_optionE.png",dpi=200); fig.savefig(f"{E}/figures/fig18_optionE.pdf")
