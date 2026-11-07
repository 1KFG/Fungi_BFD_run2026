#!/usr/bin/python3.12
"""Shape-based detection of interior dS peaks (centre <= 1.0) from the per-genome log-dS histogram,
cross-referenced with BUSCO duplication (tables/busco_genome.parquet via ASMID, joined in genome_table_final).
A peak = local maximum of the 3-bin smoothed density that has
  prominence >= PROM * peak height, centre in [0.03, 1.0], and mass within +-1 half-prominence width >= MINMASS of usable nodes.
Genomes need >= MINN usable nodes with dS>0."""
import numpy as np, pandas as pd, duckdb
from scipy.signal import find_peaks, peak_widths
E="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
EDG=np.exp(np.linspace(np.log(0.005),np.log(2.5),71)); CEN=np.sqrt(EDG[:-1]*EDG[1:])
PROM,MINMASS,MINN,CMAX,CMIN=0.35,0.08,300,1.0,0.03
H=pd.read_csv(f"{E}/data/ds_hist.tsv.gz",sep="\t",keep_default_na=False)
bc=[f"b{i:02d}" for i in range(70)]
G=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)
X=H.merge(G[["genome","ASMID","phylum","class","order","genus","genes","genes_ratio_clade","share_2gene_fam","nodes_per_1000_genes","busco_complete_pct","busco_duplicated_pct","N50_bp","contig_count","te_per_1000_genes","hybrid"]],on="genome",how="inner")
rows=[]
for d in X.to_dict('records'):
    c=np.array([float(d[b]) for b in bc]); n=c.sum()
    if n<MINN: continue
    p=c/n; s=np.convolve(np.pad(p,1,mode="edge"),np.ones(3)/3,mode="valid")
    pk,pr=find_peaks(np.pad(s,1,constant_values=0),prominence=PROM*s.max()*0.5)  # candidate maxima, padded so edge maxima are valid
    pk=pk-1; keep=[]
    for i,k in enumerate(pk):
        if k<0 or k>=70: continue
        if not (CMIN<=CEN[k]<=CMAX): continue
        prom=pr["prominences"][i]
        if prom<PROM*s[k]: continue
        w=peak_widths(np.pad(s,1,constant_values=0),[k+1],rel_height=0.5)
        lo,hi=int(max(0,np.floor(w[2][0])-1)),int(min(69,np.ceil(w[3][0])-1))
        mass=p[lo:hi+1].sum()
        if mass<MINMASS: continue
        keep.append((CEN[k],mass,s[k],prom/s[k],np.log(EDG[min(hi+1,70)]/EDG[lo])))
    base=dict(genome=d["genome"],ASMID=d["ASMID"],phylum=d["phylum"],cls=d["class"],genus=d["genus"],n_usable=int(n),
              genes=d["genes"],genes_ratio_clade=d["genes_ratio_clade"],share_2gene_fam=d["share_2gene_fam"],
              busco_dup=d["busco_duplicated_pct"],busco_complete=d["busco_complete_pct"],N50=d["N50_bp"],contigs=d["contig_count"],
              te_per_1000=d["te_per_1000_genes"],hybrid=d["hybrid"],n_peaks=len(keep))
    for j,(cen,mass,h,pm,wd) in enumerate(sorted(keep,key=lambda t:-t[1])[:2],1):
        base.update({f"peak{j}_dS":cen,f"peak{j}_mass":mass,f"peak{j}_prom":pm,f"peak{j}_logwidth":wd})
    rows.append(base)
Z=pd.DataFrame(rows)
for c in ["busco_dup","genes_ratio_clade","share_2gene_fam","N50","contigs","te_per_1000"]: Z[c]=pd.to_numeric(Z[c],errors="coerce")
Z["busco_tier"]=np.where(Z.busco_dup>=40,"high (>=40%)",np.where(Z.busco_dup>=10,"moderate (10-40%)","low (<10%)"))
Z.to_csv(f"{E}/tables/shape_peaks_all.tsv",sep="\t",index=False)
P=Z[Z.n_peaks>0].copy(); P.to_csv(f"{E}/tables/shape_peaks.tsv",sep="\t",index=False)
print("genomes evaluated",len(Z),"with >=1 interior peak <=1.0:",len(P),"; 1 peak",int((P.n_peaks==1).sum()),"; 2+",int((P.n_peaks>=2).sum()))
print(pd.crosstab(P.phylum,P.busco_tier))
print(pd.crosstab(P.cls,P.busco_tier).assign(t=lambda d:d.sum(1)).sort_values("t",ascending=False).head(15))
print(Z.groupby("busco_tier").apply(lambda d:pd.Series({"n":len(d),"frac_with_peak":(d.n_peaks>0).mean()})))
pat="Zygosaccharomyces|Verticillium_longisporum|Ogataea_polymorpha_x|Candida_orthopsilosis|Aspergillus_latus|pastorianus|Millerozyma_sorbitophila|Cryptococcus_neoformans_x"
print(Z[Z.genome.str.contains(pat)][["genome","n_peaks","peak1_dS","peak1_mass","peak2_dS","peak2_mass","busco_dup","genes_ratio_clade"]].round(3).to_string())

# ---- tiers: sharp = FWHM <= 1.2 ln units, mass >= 0.4, prominence >= 0.8
P["sharp"]=(P.peak1_logwidth<=1.2)&(P.peak1_mass>=0.4)&(P.peak1_prom>=0.8)
P["tier"]=np.where(~P.sharp,"D broad hump",np.where(P.busco_dup>=40,"A sharp + BUSCO dup>=40%",np.where(P.busco_dup>=10,"B sharp + BUSCO dup 10-40%","C sharp + BUSCO dup<10%")))
P.to_csv(f"{E}/tables/shape_peaks.tsv",sep="\t",index=False)
print(P.tier.value_counts().sort_index())
print(pd.crosstab(P.cls,P.tier).assign(t=lambda d:d.sum(axis=1)).sort_values("t",ascending=False).head(12))
print(P[P.tier.str.startswith("A")].sort_values("busco_dup",ascending=False)[["genome","cls","peak1_dS","peak1_mass","peak1_logwidth","busco_dup","genes_ratio_clade","share_2gene_fam"]].head(25).round(2).to_string())
print(P[P.tier.str.startswith("C")].sort_values("peak1_mass",ascending=False)[["genome","cls","peak1_dS","peak1_mass","peak1_logwidth","busco_dup","genes_ratio_clade","share_2gene_fam","te_per_1000"]].head(15).round(2).to_string())
