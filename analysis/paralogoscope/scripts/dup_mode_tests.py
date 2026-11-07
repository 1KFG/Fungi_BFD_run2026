#!/usr/bin/python3.12
"""Provisional duplication modes and pair-level tests (uses dup_mode.py outputs).
Informative genomes: n_eligible >= 50 and share_block_null <= 0.1 (the permutation null does not itself produce blocks).
Modes (exploratory thresholds, stated here, not validated on held-out data):
  tandem_array        share_tandem >= 0.3
  multi_contig_blocks share_block >= 0.5 and blocks on >= 5 contig pairs      (WGD / hybrid / segmental / haplotig-like; not separable here)
  single_contig_pair  share_block >= 0.5 and blocks on < 5 contig pairs       (duplicated scaffold or one segmental duplication)
  dispersed_selective share_block < 0.1, share_tandem < 0.3 and top-5 families hold >= 0.5 of recent genes
  dispersed_diffuse   share_block < 0.1, share_tandem < 0.3 and top-5 families hold < 0.5
  mixed               otherwise
Pair-level tests (eligible recent pairs, informative genomes): share of pairs in a block for pairs with a TE-domain gene, an OrsD (PF12013) gene, identical CDS, versus other pairs of the same genomes.
Outputs tables/dupmode_modes*.tsv, dupmode_pair_categories.tsv, figures/fig17_dupmode.*"""
import numpy as np, pandas as pd, duckdb, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
TE=["PF00078","PF07727","PF13456","PF00665","PF00872","PF20700","PF01609","PF13359","PF13843","PF03184","PF13358","PF02992","PF04754","PF10551","PF14372"]
G=pd.read_csv(f"{E}/tables/dupmode_genome.tsv.gz",sep="\t",low_memory=False).drop_duplicates("genome")
C=pd.read_csv(f"{E}/tables/controls.tsv",sep="\t")[["control_group","genome"]].drop_duplicates("genome"); G=G.merge(C,on="genome",how="left"); G["control_group"]=G.control_group.fillna("other")
G["informative"]=(G.n_eligible>=50)&(G.share_block_null<=0.1)
def mode(r):
    if not r.informative: return "uninformative"
    if r.share_tandem>=0.3: return "tandem_array"
    if r.share_block>=0.5: return "multi_contig_blocks" if r.n_block_contig_pairs>=5 else "single_contig_pair"
    if r.share_block<0.1 and r.share_tandem<0.3: return "dispersed_selective" if r.top5_family_gene_share>=0.5 else "dispersed_diffuse"
    return "mixed"
G["mode"]=G.apply(mode,axis=1); G.to_csv(f"{E}/tables/dupmode_genome_modes.tsv.gz",sep="\t",index=False)
inf=G[G.informative]
M=inf.groupby("mode").agg(genomes=("genome","size"),med_breadth=("breadth","median"),med_share_block=("share_block","median"),med_busco_dup=("busco_duplicated_pct","median"),med_gene_ratio=("genes_ratio_clade","median"),
     med_te_share=("te_share_of_recent_domain_hits","median"),med_identical_pair_share=("n_identical_pairs",lambda s:0),med_top5_family_share=("top5_family_gene_share","median"),med_block_dS=("block_median_dS","median")).reset_index()
M["med_identical_pair_share"]=[ (inf[inf["mode"]==m].n_identical_pairs/inf[inf["mode"]==m].n_recent_pairs).median() for m in M["mode"]]
M.to_csv(f"{E}/tables/dupmode_modes.tsv",sep="\t",index=False); print(M.round(3).to_string(),flush=True)
print("genomes:",len(G),"informative:",len(inf)); print(G["mode"].value_counts().to_dict())
cm=pd.crosstab(inf.control_group,inf["mode"]); cm.to_csv(f"{E}/tables/dupmode_modes_by_control.tsv",sep="\t"); print(cm.to_string())
bp=pd.crosstab(inf.phylum,inf["mode"]); bp.to_csv(f"{E}/tables/dupmode_modes_by_phylum.tsv",sep="\t"); print(bp.to_string(),flush=True)
# pair-level tests
con=duckdb.connect(); con.execute("pragma threads=8; pragma memory_limit='60GB'")
tl=",".join(f"'{t}'" for t in TE)
con.execute(f"create table te as select distinct protein_id from '{TB}/pfam.parquet' where split_part(pfam_acc,'.',1) in ({tl})")
con.execute(f"create table og as select distinct protein_id from '{TB}/pfam.parquet' where pfam_acc like 'PF12013%'")
con.register("infg",inf[["genome"]])
P=con.sql(f"""select r.genome, r.cls, r.status='identical' ident, (nbr_obs>=2) blk, (r.g1 in (select protein_id from te) or r.g2 in (select protein_id from te)) te, (r.g1 in (select protein_id from og) or r.g2 in (select protein_id from og)) og
   from read_parquet('{E}/data/dupmode_recent_pairs.parquet') r join infg using(genome) where r.eligible""").df()
print("eligible recent pairs in informative genomes:",len(P),flush=True)
rows=[]
for cat in ["te","og","ident"]:
    f=P[P[cat]]; o=P[~P[cat]]
    gf=f.groupby("genome").blk.agg(["size","mean"]); go=o.groupby("genome").blk.mean().rename("o"); j=gf.join(go,how="inner"); j=j[j["size"]>=5]
    d=(j["mean"]-j["o"]); p=stats.wilcoxon(d).pvalue if len(d)>=10 and (d!=0).any() else np.nan
    rows.append(dict(category={"te":"pair has a TE-domain gene","og":"pair has an OrsD gene","ident":"identical-CDS pair"}[cat],pairs=len(f),share_in_block_pooled=f.blk.mean(),share_in_block_other_pairs_pooled=o.blk.mean(),
                     genomes_compared=len(j),median_share_category=j["mean"].median(),median_share_other=j["o"].median(),frac_genomes_category_higher=(d>0).mean() if len(d) else np.nan,p_wilcoxon=p))
T=pd.DataFrame(rows); T.to_csv(f"{E}/tables/dupmode_pair_categories.tsv",sep="\t",index=False); print(T.round(4).to_string(),flush=True)
cl=P.groupby("cls").blk.agg(["size","mean"]).reset_index(); cl.to_csv(f"{E}/tables/dupmode_pair_class_block.tsv",sep="\t",index=False); print(cl)
# figure
col={"multi_contig_blocks":"#D55E00","single_contig_pair":"#E69F00","tandem_array":"#009E73","dispersed_selective":"#0072B2","dispersed_diffuse":"#56B4E9","mixed":"#999999","uninformative":"#dddddd"}
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
for m,c in col.items():
    d=G[G["mode"]==m]; ax[0].scatter(d.breadth.clip(lower=1e-4),d.share_block.fillna(-0.05),s=6,c=c,alpha=.6,label=f"{m} ({len(d)})")
for m_,mk in [("named hybrid (' x ' in name)","*"),("literature hybrid species (strain may be collapsed)","P")]:
    d=G[G.control_group==m_]; ax[0].scatter(d.breadth.clip(lower=1e-4),d.share_block.fillna(-0.05),s=60,marker=mk,facecolors="none",edgecolors="k",lw=.8)
ax[0].set_xscale("log"); ax[0].set_xlabel("breadth: share of genes in recent pairs"); ax[0].set_ylabel("share of eligible recent pairs in a collinear block"); ax[0].legend(fontsize=6,frameon=False,loc="center left"); ax[0].set_title("Modes (open symbols: named and literature hybrids)",fontsize=8)
b=G[G.informative&G.block_median_dS.notna()]; ax[1].scatter(b.n_block_contig_pairs,b.block_median_dS,s=7,c=[col[m] for m in b["mode"]],alpha=.6); ax[1].set_xscale("log"); ax[1].set_xlabel("contig pairs holding blocks"); ax[1].set_ylabel("median dS of block pairs"); ax[1].set_title("Block pairs: dS against how many contig pairs",fontsize=8)
for gname,mk in [("Verticillium_longisporum_PD589","V. longisporum"),("Puccinia_graminis_21-0","P. graminis 21-0"),("Zygosaccharomyces_bailii_ISA1307","Z. bailii")]:
    d=b[b.genome==gname]
    if len(d): ax[1].annotate(mk,(d.n_block_contig_pairs.iloc[0],d.block_median_dS.iloc[0]),fontsize=7,xytext=(3,3),textcoords="offset points")
ax[2].bar(range(len(T)),T.share_in_block_pooled,width=.35,label="category",color="#D55E00"); ax[2].bar(np.arange(len(T))+.35,T.share_in_block_other_pairs_pooled,width=.35,label="other pairs",color="#999")
ax[2].set_xticks(np.arange(len(T))+.17); ax[2].set_xticklabels(["TE-domain\ngene","OrsD\ngene","identical\nCDS"],fontsize=8); ax[2].set_ylabel("share of eligible recent pairs in a block"); ax[2].legend(frameon=False,fontsize=8); ax[2].set_title("Pair categories (informative genomes)",fontsize=8)
fig.tight_layout(); fig.savefig(f"{E}/figures/fig17_dupmode.png",dpi=200); fig.savefig(f"{E}/figures/fig17_dupmode.pdf")
