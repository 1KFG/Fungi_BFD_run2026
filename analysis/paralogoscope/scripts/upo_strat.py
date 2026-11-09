#!/usr/bin/python3.12
"""Genome-stratified checks for the PF01328 (UPO) duplication-mode result.
Question: are recent UPO-UPO pairs (dS<=0.1) enriched for (a) collinear-block membership and (b) tandem class, against other recent pairs
in the SAME genome? Mantel-Haenszel OR over genomes (RBG variance), leave-one-genome-out range, and result without the 5 genomes with most UPO pairs.
Inputs: data/dupmode_recent_pairs.parquet, data/dupmode_blocks.parquet, tables/pfam.parquet. Writes tables/upo_strat_*.tsv.
Pairs inside one genome are not independent; the genome is the stratum and the MH test treats pairs as independent within a stratum, so p-values are optimistic.
The leave-one-out and top-5-removed estimates are the more robust read."""
import numpy as np, pandas as pd, duckdb
from scipy import stats
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; T=f"{R}/tables"
c=duckdb.connect(); c.execute("pragma threads=8; pragma memory_limit='40GB'")
c.execute(f"create table dom as select distinct protein_id from '{T}/pfam.parquet' where pfam_acc like 'PF01328%'")
c.execute(f"create table rp as select * from '{E}/data/dupmode_recent_pairs.parquet'")
c.execute("create table gk as select distinct genome from rp join dom da on da.protein_id=rp.g1 join dom db on db.protein_id=rp.g2")
c.execute(f"create table bl as select * from '{E}/data/dupmode_blocks.parquet'")
c.execute("""create table pp as select rp.pid, rp.genome, rp.cls, rp.eligible, (da.protein_id is not null and db.protein_id is not null) upo, rp.c1,rp.c2,rp.i1,rp.i2
  from rp join gk using(genome) left join dom da on da.protein_id=rp.g1 left join dom db on db.protein_id=rp.g2""")
c.execute("""create table inb as select distinct pp.pid from pp join bl b on b.genome=pp.genome and
  ((b.c1=pp.c1 and b.c2=pp.c2 and pp.i1 between b.i1min and b.i1max and pp.i2 between b.i2min and b.i2max)
   or (b.c1=pp.c2 and b.c2=pp.c1 and pp.i1 between b.i2min and b.i2max and pp.i2 between b.i1min and b.i1max)) where pp.eligible""")
d=c.sql("select pp.genome, pp.cls, pp.eligible, pp.upo, (inb.pid is not null) in_block from pp left join inb using(pid)").df()
print("genomes",d.genome.nunique(),"pairs",len(d),"UPO pairs",int(d.upo.sum()),flush=True)

def mh(df,flag):
    """df: genome, upo(bool), flag(bool). Returns MH OR, 95% CI (RBG), p (MH chi2), n strata used, UPO n."""
    g=df.groupby(["genome","upo"])[flag].agg(["sum","count"]).unstack("upo").fillna(0)
    if (True not in g["count"].columns) or (False not in g["count"].columns): return (np.nan,)*5+(0,)
    a=g["sum"][True]; n1=g["count"][True]; c_=g["sum"][False]; n0=g["count"][False]
    b=n1-a; d_=n0-c_; N=n1+n0; ok=(n1>0)&(n0>0); a,b,c_,d_,N,n1,n0=[x[ok] for x in (a,b,c_,d_,N,n1,n0)]
    if len(N)==0: return (np.nan,)*5+(0,)
    R_=(a*d_/N).sum(); S=(b*c_/N).sum()
    if R_==0 or S==0: return (np.nan,)*5+(int(n1.sum()),)
    P=(a+d_)/N; Q=(b+c_)/N; r=a*d_/N; s=b*c_/N
    var=(P*r).sum()/(2*R_**2)+((P*s+Q*r).sum())/(2*R_*S)+(Q*s).sum()/(2*S**2)
    OR=R_/S; lo,hi=np.exp(np.log(OR)+np.array([-1.96,1.96])*np.sqrt(var))
    m=(a+b)*(a+c_)/N; v=(n1*n0*(a+c_)*(b+d_))/(N**2*(N-1).clip(lower=1))
    chi=(abs((a-m).sum())-0.5)**2/v.sum(); p=stats.chi2.sf(chi,1)
    return OR,lo,hi,p,int(ok.sum()),int(n1.sum())

rows=[]
for label,sub,flag in (("block_membership",d[d.eligible],"in_block"),("tandem_class",d.assign(tandem=d.cls=="tandem"),"tandem")):
    sub=sub.copy(); sub["upo"]=sub.upo.astype(bool); sub[flag]=sub[flag].astype(bool)
    top=sub[sub.upo].groupby("genome").size().sort_values(ascending=False)
    base=mh(sub,flag); rows.append((label,"all genomes",*base))
    rows.append((label,"top5 UPO genomes removed",*mh(sub[~sub.genome.isin(top.index[:5])],flag)))
    loo=[mh(sub[sub.genome!=g],flag)[0] for g in top.index[:30]]
    rows.append((label,"leave-one-out (30 largest UPO genomes) min OR",np.nanmin(loo),np.nan,np.nan,np.nan,np.nan,np.nan))
    rows.append((label,"leave-one-out max OR",np.nanmax(loo),np.nan,np.nan,np.nan,np.nan,np.nan))
    # per-genome share (genomes with >=10 UPO pairs): UPO share vs other share
    pg=sub.groupby(["genome","upo"])[flag].agg(["mean","size"]).unstack("upo")
    big=pg[(pg["size"][True]>=10)&(pg["size"][False]>=10)]
    diff=(big["mean"][True]-big["mean"][False])
    rows.append((label,f"per-genome diff (UPO-other), genomes with >=10 each: n={len(big)}, median diff, share>0",diff.median(),(diff>0).mean(),np.nan,np.nan,len(big),np.nan))
    if len(big)>=5: rows.append((label,"sign test p (share of genomes with UPO>other)",stats.binomtest(int((diff>0).sum()),len(big),0.5).pvalue,np.nan,np.nan,np.nan,len(big),np.nan))
    print(label,"top genomes by UPO pairs:",top.head(8).to_dict(),flush=True)
out=pd.DataFrame(rows,columns=["test","subset","OR_or_value","ci_lo_or_value2","ci_hi","p_mh","strata","upo_pairs"])
print(out.round(4).to_string()); out.to_csv(f"{E}/tables/upo_strat_tests.tsv",sep="\t",index=False)
