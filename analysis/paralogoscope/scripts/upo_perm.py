#!/usr/bin/python3.12
"""Permutation test for PF01328 (UPO) duplication mode of recent pairs (dS<=0.1), conditional on genome.
Statistic T = sum over genomes of (observed UPO flagged pairs - n_UPO_pairs_g * genome base rate), where the flag is (a) in a collinear block (eligible pairs) or (b) tandem class.
Null A (pair-level): shuffle the UPO label among all pairs of the genome (ignores that pairs of one family are correlated).
Null B (family-level, the one to trust): the unit is (genome, family, UPO label). Shuffle the UPO label among units within a genome, so pairs of one gene family move together.
One-sided p = (1 + #perm T >= observed) / (1 + n_perm). Also reports T for genomes with >=10 UPO pairs and with the top 5 UPO genomes removed.
Inputs as in upo_strat.py. Writes tables/upo_perm_tests.tsv."""
import numpy as np, pandas as pd, duckdb
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; T=f"{R}/tables"
NPERM=10000; rng=np.random.default_rng(20261108)
c=duckdb.connect(); c.execute("pragma threads=8; pragma memory_limit='40GB'")
c.execute(f"create table dom as select distinct protein_id from '{T}/pfam.parquet' where pfam_acc like 'PF01328%'")
c.execute(f"create table rp as select * from '{E}/data/dupmode_recent_pairs.parquet'")
c.execute("create table gk as select distinct genome from rp join dom da on da.protein_id=rp.g1 join dom db on db.protein_id=rp.g2")
c.execute(f"create table bl as select * from '{E}/data/dupmode_blocks.parquet'")
c.execute("""create table pp as select rp.pid, rp.genome, rp.family fam, rp.cls, rp.eligible, (da.protein_id is not null and db.protein_id is not null) upo, rp.c1,rp.c2,rp.i1,rp.i2
  from rp join gk using(genome) left join dom da on da.protein_id=rp.g1 left join dom db on db.protein_id=rp.g2""")
c.execute("""create table inb as select distinct pp.pid from pp join bl b on b.genome=pp.genome and
  ((b.c1=pp.c1 and b.c2=pp.c2 and pp.i1 between b.i1min and b.i1max and pp.i2 between b.i2min and b.i2max)
   or (b.c1=pp.c2 and b.c2=pp.c1 and pp.i1 between b.i2min and b.i2max and pp.i2 between b.i1min and b.i1max)) where pp.eligible""")
d=c.sql("select pp.genome, pp.fam, pp.cls, pp.eligible, pp.upo, (inb.pid is not null) in_block from pp left join inb using(pid)").df()
d["tandem"]=d.cls=="tandem"

def prep(sub,flag,level):
    """Return list of per-genome (n_units array, k array(n flagged), m array(n pairs), upo bool array)."""
    if level=="pair":
        sub=sub.assign(fam="x")  # placeholder, not used
        g=sub.groupby("genome")
        out=[]
        for _,x in g:
            out.append((x.upo.values.astype(bool),x[flag].values.astype(float),np.ones(len(x))))
        return out
    u=sub.groupby(["genome","fam","upo"]).agg(k=(flag,"sum"),m=(flag,"size")).reset_index()
    return [(x.upo.values.astype(bool),x.k.values.astype(float),x.m.values.astype(float)) for _,x in u.groupby("genome")]

def stat(units,perm):
    T=0.0
    for upo,k,m in units:
        M=m.sum()
        if M==0: continue
        lab=rng.permutation(upo) if perm else upo
        base=k.sum()/M; T+=k[lab].sum()-base*m[lab].sum()
    return T

def stat_fast(units,nperm,level):
    """Per genome: observed excess and nperm null draws with the same number of UPO units (pair level: hypergeometric; family level: random unit subsets)."""
    obs=0.0; null=np.zeros(nperm)
    for upo,k,m in units:
        M=m.sum(); n=len(upo); nu=int(upo.sum())
        if M==0 or nu==0 or nu==n: continue
        base=k.sum()/M; obs+=k[upo].sum()-base*m[upo].sum()
        if level=="pair":
            X=rng.hypergeometric(int(k.sum()),int(M-k.sum()),nu,size=nperm); null+=X-base*nu
        else:
            w=k-base*m
            idx=np.argsort(rng.random((nperm,n)),axis=1)[:,:nu]
            null+=w[idx].sum(axis=1)
    return obs,null

rows=[]
for label,flag,sub0 in (("block_membership","in_block",d[d.eligible]),("tandem_class","tandem",d)):
    sub0=sub0.copy(); sub0["upo"]=sub0.upo.astype(bool); sub0[flag]=sub0[flag].astype(int)
    top=sub0[sub0.upo].groupby("genome").size().sort_values(ascending=False)
    big=set(top[top>=10].index)
    for sname,sub in (("all genomes",sub0),("genomes with >=10 UPO pairs",sub0[sub0.genome.isin(big)]),("top5 UPO genomes removed",sub0[~sub0.genome.isin(top.index[:5])])):
        for level in ("pair","family"):
            units=prep(sub,flag,level)
            # constant-size subsets: unit-level resampling keeps the number of UPO units, not the number of UPO pairs
            obs,null=stat_fast(units,NPERM,level)
            p=(1+(null>=obs).sum())/(1+NPERM)
            rows.append((label,sname,level,len(units),obs,null.mean(),null.std(),(obs-null.mean())/null.std() if null.std()>0 else np.nan,p))
            print(rows[-1],flush=True)
out=pd.DataFrame(rows,columns=["test","subset","null_unit","genomes","T_obs","T_null_mean","T_null_sd","z","p_one_sided"])
out.to_csv(f"{E}/tables/upo_perm_tests.tsv",sep="\t",index=False); print(out.round(4).to_string())
