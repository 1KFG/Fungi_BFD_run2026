#!/usr/bin/python3.12
"""PF01328 (UPO) case study: sampling coverage, dS bins by UPO clade, duplication mode of recent pairs.
Reads (read-only): tables/pfam.parquet, tables/wgd.ks*.parquet, data/dupmode_recent_pairs.parquet, data/dupmode_blocks.parquet,
  DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full/upo_sites_bfd_all.tsv.zst (clade/catalytic annotation from another session).
Writes tables/upo_*.tsv. Recent = dS<=0.1 (as in dup_mode.py)."""
import subprocess, io, numpy as np, pandas as pd, duckdb
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; T=f"{R}/tables"
U="/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full/upo_sites_bfd_all.tsv.zst"
c=duckdb.connect(); c.execute("pragma threads=8; pragma memory_limit='40GB'")
ann=pd.read_csv(io.StringIO(subprocess.run(["zstd","-dc",U],capture_output=True,text=True,check=True).stdout),sep="\t",usecols=lambda x: x in ("protein_id","clade","clade_confidence","catalytic_ok"),low_memory=False)
c.register("ann",ann)
c.execute(f"create table dom as select distinct protein_id, split_part(protein_id,'_',1) p from '{T}/pfam.parquet' where pfam_acc like 'PF01328%'")
c.execute(f"create table sg as select species_prefix p from '{T}/wgd.ks.summary.parquet'")
c.execute(f"""create table gb as select gene, min(ds) ds from (select g1 gene, dS ds from '{T}/wgd.ks.parquet' where dS is not null and dS<=2.5
  union all select g2, dS from '{T}/wgd.ks.parquet' where dS is not null and dS<=2.5) group by 1""")
c.execute("""create table d as select dom.protein_id, dom.p, gb.ds, case when gb.ds is null then 'none' when gb.ds=0 then 'Z' when gb.ds<=0.1 then 'H1' when gb.ds<=0.3 then 'Y'
  when gb.ds<=1.5 then 'MID' else 'OLD' end bin, ann.clade, ann.catalytic_ok from dom join sg using(p) left join gb on gb.gene=dom.protein_id left join ann using(protein_id)""")
# 1 sampling coverage by genome copy number
c.execute("create table pg as select p, count(*) n, sum((ds is not null)::int) s from d group by 1")
cov=c.sql("""select case when n=1 then '1' when n<=3 then '2-3' when n<=6 then '4-6' when n<=12 then '7-12' when n<=30 then '13-30' else '>30' end copies_per_genome,
  min(n) nmin, count(*) genomes, sum(n) genes, sum(s) in_pairs, round(sum(s)/sum(n),3) frac_in_pairs from pg group by 1 order by nmin""").df()
print(cov.to_string()); cov.to_csv(f"{E}/tables/upo_sampling_coverage.tsv",sep="\t",index=False)
# 2 dS bin by clade and catalytic core (genes in a pair)
d=c.sql("select * from d where ds is not null").df(); d["clade"]=d.clade.fillna("unannotated")
t=pd.crosstab(d.clade,d.bin); t["n"]=t.sum(axis=1); t["recent_share"]=(t.get("Z",0)+t.get("H1",0))/t.n
print(t.round(3).to_string()); t.to_csv(f"{E}/tables/upo_bins_by_clade.tsv",sep="\t")
t2=pd.crosstab(d.catalytic_ok.fillna("NA"),d.bin); t2["n"]=t2.sum(axis=1); t2["recent_share"]=(t2.get("Z",0)+t2.get("H1",0))/t2.n
print(t2.round(3).to_string()); t2.to_csv(f"{E}/tables/upo_bins_by_catalytic.tsv",sep="\t")
# 3 duplication mode of recent carrier-carrier pairs vs all recent pairs
c.execute(f"create table rp as select * from '{E}/data/dupmode_recent_pairs.parquet'")
c.execute("create table up as select rp.* from rp join dom da on da.protein_id=rp.g1 join dom db on db.protein_id=rp.g2")
mode=c.sql("""select 'UPO_pairs' grp, cls, count(*) n from up group by 2 union all select 'all_recent_pairs', cls, count(*) from rp group by 2 order by 1,2""").df()
mode["share"]=mode.n/mode.groupby("grp").n.transform("sum"); print(mode.round(3).to_string()); mode.to_csv(f"{E}/tables/upo_dupmode_cls.tsv",sep="\t",index=False)
# 4 block membership of eligible UPO pairs (pair index ranges inside a block on the same contig pair)
c.execute(f"create table bl as select * from '{E}/data/dupmode_blocks.parquet'")
bm=c.sql("""select u.pid, u.genome, u.cls, u.dS, (exists(select 1 from bl b where b.genome=u.genome and ((b.c1=u.c1 and b.c2=u.c2 and u.i1 between b.i1min and b.i1max and u.i2 between b.i2min and b.i2max)
   or (b.c1=u.c2 and b.c2=u.c1 and u.i1 between b.i2min and b.i2max and u.i2 between b.i1min and b.i1max)))) in_block, u.eligible from up u where u.eligible""").df()
print("eligible UPO recent pairs",len(bm),"in block",round(bm.in_block.mean(),3) if len(bm) else None)
rp_in=c.sql("""select avg(case when exists(select 1 from bl b where b.genome=r.genome and ((b.c1=r.c1 and b.c2=r.c2 and r.i1 between b.i1min and b.i1max and r.i2 between b.i2min and b.i2max)
   or (b.c1=r.c2 and b.c2=r.c1 and r.i1 between b.i2min and b.i2max and r.i2 between b.i1min and b.i1max))) then 1.0 else 0.0 end) from (select * from rp where eligible using sample 200000 rows) r""").fetchone()[0]
print("all eligible recent pairs in block (200k sample)",round(rp_in,3))
bm.groupby("genome").agg(pairs=("pid","size"),in_block=("in_block","sum")).sort_values("pairs",ascending=False).head(30).to_csv(f"{E}/tables/upo_block_by_genome_top.tsv",sep="\t")
pd.DataFrame({"eligible_UPO_recent_pairs":[len(bm)],"in_block_share":[bm.in_block.mean() if len(bm) else np.nan],"all_eligible_in_block_share_sample":[rp_in]}).to_csv(f"{E}/tables/upo_block_summary.tsv",sep="\t",index=False)
