#!/usr/bin/python3.12
"""Test whether recent copies of a domain-carrying gene set sit near transposable-element (TE) genes.
Usage: domain_te_context.py PF12013 OrsD
TE marker: genes carrying one of 15 TE Pfam domains (same list as pfam_te_counts.py). Distance = |start difference| on the same contig (gene start positions).
Focal groups: carriers in recent pairs (dS<=0.1, bins Z/H1) and carriers in older pairs (Y/MID/OLD). Control: non-carrier genes of the same genome and protein-length class
(5 per focal gene when available). Outcome: a TE gene within 25 kb (not the gene itself). Mantel-Haenszel OR with strata = genome x length class; Robins-Breslow-Greenland CI.
Also: share of focal genes that carry a TE Pfam domain themselves, and median distance to the nearest TE gene.
Output: tables/geneset_<label>_te_context*.tsv"""
import sys, numpy as np, pandas as pd, duckdb
from scipy import stats
ACC=sys.argv[1]; LAB=sys.argv[2]
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
TE=["PF00078","PF07727","PF13456","PF00665","PF00872","PF20700","PF01609","PF13359","PF13843","PF03184","PF13358","PF02992","PF04754","PF10551","PF14372"]
con=duckdb.connect(); con.execute("pragma threads=8; pragma memory_limit='48GB'"); rng_seed=2026
con.execute(f"create table dom as select distinct protein_id from '{TB}/pfam.parquet' where pfam_acc like '{ACC}%'")
con.execute(f"""create table gb as select gene, min(ds) ds from (select g1 gene, dS ds from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5
   union all select g2, dS from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5) group by 1""")
con.execute("create table car as select d.protein_id, case when gb.ds<=0.1 then 'recent' when gb.ds is not null then 'older' else 'no_pair' end grp from dom d left join gb on gb.gene=d.protein_id")
print(con.sql("select grp,count(*) from car group by 1").fetchall(),flush=True)
# positions
con.execute(f"""create table pos as select split_part(p.protein_id,'_',1) prefix, p.protein_id, g.chrom, g.start, p.length,
   case when p.length<150 then 1 when p.length<250 then 2 when p.length<400 then 3 when p.length<700 then 4 else 5 end lc
   from '{TB}/gene_proteins.parquet' p join '{TB}/gene_info.parquet' g on g.gene_id=p.gene_id""")
tl=",".join(f"'{t}'" for t in TE)
con.execute(f"create table teg as select distinct protein_id from '{TB}/pfam.parquet' where split_part(pfam_acc,'.',1) in ({tl})")
# focal: carriers in recent / older pairs
con.execute("create table foc as select pos.*, car.grp from pos join car using(protein_id) where car.grp in ('recent','older')")
print("focal",con.sql("select grp,count(*),count(distinct prefix) from foc group by 1").fetchall(),flush=True)
# controls: non-carrier genes of focal genomes, same length class, 5 per focal gene (random within prefix x lc)
con.execute("create table need as select prefix, lc, count(*)*5 q from foc group by 1,2")
con.execute(f"""create table ctl as select * from (select pos.*, 'control' grp, row_number() over (partition by pos.prefix, pos.lc order by hash(pos.protein_id||'{rng_seed}')) rn
   from pos where pos.prefix in (select prefix from foc) and pos.protein_id not in (select protein_id from dom)) c join need n using(prefix,lc) where c.rn<=n.q""")
con.execute("create table allf as select prefix,protein_id,chrom,start,length,lc,grp from foc union all select prefix,protein_id,chrom,start,length,lc,grp from ctl")
print("control",con.sql("select count(*) from ctl").fetchone(),flush=True)
con.execute("""create table tepos as select pos.prefix, pos.protein_id, pos.chrom, pos.start from pos join teg using(protein_id) where pos.prefix in (select prefix from foc)""")
con.execute("""create table nd as select a.protein_id, min(abs(a.start-t.start)) nearest from allf a left join tepos t on t.prefix=a.prefix and t.chrom=a.chrom and t.protein_id<>a.protein_id group by 1""")
con.execute("create table res as select a.*, nd.nearest, (a.protein_id in (select protein_id from teg)) self_te from allf a join nd using(protein_id)")
d=con.sql("select * from res").df(); d["near25"]=d.nearest<=25000; d["near10"]=d.nearest<=10000
sm=d.groupby("grp").agg(n=("protein_id","size"),genomes=("prefix","nunique"),self_TE_domain=("self_te","mean"),TE_gene_within_10kb=("near10","mean"),TE_gene_within_25kb=("near25","mean"),median_nearest_TE_bp=("nearest","median")).reset_index()
print(sm.round(4).to_string(),flush=True)
def mh(f,c):
    # strata prefix x lc ; outcome near25
    rows=[]
    for (p,l),g in pd.concat([f.assign(is_f=1),c.assign(is_f=0)]).groupby(["prefix","lc"]):
        a=int(((g.is_f==1)&g.near25).sum()); b=int(((g.is_f==1)&~g.near25).sum()); cc=int(((g.is_f==0)&g.near25).sum()); dd=int(((g.is_f==0)&~g.near25).sum())
        if a+b>0 and cc+dd>0: rows.append((a,b,cc,dd))
    x=np.array(rows,float); n=x.sum(1); a,b,c_,dd=x.T; R_=a*dd/n; S=b*c_/n; P=(a+dd)/n; Q=(b+c_)/n
    OR=R_.sum()/S.sum(); se=np.sqrt((P*R_).sum()/(2*R_.sum()**2)+((P*S+Q*R_).sum())/(2*R_.sum()*S.sum())+(Q*S).sum()/(2*S.sum()**2))
    return OR,np.exp(np.log(OR)-1.96*se),np.exp(np.log(OR)+1.96*se),len(rows)
ctl=d[d.grp=="control"]; out=[]
for g in ["recent","older"]:
    OR,lo,hi,ns=mh(d[d.grp==g],ctl); out.append(dict(group=g,MH_OR_TE_within_25kb=OR,ci_lo=lo,ci_hi=hi,strata=ns)); print(g,"MH OR TE gene within 25 kb vs matched controls:",round(OR,2),round(lo,2),round(hi,2),"strata",ns,flush=True)
# recent carriers: restrict to clean genomes
G=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)
S=con.sql(f"select genome, species_prefix prefix from '{TB}/wgd.ks.summary.parquet'").df(); G=G.merge(S,on="genome"); G["busco_duplicated_pct"]=pd.to_numeric(G.busco_duplicated_pct,errors="coerce")
clean=set(G[(G.busco_duplicated_pct<20)].prefix)
OR,lo,hi,ns=mh(d[(d.grp=="recent")&d.prefix.isin(clean)],ctl[ctl.prefix.isin(clean)]); out.append(dict(group="recent, BUSCO dup<20% genomes",MH_OR_TE_within_25kb=OR,ci_lo=lo,ci_hi=hi,strata=ns)); print("recent clean",round(OR,2),round(lo,2),round(hi,2),flush=True)
sm.to_csv(f"{E}/tables/geneset_{LAB}_te_context.tsv",sep="\t",index=False); pd.DataFrame(out).to_csv(f"{E}/tables/geneset_{LAB}_te_context_MH.tsv",sep="\t",index=False)
# per genome: recent carriers' near-TE share vs controls (paired across genomes with >=3 recent carriers)
pgf=d[d.grp=="recent"].groupby("prefix").agg(n=("near25","size"),f=("near25","mean")); pgc=ctl.groupby("prefix").near25.mean().rename("c")
j=pgf.join(pgc,how="inner"); j=j[j.n>=3]; diff=(j.f-j.c); print("genomes >=3 recent carriers:",len(j),"median share near TE: carriers",round(j.f.median(),3),"controls",round(j.c.median(),3),"wilcoxon p",stats.wilcoxon(diff).pvalue,flush=True)
pd.DataFrame([dict(genomes=len(j),median_carrier=j.f.median(),median_control=j.c.median(),median_diff=diff.median(),frac_genomes_higher=(diff>0).mean(),p_wilcoxon=stats.wilcoxon(diff).pvalue)]).to_csv(f"{E}/tables/geneset_{LAB}_te_context_paired.tsv",sep="\t",index=False)
