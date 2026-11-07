#!/usr/bin/python3.12
"""Gene-level query: which paralog genes sit in young dS bins, and do they share Pfam / CAZy domains or short-protein features?
Gene set = genes in gated sampled pairs of tables/wgd.ks.parquet (identical-CDS groups: one representative per group).
A gene is assigned the MINIMUM dS over the sampled pairs it appears in (biased young when a gene has several partners; stated in the report).
Bins: Z dS==0, H1 (0,0.1], Y (0.1,0.3], MID (0.3,1.5], OLD >1.5; RECENT = Z+H1.
Test: per domain, Mantel-Haenszel odds ratio with each genome as a stratum (bin B genes vs the genome's other paralog genes);
genomes need >=5 genes on each side. Variance: Robins-Breslow-Greenland. BH-FDR across domains. Sets: all genomes / clean (assembly-suspect removed) / clean per phylum."""
import os, sys, numpy as np, pandas as pd, duckdb
from scipy import stats
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
TMP=os.environ.get("CLAUDE_JOB_DIR", E+"/data")+"/duck_tmp"; os.makedirs(TMP,exist_ok=True)
con=duckdb.connect(); con.execute(f"pragma temp_directory='{TMP}'"); con.execute("pragma memory_limit='40GB'"); con.execute("pragma threads=8")
# ---- genome sets
G=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)
G["busco_duplicated_pct"]=pd.to_numeric(G.busco_duplicated_pct,errors="coerce"); G["genes_ratio_clade"]=pd.to_numeric(G.genes_ratio_clade,errors="coerce")
SP=pd.read_csv(f"{E}/tables/shape_peaks.tsv",sep="\t"); bad=set(SP[SP.tier.str.startswith(("A","B"))].genome)
G["clean"]=~((G.busco_duplicated_pct>=20)|(G.genes_ratio_clade>=1.5)|G.genome.isin(bad)|(G.hybrid.astype(str).isin(["1","True"])))
sets={"all":G.genome, "clean":G[G.clean].genome, "clean_Ascomycota":G[G.clean&(G.phylum=="Ascomycota")].genome, "clean_Basidiomycota":G[G.clean&(G.phylum=="Basidiomycota")].genome}
print({k:len(v) for k,v in sets.items()},flush=True)
con.register("gset",pd.DataFrame({"genome":G.genome,"clean":G.clean,"phylum":G.phylum}))
# ---- gene -> min dS
con.execute(f"""create table gb as
 select k.genome, k.gene, min(k.dS) ds from (
   select genome, g1 gene, dS from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5
   union all select genome, g2, dS from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5) k
 join gset s using(genome) group by 1,2""")
con.execute("""alter table gb add column bin varchar; update gb set bin=case when ds=0 then 'Z' when ds<=0.1 then 'H1' when ds<=0.3 then 'Y' when ds<=1.5 then 'MID' else 'OLD' end""")
print(con.sql("select bin,count(*) from gb group by 1 order by 1").fetchall(),flush=True)
# ---- features
con.execute(f"create table gp as select protein_id gene, length from '{TB}/gene_proteins.parquet'")
con.execute(f"create table sp as select distinct protein_id gene from '{TB}/signalp.signal_peptide.parquet'")
con.execute("create table gf as select gb.genome,gb.gene,gb.bin,gp.length, (sp.gene is not null) sigp from gb left join gp using(gene) left join sp using(gene)")
print("length match",con.sql("select avg((length is not null)::int) from gf").fetchone(),flush=True)
# ---- domains
con.execute(f"""create table gd as select distinct gb.genome, gb.gene, 'PF:'||split_part(p.pfam_acc,'.',1) dom from gb join '{TB}/pfam.parquet' p on p.protein_id=gb.gene""")
con.execute(f"""insert into gd select distinct gb.genome, gb.gene, 'CAZy:'||regexp_replace(split_part(c.cazyme_fam,'(',1),'_.*','') from gb join '{TB}/cazy.overview.parquet' c on c.protein_id=gb.gene""")
print("domain rows",con.sql("select count(*), count(distinct dom) from gd").fetchone(),flush=True)
# name map for pfam
nm=con.sql(f"select split_part(pfam_acc,'.',1) a, any_value(pfam_id) n from '{TB}/pfam.parquet' group by 1").df(); nmap=dict(zip("PF:"+nm.a,nm.n))
# ---- per genome totals and per genome-domain counts by bin
con.execute("create table gtot as select genome, bin, count(*) n from gb group by 1,2")
con.execute("create table gtot_pf as select genome, bin, count(*) n from gb where (genome,gene) in (select genome,gene from gd where dom like 'PF:%') group by 1,2")
ann=con.sql("select bin, count(*) n, sum(((genome,gene) in (select genome,gene from gd where dom like 'PF:%'))::int) a from gb group by 1").df(); ann["frac_with_pfam"]=ann.a/ann.n; print(ann,flush=True); ann.to_csv(f"{E}/tables/gene_pfam_annotated_fraction_by_bin.tsv",sep="\t",index=False)
con.execute("create table gdc as select genome, dom, bin, count(*) a from gd join gb using(genome,gene) group by 1,2,3")
def run(B_bins,tag,prefix,totals):
    bl=",".join(f"'{b}'" for b in B_bins)
    con.execute("drop table if exists st")
    con.execute(f"""create table st as
      with t as (select genome, sum(case when bin in ({bl}) then n else 0 end) nB, sum(case when bin not in ({bl}) then n else 0 end) nR from {totals} group by 1 having nB>=5 and nR>=5),
      d as (select genome, dom, sum(case when bin in ({bl}) then a else 0 end) aB, sum(case when bin not in ({bl}) then a else 0 end) aR from gdc where dom like '{prefix}%' group by 1,2 having sum(a)>0)
      select d.genome, d.dom, d.aB, d.aR, t.nB, t.nR from d join t using(genome)""")
    out={}
    for sname,gs in sets.items():
        con.register("gs",pd.DataFrame({"genome":gs}))
        df=con.sql("""select dom, sum((aB>0)::int) n_genomes, sum(aB) aB_tot, sum(aR) aR_tot,
          sum(aB*(nR-aR)::double/(nB+nR)) num, sum((nB-aB)*aR::double/(nB+nR)) den,
          sum(((aB+nR-aR)::double/(nB+nR))*(aB*(nR-aR)::double/(nB+nR))) PR,
          sum(((aB+nR-aR)::double/(nB+nR))*((nB-aB)*aR::double/(nB+nR)) + ((nB-aB+aR)::double/(nB+nR))*(aB*(nR-aR)::double/(nB+nR))) PSQR,
          sum(((nB-aB+aR)::double/(nB+nR))*((nB-aB)*aR::double/(nB+nR))) QS
          from st join gs using(genome) group by 1 having sum((aB>0)::int)>=20 and sum(aB)>=30""").df()
        # genomes contributing nB,nR in this set (denominator context)
        tot=con.sql("select count(*) n, sum(nB) nB, sum(nR) nR from (select genome,sum(case when bin in (%s) then n else 0 end) nB,sum(case when bin not in (%s) then n else 0 end) nR from %s group by 1 having nB>=5 and nR>=5) join gs using(genome)"%(bl,bl,totals)).fetchone()
        # den==0 or tiny => unstable
        df=df[(df.den>0)&(df.num>0)].copy()
        df["OR_MH"]=df.num/df.den; df["lnOR"]=np.log(df.OR_MH)
        df["se"]=np.sqrt(df.PR/(2*df.num**2)+df.PSQR/(2*df.num*df.den)+df.QS/(2*df.den**2))
        df["z"]=df.lnOR/df.se; df["p"]=2*stats.norm.sf(abs(df.z))
        o=np.argsort(df.p.values); m=len(df); q=np.empty(m); r=df.p.values[o]*m/np.arange(1,m+1); q[o]=np.minimum.accumulate(r[::-1])[::-1]; df["q_BH"]=np.minimum(q,1)
        df["ci_lo"]=np.exp(df.lnOR-1.96*df.se); df["ci_hi"]=np.exp(df.lnOR+1.96*df.se)
        df["name"]=df.dom.map(nmap).fillna(""); df["set"]=sname; df["bins"]="+".join(B_bins); df["genomes_in_test"]=tot[0]; df["genes_B"]=tot[1]; df["genes_ref"]=tot[2]
        df=df.drop(columns=["num","den","PR","PSQR","QS"]).sort_values("p")
        df.to_csv(f"{E}/tables/gene_enrich_{tag}_{prefix.strip(chr(58)).lower()}_{sname}.tsv",sep="\t",index=False)
        sig=df[(df.q_BH<0.01)&(df.OR_MH>=2)]
        print(tag,prefix,sname,"genomes",tot[0],"genes_B",tot[1],"domains tested",len(df),"enriched q<0.01 & OR>=2:",len(sig),flush=True)
        print(sig.head(8)[["dom","name","n_genomes","aB_tot","OR_MH","ci_lo","ci_hi","q_BH"]].round(3).to_string(),flush=True)
for B,tag in [(["Z","H1"],"recent"),(["Z"],"identical"),(["Y"],"Y"),(["MID"],"MID")]:
    run(B,tag,"PF:","gtot_pf"); run(B,tag,"CAZy:","gtot")
# ---- features: per genome medians by bin; paired across genomes
F=con.sql("""select genome, bin, count(*) n, median(length) med_len, avg((length<150)::int) frac_short150, avg(sigp::int) frac_sigp from gf where length is not null group by 1,2""").df()
F=F.merge(G[["genome","phylum","clean"]],on="genome"); F.to_csv(f"{E}/data/gene_features_by_genome_bin.tsv.gz",sep="\t",index=False)
rows=[]
for sname in sets:
    gs=set(sets[sname]); f=F[F.genome.isin(gs)]
    ref=f[f.bin.isin(["OLD"])]  # compare against OLD (dS>1.5) paralogs, the most saturated reference
    for b in ["Z","H1","Y","MID"]:
        x=f[f.bin==b].set_index("genome"); r=ref.set_index("genome"); j=x.join(r,lsuffix="_B",rsuffix="_R",how="inner"); j=j[(j.n_B>=10)&(j.n_R>=10)]
        for m in ["med_len","frac_short150","frac_sigp"]:
            d=(j[m+"_B"]-j[m+"_R"]).dropna()
            if len(d)<20: continue
            rows.append(dict(set=sname,bin=b,metric=m,n_genomes=len(d),median_B=j[m+"_B"].median(),median_OLD=j[m+"_R"].median(),median_diff=d.median(),frac_genomes_higher=(d>0).mean(),p_wilcoxon=stats.wilcoxon(d).pvalue))
pd.DataFrame(rows).to_csv(f"{E}/tables/gene_features_vs_old.tsv",sep="\t",index=False)
print(pd.DataFrame(rows).round(4).to_string())
