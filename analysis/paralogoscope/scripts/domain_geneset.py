#!/usr/bin/python3.12
"""Describe the gene set carrying one Pfam domain among low-dS paralogs.
Usage: domain_geneset.py PF12013 [label]
Writes tables/geneset_<label>_*.tsv and prints the key numbers.
Reports: domain-carrying proteins per genome and clade; dS-bin distribution of the carriers against all paralog genes;
pair structure (carrier-carrier vs carrier-other); protein length and Pfam model coverage; co-occurring Pfam domains; genomic context of
recent carrier-carrier pairs (same contig, distance). A gene's dS bin is its minimum dS over sampled pairs (as in gene_recent_enrich.py)."""
import sys, numpy as np, pandas as pd, duckdb
ACC=sys.argv[1]; LAB=sys.argv[2] if len(sys.argv)>2 else ACC
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
con=duckdb.connect(); con.execute("pragma threads=8; pragma memory_limit='40GB'")
G=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)
S=con.sql(f"select genome, species_prefix from '{TB}/wgd.ks.summary.parquet'").df()
G=G.merge(S,on="genome",how="inner",suffixes=("","_s")); con.register("gt",G[["genome","species_prefix_s" if "species_prefix_s" in G else "species_prefix","phylum","class","order","genus"]].rename(columns={"species_prefix_s":"species_prefix"}))
con.execute(f"""create table dom as select protein_id, count(*) n_hits, min(domain_i_evalue) ev, min(hmm_from) hmm_from, max(hmm_to) hmm_to, max(pfam_len) pfam_len, min(ali_from) ali_from, max(ali_to) ali_to
  from '{TB}/pfam.parquet' where pfam_acc like '{ACC}%' group by 1""")
con.execute("alter table dom add column prefix varchar; update dom set prefix=split_part(protein_id,'_',1)")
n=con.sql("select count(*), count(distinct prefix) from dom").fetchone(); print("domain proteins",n,flush=True)
con.execute(f"create table pr as select protein_id, length from '{TB}/gene_proteins.parquet' where protein_id in (select protein_id from dom)")
# per-gene min dS (all genes in pairs, restricted to the focal genomes later)
con.execute(f"""create table gb as select gene, min(ds) ds from (select g1 gene, dS ds from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5
   union all select g2, dS from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5) group by 1""")
con.execute("alter table gb add column bin varchar; update gb set bin=case when ds=0 then 'Z' when ds<=0.1 then 'H1' when ds<=0.3 then 'Y' when ds<=1.5 then 'MID' else 'OLD' end")
allbin=con.sql("select bin, count(*) n from gb group by 1").df(); allbin["share"]=allbin.n/allbin.n.sum()
d=con.sql("select dom.protein_id, dom.prefix, gb.bin, gb.ds, pr.length, (dom.hmm_to-dom.hmm_from+1)::double/dom.pfam_len cov, dom.ev, dom.n_hits from dom left join gb on gb.gene=dom.protein_id left join pr using(protein_id)").df()
d["bin"]=d.bin.fillna("not in a gated pair")
dist=d.groupby("bin").size().rename("domain_genes").reset_index().merge(allbin.rename(columns={"n":"all_paralog_genes","share":"share_all"}),on="bin",how="outer").fillna(0)
dist["share_domain_in_pairs"]=dist.domain_genes/dist[dist.bin!="not in a gated pair"].domain_genes.sum()
dist.loc[dist.bin=="not in a gated pair","share_domain_in_pairs"]=np.nan
print(dist.round(4).to_string(),flush=True); dist.to_csv(f"{E}/tables/geneset_{LAB}_bins.tsv",sep="\t",index=False)
# per-genome / clade
d=d.merge(G[["species_prefix_s" if "species_prefix_s" in G else "species_prefix","genome","phylum","class","order","genus","genes","busco_duplicated_pct"]].rename(columns={"species_prefix_s":"prefix","species_prefix":"prefix"}),on="prefix",how="left")
d["recent"]=d.bin.isin(["Z","H1"])
pg=d.groupby("genome").agg(phylum=("phylum","first"),cls=("class","first"),order=("order","first"),genus=("genus","first"),genes=("genes","first"),busco_dup=("busco_duplicated_pct","first"),
        domain_proteins=("protein_id","size"),recent_genes=("recent","sum"),identical_genes=("bin",lambda x:(x=="Z").sum())).reset_index()
pg["domain_per_1000_genes"]=pg.domain_proteins/pg.genes*1000; pg["recent_share"]=pg.recent_genes/pg.domain_proteins
pg.sort_values("recent_genes",ascending=False).to_csv(f"{E}/tables/geneset_{LAB}_per_genome.tsv",sep="\t",index=False)
nall=G.genome.nunique(); print("genomes with domain",len(pg),"of",nall,"; with >=1 recent",int((pg.recent_genes>0).sum()),">=5",int((pg.recent_genes>=5).sum()),">=20",int((pg.recent_genes>=20).sum()),">=100",int((pg.recent_genes>=100).sum()))
tot=pg.recent_genes.sum(); top10=pg.recent_genes.nlargest(10).sum(); print("recent domain genes",int(tot),"top10 genomes share",round(top10/tot,3),flush=True)
cl=pg.groupby(["phylum","cls"]).agg(genomes_with=("genome","size"),domain_proteins=("domain_proteins","sum"),recent_genes=("recent_genes","sum"),median_per_1000=("domain_per_1000_genes","median")).reset_index()
tg=G.groupby(["phylum","class"]).size().rename("genomes_total").reset_index().rename(columns={"class":"cls"}); cl=cl.merge(tg,on=["phylum","cls"],how="left"); cl["frac_genomes_with"]=cl.genomes_with/cl.genomes_total
cl.sort_values("recent_genes",ascending=False).to_csv(f"{E}/tables/geneset_{LAB}_by_class.tsv",sep="\t",index=False); print(cl.sort_values("recent_genes",ascending=False).head(12).round(3).to_string(),flush=True)
print("recent_share by phylum:"); print(pg.groupby("phylum").agg(n=("genome","size"),med_domain=("domain_proteins","median"),med_recent_share=("recent_share","median")).round(3).to_string())
# protein features
f=d.groupby("bin").agg(n=("protein_id","size"),med_len=("length","median"),q25=("length",lambda x:x.quantile(.25)),q75=("length",lambda x:x.quantile(.75)),med_domain_cov=("cov","median"),frac_cov_ge_0_8=("cov",lambda x:(x>=0.8).mean()),frac_short150=("length",lambda x:(x<150).mean())).reset_index()
print(f.round(3).to_string(),flush=True); f.to_csv(f"{E}/tables/geneset_{LAB}_features.tsv",sep="\t",index=False)
# pair structure: both partners carry the domain?
con.register("dm",d[["protein_id"]])
pairs=con.sql(f"""select p.genome, p.g1, p.g2, p.dS, p.status, (p.g1 in (select protein_id from dm)) a, (p.g2 in (select protein_id from dm)) b
   from '{TB}/wgd.ks.parquet' p where p.dS is not null and p.dS<=0.1 and (p.g1 in (select protein_id from dm) or p.g2 in (select protein_id from dm))""").df()
pairs["both"]=pairs.a&pairs.b
print("pairs dS<=0.1 with >=1 domain gene:",len(pairs),"; both genes carry domain:",round(pairs.both.mean(),3),"; identical:",int((pairs.status=="identical").sum()),flush=True)
allp=con.sql(f"select count(*) from '{TB}/wgd.ks.parquet' where dS is not null and dS<=0.1").fetchone()[0]; print("all pairs dS<=0.1:",allp,"; domain pairs share",round(len(pairs)/allp,4))
# genomic context of carrier-carrier recent pairs
bp=pairs[pairs.both].copy()
con.register("bp",bp[["g1","g2"]])
loc=con.sql(f"""select bp.g1, bp.g2, a.chrom c1, a.start s1, b.chrom c2, b.start s2 from bp
  join '{TB}/gene_proteins.parquet' pa on pa.protein_id=bp.g1 join '{TB}/gene_info.parquet' a on a.gene_id=pa.gene_id
  join '{TB}/gene_proteins.parquet' pb on pb.protein_id=bp.g2 join '{TB}/gene_info.parquet' b on b.gene_id=pb.gene_id""").df()
loc["same_contig"]=loc.c1==loc.c2; loc["dist"]=(loc.s1-loc.s2).abs().where(loc.same_contig)
print("carrier-carrier recent pairs",len(loc),"same contig",round(loc.same_contig.mean(),3),"median dist (bp, same contig)",loc.dist.median(),"within 20kb:",round((loc.dist<=20000).sum()/len(loc),3),"within 100kb:",round((loc.dist<=100000).sum()/len(loc),3),flush=True)
loc.to_csv(f"{E}/tables/geneset_{LAB}_recent_pair_context.tsv.gz",sep="\t",index=False)
# co-occurring Pfam domains on carrier proteins (recent only vs all)
con.register("dr",d[d.recent][["protein_id"]]); con.register("da",d[["protein_id"]])
co=con.sql(f"""select split_part(p.pfam_acc,'.',1) acc, any_value(p.pfam_id) pname, count(distinct p.protein_id) proteins from '{TB}/pfam.parquet' p where p.protein_id in (select protein_id from dr) and p.pfam_acc not like '{ACC}%' group by 1 order by 3 desc limit 25""").df()
co["share_of_recent_carriers"]=co.proteins/len(d[d.recent]); print(co.round(3).to_string(),flush=True); co.to_csv(f"{E}/tables/geneset_{LAB}_cooccurring_pfam_recent.tsv",sep="\t",index=False)
cob=con.sql(f"select count(distinct protein_id) from '{TB}/pfam.parquet' where protein_id in (select protein_id from dr) and pfam_acc not like '{ACC}%'").fetchone()[0]; print("recent carriers with any other Pfam domain:",cob,"of",len(d[d.recent]))
