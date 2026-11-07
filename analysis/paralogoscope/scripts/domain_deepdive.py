#!/usr/bin/python3.12
"""Per-genome deep dive for a domain-carrying gene set. Usage: domain_deepdive.py PF12013 OrsD
(1) For every genome with >=5 carrier genes in a gated pair: carrier recent share against the same genome's recent share among all other paralog genes (Fisher exact), BUSCO dup, class.
    A genome whose whole paralog set is recent (assembly duplication) has no carrier-specific excess.
(2) Gene table for Sorochytrium milnesiophthora (Blastocladiomycetes, 21 recent carriers): contig, position, length, dS, partner, family, other Pfam domains.
Outputs tables/geneset_<label>_genome_enrichment.tsv, geneset_<label>_Sorochytrium_genes.tsv"""
import sys, numpy as np, pandas as pd, duckdb
from scipy import stats
ACC=sys.argv[1]; LAB=sys.argv[2]
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
con=duckdb.connect(); con.execute("pragma threads=8; pragma memory_limit='40GB'")
con.execute(f"create table dom as select distinct protein_id from '{TB}/pfam.parquet' where pfam_acc like '{ACC}%'")
con.execute(f"""create table gb as select genome, gene, min(ds) ds from (select genome, g1 gene, dS ds from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5
   union all select genome, g2, dS from '{TB}/wgd.ks.parquet' where dS is not null and dS<=2.5) group by 1,2""")
con.execute("create table gg as select genome, count(*) n, sum((ds<=0.1)::int) rec, sum((gene in (select protein_id from dom))::int) c, sum(((gene in (select protein_id from dom)) and ds<=0.1)::int) c_rec from gb group by 1")
g=con.sql("select * from gg where c>=5").df()
G=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)[["genome","class","order","genes","busco_duplicated_pct","genes_ratio_clade"]]
g=g.merge(G,on="genome",how="left")
rows=[]
for r in g.itertuples():
    a,b=r.c_rec,r.c-r.c_rec; cc,dd=(r.rec-r.c_rec),((r.n-r.c)-(r.rec-r.c_rec))
    OR,p=stats.fisher_exact([[a,b],[cc,dd]]); rows.append((r.genome,r._11 if False else None,r.n,r.rec/r.n,r.c,r.c_rec,r.c_rec/r.c,OR,p))
out=pd.DataFrame(rows,columns=["genome","_","paralog_genes","recent_share_all_paralogs","carriers_in_pairs","carriers_recent","carrier_recent_share","fisher_OR","p"]).drop(columns="_")
out=out.merge(G,on="genome",how="left"); out["q_BH"]=np.minimum(1,out.p*len(out)/out.p.rank()); 
out=out.sort_values("carriers_recent",ascending=False); out.to_csv(f"{E}/tables/geneset_{LAB}_genome_enrichment.tsv",sep="\t",index=False)
print("genomes with >=5 carriers in pairs:",len(out)); print("carrier excess (OR>1, p<0.01):",int(((out.fisher_OR>1)&(out.p<0.01)).sum()),"; no excess or depleted:",int((out.p>=0.01).sum()+((out.fisher_OR<1)&(out.p<0.01)).sum()))
hi=out[out.busco_duplicated_pct>=40]; print("genomes with BUSCO dup>=40%:",len(hi),"median recent share all paralogs",round(hi.recent_share_all_paralogs.median(),3),"; carriers",round(hi.carrier_recent_share.median(),3))
lo=out[out.busco_duplicated_pct<20]; print("BUSCO dup<20%: n",len(lo),"median recent share all paralogs",round(lo.recent_share_all_paralogs.median(),4),"; carriers",round(lo.carrier_recent_share.median(),3),"; share of genomes with excess (OR>1,p<0.01):",round(((lo.fisher_OR>1)&(lo.p<0.01)).mean(),3))
print(out.head(14)[["genome","class","busco_duplicated_pct","paralog_genes","recent_share_all_paralogs","carriers_in_pairs","carriers_recent","fisher_OR","p"]].round(4).to_string())
# Sorochytrium genes
pre=con.sql("select distinct split_part(g1,'_',1) from '%s/wgd.ks.parquet' where genome like 'Sorochytrium_milnesiophthora%%' limit 1"%TB).fetchone()[0]
con.execute(f"""create table sg as select p.protein_id, g.chrom, g.start, g.strand, p.length, b.ds from '{TB}/gene_proteins.parquet' p join '{TB}/gene_info.parquet' g on g.gene_id=p.gene_id
   left join gb b on b.gene=p.protein_id and b.genome like 'Sorochytrium_milnesiophthora%' where p.protein_id in (select protein_id from dom) and p.protein_id like '{pre}\\_%' escape '\\'""")
con.execute(f"""create table sp as select * from '{TB}/wgd.ks.parquet' where genome like 'Sorochytrium_milnesiophthora%' and (g1 in (select protein_id from sg) or g2 in (select protein_id from sg))""")
sg=con.sql("select s.*, (select count(*) from sp where (sp.g1=s.protein_id or sp.g2=s.protein_id) and sp.dS<=0.1) n_recent_pairs, (select any_value(family) from sp where sp.g1=s.protein_id or sp.g2=s.protein_id) fam, (select string_agg(distinct split_part(pf.pfam_acc,'.',1)||':'||pf.pfam_id, ',') from '"+TB+"/pfam.parquet' pf where pf.protein_id=s.protein_id and pf.pfam_acc not like '"+ACC+"%') other_pfam from sg s order by chrom, start").df()
sg.to_csv(f"{E}/tables/geneset_{LAB}_Sorochytrium_genes.tsv",sep="\t",index=False)
print("Sorochytrium carriers:",len(sg),"; contigs:",sg.chrom.nunique(),"; dS<=0.1:",int((sg.ds<=0.1).sum()),"; families:",sg.fam.nunique())
print(sg.groupby("chrom").size().sort_values(ascending=False).head(6).to_dict()); print(sg.length.describe().round(0).to_dict())
d2=sg.sort_values(["chrom","start"]); d2["gap"]=d2.groupby("chrom").start.diff(); print("median gap between neighbouring carriers on the same contig (bp):",d2.gap.median(),"; gaps<=20kb:",int((d2.gap<=20000).sum()))
print(sg[["protein_id","chrom","start","length","ds","n_recent_pairs","fam","other_pfam"]].head(25).to_string())
