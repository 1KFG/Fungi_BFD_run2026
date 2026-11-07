#!/usr/bin/python3.12
"""Duplication mode per genome: selective duplication of a gene subset (tandem arrays, dispersed TE-like copies, enriched functions)
versus block duplication (WGD, segmental, hybrid with two assembled subgenomes), from gene order in the BFD tables.
Definitions (pairs from tables/wgd.ks.parquet; gene index = rank of the gene start on its contig, from gene_info):
  recent pair      dS <= RECENT (0.1; includes identical-CDS pairs)
  class            tandem: same contig and |index difference| <= 5; cis_far: same contig, farther; trans: different contigs
  collinear nbr    another pair with dS <= NBR_DS (0.3) between the same two contigs whose two gene indices are both within K (5) genes of this pair's indices,
                   in the same direction (both increase) or inverted, and not the same pair. A pair is 'in a block' if it has >= 2 collinear neighbours (block of >= 3 pairs).
  null             the same counts after shifting gene index 2 of every pair by a pseudo-random 50-500 genes (keeps contig membership, breaks order).
  breadth          genes in recent pairs / genes in the genome.
  concentration    Jensen-Shannon divergence between the Pfam-domain composition of recent-pair genes and of all genes of the genome (log2, 0 = same).
Outputs: data/dupmode_recent_pairs.parquet, tables/dupmode_genome.tsv.gz.
Limit: wgd.ks keeps at most 5 sampled pairs per node, so neighbouring pairs of large families are under-counted; collinearity is a lower bound."""
import os, numpy as np, pandas as pd, duckdb
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
RECENT, NBR_DS, K = 0.1, 0.3, 5
TMP=os.environ.get("CLAUDE_JOB_DIR", E+"/data")+"/duck_tmp"; os.makedirs(TMP,exist_ok=True)
con=duckdb.connect(); con.execute(f"pragma temp_directory='{TMP}'; pragma threads=8; pragma memory_limit='90GB'")
TE=["PF00078","PF07727","PF13456","PF00665","PF00872","PF20700","PF01609","PF13359","PF13843","PF03184","PF13358","PF02992","PF04754","PF10551","PF14372"]
# gene index per contig
con.execute(f"""create table gi as select gene_id, chrom, start, row_number() over (partition by split_part(gene_id,'_',1), chrom order by start) idx from '{TB}/gene_info.parquet'""")
con.execute(f"create table pm as select distinct protein_id, gene_id from '{TB}/gene_proteins.parquet'")
con.execute(f"""create table pr as select k.genome, k.family, k.status, k.dS, k.g1, k.g2, a.chrom ca, a.idx ia, b.chrom cb, b.idx ib
   from '{TB}/wgd.ks.parquet' k join pm p1 on p1.protein_id=k.g1 join gi a on a.gene_id=p1.gene_id join pm p2 on p2.protein_id=k.g2 join gi b on b.gene_id=p2.gene_id
   where k.dS is not null and k.dS<={NBR_DS}""")
print("pairs dS<=0.3 with positions:",con.sql("select count(*) from pr").fetchone(),flush=True)
# normalise so (c1,i1) <= (c2,i2); add pair id and null-shifted i2
con.execute("""create table pn as select row_number() over () pid, genome, family, status, dS, g1, g2,
   case when ca<=cb then ca else cb end c1, case when ca<=cb then cb else ca end c2,
   case when ca<=cb then ia else ib end i1, case when ca<=cb then ib else ia end i2,
   50 + (hash(g1||g2) % 451)::int shift_n from pr""")
con.execute("alter table pn add column i2n bigint; update pn set i2n=i2+shift_n")
def nbrs(col2,name):
    """collinear neighbours: distinct gene on both sides (index differences non-zero), within K genes on both axes, consistent orientation; keep the larger orientation count."""
    con.execute(f"drop table if exists {name}")
    con.execute(f"""create table {name} as select pid, max(n) n from (
       select p.pid, sign((q.i1-p.i1)*(q.{col2}-p.{col2})) ori, count(distinct q.pid) n
       from pn p join pn q on q.genome=p.genome and q.c1=p.c1 and q.c2=p.c2 and q.pid<>p.pid
         and abs(q.i1-p.i1) between 1 and {K} and abs(q.{col2}-p.{col2}) between 1 and {K}
       where p.dS<={RECENT} group by 1,2) group by 1""")
nbrs("i2","nb_obs"); print("observed neighbours done",flush=True)
nbrs("i2n","nb_null"); print("null neighbours done",flush=True)
con.execute(f"""create table rp as select p.pid, p.genome, p.family, p.status, p.dS, p.g1, p.g2, p.c1, p.c2, p.i1, p.i2,
   case when p.c1=p.c2 and abs(p.i1-p.i2)<=5 then 'tandem' when p.c1=p.c2 then 'cis_far' else 'trans' end cls,
   coalesce(o.n,0) nbr_obs, coalesce(n.n,0) nbr_null from pn p left join nb_obs o using(pid) left join nb_null n using(pid) where p.dS<={RECENT}""")
con.execute(f"copy rp to '{E}/data/dupmode_recent_pairs.parquet' (format parquet, compression zstd)")
print("recent pairs:",con.sql("select count(*) from rp").fetchone(),flush=True)
# contig-pair concentration of block pairs (one contig pair holding all block pairs suggests a duplicated scaffold, many contig pairs suggest a genome-wide process)
con.execute("""create table bcp as select genome, c1, c2, count(*) n from rp where nbr_obs>=2 group by 1,2,3""")
con.execute("""create table bcm as select genome, count(*) n_block_contig_pairs, sum(n) n_block_pairs, max(n)::double/sum(n) top_contig_pair_share,
   sum(n*n)::double/(sum(n)*sum(n)) herfindahl from bcp group by 1""")
# genome-level metrics
con.execute("""create table gm as select genome, count(*) n_recent_pairs, sum((status='identical')::int) n_identical_pairs,
  avg((cls='tandem')::int) share_tandem, avg((cls='cis_far')::int) share_cis_far, avg((cls='trans')::int) share_trans,
  avg((nbr_obs>=2)::int) share_block, avg((nbr_null>=2)::int) share_block_null, avg((nbr_obs>=1)::int) share_any_nbr, avg((nbr_null>=1)::int) share_any_nbr_null from rp group by 1""")
con.execute("""create table rg as select distinct genome, g1 gene from rp union select distinct genome, g2 from rp""")
con.execute("create table fam as select genome, family, count(distinct g) ng from (select genome,family,g1 g from rp union all select genome,family,g2 from rp) group by 1,2")
con.execute("create table fm as select genome, avg((ng=2)::int) share_2gene_fam_recent, median(ng) median_family_size_recent, max(ng) max_family_size_recent, count(*) n_families_recent from fam group by 1")
con.execute("create table gcount as select genome, count(distinct g) genes_recent from (select genome,g1 g from rp union select genome,g2 from rp) group by 1")
G=con.sql("select * from gm join fm using(genome) join gcount using(genome) left join bcm using(genome)").df()
# gene totals and Pfam composition (JS divergence), TE share
S=con.sql(f"select genome, species_prefix prefix, n_pairs from '{TB}/wgd.ks.summary.parquet'").df(); G=G.merge(S,on="genome",how="left")
gt=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)[["genome","genes","phylum","class","order","genus","busco_duplicated_pct","genes_ratio_clade","hybrid","share_2gene_fam"]]
G=G.merge(gt,on="genome",how="left"); G["breadth"]=G.genes_recent/pd.to_numeric(G.genes,errors="coerce")
con.register("rgx",con.sql("select split_part(g,'_',1) prefix, g gene from (select g1 g from rp union select g2 from rp)").df())
tl=",".join(f"'{t}'" for t in TE)
pf=con.sql(f"""select split_part(p.protein_id,'_',1) prefix, split_part(p.pfam_acc,'.',1) acc, count(distinct p.protein_id) n_all,
     count(distinct case when r.gene is not null then p.protein_id end) n_rec
   from '{TB}/pfam.parquet' p left join rgx r on r.gene=p.protein_id group by 1,2""").df()
print("pfam rows",len(pf),flush=True)
def js(a,b):
    a=a/a.sum(); b=b/b.sum(); m=(a+b)/2
    kl=lambda x,y: np.sum(np.where(x>0, x*np.log2(x/np.where(y>0,y,1)),0)); return 0.5*kl(a,m)+0.5*kl(b,m)
rows=[]; rng=np.random.default_rng(2026)
for pre,g in pf.groupby("prefix"):
    nrec=int(g.n_rec.sum())
    if nrec<30: continue
    a=g.n_all.values.astype(float); o=js(a,g.n_rec.values.astype(float))
    p=a/a.sum(); nul=np.array([js(a,rng.multinomial(nrec,p).astype(float)) for _ in range(30)])   # null: nrec Pfam hits drawn at random from the proteome composition
    te=g[g.acc.isin(TE)]; rows.append((pre,o,nul.mean(),nul.std(),(o-nul.mean())/max(nul.std(),1e-6),te.n_rec.sum()/max(nrec,1),nrec))
J=pd.DataFrame(rows,columns=["prefix","js_pfam_recent_vs_proteome","js_null_mean","js_null_sd","js_z","te_share_of_recent_domain_hits","recent_domain_hits"])
G=G.merge(J,on="prefix",how="left")
G["block_excess"]=G.share_block-G.share_block_null
G.to_csv(f"{E}/tables/dupmode_genome.tsv.gz",sep="\t",index=False)
print(len(G),"genomes;",G[["breadth","share_tandem","share_cis_far","share_trans","share_block","share_block_null","js_pfam_recent_vs_proteome","js_null_mean","js_z","top_contig_pair_share"]].describe().round(3).to_string())
