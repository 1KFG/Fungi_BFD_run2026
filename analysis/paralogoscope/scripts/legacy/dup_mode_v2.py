#!/usr/bin/python3.12
"""Duplication mode per genome (v2): selective duplication of a gene subset (tandem arrays, dispersed TE-like copies, enriched functions)
versus block duplication (WGD, segmental, hybrid with two assembled subgenomes, assembly duplication). v2 follows an independent review of v1
(legacy/dup_mode_v1.py): fixes pair ordering on one contig, excludes tRNA genes, counts tandem runs and shared-gene pairs once, and uses a permutation null.
Definitions (pairs from tables/wgd.ks.parquet; gene index = rank of the gene start on its contig among protein-coding genes, ties broken by gene_id):
  recent pair   dS <= RECENT (0.1, identical-CDS pairs included). Neighbour candidates: pairs with dS <= NBR_DS (0.3).
  order         each pair is stored with (contig, index) of gene 1 <= (contig, index) of gene 2, also when both genes are on one contig.
  class         tandem: same contig, |index difference| <= 5; cis_far: same contig, farther; trans: different contigs.
  eligible      both contigs hold >= 2K+1 protein-coding genes and the pair is not a same-contig pair with |index difference| <= 2K. Block shares use eligible pairs only.
  collinear     another pair of a different family on the same two contigs, with both gene indices 1..K genes away (K = 5) and the same orientation;
                the count is min(distinct genes on side 1, distinct genes on side 2), taking the larger orientation. In a block = count >= 2.
  null          every gene's index is replaced by its rank under a pseudo-random permutation within its contig (pair membership and contig membership kept, order broken).
  block         run of block pairs on one contig pair with index gaps <= K on both sides (>= 3 pairs). Block table: pairs, genes, median dS, dS IQR, share identical, contig coverage.
  breadth       genes in recent pairs / protein-coding genes in the genome.
Limits: wgd.ks keeps at most 5 sampled pairs per tree node and collapses identical CDS, so breadth, family size and block shares are lower bounds.
Outputs: data/dupmode_recent_pairs.parquet, data/dupmode_blocks.parquet, tables/dupmode_genome.tsv.gz."""
import os, numpy as np, pandas as pd, duckdb
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; E=f"{R}/paralogoscope_pilot/eval/explore"; TB=f"{R}/tables"
RECENT, NBR_DS, K = 0.1, 0.3, 5
TMP=os.environ.get("CLAUDE_JOB_DIR", E+"/data")+"/duck_tmp"; os.makedirs(TMP,exist_ok=True)
con=duckdb.connect(); con.execute(f"pragma temp_directory='{TMP}'; pragma threads=8; pragma memory_limit='90GB'")
TE=["PF00078","PF07727","PF13456","PF00665","PF00872","PF20700","PF01609","PF13359","PF13843","PF03184","PF13358","PF02992","PF04754","PF10551","PF14372"]
con.execute(f"""create table gi as select gene_id, chrom, start, row_number() over (partition by split_part(gene_id,'_',1), chrom order by start, gene_id) idx,
   row_number() over (partition by split_part(gene_id,'_',1), chrom order by hash(gene_id||'perm2026'), gene_id) idxn,
   count(*) over (partition by split_part(gene_id,'_',1), chrom) nct from '{TB}/gene_info.parquet' where gene_type='protein_coding'""")
con.execute(f"create table pm as select distinct protein_id, gene_id from '{TB}/gene_proteins.parquet'")
con.execute(f"""create table pr as select k.genome, k.family, k.status, k.dS, k.g1, k.g2, a.chrom ca, a.idx ia, a.idxn ian, a.nct na, b.chrom cb, b.idx ib, b.idxn ibn, b.nct nb
   from '{TB}/wgd.ks.parquet' k join pm p1 on p1.protein_id=k.g1 join gi a on a.gene_id=p1.gene_id join pm p2 on p2.protein_id=k.g2 join gi b on b.gene_id=p2.gene_id
   where k.dS is not null and k.dS<={NBR_DS} and k.g1<>k.g2""")
print("pairs dS<=0.3 with positions (protein-coding only):",con.sql("select count(*) from pr").fetchone(),flush=True)
def make_pn(tag,ia,ib):
    con.execute(f"drop table if exists pn_{tag}")
    con.execute(f"""create table pn_{tag} as with s as (select *, (ca>cb or (ca=cb and {ia}>{ib})) sw from pr)
      select row_number() over () pid, genome, family, status, dS, g1, g2,
        case when sw then cb else ca end c1, case when sw then ca else cb end c2, case when sw then {ib} else {ia} end i1, case when sw then {ia} else {ib} end i2,
        case when sw then nb else na end n1, case when sw then na else nb end n2 from s""")
def nbrs(tag):
    con.execute(f"drop table if exists nb_{tag}")
    con.execute(f"""create table nb_{tag} as select pid, max(n) n from (
       select p.pid, sign((q.i1-p.i1)*(q.i2-p.i2)) ori, least(count(distinct q.i1), count(distinct q.i2)) n
       from pn_{tag} p join pn_{tag} q on q.genome=p.genome and q.c1=p.c1 and q.c2=p.c2 and q.pid<>p.pid and q.family<>p.family
         and abs(q.i1-p.i1) between 1 and {K} and abs(q.i2-p.i2) between 1 and {K}
       where p.dS<={RECENT} group by 1,2) group by 1""")
make_pn("obs","ia","ib"); nbrs("obs"); print("observed done",flush=True)
make_pn("null","ian","ibn"); nbrs("null"); print("null done",flush=True)
con.execute(f"""create table rp as select p.pid, p.genome, p.family, p.status, p.dS, p.g1, p.g2, p.c1, p.c2, p.i1, p.i2, p.n1, p.n2,
   case when p.c1=p.c2 and abs(p.i1-p.i2)<=5 then 'tandem' when p.c1=p.c2 then 'cis_far' else 'trans' end cls,
   (p.n1>={2*K+1} and p.n2>={2*K+1} and not (p.c1=p.c2 and abs(p.i1-p.i2)<={2*K})) eligible,
   coalesce(o.n,0) nbr_obs from pn_obs p left join nb_obs o using(pid) where p.dS<={RECENT}""")
# the null table has its own pairs/order; attach null counts by pair identity (g1,g2) to keep a single pair-level table
con.execute("""create table nn as select n.g1, n.g2, coalesce(b.n,0) nbr_null, (n.n1>=11 and n.n2>=11 and not (n.c1=n.c2 and abs(n.i1-n.i2)<=10)) elig_null from pn_null n left join nb_null b using(pid) where n.dS<=0.1""")
con.execute("create table rp2 as select r.*, coalesce(m.nbr_null,0) nbr_null, coalesce(m.elig_null,false) elig_null from rp r left join nn m on m.g1=r.g1 and m.g2=r.g2")
con.execute(f"copy rp2 to '{E}/data/dupmode_recent_pairs.parquet' (format parquet, compression zstd)")
print("recent pairs:",con.sql("select count(*) from rp2").fetchone(),flush=True)
# genome metrics
con.execute("""create table gm as select genome, count(*) n_recent_pairs, sum((status='identical')::int) n_identical_pairs, avg((cls='tandem')::int) share_tandem,
  avg((cls='cis_far')::int) share_cis_far, avg((cls='trans')::int) share_trans, sum(eligible::int) n_eligible, avg(eligible::int) share_eligible,
  avg(case when eligible then (nbr_obs>=2)::int end) share_block, avg(case when elig_null then (nbr_null>=2)::int end) share_block_null from rp2 group by 1""")
con.execute("create table fam as select genome, family, count(distinct g) ng from (select genome,family,g1 g from rp2 union all select genome,family,g2 from rp2) group by 1,2")
con.execute("""create table fm as select genome, avg((ng=2)::int) share_2gene_fam_recent, median(ng) median_family_size_recent, max(ng) max_family_size_recent, count(*) n_families_recent,
   sum(ng) sum_fam_genes from fam group by 1""")
con.execute("""create table fr as select genome, sum(ng) filter (where rk<=5)::double/sum(ng) top5_family_gene_share from (select genome, ng, row_number() over (partition by genome order by ng desc) rk from fam) group by 1""")
con.execute("create table gcount as select genome, count(distinct g) genes_recent from (select genome,g1 g from rp2 union select genome,g2 from rp2) group by 1")
G=con.sql("select * from gm join fm using(genome) join fr using(genome) join gcount using(genome)").df()
# blocks as units (eligible block pairs; run = same contig pair, index gaps <= K on both sides)
B=con.sql("select genome,c1,c2,i1,i2,dS,status,n1,n2,g1,g2 from rp2 where eligible and nbr_obs>=2 order by genome,c1,c2,i1,i2").df()
if len(B):
    g=B.groupby(["genome","c1","c2"],sort=False)
    gap1=g.i1.diff().abs(); gap2=g.i2.diff().abs(); new=(gap1.isna())|(gap1>K)|(gap2>K)
    B["blk"]=new.cumsum()
    blocks=B.groupby("blk").agg(genome=("genome","first"),c1=("c1","first"),c2=("c2","first"),n_pairs=("dS","size"),i1min=("i1","min"),i1max=("i1","max"),i2min=("i2","min"),i2max=("i2","max"),
          median_dS=("dS","median"),dS_q25=("dS",lambda x:x.quantile(.25)),dS_q75=("dS",lambda x:x.quantile(.75)),share_identical=("status",lambda x:(x=="identical").mean()),n1=("n1","first"),n2=("n2","first")).reset_index()
    blocks=blocks[blocks.n_pairs>=3].copy()
    blocks["span1"]=blocks.i1max-blocks.i1min+1; blocks["span2"]=blocks.i2max-blocks.i2min+1
    blocks["retention"]=blocks.n_pairs/np.minimum(blocks.span1,blocks.span2)          # pairs per gene of the shorter span
    blocks["coverage_smaller_contig"]=np.minimum(blocks.span1/blocks.n1, blocks.span2/blocks.n2)  # span of the block over the contig it sits on, the smaller fraction
    blocks.to_parquet(f"{E}/data/dupmode_blocks.parquet",index=False)
    bg=blocks.groupby("genome").agg(n_blocks=("blk","size"),n_block_pairs=("n_pairs","sum"),n_block_contig_pairs=("c1",lambda s:0),block_median_dS=("median_dS","median"),
        block_dS_mad=("median_dS",lambda x:float(np.median(np.abs(x-np.median(x))))),median_block_pairs=("n_pairs","median"),max_block_pairs=("n_pairs","max"),
        block_retention_median=("retention","median"),block_coverage_median=("coverage_smaller_contig","median"),block_share_identical=("share_identical","mean")).reset_index()
    cp=blocks.groupby("genome").apply(lambda d:d[["c1","c2"]].drop_duplicates().shape[0],include_groups=False).rename("n_block_contig_pairs2").reset_index()
    bg=bg.drop(columns="n_block_contig_pairs").merge(cp.rename(columns={"n_block_contig_pairs2":"n_block_contig_pairs"}),on="genome",how="left")
    bg["blocks_per_contig_pair"]=bg.n_blocks/bg.n_block_contig_pairs
    G=G.merge(bg,on="genome",how="left")
print("blocks:",len(blocks) if len(B) else 0,flush=True)
S=con.sql(f"select genome, species_prefix prefix from '{TB}/wgd.ks.summary.parquet'").df(); G=G.merge(S,on="genome",how="left")
gt=pd.read_csv(f"{E}/data/genome_table_final.tsv.gz",sep="\t",keep_default_na=False,na_values=[""],low_memory=False)[["genome","genes","phylum","class","order","genus","busco_duplicated_pct","genes_ratio_clade","hybrid","share_2gene_fam"]]
G=G.merge(gt,on="genome",how="left"); G["breadth"]=G.genes_recent/pd.to_numeric(G.genes,errors="coerce")
# Pfam composition (descriptive only: no significance claim, hits are not independent) and TE share among recent-pair genes
con.register("rgx",con.sql("select split_part(g,'_',1) prefix, g gene from (select g1 g from rp2 union select g2 from rp2)").df())
pf=con.sql(f"""select split_part(p.protein_id,'_',1) prefix, split_part(p.pfam_acc,'.',1) acc, count(distinct p.protein_id) n_all, count(distinct case when r.gene is not null then p.protein_id end) n_rec
   from '{TB}/pfam.parquet' p left join rgx r on r.gene=p.protein_id group by 1,2""").df()
def js(a,b):
    a=a/a.sum(); b=b/b.sum(); m=(a+b)/2
    f=lambda x,y: np.sum(np.where(x>0, x*np.log2(np.where(x>0,x,1)/np.where(y>0,y,1)),0)); return 0.5*f(a,m)+0.5*f(b,m)
rows=[]
for pre,g in pf.groupby("prefix"):
    nrec=int(g.n_rec.sum())
    if nrec<30: continue
    te=g[g.acc.isin(TE)]; rows.append((pre,js(g.n_all.values.astype(float),g.n_rec.values.astype(float)),te.n_rec.sum()/max(nrec,1),nrec))
J=pd.DataFrame(rows,columns=["prefix","js_pfam_recent_vs_proteome_descriptive","te_share_of_recent_domain_hits","recent_domain_hits"]); G=G.merge(J,on="prefix",how="left")
G["block_excess"]=G.share_block-G.share_block_null
G.to_csv(f"{E}/tables/dupmode_genome.tsv.gz",sep="\t",index=False)
print(len(G),"genomes;",G[["breadth","share_tandem","share_cis_far","share_trans","share_eligible","share_block","share_block_null","top5_family_gene_share","n_blocks","block_median_dS"]].describe().round(3).to_string())
