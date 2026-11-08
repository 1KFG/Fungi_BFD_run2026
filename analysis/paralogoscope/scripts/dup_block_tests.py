#!/usr/bin/python3.12
"""Block-level descriptors by provisional mode and control group (uses dup_mode.py block table and dup_mode_tests.py modes).
Per block: retention (share of genes in each span that have a sampled partner; the smaller of the two sides), strand concordance (direct blocks: gene strands agree; inverted: disagree),
copies (1 + number of other blocks overlapping the busiest region by >= 50%; 2 = a two-copy block), coverage of the smaller contig, end-to-end flag (coverage >= 0.8), median dS.
Block classes: scaffold_like (end-to-end and median dS <= 0.01), multi_copy (copies > 2), two_copy (otherwise).
Expectations (from the review): WGD or hybrid -> copies 2, one dS mode across blocks; segmental repeats -> copies > 2; duplicated scaffold or haplotig -> end-to-end, dS ~ 0.
Usage: DUP_TAG=_w03 dup_block_tests.py   Outputs tables/dupmode_block_summary{TAG}.tsv, dupmode_block_controls{TAG}.tsv"""
import os, numpy as np, pandas as pd
E="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"; TAG=os.environ.get("DUP_TAG","")
B=pd.read_parquet(f"{E}/data/dupmode_blocks{TAG}.parquet"); G=pd.read_csv(f"{E}/tables/dupmode_genome_modes{TAG}.tsv.gz",sep="\t",low_memory=False).drop_duplicates("genome")
B=B.merge(G[["genome","mode","control_group","informative","busco_duplicated_pct"]],on="genome",how="left")
B["cls"]=np.where(B.end_to_end&(B.median_dS<=0.01),"scaffold_like",np.where(B.copies>2,"multi_copy","two_copy"))
print("blocks:",len(B),"in informative genomes:",int(B.informative.fillna(False).sum()))
def summ(d):
    return pd.Series(dict(genomes=(d["genome"].nunique() if "genome" in d else 1),blocks=len(d),median_pairs=d.n_pairs.median(),retention=d.retention.median(),retention_asym=d.retention_asym.median(),strand_concordance=d.strand_concordance.median(),
        share_strand_concordance_ge_0_9=(d.strand_concordance>=0.9).mean(),copies_median=d.copies.median(),share_multi_copy=(d.cls=="multi_copy").mean(),share_scaffold_like=(d.cls=="scaffold_like").mean(),
        coverage_smaller_contig=d.coverage_smaller_contig.median(),median_dS=d.median_dS.median()))
I=B[B.informative.fillna(False)]
S=I.groupby("mode").apply(summ,include_groups=False).reset_index(); S.to_csv(f"{E}/tables/dupmode_block_summary{TAG}.tsv",sep="\t",index=False); print(S.round(3).to_string())
Cc=B[B.control_group.notna() & (B.control_group!="other")].groupby("control_group").apply(summ,include_groups=False).reset_index(); Cc.to_csv(f"{E}/tables/dupmode_block_controls{TAG}.tsv",sep="\t",index=False); print(Cc.round(3).to_string())
names="Verticillium_longisporum|Zygosaccharomyces_bailii|Zygosaccharomyces_parabailii|Ogataea_polymorpha_x|Puccinia_graminis_21|Aspergillus_latus|Cryptococcus_neoformans_x|Saccharomyces_pastorianus"
K=B[B.genome.str.contains(names)].groupby("genome").apply(summ,include_groups=False).reset_index(); K.to_csv(f"{E}/tables/dupmode_block_named{TAG}.tsv",sep="\t",index=False); print(K.round(3).to_string())
# genome-level class shares
GS=B.groupby("genome").apply(lambda d:pd.Series(dict(blocks=len(d),share_multi_copy=(d.cls=="multi_copy").mean(),share_scaffold_like=(d.cls=="scaffold_like").mean(),share_two_copy=(d.cls=="two_copy").mean())),include_groups=False).reset_index().merge(G[["genome","mode","phylum","class","busco_duplicated_pct","informative"]],on="genome")
GS.to_csv(f"{E}/tables/dupmode_block_genome_classes{TAG}.tsv.gz",sep="\t",index=False)
x=GS[GS.informative.fillna(False)]; print("informative genomes with blocks:",len(x),"| mostly scaffold_like (>=0.5):",int((x.share_scaffold_like>=.5).sum()),"| mostly multi_copy:",int((x.share_multi_copy>=.5).sum()),"| mostly two_copy:",int((x.share_two_copy>=.5).sum()))
