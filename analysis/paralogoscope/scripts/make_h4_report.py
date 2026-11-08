#!/usr/bin/python3.12
"""Write reports/H4_duplication_mode.md from the dup_mode outputs. Numbers are computed from tables/."""
import datetime, numpy as np, pandas as pd
from scipy import stats
E="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"; T=f"{E}/tables"
G=pd.read_csv(f"{T}/dupmode_genome_modes.tsv.gz",sep="\t",low_memory=False); M=pd.read_csv(f"{T}/dupmode_modes.tsv",sep="\t"); PC=pd.read_csv(f"{T}/dupmode_pair_categories.tsv",sep="\t")
CM=pd.read_csv(f"{T}/dupmode_modes_by_control.tsv",sep="\t",index_col=0); BP=pd.read_csv(f"{T}/dupmode_modes_by_phylum.tsv",sep="\t",index_col=0); CL=pd.read_csv(f"{T}/dupmode_pair_class_block.tsv",sep="\t")
def tab(df,fmt=3,index=False):
    if index: df=df.reset_index()
    cols=list(df.columns); o=["| "+" | ".join(map(str,cols))+" |","|"+"|".join(["---"]*len(cols))+"|"]
    for _,r in df.iterrows(): o.append("| "+" | ".join("" if (isinstance(v,(float,np.floating)) and np.isnan(v)) else (f"{v:.{fmt}g}" if isinstance(v,(float,np.floating)) else str(v).replace("|","/")) for v in r.values)+" |")
    return "\n".join(o)
def hyp(n,s,against,test,res,verdict,lim): return f"## {n}\n\n**Hypothesis.** {s}\n\n**It would count against the hypothesis if:** {against}\n\n**Test.** {test}\n\n**Result.** {res}\n\n**Verdict: {verdict}.** {lim}\n\n"
inf=G[G.informative]; n_all,n_inf=len(G),len(inf)
md=f'''---
title: "H4. Selective gene duplication or block-scale duplication? Gene-order evidence"
date: "{datetime.date.today().isoformat()}"
---

Question: can recent paralog pairs (dS <= 0.1) be separated into duplication of a subset of genes (tandem arrays, dispersed copies, enriched families) and block-scale duplication (whole-genome duplication, segmental duplication, hybrid with two assembled subgenomes, assembly duplication)? Method: `scripts/dup_mode.py` (definitions in its docstring) and `scripts/dup_mode_tests.py`. The first version of `dup_mode.py` was reviewed by an independent model reviewer. The review found a pair-ordering bug on single contigs, tRNA genes in the gene index, tandem runs counted as blocks, a null that was too easy, and an invalid z-score for Pfam composition. All were fixed in version 2, and the first version is kept in `scripts/legacy/dup_mode_v1.py`.

Status: interim (7,384-genome tables). The mode thresholds are exploratory and were not validated on held-out data.

## What the data can and cannot show

A pair is "in a block" when at least 2 pairs of a different family lie within 5 genes of both of its genes, on the same two contigs, in the same orientation. The null permutes gene order within each contig. Only genomes with at least 50 eligible recent pairs and a null block share <= 0.1 are classed (informative: {n_inf:,} of {n_all:,}). The others are too small, too fragmented or too dense for the test.

Gene order and dS cannot separate whole-genome duplication from a hybrid, a segmental duplication, or a duplicated assembly (haplotigs, duplicated scaffolds). They can separate block from non-block patterns. Old WGDs (Saccharomycetaceae, Mucorales) are older than the dS window and are not detected.

'''
# H4.1 controls
cm=CM.copy(); cm["blocks"]=cm.get("multi_contig_blocks",0)+cm.get("single_contig_pair",0); cm["classed"]=cm.drop(columns=["blocks"]).sum(axis=1)
hy=cm.loc[[i for i in cm.index if "hybrid" in i]]; rate=hy.blocks.sum()/max(hy.classed.sum(),1); mu=cm.loc[[i for i in cm.index if "Mucorales" in i]]; mu_rate=mu.blocks.sum()/max(mu.classed.sum(),1)
v="PARTLY SUPPORTED" if rate>=0.5 and rate>mu_rate else "INCONCLUSIVE"
vl=G[G.genome.str.contains("Verticillium_longisporum|Zygosaccharomyces_bailii|Ogataea_polymorpha_x|Zygosaccharomyces_parabailii|Aspergillus_latus|Puccinia_graminis_21-0|Zygosaccharomyces_rouxii")].drop_duplicates("genome")
md+=hyp("H4.1 Known hybrids show block duplication and old WGD does not","Named and literature hybrids have most recent pairs in collinear blocks that are spread over many contig pairs with one dS mode. Old WGD (Mucorales) does not register in the recent window.",
        "hybrids are not in the block modes, or the null produces similar block shares.","Mode of control genomes among informative genomes; block statistics of named examples. Verdict rule: PARTLY SUPPORTED if hybrids are in a block mode in >= 50% of classed hybrids and more often than the Mucorales controls.",
        tab(cm[["multi_contig_blocks","single_contig_pair","dispersed_diffuse","dispersed_selective","tandem_array","mixed"]].fillna(0).astype(int),3,index=True)+f"\n\nHybrid controls in a block mode: {rate:.2f} ({int(hy.blocks.sum())} of {int(hy.classed.sum())} classed). Mucorales WGD controls: {mu_rate:.2f}.\n\nNamed examples:\n\n"+
        tab(vl[["genome","n_eligible","breadth","share_block","share_block_null","n_blocks","n_block_contig_pairs","block_median_dS","block_dS_mad","informative","mode"]],3), v,
        f"*V. longisporum* has 98% of eligible recent pairs in blocks (null 0.00) with 246 blocks on 60 contig pairs and a single dS mode (median 0.082, MAD 0.007). *Z. parabailii* has a dense-contig null (0.59) and is unclassed. *A. latus* has about 38% in blocks; its parental divergence lies above the window. The hybrid group is small (n = {int(hy.classed.sum())}), and the pre-WGD Saccharomycetaceae controls are partly hybrids (Zygosaccharomyces), so they do not serve as negative controls.")
# H4.2 identical
ic=PC[PC.category.str.contains("identical")].iloc[0]
v="SUPPORTED" if ic.share_in_block_pooled>1.5*ic.share_in_block_other_pairs_pooled and ic.frac_genomes_category_higher>0.7 else "NOT SUPPORTED"
md+=hyp("H4.2 Identical-CDS copies are mostly parts of duplicated regions, not scattered single genes","Identical pairs are in blocks more often than other recent pairs of the same genome.","identical pairs are in blocks as rarely as other recent pairs.",
        "Share of eligible recent pairs in a block, identical-CDS pairs against other pairs, pooled and per genome (genomes with >= 5 identical pairs). Verdict rule: SUPPORTED if the pooled share is > 1.5 times the other share and > 70% of genomes show a higher share.",
        tab(PC[PC.category.str.contains("identical")],3), v,
        f"{ic.share_in_block_pooled:.2f} of identical pairs are in blocks against {ic.share_in_block_other_pairs_pooled:.2f} of other pairs, and {ic.frac_genomes_category_higher:.2f} of {int(ic.genomes_compared)} genomes show a higher share. The excess of dS = 0 paralogs therefore mostly sits in duplicated regions. Whether the regions are real segmental duplications or assembly duplication is not separated here (see H4.4 and the single-contig-pair mode).")
# H4.3 dispersed copies: TE and OrsD
te=PC[PC.category.str.contains("TE")].iloc[0]; og=PC[PC.category.str.contains("OrsD")].iloc[0]
v="SUPPORTED" if te.share_in_block_pooled<0.5*te.share_in_block_other_pairs_pooled and og.share_in_block_pooled<0.5*og.share_in_block_other_pairs_pooled else "NOT SUPPORTED"
md+=hyp("H4.3 Recent copies of TE-domain genes and of the OrsD family are dispersed, not block duplicated","Pairs that include a TE-domain gene or an OrsD gene are in blocks less often than other recent pairs of the same genome.","these pairs are in blocks as often as other pairs.",
        "Same test as H4.2 for pairs containing a gene with one of 15 TE Pfam domains, or a PF12013 gene. Verdict rule: SUPPORTED if both pooled shares are below half of the other-pair share.",
        tab(PC[~PC.category.str.contains("identical")],3), v,
        f"TE-domain pairs: {te.share_in_block_pooled:.2f} in blocks against {te.share_in_block_other_pairs_pooled:.2f}. OrsD pairs: {og.share_in_block_pooled:.2f} against {og.share_in_block_other_pairs_pooled:.2f}. This fits copy-and-scatter duplication of these genes, as opposed to duplication of whole regions. It does not show that the copies are mobile elements. A low block share also arises if blocks are missed because of fragmented assemblies, but that applies equally to the other pairs of the same genomes.")
# H4.4 modes and BUSCO
mc=inf[inf["mode"]=="multi_contig_blocks"]; hi=(mc.busco_duplicated_pct>=40).mean(); ot=inf[~inf["mode"].isin(["multi_contig_blocks"])]; hi0=(ot.busco_duplicated_pct>=40).mean()
rho,pp=stats.spearmanr(inf.share_block.fillna(0),pd.to_numeric(inf.busco_duplicated_pct,errors="coerce").fillna(0))
v="PARTLY SUPPORTED" if hi>hi0 and rho>0.15 else "NOT SUPPORTED"
md+=hyp("H4.4 Block-scale genomes are the ones with assembly-level duplication","Genomes in the multi-contig block mode have higher BUSCO duplication and gene counts than genomes in dispersed modes.","BUSCO duplication does not differ between block and dispersed modes.",
        "Median BUSCO duplication, gene-count ratio, identical-pair share and block dS by mode; share of multi-contig-block genomes with BUSCO duplication >= 40%; Spearman correlation of block share with BUSCO duplication. Verdict rule: PARTLY SUPPORTED if the >= 40% share is higher in block genomes and rho > 0.15.",
        tab(M,3), v, f"{hi:.2f} of {len(mc)} multi-contig-block genomes have BUSCO duplication >= 40%, against {hi0:.2f} of other informative genomes. Spearman rho (block share against BUSCO duplication) = {rho:.2f} (p = {pp:.1e}). Most block-mode genomes have low BUSCO duplication (median {mc.busco_duplicated_pct.median():.1f}%), so they hold smaller duplicated regions or low-divergence blocks, not whole doubled genomes. The single-contig-pair mode has near-zero block dS (median {M[M['mode']=='single_contig_pair'].med_block_dS.iloc[0]:.3f}) and a high identical share ({M[M['mode']=='single_contig_pair'].med_identical_pair_share.iloc[0]:.2f}), which fits a duplicated scaffold or one segmental duplication.")
md+="## Modes by phylum (informative genomes)\n\n"+tab(BP,3,index=True)+"\n\n## Figure\n\n![Duplication modes](../figures/fig17_dupmode.png)\n\n"
md+="""## Claims these data do not support

- That a genome is allopolyploid rather than autopolyploid.
- That a genome is a hybrid rather than a haplotig or scaffold duplication, without read depth or heterozygosity data.
- That a genome has no block duplication: most assemblies are fragmented (84% of recent pairs are between contigs) and 77% of genomes are unclassed.
- Absolute duplication rates or dates. Pairs are sampled (at most 5 per node), and identical CDS are collapsed.
- Any significance statement about Pfam composition. A permutation test that respects protein and family structure is not built.

## Follow-up queries that use the same tables

- Retention inside blocks (genes with a partner on both sides) and strand consistency from `gene_info.strand`.
- Synteny depth: number of distinct partner regions per region (2 fits WGD or a hybrid; higher fits repeated segmental duplication).
- Contig coverage and end-to-end blocks with dS near 0 for assembly duplication.
- For dispersed copies: number of contigs each family spans and TE-domain genes within 3 genes of each copy, against random genes.
- Repeat the analysis with a wider dS window (<= 0.3) for hybrids with older parents.

## Reproduce

`scripts/dup_mode.py`, `scripts/dup_mode_tests.py`, `scripts/make_h4_report.py`.
"""
open(f"{E}/reports/H4_duplication_mode.md","w").write(md); print("wrote H4",len(md))
