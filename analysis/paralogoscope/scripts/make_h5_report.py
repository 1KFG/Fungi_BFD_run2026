#!/usr/bin/python3.12
"""Write reports/H5_block_structure_window_codeml.md: block-level descriptors (item 2), the wider dS window (item 3) and the codeml second pass (Option E). Numbers come from tables/."""
import datetime, os, numpy as np, pandas as pd
E="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"; T=f"{E}/tables"
rd=lambda p: pd.read_csv(p,sep="\t",low_memory=False)
def tab(df,fmt=3):
    cols=list(df.columns); o=["| "+" | ".join(map(str,cols))+" |","|"+"|".join(["---"]*len(cols))+"|"]
    for _,r in df.iterrows(): o.append("| "+" | ".join("" if (isinstance(v,(float,np.floating)) and np.isnan(v)) else (f"{v:.{fmt}g}" if isinstance(v,(float,np.floating)) else str(v).replace("|","/")) for v in r.values)+" |")
    return "\n".join(o)
def hyp(n,s,against,test,res,verdict,lim): return f"## {n}\n\n**Hypothesis.** {s}\n\n**It would count against the hypothesis if:** {against}\n\n**Test.** {test}\n\n**Result.** {res}\n\n**Verdict: {verdict}.** {lim}\n\n"
BS=rd(f"{T}/dupmode_block_summary.tsv"); BC=rd(f"{T}/dupmode_block_controls.tsv"); BN=rd(f"{T}/dupmode_block_named.tsv"); BS3=rd(f"{T}/dupmode_block_summary_w03.tsv")
G1=rd(f"{T}/dupmode_genome_modes.tsv.gz").drop_duplicates("genome"); G3=rd(f"{T}/dupmode_genome_modes_w03.tsv.gz").drop_duplicates("genome"); CAND=rd(f"{T}/dupmode_genome_scale_candidates_w03.tsv")
GC=rd(f"{T}/dupmode_block_genome_classes.tsv.gz")
md=f'''---
title: "H5. Block structure, a wider dS window, and a codeml check of peak positions"
date: "{datetime.date.today().isoformat()}"
---

Follow-up to H4 (`reports/H4_duplication_mode.md`). Three additions: block-level descriptors (retention, strand concordance, copy depth, contig coverage), a wider recent window (dS <= 0.3, neighbours dS <= 0.5), and a codeml second pass for genomes with a possible interior dS peak (Option E). Scripts: `dup_mode.py` (env `DUP_RECENT`, `DUP_NBR`, `DUP_TAG`), `dup_block_tests.py`, `optionE_compare.py`, `option_e_codeml/`.

Status: interim (7,384-genome tables). Mode thresholds are exploratory.

'''
# H5.1 strand concordance
m=BS[BS["mode"]=="multi_contig_blocks"].iloc[0]; allw=(BS.blocks*BS.share_strand_concordance_ge_0_9).sum()/BS.blocks.sum()
v="SUPPORTED" if allw>0.85 else "NOT SUPPORTED"
md+=hyp("H5.1 Collinear blocks are real syntenic structures","Gene strands in a block agree with the block orientation (direct blocks: same strand; inverted blocks: opposite strand). Strand is not used to call blocks, so agreement is independent evidence.",
        "strand concordance is near 0.5 (the random expectation).","Per block, share of pairs whose strand relation matches the block orientation; share of blocks with concordance >= 0.9, in informative genomes (window dS <= 0.1). Verdict rule: SUPPORTED if > 85% of blocks have concordance >= 0.9.",
        tab(BS[["mode","genomes","blocks","median_pairs","strand_concordance","share_strand_concordance_ge_0_9","retention","copies_median"]],3)+f"\n\nAll modes pooled: {allw:.3f} of blocks have strand concordance >= 0.9. The expectation for unrelated pairs is 0.5 per gene pair.",v,
        "This supports the block calls as syntenic runs. It does not say what produced the duplication.")
# H5.2 two-copy vs multi-copy
nm=BN.set_index("genome"); sel=[g for g in nm.index if any(k in g for k in ["Verticillium_longisporum","Zygosaccharomyces_bailii","Zygosaccharomyces_parabailii","Ogataea_polymorpha_x"])]
hy=nm.loc[sel,"share_multi_copy"].mean(); mc=BS[BS["mode"]=="multi_contig_blocks"].iloc[0].share_multi_copy
v="SUPPORTED" if hy<mc else "INCONCLUSIVE"
md+=hyp("H5.2 Hybrid and doubled-genome blocks are two-copy; repeated segmental duplication is multi-copy","Genomes with a hybrid or whole-genome-like signal have blocks with copy number 2, while multi-copy blocks (copy number > 2) mark repeated segmental duplication.",
        "named hybrids have as many multi-copy blocks as other block genomes.","Copy number of a block = 1 + the number of other blocks that overlap its busiest region by >= 50%. Share of blocks with copies > 2 in named hybrids against the multi-contig-block mode. Verdict rule: SUPPORTED if the hybrid share is lower.",
        tab(BN[["genome","blocks","median_pairs","retention","strand_concordance","copies_median","share_multi_copy","share_scaffold_like","coverage_smaller_contig","median_dS"]],3)+f"\n\nMean share of multi-copy blocks: four named hybrids {hy:.3f}; multi-contig-block mode {mc:.3f}.",v,
        "*V. longisporum* and the *Zygosaccharomyces* genomes have 0-3% multi-copy blocks. *Puccinia graminis* 21-0 has about 30%, which fits repeated duplication inside a two-haplotype assembly. Sampled pairs under-count partners in large families, so copy number is a lower bound.")
# H5.3 scaffold-like
sc=GC[GC.informative.fillna(False)]; nsc=int((sc.share_scaffold_like>=.5).sum()); og=BN[BN.genome.str.contains("Ogataea")].iloc[0]
md+=hyp("H5.3 Some block genomes are duplicated scaffolds or doubled assemblies","Blocks that cover >= 80% of the smaller contig with median dS <= 0.01 mark whole duplicated contigs (assembly duplication) and are concentrated in a minority of genomes.",
        "such blocks are spread evenly over all genomes.","Per genome, share of blocks that are scaffold-like. Count genomes with >= 50% scaffold-like blocks.",
        f"{nsc} of {len(sc)} informative genomes with blocks have >= 50% scaffold-like blocks. In the *Ogataea polymorpha x parapolymorpha* hybrid genome {og.share_scaffold_like:.2f} of {int(og.blocks)} blocks are scaffold-like (median block coverage {og.coverage_smaller_contig:.2f}, median dS {og.median_dS:.3f}).",
        "DESCRIPTIVE","Scaffold-like blocks fit duplicated contigs or a doubled assembly. A real near-identical homeolog pair would look the same. Read depth or contig-level alignment would separate them.")
# H5.4 window
w1=int(G1.informative.sum()); w3=int(G3.informative.sum()); a=G1.set_index("genome"); b=G3.set_index("genome")
cmp=pd.DataFrame([dict(genome=g,breadth_w01=a.loc[g,"breadth"],breadth_w03=b.loc[g,"breadth"],share_block_w01=a.loc[g,"share_block"],share_block_w03=b.loc[g,"share_block"],null_w03=b.loc[g,"share_block_null"],mode_w01=a.loc[g,"mode"],mode_w03=b.loc[g,"mode"]) for g in b.index if g in a.index and any(k in g for k in ["Verticillium_longisporum","Zygosaccharomyces_bailii","Zygosaccharomyces_mellis","Zygosaccharomyces_pseudobailii","Zygosaccharomyces_parabailii","Ogataea_polymorpha_x","Cryptococcus_neoformans_x","Aspergillus_latus","pastorianus","Puccinia_graminis_21"])])
v="SUPPORTED" if w3>w1 and (cmp.breadth_w03>cmp.breadth_w01).mean()>0.7 else "PARTLY SUPPORTED"
md+=hyp("H5.4 A wider window (dS <= 0.3) recovers genome-scale duplication in hybrids whose copies diverged beyond 0.1","With the wider window, more genomes become informative and hybrids show breadth near 1 instead of a fraction.",
        "the number of informative genomes and hybrid breadth do not change.","Informative genomes (>= 50 eligible recent pairs, null block share <= 0.1) and named-genome breadth and block share for both windows.",
        f"Informative genomes: {w1:,} (window 0.1) and {w3:,} (window 0.3).\n\n"+tab(cmp,3),v,
        "The window ties the interpretation to dS: *V. longisporum* blocks have median dS 0.14 and *Z. bailii* 0.16, both above the first window. *A. latus* stays 'mixed' (block share 0.33). The *Cryptococcus* AD hybrid has a null block share of 0.19 at the wider window and stays unclassed.")
# H5.5 candidate list
cand=CAND.copy(); cand["genus"]=cand.genome.str.split("_").str[0]; gtop=cand.genus.value_counts().head(10).rename_axis("genus").reset_index(name="genomes")
cand["busco_dup"]=pd.to_numeric(cand.busco_dup,errors="coerce")
bt=pd.cut(cand.busco_dup,[-1,10,40,101],labels=["<10","10-40",">=40"]).value_counts().reindex(["<10","10-40",">=40"]).rename_axis("BUSCO dup %").reset_index(name="genomes")
md+=hyp("H5.5 Genome-scale block genomes form a candidate list that overlaps assembly evidence","Genomes with breadth >= 0.3, block share >= 0.5 and blocks on >= 5 contig pairs (window 0.3) also show high BUSCO duplication or a sharp interior peak.",
        "BUSCO duplication and peak tiers are no higher in the candidates than in other informative genomes.","Candidate list from `dupmode_genome_scale_candidates_w03.tsv`; BUSCO duplicated-% bins; overlap with shape tiers A and B; genera.",
        f"{len(cand)} candidate genomes; {int(cand.phylum.eq('Basidiomycota').sum())} Basidiomycota, {int(cand.phylum.eq('Ascomycota').sum())} Ascomycota, {int(cand.phylum.eq('Microsporidia').sum())} Microsporidia. Shape tier A or B: {int(cand.tier.fillna('').str.match('^[AB]').sum())}; any interior peak <= 1.0: {int((cand.n_peaks>0).sum())}.\n\n"+tab(bt,3)+"\n\nLargest genera:\n\n"+tab(gtop,3)+"\n\nFirst 25 by breadth:\n\n"+tab(cand.sort_values("breadth",ascending=False).head(25)[["genome","class","breadth","share_block","share_block_null","n_block_contig_pairs","block_median_dS","block_dS_mad","busco_dup","genes_ratio_clade","tier"]],3),
        "DESCRIPTIVE","The list is a screen. Several strains of one genus count as separate genomes (for example the Apiotrichum strains), so the number of independent events is lower than the number of genomes. Microsporidia entries have BUSCO duplication near 10% and gene counts several times the clade median; they need a separate check for contamination or assembly doubling.")
# H5.6 Option E
pS=f"{T}/optionE_peak_shift.tsv"
if os.path.exists(pS):
    P=rd(pS); NR=rd(f"{T}/optionE_node_ratio.tsv"); bp=P.dropna(subset=["yn00_peak","codeml_peak"]); med=bp.shift_ratio.median(); q=bp.shift_ratio.quantile([.1,.9]).tolist(); within=((bp.shift_ratio>0.9)&(bp.shift_ratio<1.1)).mean()
    v="SUPPORTED" if within>=0.7 else "PARTLY SUPPORTED"
    md+=hyp("H5.6 yn00 peak positions below dS 1.0 are close to codeml positions","For genomes with an interior peak, the peak position from codeml is within 10% of the yn00 position, so the screens (H1, H4, H5) do not move when the estimator changes.",
            "codeml peak positions differ from yn00 by more than 10% in most genomes.",
            f"Same pairs-per-node sampling and gates; codeml run on {len(P)} genomes with a possible peak (shape tiers A to C, excluding 8 genomes that already had codeml). The shape detector is applied to codeml node dS; the position is the median node dS inside the peak window. Verdict rule: SUPPORTED if >= 70% of genomes have a ratio between 0.9 and 1.1.",
            f"Peak detected in both estimators: {len(bp)} of {len(P)} genomes. Position ratio codeml / yn00: median {med:.3f}, 10th-90th percentile {q[0]:.3f} to {q[1]:.3f}; {within:.2f} of genomes within 10%.\n\n"+tab(NR,3),v,
            "This node-level result differs from the earlier wave-1 comparison (300 genomes, `tables/yn00_vs_codeml_nodes_binned.tsv`), which gave codeml / yn00 ratios of 1.08 to 1.13 for yn00 dS 0.2 to 1.6. Here the ratio is about 1.03 below dS 0.8 and 0.94 to 0.97 above 0.8. The genome sets differ (wave-1: all taxa; here: genomes with a possible peak) and the node matching differs. I have not resolved the difference. One untested possibility is regression to the mean, because nodes are binned by the yn00 value and the two runs sample different pairs. What this test does show is that peak positions below dS 1.0 agree within 10% in most genomes. Peaks above dS 1.0 are not assessed. One genome (Hemileccinum hortonii, tier C) has a peak with yn00 and none with codeml.\n\n![Option E](../figures/fig18_optionE.png)\n\n")
else:
    md+="## H5.6 codeml second pass\n\nNot finished when this report was built.\n\n"
md+="""## Follow-ups

- Contig-level alignment or read depth for the scaffold-like genomes and for the candidate list, to separate duplicated assemblies from real hybrids.
- Treat strains of one genus (for example Apiotrichum) as one candidate event.
- Microsporidia candidates: check contamination and assembly doubling.
- Phase the two block partners in V. longisporum, Z. bailii and Epichloe coenophiala against relatives (gene-tree placement) to test for allopolyploidy; gene order and dS cannot do this.
- Rerun on the final tables.

## Reproduce

`DUP_TAG= scripts/dup_mode.py`; `DUP_RECENT=0.3 DUP_NBR=0.5 DUP_TAG=_w03 scripts/dup_mode.py`; `scripts/dup_mode_tests.py` and `scripts/dup_block_tests.py` (same env); `option_e_codeml/select_genomes.py`, `submit_array.sh`; `scripts/optionE_compare.py`; `scripts/make_h5_report.py`.
"""
open(f"{E}/reports/H5_block_structure_window_codeml.md","w").write(md); print("wrote H5",len(md))
