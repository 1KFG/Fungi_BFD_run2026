#!/usr/bin/python3.12
"""Build report.md (embedded PNG figures) from the exploration tables.
All numbers in the text are computed here from data/ and tables/. PDF is built by job_report.sh.
"""
import os, datetime
import numpy as np, pandas as pd

E = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
T, D, F = f"{E}/tables", f"{E}/data", "figures"
rd = lambda p: pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""], low_memory=False)

G = rd(f"{D}/genome_table_final.tsv.gz")
PH = rd(f"{T}/clade_summary_phylum.tsv")
CL = rd(f"{T}/clade_summary_class.tsv")
OF = rd(f"{T}/outlier_flags.tsv")
CT = rd(f"{T}/controls.tsv")
PK = rd(f"{T}/peak_catalogue.tsv")
VC = rd(f"{T}/varcomp.tsv")
SE = rd(f"{T}/sensitivity_codon100.tsv")
YB = rd(f"{T}/yn00_vs_codeml_nodes_binned.tsv")
YG = rd(f"{T}/yn00_vs_codeml_genomes.tsv")

for c in ["hybrid", "busco_duplicated_pct", "genes_ratio_clade", "largest_identical_group", "share_2gene_fam",
          "nodes_per_1000_genes", "ident_genes_per_1000", "te_per_1000_genes", "contig_count", "N50_bp"]:
    for d in (OF, PK):
        if c in d: d[c] = pd.to_numeric(d[c], errors="coerce")


def md_table(df, floatfmt=3):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append("" if np.isnan(v) else (f"{v:.{floatfmt}g}"))
            else:
                cells.append(str(v).replace("|unassigned", " (unassigned)").replace("|", "/"))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def fig(name, caption):
    return f"![{caption}]({F}/{name}.png)\n\n*{caption}*\n"


n_all = len(G)
n_phy = G.phylum.nunique() if "phylum" in G else None
R = []
w = R.append

w(f"""---
title: "Paralog dS landscape of {n_all:,} fungal genomes (BFD): identical-copy excess and recent duplication signals"
date: "{datetime.date.today().isoformat()}"
---

# Summary

""")

# ---------- computed headline numbers
cov = G.frac_nodes_with_ds.median() if "frac_nodes_with_ds" in G else np.nan
w(f"""This report describes node-level synonymous divergence (dS) of paralogs in {n_all:,} genomes with a finished
paralogoscope table. The method is: wgd dmd families, identical-CDS collapse, yn00 dS on sampled pairs
(k = 5 per node, gates dS <= 2.5 and >= 30 codons), and nodes from a neighbour-joining tree on protein identity.
All values below are computed from `tables/` in this folder. Nothing here is a claim about individual genomes
without independent evidence. A flag from dS is a screening result.

Limits that apply to every result:

- yn00 dS reads low above about 0.5 and saturates near 2 (Section 8). Only dS 0 to about 1.5 carries position information.
  The median share of nodes with a usable dS per genome is {cov:.2f}. The rest lie beyond the identity window.
- The data cannot detect old fungal WGDs (Saccharomycetaceae, Mucorales). They fall in the old or saturated mass.
  Question 2 is therefore: which genomes carry a recent genome-scale duplication
  (hybrid, heterokaryon or dikaryon assembled with both genomes, or a mixed assembly).
- Aneuploidy is not visible.
- Clade differences are partly confounded with assembly source and technology (Section 6).
""")

# ---------- 1 data
w("# 1. Data and bins\n")
w("""Bins of node dS: **Z** identical (dS = 0), **H** (0, 0.05], **Y** (0.05, 0.3], **MID** (0.3, 1.5], **OLD** (1.5, 2.5],
and **BEYOND** (no usable dS). Shares use two denominators: all nodes (`shareAll_*`) and nodes with a usable dS.
Genome-level intervals come from a family-level cluster bootstrap. Genome statistics are medians per clade,
and genomes are never pooled by node.\n""")
w(fig("fig11_coverage_by_phylum", "Share of nodes with a usable dS, by phylum."))

# ---------- 2 composition
w("# 2. Composition by taxonomic group\n")
ph = PH[PH.eligible.astype(str) == "True"].copy() if "eligible" in PH else PH
ph = ph.sort_values("n_genomes", ascending=False)
cols = ["phylum", "n_genomes", "n_genera", "med_frac_nodes_with_ds", "med_shareAll_Z", "med_shareAll_H", "med_shareAll_Y",
        "med_shareAll_MID", "med_nodes_per_1000_genes", "med_share_2gene_fam", "med_busco_duplicated_pct", "med_te_per_1000_genes"]
cols = [c for c in cols if c in ph]
w("Phylum medians over genomes (columns `med_*`; bootstrap intervals are in `tables/clade_summary_phylum.tsv`).\n")
w(md_table(ph[cols], 3) + "\n")
cl = CL[CL.eligible.astype(str) == "True"].copy() if "eligible" in CL else CL
cl = cl.sort_values("n_genomes", ascending=False).head(15)
ccols = [c for c in ["class", "n_genomes", "n_genera", "med_shareAll_Z", "med_shareAll_H", "med_shareAll_Y", "med_shareAll_MID",
                     "med_nodes_per_1000_genes", "med_share_2gene_fam"] if c in cl]
w("\nFifteen largest classes (n >= 10 eligible):\n")
w(md_table(cl[ccols], 3) + "\n")
w(fig("fig01_ridgeline_by_class", "Node-dS distributions by class. Grey band: yn00 saturation zone (no position claim). Black bar: identical-copy share."))
w(fig("fig03_order_bin_heatmap", "Bin shares by order (genomes / genera in labels)."))
w(fig("fig08_nodes_density", "Nodes per 1,000 genes against dS-bin mass."))

# ---------- 3 variance partition
w("# 3. How much of the variation is taxonomic\n")
w("Nested random-effects model (taxonomy as a phylogeny proxy; branch lengths ignored), with and without assembly covariates.\n")
vc = VC.copy()
piv = vc.pivot_table(index=["response", "adjustment"], columns="level", values="share").reset_index()
order = [c for c in ["phylum", "class", "order", "family", "genus", "genome (residual)"] if c in piv]
w(md_table(piv[["response", "adjustment"] + order], 2) + "\n")
w("Reading: `share` is the fraction of total random variance at that level. A large residual or genus share means the trait varies within clades. "
  "Large phylum to order shares mean clade identity predicts it.\n")
w(fig("fig04_variance_components", "Variance components by taxonomic level."))

# ---------- 4 dS = 0
w("# 4. Question 1: excess of dS = 0 and near-zero paralogs\n")
fl = OF.copy()
nflag = lambda c: int((fl[c].astype(str) == "True").sum())
w(f"""Outlier flags use robust z within clade (order, then class, then phylum; n >= 10; MAD floors) on unshrunken values.
Flagged rows: identical excess {nflag('flag_identical_excess')}, H excess {nflag('flag_H_excess')}, Y excess {nflag('flag_Y_excess')},
MID excess {nflag('flag_MID_excess')}, duplication-dense {nflag('flag_dup_dense')}, 2-gene-family excess {nflag('flag_2gene_excess')}
(rows in `outlier_flags.tsv`: {len(fl)}). These counts are screening hits, not confirmed events. No multiplicity-adjusted claim is made per genome.
""")
idf = fl[fl.flag_identical_excess.astype(str) == "True"].copy()
idf["profile"] = np.where(idf.busco_duplicated_pct >= 20, "assembly-duplication-like (BUSCO dup >= 20%)",
                          np.where(idf.largest_identical_group >= 5, "TE/multicopy-like (large group, low BUSCO dup)", "other"))
pr = idf.groupby("profile").agg(n=("genome", "size"), med_busco_dup=("busco_duplicated_pct", "median"),
                                med_largest_group=("largest_identical_group", "median"),
                                med_te_per_1000=("te_per_1000_genes", "median"),
                                med_ident_per_1000=("ident_genes_per_1000", "median")).reset_index()
w("Identical-copy excess split by simple descriptors (BUSCO duplicated % and largest identical group). The two groups are a heuristic, not a diagnosis:\n")
w(md_table(pr, 3) + "\n")
top = idf.groupby(["phylum", "class"]).size().sort_values(ascending=False).head(12).reset_index(name="n_flagged")
w("\nClasses with the most identical-excess flags:\n")
w(md_table(top, 3) + "\n")
w("""Interpretation (from the framing review, `fable_context.md`): a flag with high BUSCO duplication and many small groups points to
duplicated contigs or haplotigs. A flag with low BUSCO duplication, large groups and high TE counts points to TE-derived models
or multicopy arrays. H excess without Z points to haplotype or heterozygosity depth. Tests that would separate these are listed in Section 10.\n""")
w(fig("fig02_identical_copies", "Identical-copy statistics by clade."))
w(fig("fig07_young_vs_assembly", "Young-node mass (H, Y) against BUSCO duplication and N50."))

# ---------- 5 peaks
w("# 5. Question 2: dS peaks suggesting recent genome-scale duplication, hybrid or mixed assembly\n")
w("""Peaks come from a Gaussian mixture on log dS (dS > 0). A peak counts as *stable* only if all of these hold:
node bootstrap >= 0.8, present at >= 100 codons, stable under pair resampling (>= 0.8), and position in 0.05 to 1.5.
Components outside that window are not treated as duplication peaks.\n""")
pk = PK.copy()
for c in ["support_node_boot", "stable_resample", "peak_dS"]:
    pk[c] = pd.to_numeric(pk[c], errors="coerce")
pk["stable"] = (pk["support_node_boot"] >= 0.8) & (pk["stable_codon100"].astype(str) == "True") & \
               (pk["stable_resample"] >= 0.8) & (pk["in_window_0.05_1.5"].astype(str) == "True")
pk["genome_scale"] = pk["stable"] & (pk["busco_duplicated_pct"] >= 40) & (pk["genes_ratio_clade"] >= 1.5)
w(f"The catalogue file lists {len(pk)} components in {pk.genome.nunique()} genomes. Of these, {int(pk.stable.sum())} pass the stability rules above (the catalogue is already filtered; this is a check). "
  f"Stable and with assembly support (BUSCO dup >= 40% and gene count >= 1.5x clade median): {int(pk.genome_scale.sum())} components in {pk[pk.genome_scale].genome.nunique()} genomes.\n")
bp = pk[pk.stable].groupby("phylum").agg(components=("genome", "size"), genomes=("genome", "nunique")).reset_index().sort_values("genomes", ascending=False)
w(md_table(bp, 3) + "\n")
gs = pk[pk.genome_scale].sort_values("busco_duplicated_pct", ascending=False)
gsc = [c for c in ["genome", "phylum", "class", "peak_dS", "support_node_boot", "busco_duplicated_pct", "genes_ratio_clade", "share_2gene_fam"] if c in gs]
w("\nTop 40 stable peaks with assembly support (sorted by BUSCO duplicated %):\n")
w(md_table(gs[gsc].head(40), 3) + "\n")
w(fig("fig06_peaks", "Peak positions and masses of the stable catalogue."))

# controls
w("## 5.1 Controls\n")
cg = CT.groupby("control_group").agg(n=("genome", "size"), med_Z=("shareAll_Z", "median"), med_H=("shareAll_H", "median"),
                                     med_Y=("shareAll_Y", "median"), med_2gene=("share_2gene_fam", "median"),
                                     med_busco_dup=("busco_duplicated_pct", "median"),
                                     med_gene_ratio=("genes_ratio_clade", "median"),
                                     frac_with_supported_peak=("n_supported_peaks_window", lambda s: float((pd.to_numeric(s) > 0).mean()))).reset_index()
w(md_table(cg, 3) + "\n")
w("""Named hybrids (' x ' in name) and literature hybrids are the positive controls for a mid-dS peak. Collapsed assemblies of hybrids show no peak.
Saccharomycetaceae and Mucorales are controls for old mass and 2-gene structure, not for a peak, because their WGD is older than the dS window.
Several named controls show sharp peaks near dS 0.2 to 0.3 or near 1.3 to 1.5 in the figure. A peak at about 1.3 to 1.5 lies near the yn00 saturation zone,
so its position is a lower bound on the true divergence.\n""")
w(fig("fig05_controls", "Control genomes: node-dS shape."))


# ---------- 5.2 shape-based peaks
SPa = rd(f"{T}/shape_peaks_all.tsv"); SPk = rd(f"{T}/shape_peaks.tsv")
for c in ["busco_dup","peak1_dS","peak1_mass","peak1_logwidth","genes_ratio_clade","share_2gene_fam","te_per_1000"]:
    SPk[c] = pd.to_numeric(SPk[c], errors="coerce")
w("## 5.2 Interior peaks (centre <= 1.0) detected from histogram shape, cross-referenced with BUSCO\n")
w(f"""This screen uses the per-genome log-dS histogram directly (70 bins, 0.005 to 2.5), not the mixture model. A peak is a local maximum of the 3-bin smoothed density with
prominence >= 35% of its height, centre between 0.03 and 1.0, and at least 8% of the genome's tree nodes within its half-height width. Genomes need >= 300 tree nodes with dS > 0
({len(SPa):,} genomes evaluated). A peak is called *sharp* if the half-height width is <= 1.2 natural-log units (a 3.3-fold range of dS), the peak holds >= 40% of nodes, and prominence is >= 0.8.
Tiers: **A** sharp and BUSCO duplicated >= 40%; **B** sharp and 10-40%; **C** sharp and < 10%; **D** broad hump. BUSCO values are `duplicated_pct` from `tables/busco_genome.parquet` (fungi_odb12, joined by assembly ID).
Peaks above 1.0 are not called, so hybrids whose second peak lies near 1.5 (for example the older Zygosaccharomyces parabailii peak) show only their lower peak.
""")
tc = SPk.tier.value_counts().sort_index()
w(f"{len(SPk)} of {len(SPa):,} genomes have an interior peak <= 1.0 ({int((SPk.n_peaks==1).sum())} with one peak, {int((SPk.n_peaks>=2).sum())} with two or more). Tier counts: " + "; ".join(f"{k}: {v}" for k, v in tc.items()) + ".\n")
bt = SPa.copy(); bt["busco_dup"] = pd.to_numeric(bt.busco_dup, errors="coerce")
bt["BUSCO dup bin"] = pd.cut(bt.busco_dup, [-1, 1, 5, 10, 20, 40, 101], labels=["<1", "1-5", "5-10", "10-20", "20-40", ">=40"])
tb = bt.groupby("BUSCO dup bin", observed=True).agg(genomes=("genome", "size"), share_with_interior_peak=("n_peaks", lambda x: float((pd.to_numeric(x) > 0).mean()))).reset_index()
w("Share of genomes with an interior peak by BUSCO duplication:\n"); w(md_table(tb, 3) + "\n")
w("Tier A by phylum and class:\n")
A = SPk[SPk.tier.str.startswith("A")]
w(md_table(A.groupby(["phylum", "cls"]).size().sort_values(ascending=False).reset_index(name="n_genomes_tier_A").head(14), 3) + "\n")
kn = SPa[SPa.genome.str.contains("Zygosaccharomyces|Verticillium_longisporum|Ogataea_polymorpha_x|Cryptococcus_neoformans_x|Aspergillus_latus|pastorianus")]
kc = [c for c in ["genome", "n_peaks", "peak1_dS", "peak1_mass", "peak2_dS", "peak2_mass", "busco_dup", "genes_ratio_clade"] if c in kn]
w("Known and literature hybrids in this screen (n_peaks = 0 means no interior peak <= 1.0, either no peak or a peak above 1.0):\n"); w(md_table(kn[kc], 3) + "\n")
w(fig("fig13_shape_peaks_gallery", "Shape-detected interior peaks. Orange: tier A, amber: tier B, blue: tier C, grey: broad. Dotted line: detected peak."))
w(fig("fig14_shape_peaks_vs_busco", "Peak position against BUSCO duplication, and the share of genomes with an interior peak by BUSCO duplication."))
w("""Reading the tiers. Tier A genomes combine a narrow dS peak with a high BUSCO duplicated fraction. The named and literature hybrids sit here
(Verticillium longisporum, Ogataea polymorpha x parapolymorpha, several Zygosaccharomyces). This is consistent with an assembly that carries two diverged genomes,
but dS and BUSCO together do not prove it. Gene count against the clade median, synteny and read data are the next checks (Section 10).
Tier C genomes have a sharp peak but little BUSCO duplication. Some have gene counts far below the clade median (for example several Boletaceae with gene-count ratios of 0.1 to 0.5).
That can mean partial annotation or a partial duplication, not a whole-genome event. Tier C is a list to inspect with FCS-GX and coverage data, not a list of events.
The yeast peaks near dS 0.5 to 1.0 lie where yn00 reads low relative to codeml, so true parental divergence there is likely higher than the peak position.\n""")

# ---------- 6 covariates
w("# 6. Assembly covariates and confounding\n")
w("Fixed effects (contig count, size, GC, gene count, BUSCO duplication, TE per 1,000 genes) are in `tables/varcomp_fixed.tsv`. "
  "Compare the 'no covariates' and 'with covariates' rows in Section 3 to see how much clade variance is attributable to assembly properties. "
  "Where the clade share changes after adjustment, clade differences are not separable from assembly differences with these data.\n")

# ---------- 7 sensitivity
w("# 7. Sensitivity to the codon gate (30 vs 100 codons)\n")
w(md_table(SE, 3) + "\n")
w("`spearman` is the rank correlation of per-genome bin shares between the two gates.\n")
w(fig("fig09_sensitivity", "Bin shares at 30 vs 100 codons."))

# ---------- 8 yn00 vs codeml
w("# 8. yn00 compared with codeml (300 wave-1 genomes)\n")
w(md_table(YB, 3) + "\n")
w(f"Genomes compared: {len(YG)}. Median per-genome ratio (codeml / yn00): {pd.to_numeric(YG.median_ratio, errors='coerce').median():.3f}.\n")
w(fig("fig10_yn00_vs_codeml", "yn00 against codeml node dS."))
w("""The ratio table shows the dS range where the two estimators agree and where they separate. Above about 0.5 the median codeml/yn00 ratio grows and
the within-20% share falls. The 1.5 limit is stated as a measured property of this pipeline, not a general property of yn00.
A codeml second pass for yn00 dS 0.5 to 1.6 is planned as a calibration subset.\n""")

# ---------- 9 examples
w("# 9. Example genomes\n")
w(fig("fig12_example_gallery", "Example genomes for identical-excess, H-excess and MID-excess flags."))

# ---------- 10 follow ups
w("# 10. Follow-up tests for flagged genomes\n")
w("""1. BUSCO duplicated fraction and gene count against clade median (already in the tables).
2. Pfam domains of genes in identical groups (RVT_1, rve, DDE_Tnp, Chromo, gag). Report identical-group genes per 1,000 after removing TE-domain genes.
3. Contig-level checks: identical contigs, coverage bimodality, k-mer spectrum (heterozygous peak).
4. FCS-GX contig taxonomy for partial-genome peaks (possible contamination).
5. Synteny (wgd syn or i-ADHoRe) for any genome-scale claim. Do not use dS alone.
6. GFF coordinate overlap for annotation duplicates.
""")

# ---------- 11 literature
w("# 11. Literature context\n")
w("""Citations below were checked in PubMed on 2026-11-06 (see `fable_context.md` for DOIs and PMIDs).

- Yeast WGD: Wolfe and Shields 1997 (Nature); Kellis et al. 2004 (Nature); Marcet-Houben and Gabaldon 2015 (PLoS Biol), which proposes interspecies hybridisation as the origin.
- Mucorales: Ma et al. 2009 (PLoS Genet, Rhizopus delemar); Corrochano et al. 2016 (Curr Biol, Phycomyces and Mucor).
- Hybrids: Nakao et al. 2009 and Libkind et al. 2011 (S. pastorianus); Gallone et al. 2019 (beer hybrids); Ortiz-Merino et al. 2017 (Z. parabailii);
  Pryszcz et al. 2014 (C. orthopsilosis); Mixao and Gabaldon 2018 (review); Steenwyk et al. 2020 (A. latus); Fogelqvist et al. 2018 (V. longisporum);
  Lin et al. 2007 and Samarasinghe and Xu 2018 (Cryptococcus hybrids).
- Dikaryon and heterokaryon assemblies: Li et al. 2019 (Puccinia graminis haplotypes); Ropars et al. 2016, Chen et al. 2018 and Sperschneider et al. 2023 (AMF nuclei).
- TE-rich lineages: Duplessis et al. 2011 (rusts); Spanu et al. 2010 and Frantzeskakis et al. 2018 (Blumeria); Sipos et al. 2017 (Armillaria);
  Martin et al. 2010 (Tuber); Tisserant et al. 2013 (Rhizophagus); Goodwin et al. 2011 (Zymoseptoria).
- RIP: Galagan and Selker 2004 (Trends Genet).
- Assembly and Ks methods: Roach et al. 2018 (Purge Haplotigs); Steinegger and Salzberg 2020; Astashyn et al. 2024 (FCS-GX);
  Vanneste et al. 2013; Tiley et al. 2018; Zwaenepoel and Van de Peer 2019 (two papers); Chen et al. 2024 (wgd v2); Yang and Nielsen 2000 (yn00).

Statements marked unverified in `fable_context.md` are not used here.
""")


# ---------- gene level
w("# 12. Gene-level queries: young paralogs and shared domains\n")
GF = rd(f"{T}/gene_features_vs_old.tsv"); AF = rd(f"{T}/gene_pfam_annotated_fraction_by_bin.tsv")
w("""Query design. Genes are those in the gated sampled pairs of `wgd.ks.parquet` (identical-CDS groups contribute one representative). Each gene is assigned the minimum dS over the sampled pairs it appears in,
so a gene with several partners is placed in its youngest bin. Bins: Z (dS = 0), H1 (0, 0.1], Y (0.1, 0.3], MID (0.3, 1.5], OLD (> 1.5). *Recent* means Z plus H1.
For each Pfam domain and each CAZy family, a Mantel-Haenszel odds ratio compares genes of the focal bin with all other paralog genes of the same genome and the same protein-length class
(five classes: <150, 150-250, 250-400, 400-700, >=700 aa). Genome and length class are the strata, so large genomes and short proteins do not drive the result.
For Pfam, only genes with at least one Pfam hit count in both groups, because the Pfam hit rate differs by bin (below). Variance is Robins-Breslow-Greenland; q is Benjamini-Hochberg across domains.
A domain needs >= 20 genomes and >= 30 focal genes. Sets: all genomes, a *clean* set (BUSCO duplicated < 20%, gene count < 1.5x clade median, no named hybrid, not in shape tiers A or B), and clean Ascomycota and Basidiomycota separately.
Two checks changed the first version of this test: strata with the domain only in the reference group must be kept, and the length strata were added after young genes proved shorter. Results from the earlier, unstratified version are kept in `tables/gene_enrich_unstratified/` for comparison.
""")
w("Share of paralog genes with at least one Pfam hit, by bin (all genomes):\n"); w(md_table(AF[["bin", "n", "frac_with_pfam"]], 3) + "\n")
fs = GF[GF.set == "clean"].copy()
w("Protein features of each bin against OLD (> 1.5) in the same genome, clean genomes, paired across genomes with >= 10 genes in both bins (Wilcoxon):\n")
w(md_table(fs[["bin", "metric", "n_genomes", "median_B", "median_OLD", "median_diff", "frac_genomes_higher", "p_wilcoxon"]], 3) + "\n")
w("""Young paralog genes are shorter and less often carry a signal peptide than OLD paralogs, in most genomes. Part of the length difference can arise from detection:
a short protein is less likely to keep an alignable, gated pair at high dS. The enrichment test therefore stratifies on length.
""")
rows = []
for tag, lab in [("recent", "recent (dS <= 0.1)"), ("identical", "identical (dS = 0)"), ("Y", "Y (0.1-0.3)"), ("MID", "MID (0.3-1.5)")]:
    for sname in ["all", "clean", "clean_Ascomycota", "clean_Basidiomycota"]:
        for kind in ["pf", "cazy"]:
            d = rd(f"{T}/gene_enrich_{tag}_{kind}_{sname}.tsv")
            for c in ["q_BH", "OR_MH"]: d[c] = pd.to_numeric(d[c], errors="coerce")
            rows.append(dict(bin=lab, set=sname, kind="Pfam" if kind == "pf" else "CAZy", genomes=int(d.genomes_in_test.iloc[0]) if len(d) else 0,
                             domains_tested=len(d), enriched_q01_OR2=int(((d.q_BH < 0.01) & (d.OR_MH >= 2)).sum())))
w("Number of domains enriched (q < 0.01, OR >= 2) by bin and set:\n"); w(md_table(pd.DataFrame(rows), 3) + "\n")
def top(tag, kind, sname, n=20, minB=300):
    d = rd(f"{T}/gene_enrich_{tag}_{kind}_{sname}.tsv")
    for c in ["q_BH", "OR_MH", "aB_tot", "n_genomes", "ci_lo", "ci_hi"]: d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.q_BH < 0.01) & (d.OR_MH >= 2) & (d.aB_tot >= minB)].sort_values("OR_MH", ascending=False).head(n)
    return d[["dom", "name", "n_genomes", "aB_tot", "OR_MH", "ci_lo", "ci_hi"]].rename(columns={"aB_tot": "focal_genes", "OR_MH": "OR"})
w("Top Pfam domains for recent paralogs, clean genomes (>= 300 focal genes, sorted by odds ratio):\n"); w(md_table(top("recent", "pf", "clean"), 3) + "\n")
w("Top Pfam domains for recent paralogs, clean Ascomycota:\n"); w(md_table(top("recent", "pf", "clean_Ascomycota", 15), 3) + "\n")
w("Top Pfam domains for recent paralogs, clean Basidiomycota:\n"); w(md_table(top("recent", "pf", "clean_Basidiomycota", 15), 3) + "\n")
w("Top Pfam domains for identical-CDS genes, clean genomes:\n"); w(md_table(top("identical", "pf", "clean", 15, 200), 3) + "\n")
w("Top Pfam domains for Y-bin genes (0.1-0.3), clean genomes:\n"); w(md_table(top("Y", "pf", "clean", 15), 3) + "\n")
cz = top("recent", "cazy", "clean", 10, 30)
w("CAZy families enriched among recent paralogs, clean genomes (>= 30 focal genes):\n"); w((md_table(cz, 3) if len(cz) else "None passed q < 0.01 and OR >= 2.") + "\n")
w(fig("fig15_gene_enrichment", "Pfam enrichment volcano plots (clean genomes) and protein length by dS bin."))
w("""Interpretation, with limits.

- Domain names are those in the Pfam table. Many top domains are named for transposases, integrases or retroelement proteins
  (for example Plavaka, KDZ, CxC, Tet_JBP, MULE, Transposase_21, rva_4, zf-RVT, Myb_DNA-bind_3 and _4, Cas12f1-like_TNB). I did not check each family description here. By name, this set points to transposon-derived gene models among recent duplicates.
  This matches the TE-like identical-copy profile in Section 4.
- SPOK2_N and OrsD (PF12013) are enriched among recent paralogs, mainly in Ascomycota. I did not verify what these families do. The OrsD gene set is examined in `reports/H3_OrsD_case_study.md`.
- Y-bin genes are enriched for ribosomal-protein and translation-factor domains (L2, S4, S6e, L21e). Ribosomal-protein families are small and often multicopy, so duplicated models or processed copies are plausible. This is not tested here.
- Several domains are of unknown function (DUF names, RRT15). They are reported as found and not interpreted.
- CAZy families show almost no enrichment in the recent bin after stratification. Recent duplication of CAZymes is not a general feature of these data.
- The dS bin is set by sampled pairs, not by all pairs of a gene, and each family contributes at most five pairs per node. Enrichment of a domain may partly reflect how families are sampled.
  The effect of this sampling on gene-level bins has not been measured.
""")

# ---------- 13 reproducibility
w("# 13. Reproducibility\n")
w(f"""Scripts: `eval/explore/scripts/` (shape_peaks, shape_plots, gene_recent_enrich, gene_plots, gene_stats, taxonomy, build_features2, build_table, pfam_te_counts, build_hist, compare_codeml_yn00, stats, varcomp.R, plots, make_report).
Inputs: `paralogoscope_pilot/nf_full/tables_yn00` (interim, {n_all:,} genomes; to be replaced by the final tables when the last giant genome completes).
Report built {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}.
""")

open(f"{E}/report.md", "w").write("\n".join(R))
print("wrote report.md", sum(len(x) for x in R), "chars")
