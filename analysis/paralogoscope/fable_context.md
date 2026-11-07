# Paralog dS survey of 7,384 fungal genomes: scientific framing review

Author: peer review by Claude (Fable 5.1), 2026-11-06.
Scope: framing and statistics only. No code. Citations marked "verified" were checked in PubMed on 2026-11-06 (author, year, journal, DOI/PMID). Everything marked "unverified" is from memory and must be checked before use.

Data assumptions used here are exactly those in the brief: wgd dmd families (DIAMOND E<1e-5, MCL I=1.5), byte-identical CDS collapsed into groups, yn00 dS on the rest, pairs dropped at dS>2.5 or <30 codons, NJ-on-protein-identity node assignment, up to 5 sampled gated pairs per node, per-genome summaries, taxonomy to species, a " x " name flag. Five huge-family genomes and Pachysolen are absent.

---

## 0. Two framing corrections before anything else

1. **This data cannot detect most fungal WGDs.** The known fungal WGDs (Saccharomycetaceae, Mucorales/Rhizopus, Phycomyces) are old. Their retained paralog pairs sit at or past yn00 saturation (dS ≥ 2), so they fall into the "old/saturated" mass or are dropped by the dS>2.5 gate. The yeast WGD was recognised from synteny and ohnolog pairs, not from a Ks peak (Wolfe & Shields 1997; Kellis et al. 2004, verified). The usable dS window (0 to roughly 1.5) means question 2 is really: **"which genomes carry a recent genome-scale duplication, i.e. an allopolyploid/hybrid, a dikaryotic or heterokaryotic assembly, or a mixed/contaminated assembly."** Write the question that way. Rhizopus/Mucorales and Saccharomyces are *not* positive controls for a mid-dS peak; they are positive controls for the *old* mass and for ohnolog-like 2-gene family structure. Use the " x " hybrids and known allodiploids as the peak-positive controls.

2. **Aneuploidy is invisible.** Extra chromosome copies collapse into one contig in any standard assembly and produce no paralog pair. Drop "aneuploid" from the list of things this data can show. Segmental duplication is detectable but is small-mass and dS-mixed.

---

## A. Prioritised hypotheses with expected signature, confounders and discriminating tests

Signatures are stated for node-averaged yn00 dS in the 0-2 window, per genome.

### Question 1: excess of dS = 0 / near-zero paralogs

Separate three bins before interpretation. They have different causes:

- **Bin Z (byte-identical CDS groups).** Copies that have had no time to diverge, or copies that are the same DNA assembled twice.
- **Bin H (dS ≈ 0.002-0.05, nonzero).** The size of synonymous heterozygosity in a diploid/dikaryon. Allelic haplotypes land here, not in Z, unless the strain is homozygous.
- **Bin Y (dS 0.05-0.3).** Real young duplicates, and RIP-mutated duplicates in Pezizomycotina.

| Rank | Hypothesis | Expected signature in this data | Main confounders | Discriminating test / extra data |
|---|---|---|---|---|
| 1 | **Haplotype duplication in the assembly** (unpurged haplotigs; dikaryotic Basidiomycota such as rusts and fruiting-body Agaricomycetes assembled with both nuclei; AMF heterokaryons with two nuclear genomes) | Large mass in bin H (not Z); 2-gene families dominate; nodes per 1,000 genes close to 500 (i.e. nearly every gene has one partner); peak is narrow because it is one coalescent depth; genome gene count ~2x clade median. Rusts: Li et al. 2019 (verified) assembled the two haplotypes of Puccinia graminis separately. AMF: Sperschneider et al. 2023 (verified) assembled two co-existing nuclear genomes in Rhizophagus heterokaryons. | A true young allopolyploid gives a similar shape but at larger dS (parental divergence). | BUSCO duplicated fraction (>40-50% is diagnostic); gene count vs clade median; assembly contig count / N50 and whether a purge step (purge_haplotigs, Roach et al. 2018 verified; purge_dups) was run; read-depth bimodality (half-depth contigs); k-mer spectrum (GenomeScope-style heterozygous peak); contig-level self-synteny with near-identical collinear blocks. |
| 2 | **TE-derived gene models** (annotated transposase/RT/gag-pol ORFs) | Bin Z dominated by a few *large* identical groups (largest group size large); families of size >>2; such groups recur across genomes of the same clade. Expected in TE-rich lineages: Pucciniales (Duplessis et al. 2011 verified), Blumeria (Spanu et al. 2010; Frantzeskakis et al. 2018 verified: recent TE burst), Armillaria (Sipos et al. 2017 verified), Tuber (Martin et al. 2010 verified: genome ~58% TE, number from memory, unverified), Glomeromycota (Tisserant et al. 2013 verified). | Real multicopy gene arrays (histones, rDNA-linked, effector arrays) also give large identical groups. | Pfam/InterPro domains of the genes in identical groups (RVT_1, rve, DDE_Tnp, Chromo, gag); overlap with RepeatMasker/EDTA annotation; protein length distribution; fraction of the group lacking any non-TE domain. Report "identical-group genes per 1,000 genes after removing TE-domain genes" as the cleaned statistic. |
| 3 | **Duplicated contigs from merged/hybrid assemblies or re-included alt scaffolds** | Bin Z with many *small* (size 2) identical groups spread over the genome, no TE domains, not confined to one family; often paired with BUSCO duplication; sometimes also "N50 high but contig count anomalous". | Hard to separate from rank 1 if the strain is homozygous (haplotigs then are byte-identical). | Contig pair identity over whole contig length (identical contigs, not just identical CDS); NCBI assembly record (alt loci, "diploid assembly" flag); read-depth per contig. |
| 4 | **Recent tandem / segmental duplication (biological)** | Bin Z and bin Y, small mass, families of size 2-5, mixed dS inside a family; positions adjacent on a contig. | Assembly collapse hides some; RIP erases some. | Gene adjacency (same contig, distance), orientation; nodes with mixed dS inside family; GC3 and codon usage normal for the genome. |
| 5 | **Multicopy conserved genes** (histones, ubiquitin, tubulin isoforms, rRNA-adjacent, effector arrays in pathogens) | Few families, each with several identical or near-identical members; recur across all genomes of a clade (so the clade *mean* is high, not an outlier). | TE models. | Annotation/Pfam of the family; cross-genome consistency of the family's identical-copy count. |
| 6 | **Mixed culture / contamination by a close relative** | Mass at dS equal to host-contaminant divergence (could be bin H if same species, bin Y if sibling species, or mid peak if different species); mass is a *fraction* of the genome, uneven; contigs differ in coverage and often in GC. | Allopolyploid (full-genome, synteny preserved). | FCS-GX (Astashyn et al. 2024 verified) or contig-level taxonomic assignment; coverage bimodality; GC3 bimodality of the pair subset vs genome; BUSCO duplication *below* what a hybrid would give. |
| 7 | **RIP-active lineages show a deficit of bin Z/H and an excess in bin Y with a distinct mutation spectrum** (Pezizomycotina) | Few identical copies; excess pairs at dS 0.1-0.5 with dN/dS near 1 and extreme transition bias. RIP mutates duplicated sequence >~400 bp and >~80% identity (Galagan & Selker 2004 verified; thresholds from memory, unverified). | True young duplicates under relaxed selection. | Transition/transversion ratio and CpA->TpA bias (RIP indices) in the pair alignments; dN/dS near 1; such pairs often sit near TEs. |
| 8 | **Annotation transfer duplicates** (same gene predicted twice on overlapping models, or liftover from a reference with duplicate models) | Bin Z pairs that overlap on the same contig coordinates or are adjacent with identical CDS. | Tandem arrays. | Coordinate overlap check in the GFF; CDS coordinates identical or nested. |

Priority for a report: ranks 1-3 explain most outliers in public fungal assemblies; rank 7 explains most *clade-level* differences between Pezizomycotina and the rest; ranks 4-5 are the biology left after removing 1-3 and 7-8.

### Question 2: distinct mid-dS peaks (recent genome-scale duplication, hybrid, contamination)

| Rank | Hypothesis | Expected signature | Confounders | Discriminating test / extra data |
|---|---|---|---|---|
| 1 | **Allopolyploid / interspecific hybrid with both subgenomes assembled** | One narrow peak at the parental divergence (position depends on parents; width is set by sampling noise and rate variation, not by biology); 2-gene families dominate the peak; nodes per 1,000 genes high; gene count ~2x; BUSCO duplication high (>50%); long collinear blocks between the subgenomes. Known cases: Saccharomyces pastorianus (S. cerevisiae x S. eubayanus; Nakao et al. 2009, Libkind et al. 2011 verified); other beer hybrids (Gallone et al. 2019 verified); Zygosaccharomyces parabailii (Ortiz-Merino et al. 2017 verified); Candida orthopsilosis hybrids (Pryszcz et al. 2014 verified) and other Candida hybrids reviewed by Mixão & Gabaldón (2018, verified); Aspergillus latus allodiploid (Steenwyk et al. 2020 verified); Verticillium longisporum (Fogelqvist et al. 2018 verified); Cryptococcus AD hybrids (Lin et al. 2007 verified; Samarasinghe & Xu 2018 verified). | Dikaryon/haplotig duplication (peak at heterozygosity level, i.e. much lower dS); contamination. Many hybrid assemblies are *collapsed* short-read assemblies; then the homeologs are missing and the genome looks normal. Loss of heterozygosity regions reduce peak mass. | BUSCO duplication; gene count ratio; self-synteny block lengths; " x " names and literature list as positive controls; for yeasts, compare peak position to the dS between the two parent species' orthologs (should match). |
| 2 | **Contamination / mixed isolate from a related species** | Peak at the divergence, but mass is partial and uneven; contigs with the paralog partner show distinct coverage/GC; BUSCO duplication moderate. | Hybrid with heavy subgenome loss. | FCS-GX / contig taxonomy; coverage; GC3 of the paralog subset. Public fungal assemblies contain contamination at non-trivial rates (Steinegger & Salzberg 2020 verified for GenBank broadly; a fungal-specific survey number is unverified). |
| 3 | **Haplotype duplication (ranks 1 of Q1) misread as a peak** | Peak at dS ≈ heterozygosity (usually <0.05). | | Same tests as Q1 rank 1. Rule: a peak at dS<0.05 is "young burst", not "WGD". |
| 4 | **Lineage-specific burst of segmental/tandem duplication or TE-driven gene amplification** | Broad hump, mixed family sizes, not 2-gene dominated, BUSCO duplication low. | | Family-size structure; TE domain content; adjacency. |
| 5 | **Old WGD (Saccharomycetaceae, Mucorales, Phycomyces)** | No mid peak. Elevated old/saturated mass and elevated count of 2-gene families relative to sister clades. Rhizopus delemar WGD was inferred from paralogous syntenic blocks (Ma et al. 2009 verified); Phycomyces and Mucor show extensive duplication (Corrochano et al. 2016 verified). Yeast WGD: ~450-550 ohnolog pairs in S. cerevisiae (Kellis 2004 reported 457, from memory, unverified). | Saturation pile-up near dS 2 looks like a "peak" in mixture models. | Treat as old mass only. If a WGD claim is needed, use synteny (wgd syn / i-ADHoRe), not dS. |

**Signatures summary for a one-line rule set:**
- Position: <0.05 = heterozygosity/haplotigs; 0.05-1.5 = candidate recent event; >1.5 = old/saturated, no position claim.
- Mass and family structure: genome-scale events give 2-gene family dominance and nodes/1,000 genes approaching 500; local events do not.
- Fraction of genome: hybrid ~ whole; contamination ~ partial; TE burst ~ partial but family-concentrated.
- Composition: contamination gives GC3/codon-usage bimodality of the paralog subset; hybrids of close parents do not; hybrids of distant parents may.
- Independent data: BUSCO duplication, read depth, k-mer spectrum, FCS-GX, synteny. Without at least one of these, the dS peak alone supports "candidate", not a conclusion.

---

## B. Critique of the statistical plan

### B1. Units and pseudo-replication
- The node is not the unit of inference for genome-level claims; the genome is. Within a genome, nodes inside one family share history and a few big families (TEs) can contribute hundreds of nodes. Wilson intervals on node proportions assume independent nodes and will be too narrow. Use a family-level cluster bootstrap for per-genome proportions, or report Wilson intervals as optimistic and say so.
- Across genomes, species within a genus and genera within a family are not independent. With 1,860 species representatives plus 5,530 singletons the species level is handled, but genus/family/order structure is not. Kruskal-Wallis across clades treats genomes as exchangeable; a clade of 300 Aspergillus-like genomes will swamp a clade of 5. The p-values will be tiny and uninformative.
- "Share of gene pairs that are identical" is dominated by large families (pairs grow as n^2). Prefer genes-in-identical-groups per 1,000 genes and the largest identical group size.

### B2. Phylogenetic non-independence without a tree
Defensible options, in order of preference:
1. **Nested random-effects model on the taxonomy**: genome-level response (e.g. logit of young-burst mass, or log nodes per 1,000 genes) with random intercepts phylum/class/order/family/genus, plus fixed effects for assembly covariates (see B6). Report the variance partition (how much variation sits at order vs genus vs genome) and the clade BLUPs with CIs. This is the standard taxonomy-as-proxy approach when no tree exists (a formal reference for taxonomic random effects as a phylogeny proxy: unverified, but the practice is common in comparative ecology). Caveat: it treats the taxonomy as a star-within-star phylogeny; branch lengths are ignored.
2. **Within-order comparisons only**, with a minimum n per order (e.g. >=10 genomes, >=3 genera), and a sensitivity run with one random genome per genus repeated 100-1,000 times. Report the distribution of clade medians over resamples.
3. **Phylogenetic ANOVA** (Garland et al. 1993 style; citation unverified) needs a tree. A tree is obtainable from BUSCO single-copy orthologs on 7,384 genomes but that is a separate project. Do not promise it.
Do not use Kruskal-Wallis across phyla as the headline test. Use it only as a screening statistic.

### B3. The dS window and gates
- The <30 codon gate removes short genes and short alignments; the dS>2.5 gate removes old pairs. Both truncate the denominator differently per genome (TE-rich and fast-evolving genomes lose more pairs). Report the number and fraction of dropped pairs and nodes per genome as a column, and compute "young mass" with two denominators (all nodes; nodes with usable dS). If a result changes between them, it is a gate effect.
- 30 codons is too short for stable yn00 dS. yn00 variance is very high below ~100 codons and dS values of 0 or very large are common there. Run the main analysis at >=100 (or >=150) aligned codons as a sensitivity check and report both. If the young-burst mass differs materially, the short-alignment pairs are driving it.
- Near-zero dS is quantised: one synonymous difference over ~60 synonymous sites gives dS ≈ 0.017. Bins below ~0.02 are therefore coarse and length-dependent. Do not split bins finer than 0.02 below dS 0.1.
- Genetic code: confirm that CTG-clade yeasts (code 12) and any other non-standard-code genomes were run with the right code. The wrong code changes dN far more than dS, so dS results survive, but dN/dS-based RIP tests (A, rank 7) will not.
- The 1.5 cutoff for "unreliable" is reasonable given the measured yn00 behaviour (reads low above 0.5, saturates near 2). State it as a measured property of this pipeline, not as a universal yn00 fact.

### B4. Mixture models and peak detection
- Gaussian mixtures on log dS: the identical groups (dS = 0) cannot be logged. Model them as a separate point mass and fit the mixture to dS>0 only. Add a floor (e.g. 0.005) for nonzero but tiny values and state it.
- BIC on log-dS GMMs over-splits skewed, truncated distributions and routinely places components at both gate boundaries (near the floor and near 2). Vanneste et al. 2013 (verified) and Tiley et al. 2018 (verified) document that Ks mixture components are unstable and that the lowest-Ks and saturation regions produce spurious components. Zwaenepoel & Van de Peer 2019 (verified, MBE) showed that a background of small-scale duplication under a duplication-loss process produces a monotone decaying distribution that mixture models readily "split" into peaks. Any component with mean >1.5 or <0.05 is not a WGD/hybrid peak by definition; say so in the method.
- Require a peak to be (a) present in the KDE with a fixed bandwidth on log scale, (b) stable under node bootstrap, (c) stable under re-sampling of the 5 pairs per node (the per-node mean has 11-17% relative error at the 90th percentile; simulate it), and (d) stable under the >=100-codon gate. Report peak position with the bootstrap interval and peak mass (fraction of usable nodes within the component).
- Minimum nodes >=200 for mixture fitting is fine, but it biases the "peak catalogue" toward large, duplication-rich genomes. Report how many genomes per clade qualify.
- wgd v2 (Chen, Zwaenepoel & Van de Peer 2024, verified) implements node-weighted Ks distributions, mixture fitting and synteny; its documentation is the right place to align conventions.

### B5. Outlier flagging and shrinkage
- Empirical-Bayes shrinkage toward the clade mean assumes exchangeability. Hybrids, contaminated and haplotig-duplicated genomes are exactly the non-exchangeable cases; shrinkage will pull them toward the clade and hide them. Use shrunken estimates for clade summaries and *unshrunken* values for outlier flags. Say which is used in each table.
- Robust z (median/MAD on logit) within class: with 7,384 genomes some classes have <5 members and the MAD is undefined or zero. Set a minimum (e.g. n>=10) and fall back to phylum. Report the reference set used per genome.
- Multiplicity on 7,384 genomes x several statistics: BH-FDR is fine for the screen. Then treat flags as a screening list and confirm the top candidates with independent evidence (BUSCO duplication, FCS-GX, coverage). Do not report "N contaminated genomes" from dS alone.

### B6. Missing covariates
Assembly and annotation properties will explain a large part of the variation and are confounded with taxonomy (whole clades were sequenced by one centre with one technology). At minimum obtain: assembly N50, contig count, assembly size, BUSCO complete/duplicated/fragmented, sequencing technology/assembler if available, annotation source (funannotate wave vs external), gene count, and whether a haplotig purge was run. Add them as fixed effects in the model and show the marginal effect of taxonomy after adjustment. If these are unavailable, the report must say that clade differences are not separable from assembly-technology differences.

### B7. Effect sizes to report
- For clade contrasts: difference in medians and ratio of medians with cluster-bootstrap (by genus) CIs; Cliff's delta is acceptable as a scale-free complement. Report n genomes and n genera per clade next to every effect.
- For genome flags: the raw statistic, its within-clade robust z, and the independent-evidence columns. No p-values needed for individual genomes.
- For the mixed model: variance components and the share of total variance at each taxonomic level. This is the most honest single answer to "do clades differ".

### B8. Controls
- " x " name flag: good positive control, but incomplete. Many hybrids are not named with " x " (S. pastorianus, C. orthopsilosis hybrid strains, V. longisporum, A. latus, Z. parabailii, Cryptococcus AD hybrids, Pichia/Millerozyma sorbitophila). Build an explicit list from the literature (section C) and treat "named hybrid" and "literature hybrid" as two control sets. Also expect that some named hybrids were assembled collapsed and will show no peak; that is a result about assemblies, not a failure of the method.
- Saccharomyces/Lachancea/Kluyveromyces and Rhizopus: controls for old mass and 2-gene family structure, not for a mid peak (section 0).
- Negative controls for haplotig duplication: haploid/homokaryotic strains with long-read assemblies that were purged (e.g. well-curated model strains). Negative controls for contamination: reference genomes passed through FCS-GX.

---

## C. Published and evolutionary context

Only items with "verified" were checked in PubMed. Numbers in this section that are not in the verified paper titles/abstracts are from memory and marked unverified.

### C1. Whole-genome duplications
- **Saccharomycetaceae post-WGD clade.** Wolfe & Shields 1997, Nature, 10.1038/42711, PMID 9192896 (verified). Kellis, Birren & Lander 2004, Nature, 10.1038/nature02424, PMID 15004568 (verified; 457 ohnolog pairs, number unverified). Marcet-Houben & Gabaldón 2015, PLoS Biol, 10.1371/journal.pbio.1002220, PMID 26252497 (verified): the WGD is better explained as an ancient interspecies hybridisation (allopolyploidy) between two lineages related to Kluyveromyces/Lachancea/Eremothecium. Ohnolog dS is saturated; expect only old-mass signal. Post-WGD genera: Saccharomyces, Nakaseomyces (C. glabrata), Kazachstania, Naumovozyma, Tetrapisispora, Vanderwaltozyma. Pre-WGD: Kluyveromyces, Lachancea, Eremothecium, Zygosaccharomyces, Torulaspora. Broader context: Shen et al. 2018, Cell, 10.1016/j.cell.2018.10.023, PMID 30415838 (verified); Opulente et al. 2024, Science, 10.1126/science.adj4503, PMID 38662846 (verified; 1,154 Saccharomycotina genomes, number from memory).
- **Rhizopus delemar/oryzae WGD.** Ma et al. 2009, PLoS Genet, 10.1371/journal.pgen.1000549, PMID 19578406 (verified): WGD inferred from paralogous syntenic blocks and gene-family expansion. Ks of the WGD pairs was reported but I do not recall the value; unverified. Expect old mass, not a peak.
- **Mucorales/Phycomyces.** Corrochano et al. 2016, Curr Biol, 10.1016/j.cub.2016.04.038, PMID 27238284 (verified): extensive genome duplication in Phycomyces blakesleeanus and Mucor circinelloides, with expansion of signalling genes. Same caveat as Rhizopus.
- **Zygosaccharomyces parabailii.** Ortiz-Merino et al. 2017, PLoS Biol, 10.1371/journal.pbio.2002128, PMID 28510588 (verified): interspecies hybrid whose fertility was restored by WGD after a failed mating-type switch. Two subgenomes at moderate divergence; this is a mid-peak positive control if both subgenomes are in the assembly.
- Microsporidia, Chytridiomycota, Blastocladiomycota: I am not aware of a verified WGD report usable here. Microsporidia have compact genomes and fast evolution; dS saturates quickly and paralog counts are low. Treat them as "few nodes, mostly old" and do not over-interpret.

### C2. Hybridisation / allopolyploidy (recent; mid-peak positive controls if assembled uncollapsed)
- **Saccharomyces pastorianus** (S. cerevisiae x S. eubayanus). Nakao et al. 2009, DNA Res, 10.1093/dnares/dsp003, PMID 19261625 (verified). Libkind et al. 2011, PNAS, 10.1073/pnas.1105430108, PMID 21873232 (verified; S. eubayanus is the non-cerevisiae parent). Other beer hybrids (S. cerevisiae x S. kudriavzevii, x S. uvarum): Gallone et al. 2019, Nat Ecol Evol, 10.1038/s41559-019-0997-9, PMID 31636425 (verified). Expected homeolog dS: equal to the dS between the parent species; measure it directly from parent-genome orthologs in the dataset rather than quoting a number.
- **Candida orthopsilosis** hybrid clinical strains: Pryszcz et al. 2014, GBE, 10.1093/gbe/evu082, PMID 24747362 (verified). **C. metapsilosis**, **C. inconspicua**, **C. albicans** parental divergence hypothesis and others: reviewed in Mixão & Gabaldón 2018, Yeast, 10.1002/yea.3242, PMID 28681409 (verified). Note: C. albicans is not post-WGD; it is a CTG-clade (code 12) yeast. Nakaseomyces (C. glabrata) is post-WGD.
- **Aspergillus latus** allodiploid (A. spinulosporus x a relative of A. quadrilineatus): Steenwyk et al. 2020, Curr Biol, 10.1016/j.cub.2020.04.071, PMID 32502407 (verified). A Pezizomycotina positive control.
- **Verticillium longisporum** (allodiploid hybrid, V. dahliae-like parents): Fogelqvist et al. 2018, BMC Genomics, 10.1186/s12864-017-4407-x, PMID 29298673 (verified).
- **Cryptococcus AD hybrids** (C. neoformans x C. deneoformans) and gattii-complex hybrids: Lin et al. 2007, PLoS Genet, 10.1371/journal.pgen.0030186, PMID 17953489 (verified); Samarasinghe & Xu 2018, Infect Genet Evol, 10.1016/j.meegid.2018.10.011, PMID 30342094 (verified); Samarasinghe et al. 2020, Genes, 10.3390/genes11010101, PMID 31963231 (verified). Parent divergence is large enough that homeolog dS may be near or beyond 1; check.
- **Puccinia graminis Ug99** arose by somatic hybridisation (nuclear exchange), with haplotypes assembled separately: Li et al. 2019, Nat Commun, 10.1038/s41467-019-12927-7, PMID 31699975 (verified). Rust assemblies with two haplotypes are the canonical "whole genome at heterozygosity dS" case.
- **Arbuscular mycorrhizal heterokaryons** carry two nuclear genomes: Ropars et al. 2016, Nat Microbiol, 10.1038/nmicrobiol.2016.33, PMID 27572831 (verified); Chen et al. 2018, eLife, 10.7554/eLife.39813, PMID 30516133 (verified); Sperschneider et al. 2023, Nat Microbiol, 10.1038/s41564-023-01495-8, PMID 37884816 (verified; both genomes assembled with HiFi + Hi-C). Glomeromycota genomes in this dataset may therefore show a dS-low whole-genome doubling if a heterokaryon was assembled unphased; the five excluded huge-family genomes are also from this group.
- Fusarium: I know of no verified allopolyploid Fusarium genome. Unverified; do not list Fusarium as a hybrid control unless a source is found.
- Pichia/Millerozyma sorbitophila hybrid (Leh-Louis et al. 2012, unverified) and Zygosaccharomyces rouxii hybrid strains (unverified).

### C3. Excess identical/near-identical paralogs expected from TE content or multicopy arrays
- Pucciniales: Duplessis et al. 2011, PNAS, 10.1073/pnas.1019315108, PMID 21536894 (verified; large TE-rich genomes, Melampsora and Puccinia).
- Blumeria graminis: Spanu et al. 2010, Science, 10.1126/science.1194573, PMID 21148392 (verified); Frantzeskakis et al. 2018, BMC Genomics, 10.1186/s12864-018-4750-6, PMID 29788921 (verified; recent TE burst, "one-speed" genome).
- Armillaria: Sipos et al. 2017, Nat Ecol Evol, 10.1038/s41559-017-0347-8, PMID 29085064 (verified; genome expansion).
- Tuber melanosporum: Martin et al. 2010, Nature, 10.1038/nature08867, PMID 20348908 (verified; TE-rich genome).
- Glomeromycota (Rhizophagus irregularis): Tisserant et al. 2013, PNAS, 10.1073/pnas.1313452110, PMID 24277808 (verified; TE-rich, expanded families).
- Zymoseptoria tritici accessory chromosomes: Goodwin et al. 2011, PLoS Genet, 10.1371/journal.pgen.1002070, PMID 21695235 (verified; dispensome). Accessory chromosomes carry repeats and duplicated/degenerate genes; a Zymoseptoria dS peak claim should be checked against chromosome of origin.
- Effector arrays and histone clusters: known from many pathogen genomes; no single citation verified here.

### C4. RIP (Pezizomycotina) and the dS=0 deficit
- Galagan & Selker 2004, Trends Genet, 10.1016/j.tig.2004.07.007, PMID 15313550 (verified): RIP removes or mutates duplicated sequence; it is a reason RIP-active Pezizomycotina have few recent identical paralogs and why their young duplicates show inflated dS with a transition bias. Expect Pezizomycotina to sit low on identical-copy statistics for a mechanistic reason, not because duplication is rarer.

### C5. Assembly and contamination pitfalls
- Haplotig duplication in diploid/dikaryotic long-read assemblies: Roach et al. 2018, BMC Bioinformatics, 10.1186/s12859-018-2485-7, PMID 30497373 (verified; Purge Haplotigs).
- Database contamination: Steinegger & Salzberg 2020, Genome Biol, 10.1186/s13059-020-02023-1, PMID 32398145 (verified; >2 million contaminated GenBank entries). Astashyn et al. 2024, Genome Biol, 10.1186/s13059-024-03198-7, PMID 38409096 (verified; FCS-GX). A fungal-specific contamination rate for RefSeq/GenBank assemblies: unverified.
- Ks-plot methodology: Vanneste et al. 2013, MBE, 10.1093/molbev/mss214, PMID 22936721 (verified); Tiley et al. 2018, GBE, 10.1093/gbe/evy200, PMID 30239709 (verified); Zwaenepoel & Van de Peer 2019, Bioinformatics, 10.1093/bioinformatics/bty915, PMID 30398564 (verified; wgd); Zwaenepoel & Van de Peer 2019, MBE, 10.1093/molbev/msz088, PMID 31004147 (verified); Chen et al. 2024, Bioinformatics, 10.1093/bioinformatics/btae272, PMID 38632086 (verified; wgd v2). Yang & Nielsen 2000, MBE, 10.1093/oxfordjournals.molbev.a026236, PMID 10666704 (verified; yn00).

---

## D. Ten most informative figures/tables

1. **Node-dS distribution by class (ridgeline, each genome weighted equally, dS=0 shown as a separate bar).** Reader sees the clade-level shape (young/mid/old) and that most clades are monotone-decaying with no mid peak.
2. **Scatter: identical-group genes per 1,000 genes vs largest identical group size, coloured by phylum, TE-domain genes removed in a second panel.** Separates TE amplification (few large groups) from assembly duplication (many small groups) and shows how much of the signal disappears after TE filtering.
3. **Heatmap: order (n>=10 genomes, n>=3 genera) x dS-bin mass (Z, H, Y, mid 0.3-1.5, old >1.5, dropped).** The "composition and duplication potential" summary; include n genomes and n genera as marginal columns.
4. **Variance-components plot from the nested random-effects model** for young mass, mid mass, nodes/1,000 genes, with assembly covariates in and out. Answers "do clades differ, or do technologies differ".
5. **Positive-control panel: dS histograms for literature hybrids/allodiploids and " x "-named genomes next to their parent-genus genomes**, annotated with BUSCO duplication and gene-count ratio. Shows what a real homeolog peak looks like in this pipeline and which hybrids were assembled collapsed.
6. **Peak catalogue table**: every genome with a supported peak in 0.05-1.5 (bootstrap-stable, >=100-codon-stable): peak position with CI, mass, 2-gene family fraction, nodes/1,000 genes, BUSCO duplicated %, gene count / clade median, N50, FCS-GX result if available, literature status. This is the deliverable for question 2.
7. **Young-burst mass vs BUSCO duplicated % (and vs contig N50)**, coloured by sequencing technology if known. Shows how much of the "young burst" is haplotype duplication.
8. **Nodes per 1,000 genes vs gene count, by phylum, with Saccharomycetaceae post-WGD and Mucorales highlighted.** Duplication potential vs genome size; old-WGD clades should show elevated 2-gene family counts without a mid peak.
9. **RIP panel for Pezizomycotina vs other phyla: Ts/Tv (or RIP index) of paralog pairs by dS bin.** Explains the Pezizomycotina deficit of identical copies mechanistically.
10. **Sensitivity figure: peak detection and bin masses under codon gate 30 vs 100 vs 150, and under re-sampling of the 5 pairs per node.** Shows the results that are robust to the pipeline choices.

---

## E. What not to claim from this data

- Do not claim WGD from a dS peak. Claim "recent genome-scale duplication candidate"; WGD vs hybrid vs contamination vs haplotype duplication needs synteny, BUSCO, depth or contig taxonomy.
- Do not claim aneuploidy or ploidy level. Assemblies collapse copies.
- Do not report absolute ages or rates. yn00 dS in this window is not calibrated, rates vary by lineage, and the per-node mean carries 11-17% error.
- Do not interpret a peak at dS <0.05 as a duplication event. It is at the heterozygosity scale and is assembly-dependent.
- Do not interpret byte-identical copies as biology until TE-domain genes, overlapping models and duplicated contigs are removed.
- Do not say "Pezizomycotina duplicate less" from a deficit of dS=0 copies; RIP and deletion remove the evidence.
- Do not compare raw clade proportions across phyla without assembly covariates; whole clades differ in technology and curation.
- Do not report a count of contaminated genomes from dS alone; report a screened candidate list with the independent evidence columns.
- Do not generalise from orders with one or two genomes; report them in the tables but exclude them from clade contrasts.
- Do not treat Saccharomyces or Rhizopus as evidence that the method detects WGD; they cannot show a mid peak in this window.
- Do not cite numbers from this document that are marked unverified (ohnolog pair counts, TE percentages, RIP thresholds, contamination rates) without checking the source.
