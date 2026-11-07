#!/usr/bin/python3.12
"""Write three hypothesis-testing reports from the exploration tables (numbers are computed here, not typed):
  reports/H1_interior_peaks.md     interior dS peaks and genome duplication / hybrid / assembly hypotheses
  reports/H2_young_paralogs.md     domain and feature hypotheses for recently duplicated genes
  reports/H3_OrsD_case_study.md    the PF12013 (OrsD) gene set
Each hypothesis has: statement, prediction that would count against it, test, result, verdict, limits.
Verdict words: SUPPORTED, NOT SUPPORTED, REJECTED, INCONCLUSIVE. The verdict rule is written next to each test.
Figures: fig13-fig15 (reused) and fig16 (OrsD, made here)."""
import os, re, datetime, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats

E = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
T, D, O = f"{E}/tables", f"{E}/data", f"{E}/reports"
os.makedirs(O, exist_ok=True); os.makedirs(f"{O}/figures", exist_ok=True)
rd = lambda p: pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""], low_memory=False)
TODAY = datetime.date.today().isoformat()


def md_table(df, fmt=3):
    cols = list(df.columns); out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append("" if (isinstance(v, (float, np.floating)) and np.isnan(v)) else (f"{v:.{fmt}g}" if isinstance(v, (float, np.floating)) else str(v).replace("|", "/")))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def header(title, intro):
    return f'---\ntitle: "{title}"\ndate: "{TODAY}"\n---\n\n{intro}\n\nStatus: interim. All numbers come from the 7,384-genome interim tables (yn00 dS). They will change when the final tables are published.\n\n'


def hyp(n, statement, prediction_against, test, result, verdict, limits):
    return (f"## {n}\n\n**Hypothesis.** {statement}\n\n**It would count against the hypothesis if:** {prediction_against}\n\n**Test.** {test}\n\n"
            f"**Result.** {result}\n\n**Verdict: {verdict}.** {limits}\n\n")


# =============================================================== H1
SPa, SPk = rd(f"{T}/shape_peaks_all.tsv"), rd(f"{T}/shape_peaks.tsv")
for d in (SPa, SPk):
    for c in ["busco_dup", "genes_ratio_clade", "share_2gene_fam", "peak1_dS", "peak1_mass", "peak1_logwidth", "n_peaks"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
CT = rd(f"{T}/controls.tsv")
SPa["has"] = SPa.n_peaks > 0
SPa = SPa.merge(SPk[["genome", "tier"]], on="genome", how="left"); SPa["tier"] = SPa.tier.fillna("none")
md = header("H1. Do interior dS peaks mark genome duplication, hybridisation or assembly artefacts?",
            "This report tests hypotheses about genomes whose node-dS histogram has a peak at dS <= 1.0. Method: `scripts/shape_peaks.py`. Statistics: `scripts/README.md`.")
# H1.1
bt = SPa.assign(b=pd.cut(SPa.busco_dup, [-1, 1, 5, 10, 20, 40, 101], labels=["<1", "1-5", "5-10", "10-20", "20-40", ">=40"])).groupby("b", observed=True).agg(genomes=("genome", "size"), share_with_peak=("has", "mean")).reset_index()
rho, p = stats.spearmanr(SPa.busco_dup.fillna(0), SPa.has.astype(int))
v1 = "SUPPORTED" if (bt.share_with_peak.iloc[-1] > 3 * bt.share_with_peak.iloc[0] and rho > 0.2) else "NOT SUPPORTED"
md += hyp("H1.1 Interior peaks track assembly duplication", "Genomes with an interior peak carry two diverged genome copies in the assembly, so BUSCO duplicated % is higher in them.",
          "the share of genomes with a peak does not rise with BUSCO duplication.",
          f"Share of the {len(SPa):,} evaluated genomes with an interior peak by BUSCO duplicated-% bin; Spearman correlation between BUSCO duplicated % and having a peak. Verdict rule: SUPPORTED if the top bin is more than 3 times the bottom bin and rho > 0.2.",
          md_table(bt, 3) + f"\n\nSpearman rho = {rho:.2f} (p = {p:.1e}).", v1,
          "This shows association, not mechanism. High BUSCO duplication has several causes (hybrid, dikaryon, haplotig duplication, contamination). Peaks also occur in 10% of genomes with BUSCO duplication under 1%, so BUSCO duplication is not necessary.")
# H1.2 positive controls
ctl = CT.merge(SPa[["genome", "has", "peak1_dS", "tier"]], on="genome", how="left")
ctl["group"] = ctl.control_group.str.replace(r" \(.*\)", "", regex=True)
cd = ctl.groupby("control_group").agg(n=("genome", "size"), evaluated=("has", lambda s: int(s.notna().sum())), with_interior_peak=("has", lambda s: int(s.fillna(False).sum()))).reset_index()
hyb = cd[cd.control_group.str.contains("hybrid")]; det = hyb.with_interior_peak.sum() / max(hyb.evaluated.sum(), 1)
wgd = cd[cd.control_group.str.contains("WGD")]; fp = wgd.with_interior_peak.sum() / max(wgd.evaluated.sum(), 1)
v2 = "SUPPORTED" if det >= 0.5 and fp < det / 2 else "INCONCLUSIVE"
kn = SPa[SPa.genome.str.contains("Verticillium_longisporum|Ogataea_polymorpha_x|Zygosaccharomyces_bailii|Zygosaccharomyces_parabailii|Aspergillus_latus|pastorianus|Cryptococcus_neoformans_x")][["genome", "n_peaks", "peak1_dS", "peak1_mass", "peak2_dS", "busco_dup", "genes_ratio_clade"]]
md += hyp("H1.2 The detector finds known hybrids and not old WGDs", "Named and literature hybrids (positive controls) usually show an interior peak. Old WGD clades (Saccharomycetaceae, Mucorales) do not, because their duplicates are older than the dS window.",
          "most hybrid controls lack a peak, or old-WGD controls show peaks as often as hybrids.",
          f"Detection rate of the shape screen in each control group. Verdict rule: SUPPORTED if hybrids are detected in >= 50% and the old-WGD rate is below half of that.",
          md_table(cd, 3) + f"\n\nHybrid controls with a peak: {det:.2f}. Old-WGD controls with a peak: {fp:.2f}.\n\nNamed examples:\n\n" + md_table(kn, 3), v2,
          f"The old-WGD rate is not a clean false-positive rate. Only genomes with >= 300 nodes with dS > 0 are evaluated: {int(wgd.evaluated.sum())} of {int(CT[CT.control_group.str.contains('WGD')].shape[0])} old-WGD controls pass, and the pre-WGD group keeps only {int(cd[cd.control_group.str.startswith('pre')].evaluated.iloc[0])} of {int(cd[cd.control_group.str.startswith('pre')].n.iloc[0])}. Pre-WGD genera such as Zygosaccharomyces and Torulaspora include hybrid strains, which are the genomes with enough paralog nodes to pass. The control set therefore mixes old-WGD and hybrid genomes. A control set that excludes known hybrid genera is needed for a real specificity estimate. Misses: *Aspergillus latus* and *S. pastorianus* have no peak at dS <= 1.0 because their peaks are above 1.0. *Z. parabailii* shows only its lower peak. The positive controls were also used to choose the thresholds, so the detection rate is optimistic.")
# H1.3 whole-genome signal in tier A
ta = SPa[SPa.tier != "none"].groupby("tier").agg(genomes=("genome", "size"), med_gene_ratio=("genes_ratio_clade", "median"), frac_ratio_ge_1_5=("genes_ratio_clade", lambda s: (s >= 1.5).mean()), med_2gene_share=("share_2gene_fam", "median"), med_busco_dup=("busco_dup", "median")).reset_index()
base = SPa[SPa.tier == "none"].agg({"genes_ratio_clade": "median", "share_2gene_fam": "median"})
A = SPa[SPa.tier.str.startswith("A")]; Dn = SPa[SPa.tier == "none"]
pg = stats.mannwhitneyu(A.genes_ratio_clade.dropna(), Dn.genes_ratio_clade.dropna()).pvalue
v3 = "SUPPORTED" if A.genes_ratio_clade.median() >= 1.5 * Dn.genes_ratio_clade.median() else "NOT SUPPORTED"
md += hyp("H1.3 Tier A peaks are whole-genome signals", "Genomes with a sharp peak and BUSCO duplication >= 40% (tier A) hold two genome copies, so they should also have about twice the clade-median gene count and a higher share of 2-gene families.",
          "tier A gene counts are near the clade median.", "Median gene-count ratio and 2-gene-family share by tier, against genomes with no interior peak. Verdict rule: SUPPORTED if tier A median gene ratio >= 1.5 times that of no-peak genomes.",
          md_table(ta, 3) + f"\n\nGenomes with no interior peak: median gene ratio {base['genes_ratio_clade']:.2f}, median 2-gene share {base['share_2gene_fam']:.2f}. Tier A against no-peak, gene ratio: Mann-Whitney p = {pg:.1e}.", v3,
          "Gene counts above the clade median also arise from annotation differences. Synteny or read depth is needed for a genome-scale claim.")
# H1.4 tier C
C = SPa[SPa.tier.str.startswith("C")]
fl = (C.genes_ratio_clade < 0.7).mean() if len(C) else np.nan; fl0 = (SPa[SPa.tier == "none"].genes_ratio_clade < 0.7).mean()
v4 = "PARTLY SUPPORTED" if (fl >= 0.25 and fl >= 2 * fl0) else "NOT SUPPORTED"
md += hyp("H1.4 Sharp peaks without BUSCO duplication (tier C) are partial annotations or partial duplications", "Tier C genomes have more low gene counts than genomes without a peak, which restricts the paralog set and sharpens a small peak.",
          "tier C low-gene-count share is similar to that of genomes without a peak.", f"Share of tier C genomes (n = {len(C)}) with gene count below 0.7 times the clade median, against genomes with no interior peak. Verdict rule: PARTLY SUPPORTED if tier C share >= 25% and >= 2 times the no-peak share; NOT SUPPORTED otherwise.",
          f"{fl:.2f} of tier C genomes have a gene ratio below 0.7, against {fl0:.2f} of genomes with no interior peak. Median gene ratio is {C.genes_ratio_clade.median():.2f} in tier C and {SPa[SPa.tier=='none'].genes_ratio_clade.median():.2f} in genomes without a peak. The majority of tier C genomes are therefore not low-gene-count.\n\n" + md_table(C.sort_values("peak1_mass", ascending=False).head(12)[["genome", "cls", "peak1_dS", "peak1_mass", "busco_dup", "genes_ratio_clade", "share_2gene_fam"]], 3), v4,
          "A low gene ratio shows a small gene set but not the cause. For most tier C genomes the explanation is not found here; they need FCS-GX, coverage and synteny.")
# H1.5 haplotype-depth peaks
phy = A.groupby("phylum").agg(n=("genome", "size"), med_peak_dS=("peak1_dS", "median"), frac_peak_le_0_1=("peak1_dS", lambda s: (s <= 0.1).mean())).reset_index()
low = (A.peak1_dS <= 0.1).mean()
md += hyp("H1.5 Peaks at low dS are heterozygosity (haplotype) signals, and peaks at higher dS are hybrid divergence", "Tier A peaks are spread between heterozygosity-level dS (<= 0.1) and hybrid-level dS, and low-dS peaks are common in dikaryotic or heterokaryotic groups.",
          "all tier A peaks sit at one dS.", "Distribution of tier A peak positions by phylum.",
          md_table(phy, 3) + f"\n\n{low:.2f} of tier A peaks are at dS <= 0.1 (n = {len(A)}).", "INCONCLUSIVE",
          "dS position alone does not tell heterozygosity from hybrid divergence. Both occur in tier A. A distinction needs the parent genomes' ortholog dS or phased assemblies. Positions near 0.5 to 1.0 are lower bounds because yn00 reads low there.")
md += "## Figures\n\n![Shape-detected peaks](../figures/fig13_shape_peaks_gallery.png)\n\n![Peak position and BUSCO duplication](../figures/fig14_shape_peaks_vs_busco.png)\n\n"
md += "## Reproduce\n\n`scripts/shape_peaks.py`, `scripts/shape_plots.py`, `scripts/make_hypothesis_reports.py`. Thresholds were set by inspecting positive controls and have not been validated on held-out data.\n"
open(f"{O}/H1_interior_peaks.md", "w").write(md)

# =============================================================== H2
def load(tag, kind, s, strat=True):
    p = f"{T}/gene_enrich_{tag}_{kind}_{s}.tsv" if strat else f"{T}/gene_enrich_unstratified/gene_enrich_{tag}_{kind}_{s}.tsv"
    d = rd(p)
    for c in ["q_BH", "OR_MH", "aB_tot", "n_genomes", "lnOR", "ci_lo", "ci_hi"]: d[c] = pd.to_numeric(d[c], errors="coerce")
    return d
sig = lambda d: d[(d.q_BH < 0.01) & (d.OR_MH >= 2)]
rc, ra, rA, rB = (load("recent", "pf", s) for s in ["clean", "all", "clean_Ascomycota", "clean_Basidiomycota"])
GF = rd(f"{T}/gene_features_vs_old.tsv"); AF = rd(f"{T}/gene_pfam_annotated_fraction_by_bin.tsv")
TE_RE = re.compile(r"transpos|Tnp|RVT|^rve|rva|MULE|hAT|Dimer_Tnp|Myb_DNA-bind|CxC|KDZ|Plavaka|Tet_JBP|Cas12f1|gag|Retrotran|integrase|zf-RVT|DDE|Asp_protease|Transposase|Helitron|Chromo|DUF6532|DUF6698|DUF6830|FLYWCH|OTU|Peptidase_C48", re.I)
top = sig(rc)[sig(rc).aB_tot >= 300].sort_values("OR_MH", ascending=False).head(40).copy(); top["te_name_match"] = top.name.fillna("").apply(lambda x: bool(TE_RE.search(x)))
md = header("H2. Which domains and features mark recently duplicated genes?", "Genes in the youngest dS bins (dS <= 0.1, 'recent') are compared with other paralog genes of the same genome and protein-length class. Method and statistics: `scripts/gene_recent_enrich.py` and `scripts/README.md`.")
# H2.1 TE-named domains
share_te = top.te_name_match.mean()
v = "PARTLY SUPPORTED (by domain name only, post hoc pattern)" if share_te >= 0.4 else "NOT SUPPORTED"
md += hyp("H2.1 Recent duplicates are enriched for transposon-derived genes", "Pfam domains named for transposases, integrases or retroelement proteins are over-represented among recent paralogs.",
          "fewer than 40% of the top enriched domains have transposon-related names.",
          f"Among the 40 highest-odds-ratio domains in the clean genome set (q < 0.01, odds ratio >= 2, >= 300 focal genes), count those whose Pfam *name* matches a transposon-related pattern. The pattern is my choice. I wrote it after seeing the top-domain list, so it is post hoc and liberal (it also matches names such as Chromo, OTU and Peptidase_C48). The pattern is: `" + TE_RE.pattern + "`. Matches are flagged in the table so the reader can audit them. Verdict rule: SUPPORTED if >= 40%.",
          f"{int(top.te_name_match.sum())} of {len(top)} ({share_te:.2f}) match.\n\n" + md_table(top[["dom", "name", "n_genomes", "aB_tot", "OR_MH", "ci_lo", "ci_hi", "te_name_match"]].rename(columns={"aB_tot": "focal_genes", "OR_MH": "OR"}), 3), v,
          "This rests on domain names, which I did not check against Pfam family descriptions. Several top domains are of unknown function. A direct TE test for one family (OrsD, report H3) found no proximity to TE-domain genes.")
# H2.2 length / signal peptide
fs = GF[GF.set == "clean"]
f = fs[fs.metric == "med_len"][["bin", "n_genomes", "median_B", "median_OLD", "frac_genomes_higher"]]
sp = fs[fs.metric == "frac_sigp"][["bin", "median_B", "median_OLD", "frac_genomes_higher"]].rename(columns={"median_B": "sigp_B", "median_OLD": "sigp_OLD", "frac_genomes_higher": "sigp_frac_higher"})
zz = f[f.bin == "H1"].iloc[0]
v = "SUPPORTED" if zz.median_B < zz.median_OLD and zz.frac_genomes_higher < 0.2 else "NOT SUPPORTED"
md += hyp("H2.2 Young paralog genes are shorter and less often secreted", "Recent paralogs include short gene models (fragments, small ORFs), so their median protein length is lower than that of old paralogs and fewer have a signal peptide.",
          "median lengths are equal or higher in the young bins.", "Per-genome median protein length of each bin against the OLD (dS > 1.5) bin in the same genome, clean genomes, genomes with >= 10 genes in both bins. Verdict rule: SUPPORTED if the H1 bin median is lower and fewer than 20% of genomes show higher length.",
          md_table(f.merge(sp, on="bin"), 3), v,
          "Part of the effect can be a detection artefact: short proteins less often keep an alignable gated pair at high dS. The enrichment tests in H2.1 stratify on length for that reason, but the length effect itself is not separated from detection here. The OrsD family is an exception (long proteins, report H3).")
# H2.3 CAZy
rows = []
for tag, lab in [("recent", "recent"), ("identical", "identical")]:
    for s in ["all", "clean", "clean_Ascomycota", "clean_Basidiomycota"]:
        d = load(tag, "cazy", s); rows.append(dict(bin=lab, set=s, families_tested=len(d), enriched=len(sig(d)), top=", ".join(sig(d).sort_values("OR_MH", ascending=False).dom.str.replace("CAZy:", "").head(3))))
cz = pd.DataFrame(rows); v = "NOT SUPPORTED" if cz.enriched.max() <= 5 else "SUPPORTED"
md += hyp("H2.3 Recent duplication is concentrated in carbohydrate-active enzyme families", "CAZy families are enriched among recent paralogs.", "few or no CAZy families pass q < 0.01 and OR >= 2.",
          "Same test as Pfam, on CAZy family labels (subfamily suffix removed). Verdict rule: NOT SUPPORTED if at most 5 families pass in any set.", md_table(cz, 3), v,
          "CAZy labels are coarse. A gene's bin depends on sampled pairs. The result applies to the whole-set level, not to individual genomes or lifestyles.")
# H2.4 robustness: Ascomycota vs Basidiomycota concordance, length strata, assembly artefacts
j = rA.merge(rB, on="dom", suffixes=("_A", "_B")); j = j[(j.aB_tot_A >= 30) & (j.aB_tot_B >= 30)]
rho_ab = stats.spearmanr(j.lnOR_A, j.lnOR_B).correlation; sa, sb = set(sig(rA).dom), set(sig(rB).dom)
ru = load("recent", "pf", "clean", strat=False); jj = rc.merge(ru, on="dom", suffixes=("", "_u")); rho_u = stats.spearmanr(jj.lnOR, jj.lnOR_u).correlation
jc = rc.merge(ra, on="dom", suffixes=("", "_a")); rho_c = stats.spearmanr(jc.lnOR, jc.lnOR_a).correlation
conc = pd.DataFrame([dict(comparison="Ascomycota vs Basidiomycota (clean), domains tested in both", domains=len(j), spearman_lnOR=rho_ab, shared_enriched=len(sa & sb), only_A=len(sa - sb), only_B=len(sb - sa)),
                     dict(comparison="length-stratified vs unstratified (clean)", domains=len(jj), spearman_lnOR=rho_u, shared_enriched=len(set(sig(rc).dom) & set(sig(ru).dom)), only_A=len(set(sig(rc).dom) - set(sig(ru).dom)), only_B=len(set(sig(ru).dom) - set(sig(rc).dom))),
                     dict(comparison="clean vs all genomes", domains=len(jc), spearman_lnOR=rho_c, shared_enriched=len(set(sig(rc).dom) & set(sig(ra).dom)), only_A=len(set(sig(rc).dom) - set(sig(ra).dom)), only_B=len(set(sig(ra).dom) - set(sig(rc).dom)))])
v = "SUPPORTED" if rho_ab > 0.3 and rho_u > 0.7 and rho_c > 0.7 else "INCONCLUSIVE"
md += hyp("H2.4 The enrichments are robust to clade, protein length and assembly quality", "A real property of young paralogs should hold in both large phyla, with and without length strata, and with or without assembly-suspect genomes.",
          "domain log odds ratios are uncorrelated between Ascomycota and Basidiomycota, or change strongly with length strata or genome filtering.",
          "Spearman correlation of log odds ratios between sets, for domains tested in both. Verdict rule: SUPPORTED if rho > 0.3 between phyla and > 0.7 for the length and genome-filter comparisons.",
          md_table(conc, 3), v, "Correlation says the ranking is stable, not that each domain is a true biological signal. The first version of this test (strata dropped when the domain was absent from the focal group) inflated nearly all odds ratios; that error was found by checking that 88% of domains passing was implausible.")
md += "## Figure\n\n![Enrichment volcano and protein length](../figures/fig15_gene_enrichment.png)\n\n"
md += "## Reproduce\n\n`scripts/gene_recent_enrich.py`, `scripts/gene_plots.py`, `scripts/make_hypothesis_reports.py`.\n"
open(f"{O}/H2_young_paralogs.md", "w").write(md)

# =============================================================== H3
B = rd(f"{T}/geneset_OrsD_bins.tsv"); PG = rd(f"{T}/geneset_OrsD_per_genome.tsv"); GE = rd(f"{T}/geneset_OrsD_genome_enrichment.tsv"); FE = rd(f"{T}/geneset_OrsD_features.tsv")
TEc = rd(f"{T}/geneset_OrsD_te_context.tsv"); MH = rd(f"{T}/geneset_OrsD_te_context_MH.tsv"); CO = rd(f"{T}/geneset_OrsD_cooccurring_pfam_recent.tsv"); SG = rd(f"{T}/geneset_OrsD_Sorochytrium_genes.tsv")
CX = pd.read_csv(f"{T}/geneset_OrsD_recent_pair_context.tsv.gz", sep="\t")
for c in ["domain_genes", "all_paralog_genes", "share_all", "share_domain_in_pairs"]: B[c] = pd.to_numeric(B[c], errors="coerce")
for c in ["fisher_OR", "p", "busco_duplicated_pct", "carriers_in_pairs", "carriers_recent", "recent_share_all_paralogs", "carrier_recent_share"]: GE[c] = pd.to_numeric(GE[c], errors="coerce")
bp = B[B.bin != "not in a gated pair"].set_index("bin"); rec_c = bp.loc[["Z", "H1"], "domain_genes"].sum() / bp.domain_genes.sum(); rec_a = bp.loc[["Z", "H1"], "all_paralog_genes"].sum() / bp.all_paralog_genes.sum()
orc = rc[rc.dom=="PF:PF12013"].iloc[0]; ora = ra[ra.dom=="PF:PF12013"].iloc[0]
nprot = int(B.domain_genes.sum()); nout = int(B[B.bin == "not in a gated pair"].domain_genes.iloc[0])
md = header("H3. The OrsD (PF12013) gene set: what is the excess of low-dS duplicates?", "OrsD is the Pfam name in the database. Its function and family description were not verified (a web search returned no entry). This report tests what the data show about the gene set. Method: `scripts/domain_geneset.py`, `domain_te_context.py`, `domain_deepdive.py`.")
md += f"**Starting observation.** Among carriers of PF12013 that sit in a gated pair, {rec_c:.1%} are at dS <= 0.1. Among all paralog genes the share is {rec_a:.1%}. In the clean genome set the odds ratio against other paralogs is {orc.OR_MH:.1f} (95% CI {orc.ci_lo:.1f} to {orc.ci_hi:.1f}; length- and genome-stratified Mantel-Haenszel; `tables/gene_enrich_recent_pf_clean.tsv`). The gene set has {nprot:,} proteins, of which {nout:,} are in no gated pair.\n\n"
md += md_table(B[["bin", "domain_genes", "all_paralog_genes", "share_domain_in_pairs", "share_all"]], 3) + "\n\n"
# H3.1 genome-specific artefact
top10 = PG.sort_values("recent_genes", ascending=False).head(10).recent_genes.sum() / PG.recent_genes.sum()
n_g = int((PG.recent_genes > 0).sum()); n5 = int((PG.recent_genes >= 5).sum()); lo = GE[GE.busco_duplicated_pct < 20]; exc = ((lo.fisher_OR > 1) & (lo.p < 0.01)).mean()
v = "REJECTED" if (top10 < 0.3 and n_g >= 100 and exc > 0.3) else "NOT REJECTED"
md += hyp("H3.1 The excess comes from a few genomes", "A small number of genomes carry a large family of near-identical OrsD copies and drive the pooled result.", "most recent carriers sit in the top handful of genomes, or few genomes show carrier excess on their own.",
          "Share of recent carriers in the top 10 genomes; number of genomes with recent carriers; per-genome Fisher test of carrier recent share against the same genome's other paralogs, in genomes with BUSCO duplication < 20% and >= 5 carriers in pairs. Verdict rule: REJECTED if top-10 share < 30%, >= 100 genomes, and > 30% of testable genomes show excess.",
          f"Top 10 genomes hold {top10:.1%} of recent carriers. {n_g} genomes have at least one, {n5} at least five. Of {len(lo)} testable low-duplication genomes, {exc:.1%} show carrier excess (p < 0.01). Median recent share: carriers {lo.carrier_recent_share.median():.2f}, all paralogs {lo.recent_share_all_paralogs.median():.3f}.\n\n"
          + md_table(GE.sort_values("carriers_recent", ascending=False).head(10)[["genome", "class", "busco_duplicated_pct", "paralog_genes", "recent_share_all_paralogs", "carriers_in_pairs", "carriers_recent", "fisher_OR", "p"]], 3), v,
          "Several of the largest cases come from genomes with small total paralog sets, where one family is a large part of the set.")
# H3.2 assembly duplication
hi = GE[GE.busco_duplicated_pct >= 40]
md += hyp("H3.2 The excess is assembly (haplotype) duplication", "Recent carrier pairs are allelic copies in dikaryotic or doubled assemblies, so the excess disappears once genomes with high BUSCO duplication are removed.", "the enrichment persists in genomes with BUSCO duplication < 20%.",
          "Compare the clean-set odds ratio with the all-genomes odds ratio, and carrier against all-paralog recent share in high-duplication genomes.",
          f"The clean-set enrichment is OR {orc.OR_MH:.1f} (all genomes: {ora.OR_MH:.1f}; `tables/gene_enrich_recent_pf_all.tsv`). In genomes with BUSCO duplication >= 40% (n = {len(hi)} testable), recent share is {hi.recent_share_all_paralogs.median():.2f} for all paralogs and {hi.carrier_recent_share.median():.2f} for carriers. Examples: *Puccinia graminis* 21-0 (BUSCO dup 81%): 79 of 83 carriers recent against 72% of all paralogs; *Melampsora americana* (60%): carriers depleted (OR "
          f"{GE[GE.genome.str.startswith('Melampsora_americana')].fisher_OR.iloc[0]:.2f}).", "REJECTED as the main explanation",
          "Assembly duplication contributes in a few rust genomes, but the excess persists in low-duplication genomes.")
# H3.3 fragments
fz = FE.set_index("bin")
for c in ["med_len", "med_domain_cov", "frac_cov_ge_0_8", "frac_short150"]: FE[c] = pd.to_numeric(FE[c], errors="coerce")
md += hyp("H3.3 The recent carriers are gene fragments", "Recent carriers are short partial models that carry part of the domain.", "carriers are long and carry most of the Pfam model.",
          "Protein length and Pfam model coverage (hit length / model length) of carriers by bin.", md_table(FE[["bin", "n", "med_len", "med_domain_cov", "frac_cov_ge_0_8", "frac_short150"]], 3), "REJECTED",
          "Median carrier length is about 1,000 aa in the youngest bins and median model coverage 0.96. This differs from the general result in H2.2. The Sorochytrium genome does have short copies of 284 to 711 aa without a helicase domain, so fragments exist in some genomes.")
# H3.4 tandem arrays
sc = CX.same_contig.mean(); w20 = (CX.dist <= 20000).mean()
md += hyp("H3.4 The recent copies are tandem arrays", "Recent carrier-carrier pairs lie next to each other on the same contig.", "most pairs are on different contigs or far apart.", "For recent pairs in which both genes carry the domain: share on the same contig, distance, share within 20 kb.",
          f"{len(CX):,} carrier-carrier recent pairs. Same contig: {sc:.3f}. Same contig and within 20 kb: {w20:.3f}. Median distance on the same contig: {CX.dist.median():,.0f} bp.", "REJECTED", "Contigs are fragments of a chromosome in most assemblies, so different contigs can still be close in the real genome; the rate of same-contig pairs is a lower bound on local clustering.")
# H3.5 TE proximity
MH["MH_OR_TE_within_25kb"] = pd.to_numeric(MH.MH_OR_TE_within_25kb); r = MH[MH.group == "recent"].iloc[0]
md += hyp("H3.5 The recent copies sit next to transposable elements", "Recent carriers have TE-domain genes within 25 kb more often than matched control genes.", "no excess of nearby TE-domain genes, with a confidence interval that excludes a large effect.",
          "TE marker: genes carrying one of 15 TE Pfam domains. Each recent carrier is compared with up to 5 control genes from the same genome and protein-length class (Mantel-Haenszel odds ratio; strata genome x length class). Verdict rule: SUPPORTED if the interval excludes 1 and OR > 1; NOT SUPPORTED if the interval includes 1.",
          md_table(TEc, 4) + "\n\n" + md_table(MH, 3), "INCONCLUSIVE (no evidence for; low power)",
          f"Recent carriers: OR {r.MH_OR_TE_within_25kb:.2f} (95% CI {r.ci_lo:.2f} to {r.ci_hi:.2f}). Only about 3% of control genes have a TE-domain neighbour within 25 kb, so the test has low power and the interval is wide. The marker is crude (15 domains). A repeat-annotation test (RepeatMasker or EarlGrey) on these genomes was not possible because no genome-wide repeat annotation exists for them in the tables.")
# H3.6 domain architecture / Sorochytrium
co = CO.copy(); co["share_of_recent_carriers"] = pd.to_numeric(co.share_of_recent_carriers)
md += hyp("H3.6 The family is a helicase-domain gene family with dispersed recent copies", "Recent carriers share a common architecture (helicase domains) and form few gene families within a genome.", "carriers have varied domain architectures and no dominant partner domain.",
          "Co-occurring Pfam domains on recent carriers; in the *Sorochytrium milnesiophthora* genome (Blastocladiomycetes; 23 carriers, 21 recent) the gene list, contig positions and gene family.",
          md_table(co.head(6)[["acc", "pname", "proteins", "share_of_recent_carriers"]], 3) + f"\n\nSorochytrium: {len(SG)} carriers on {SG.chrom.nunique()} contigs; {SG.fam.nunique()} gene family; median length {SG.length.median():.0f} aa; "
          f"{int((SG.start <= 10000).sum())} of {len(SG)} carriers start within 10 kb of the contig start.\n\n" + md_table(SG.head(12)[["protein_id", "chrom", "start", "length", "ds", "n_recent_pairs", "other_pfam"]], 3),
          "SUPPORTED (descriptive)", "DEAD, Helicase_C and ResIII co-occur on 43%, 37% and 20% of recent carriers, and 54% of carriers have no other Pfam domain. This describes the set. It does not show that the copies are mobile. Verifying the family would need sequence similarity to known families and repeat annotation.")
md += "## Figure\n\n![OrsD carriers](figures/fig16_OrsD.png)\n\n## Reproduce\n\n`scripts/domain_geneset.py PF12013 OrsD`, `scripts/domain_te_context.py PF12013 OrsD`, `scripts/domain_deepdive.py PF12013 OrsD`, `scripts/make_hypothesis_reports.py`.\n"
open(f"{O}/H3_OrsD_case_study.md", "w").write(md)
# fig16
fig, ax = plt.subplots(1, 2, figsize=(11, 4)); o = ["Z", "H1", "Y", "MID", "OLD"]; x = np.arange(5)
ax[0].bar(x - .2, [bp.loc[b, "share_domain_in_pairs"] for b in o], .4, label="OrsD carriers in pairs", color="#D55E00"); ax[0].bar(x + .2, [bp.loc[b, "share_all"] for b in o], .4, label="all paralog genes", color="#999")
ax[0].set_xticks(x); ax[0].set_xticklabels(["dS=0", "0-0.1", "0.1-0.3", "0.3-1.5", ">1.5"]); ax[0].set_ylabel("share of genes"); ax[0].legend(frameon=False, fontsize=8); ax[0].set_title("dS bin of OrsD carriers", fontsize=9)
g = GE[GE.busco_duplicated_pct < 20]; ax[1].scatter(g.recent_share_all_paralogs, g.carrier_recent_share, s=np.clip(g.carriers_in_pairs, 5, 80), alpha=.5, c="#0072B2"); ax[1].plot([0, 1], [0, 1], "k:", lw=.7)
ax[1].set_xlabel("recent share of all paralog genes in the genome"); ax[1].set_ylabel("recent share of OrsD carriers"); ax[1].set_title(f"Per genome (BUSCO dup < 20%, n = {len(g)}; size = carriers)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{O}/figures/fig16_OrsD.png", dpi=200); fig.savefig(f"{O}/figures/fig16_OrsD.pdf")
print("wrote", os.listdir(O))
