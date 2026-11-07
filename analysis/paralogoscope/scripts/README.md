# Paralog dS exploration: scripts, statistics and how to rerun

Folder: `paralogoscope_pilot/eval/explore/`. Scripts are in `scripts/`. Outputs go to `data/` (large intermediates), `tables/` (small tables) and `figures/`.
The report is `report.md` and `report.pdf`. All numbers in the report text are computed by `make_report.py` from `tables/` and `data/`.

Run everything with SLURM dependencies: `bash scripts/run_all.sh`. Each `job_*.sh` runs one script. Use `/usr/bin/python3.12` (the default `python3` on the login node is 3.9).
Inputs are the paralogoscope tables (`nf_full/tables_yn00`, later `nf_full/tables_yn00_final`), `tables/busco_genome.parquet`, `tables/pfam.parquet`, `tables/cazy.overview.parquet`,
`tables/gene_proteins.parquet`, `tables/signalp.signal_peptide.parquet`, `tables/asm_stats.parquet` and the run sample sheet.
When the final tables are published, change `T=` in `job_features2.sh` (and the node-table paths in `job_hist.sh` and `job_compare.sh`), then rerun `run_all.sh`.

## Pipeline order

| Step | Script (job) | Input | Output | What it computes |
|---|---|---|---|---|
| 1 | `gene_stats.py` (`job_gene_stats.sh`) | CDS FASTA files | `data/gene_stats.tsv` | Genes, mean CDS length, GC, GC3 per genome. |
| 2 | `build_features2.py` (`job_features2.sh`) | `wgd.ks.nodes/pairs/summary.parquet` | `data/features2.tsv.gz`, `gmm_components2.tsv.gz`, `peaks2.tsv.gz` | Bin shares, cluster-bootstrap CIs, Gaussian mixture on log dS, peak stability flags. |
| 3 | `build_hist.py` (`job_hist.sh`) | node table | `data/ds_hist.tsv.gz` | 70 log-spaced bins, 0.005 to 2.5, tree nodes with dS > 0, plus counts of dS = 0 nodes, identical groups and nodes without dS. |
| 4 | `pfam_te_counts.py` (`job_pfam_te.sh`) | `tables/pfam.parquet` (DuckDB) | `data/pfam_te_counts.tsv` | Proteins with any Pfam hit and with one of 15 TE Pfam domains per genome. |
| 5 | `compare_codeml_yn00.py` (`job_compare.sh`) | wave-1 node tables from both estimators | `tables/yn00_vs_codeml_*.tsv`, `data/yn00_vs_codeml_pairs.tsv.gz` | Agreement of yn00 and codeml node dS by dS bin. |
| 6 | `taxonomy.py`, `build_table.py` | steps 1, 2, 4, sample sheet, asm_stats, BUSCO | `data/taxonomy.tsv`, `data/genome_table.tsv.gz` | One row per genome: features, taxonomy, assembly covariates. Derived columns: `nodes_per_1000_genes`, `share_2gene_fam`, `genes_ratio_clade`, `te_per_1000_genes`. |
| 7 | `stats.py` (`job_stats.sh`) | steps 2, 3, 6 | `tables/clade_summary_*`, `class_effects_vs_phylum`, `class_medians_one_per_genus`, `outlier_flags`, `controls`, `peak_catalogue`, `sensitivity_codon100`; `data/genome_table_final.tsv.gz` | Clade tables, effect sizes, outlier flags, controls, stable-peak catalogue, gate sensitivity. |
| 8 | `varcomp.R` (`job_varcomp.sh`) | `data/genome_table_final.tsv.gz` | `tables/varcomp.tsv`, `varcomp_fixed.tsv` | Nested random-effects model on taxonomy (lme4). |
| 9 | `shape_peaks.py`, `shape_plots.py` (`job_shape.sh`) | `data/ds_hist.tsv.gz`, genome table (BUSCO joined by assembly ID) | `tables/shape_peaks.tsv`, `shape_peaks_all.tsv`; fig13, fig14 | Interior peak detection (centre <= 1.0) from histogram shape, tiers by BUSCO duplication. |
| 10 | `gene_recent_enrich.py` (`job_gene_enrich.sh`) | wgd.ks pairs, Pfam, CAZy, proteins, SignalP, shape tiers | `tables/gene_enrich_<bin>_<pf\|cazy>_<set>.tsv`, `gene_features_vs_old.tsv`, `gene_pfam_annotated_fraction_by_bin.tsv`, `data/gene_features_by_genome_bin.tsv.gz` | Gene-level bins, Pfam/CAZy enrichment, protein features. |
| 11 | `gene_plots.py` (`job_gene_plots.sh`), `plots.py` (`job_plots.sh`) | tables and data | `figures/fig01 ... fig15` (PNG and PDF) | Figures. |
| 12 | `make_report.py`, `job_report.sh` | tables, figures | `report.md`, `report.pdf` | Report. PDF: system texlive xelatex with DejaVu Sans. |

`legacy/gene_recent_enrich_unstratified.py` is the first version of step 10 without protein-length strata. Its tables are in `tables/gene_enrich_unstratified/`. It is kept for comparison only.

## Statistical methods (what each number means)

**Bins of node dS.** Z: identical CDS (dS = 0). H: (0, 0.05]. Y: (0.05, 0.3]. MID: (0.3, 1.5]. OLD: (1.5, 2.5]. BEYOND: no usable dS (pair outside the identity window or failed the gates).
Shares use two denominators: all nodes (`shareAll_*`) and nodes with a usable dS. The gates are dS <= 2.5 and >= 30 codons. A 100-codon sensitivity run is in `sensitivity_codon100.tsv`.

**Per-genome intervals.** Family-level cluster bootstrap (200 draws, `build_features2.py`). Nodes within a family are not independent, so node-level Wilson intervals would be too narrow.

**Clade summaries (`stats.py`).** Medians over genomes. Intervals from a genus-cluster bootstrap (500 draws). A clade needs >= 10 genomes and >= 3 genera to be eligible.
Class effects against the rest of the phylum use the ratio of medians with a genus-cluster bootstrap CI and Cliff's delta. One random genome per genus, 100 draws, tests whether dense sampling of one genus drives a class median.

**Outlier flags (`stats.py`).** Robust z of the raw (unshrunken) value within the most specific clade with >= 10 genomes (order, then class, then phylum). The scale is MAD with a floor.
Flags are screening hits. They carry the evidence columns (BUSCO duplicated, gene-count ratio, 2-gene family share, contig count, N50, TE per 1,000 genes). No per-genome p-value is reported.

**Peak catalogue (`build_features2.py`, `stats.py`).** Gaussian mixture on log dS for genomes with >= 200 positive-dS nodes. A peak is stable if the node bootstrap support is >= 0.8,
it is present at >= 100 codons, it is stable under pair resampling (>= 0.8) and its position is in 0.05 to 1.5. The catalogue file contains only components that pass.

**Variance partition (`varcomp.R`).** Linear mixed models with nested random intercepts (phylum, class, order, family, genus) and a genome residual. Fixed effects in the covariate model: contig count, size, GC, gene count, BUSCO duplicated, TE per 1,000 genes.
The `share` column is the fraction of total random variance at each level. The taxonomy stands in for a phylogeny and ignores branch lengths.

**Shape-based peaks (`shape_peaks.py`).**
- Smoothed density: 3-bin moving average of the node-count histogram divided by the total.
- Peak: local maximum with prominence >= 0.35 x its height, centre in [0.03, 1.0], and >= 8% of tree nodes inside the half-height width. Genomes need >= 300 tree nodes with dS > 0.
- Sharp: half-height width <= 1.2 ln-units, peak mass >= 0.4 and prominence >= 0.8.
- Tiers: A sharp and BUSCO duplicated >= 40%; B sharp and 10-40%; C sharp and < 10%; D broad.
- Peaks above 1.0 are not called. These thresholds were chosen by inspection of named hybrids and Zygosaccharomyces genomes. They have not been tuned or validated on held-out data.

**Gene-level enrichment (`gene_recent_enrich.py`).**
- Gene set: genes in gated sampled pairs of `wgd.ks.parquet`. Identical groups contribute one representative gene. Each gene takes its minimum dS over the sampled pairs it is in.
- Focal bin B (recent = Z + H1 with dS <= 0.1, or Z, Y, MID alone). Reference: all other paralog genes of the same genome.
- Strata: genome x protein-length class (< 150, 150-250, 250-400, 400-700, >= 700 aa). A stratum needs >= 5 genes on each side.
- Pfam: only genes with >= 1 Pfam hit count in both groups, because the hit rate differs by bin (about 0.46 to 0.53 for young bins and 0.83 for OLD). CAZy: all paralog genes are the denominator.
- Mantel-Haenszel odds ratio per domain: OR = sum(a*d/n) / sum(b*c/n) over strata, where a = focal genes with the domain, b = focal genes without, c = reference genes with, d = reference genes without.
  Variance of ln(OR): Robins-Breslow-Greenland. Two-sided z test, Benjamini-Hochberg q across the domains tested in each set.
- A domain is tested if it occurs in >= 20 genomes among focal genes and in >= 30 focal genes. The test table reports "enriched" as q < 0.01 and OR >= 2.
- Sets: `all`; `clean` (BUSCO duplicated < 20%, gene count < 1.5x clade median, no named hybrid, not in shape tiers A or B); `clean_Ascomycota`; `clean_Basidiomycota`.
- Protein features: per-genome median protein length, share of proteins < 150 aa and share with a signal peptide, bin versus OLD, paired Wilcoxon across genomes with >= 10 genes in both bins.

## Known limits and bugs fixed

- The gene-level bin comes from sampled pairs (at most 5 per node), not all pairs of a gene. The effect of this sampling on gene-level bins is not measured.
- A first version of the enrichment test kept only strata where the focal group had the domain. That inflated nearly every odds ratio (about 88% of domains passed). The fix keeps every stratum where the domain occurs in either group.
  Adding length strata then reduced the "recent" Pfam hits further (see `tables/gene_enrich_unstratified/` against `tables/gene_enrich_*`).
- Domain function is read from the Pfam names in the database. Functional statements in the report are by name and need checking against Pfam and the literature before publication.
- Shape peaks near dS 0.5 to 1.5 are positions as read by yn00. yn00 reads lower than codeml above about 0.5, so true divergence is likely higher.
- The current tables are interim (7,384 genomes). One giant-family genome (*Jimgerdemannia lactiflua*) and *Pachysolen* are not yet in them.

## Gene-set and hypothesis scripts (added 2026-11-07)

| Script (job) | Output | What it does |
|---|---|---|
| `domain_geneset.py PF12013 OrsD` (`job_orsd.sh`) | `tables/geneset_OrsD_*.tsv` | Describes the gene set carrying one Pfam domain: dS-bin shares against all paralog genes, per-genome and per-class counts, protein length and Pfam model coverage, co-occurring domains, contig distance of recent pairs. |
| `domain_te_context.py PF12013 OrsD` (`job_orsd_te.sh`) | `tables/geneset_OrsD_te_context*.tsv` | Tests whether recent carriers sit near genes with one of 15 TE Pfam domains. Controls: up to 5 non-carrier genes per focal gene from the same genome and protein-length class. Mantel-Haenszel odds ratio, strata genome x length class. |
| `domain_deepdive.py PF12013 OrsD` (`job_orsd_deep.sh`) | `tables/geneset_OrsD_genome_enrichment.tsv`, `geneset_OrsD_Sorochytrium_genes.tsv` | Per-genome Fisher test of carrier recent share against the genome's other paralogs; gene table of one genome. |
| `make_hypothesis_reports.py` (`job_hypothesis_reports.sh`) | `reports/H1_*.md/.pdf`, `H2_*`, `H3_*`, `reports/figures/fig16_OrsD.*` | Writes three hypothesis-testing reports. Every number is computed from the tables. Verdict rules are stated in each test. |

Limits that apply to the hypothesis reports: interim tables (7,384 genomes); detector thresholds were chosen by inspecting positive controls; the TE-name pattern in H2.1 was written after seeing the result list; the TE marker in H3.5 is 15 Pfam domains and has low power.
## Repository copy (interim)

This folder was committed to `Fungi_BFD_run2026` as `analysis/paralogoscope/` from the working copy in `paralogoscope_pilot/eval/explore/`. The numbers are from the interim tables (7,384 genomes).
In the repository, tables larger than 200 KB are stored gzipped (`*.tsv.gz`); pandas reads them directly. The scripts still read and write the uncompressed working copy under `/bigdata/.../paralogoscope_pilot/eval/explore/`, so rerun them there.
`data/` (large intermediates) and `tables/gene_enrich_unstratified/` (first, unstratified version of the enrichment) are not committed. `run_all.sh` regenerates `data/`.
