#!/usr/bin/python3.12
"""Join the per-genome paralog features with taxonomy, gene statistics and assembly covariates into data/genome_table.tsv.
Derived columns: shares of the dS bins with both denominators, nodes per 1,000 genes, paralog gene fraction, identical-group genes per 1,000 genes,
two-gene family share, gene-count ratio to the clade median (order, else class, else phylum), TE-domain protein counts where Pfam counts exist."""
import numpy as np, pandas as pd
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"; D = f"{R}/paralogoscope_pilot/eval/explore/data"
F = pd.read_csv(f"{D}/features2.tsv.gz", sep="\t")
T = pd.read_csv(f"{D}/taxonomy.tsv", sep="\t").drop_duplicates("genome")
G = pd.read_csv(f"{D}/gene_stats.tsv", sep="\t")
# covariates from the current tables/ (assembly ID -> locus tag through tables/species)
sp = pd.read_parquet(f"{R}/tables/species.parquet", columns=["ASMID", "LOCUSTAG"]).rename(columns={"LOCUSTAG": "species_prefix"})
asm = pd.read_parquet(f"{R}/tables/asm_stats.parquet", columns=["ASMID", "contig_count", "total_length_bp", "N50_bp", "L50", "gc_pct", "masked_pct", "t2t_scaffolds"])
bus = pd.read_parquet(f"{R}/tables/busco_genome.parquet", columns=["ASMID", "complete_pct", "duplicated_pct"]).rename(columns={"complete_pct": "busco_complete_pct", "duplicated_pct": "busco_duplicated_pct"})
gp = sp.merge(asm, on="ASMID", how="left").merge(bus, on="ASMID", how="left").drop(columns="ASMID").drop_duplicates("species_prefix")
X = F.merge(T, on="genome", how="left").merge(G, on="genome", how="left").merge(gp, on="species_prefix", how="left")
X = X.rename(columns={"PHYLUM": "phylum", "SUBPHYLUM": "subphylum", "CLASS": "class", "ORDER": "order", "FAMILY": "family_tax", "GENUS": "genus"})
# missing ranks: nested labels (parent|unassigned), a missing genus makes the genome its own group, so unrelated genomes are not pooled
X["phylum"] = X.phylum.fillna("unassigned")
X["class"] = X["class"].fillna(X.phylum + "|unassigned")
X["order"] = X["order"].fillna(X["class"] + "|unassigned")
X["family_tax"] = X.family_tax.fillna(X["order"] + "|unassigned")
X["genus"] = X.genus.fillna("g:" + X.genome)
# shares (usable denominator: nodes with a dS + identical groups; all-node denominator: every tree node + identical groups)
for nm in ["Z", "H", "Y", "MID", "OLD"]:
    X["share_" + nm] = X["n_" + nm] / X.n_nodes_usable.replace(0, np.nan)
    X["shareAll_" + nm] = X["n_" + nm] / X.n_nodes_all.replace(0, np.nan)
X["frac_nodes_with_ds"] = X.n_nodes_usable / X.n_nodes_all.replace(0, np.nan)
X["nodes_per_1000_genes"] = 1000 * X.n_nodes_all / X.genes
X["paralog_gene_frac"] = X.genes_in_families / X.genes
X["ident_genes_per_1000"] = 1000 * X.genes_in_identical_groups / X.genes
X["share_2gene_fam"] = X.f2 / X.n_families.replace(0, np.nan)
X["frac_genes_in_fam_gt50"] = X.genes_in_big_families_gt50 / X.genes_in_families.replace(0, np.nan)
def clade_med(col):
    out = pd.Series(np.nan, index=X.index)
    for lvl in ["phylum", "class", "order"]:   # most specific level with >= 5 genomes wins
        m = X.groupby(lvl)[col].transform("median"); n = X.groupby(lvl)[col].transform("size")
        out = np.where(n >= 5, m, out)
    return out
X["shareAll_BEYOND"] = 1 - X.frac_nodes_with_ds
X["genes_ratio_clade"] = X.genes / clade_med("genes")
# TE-domain proteins (Pfam search of tables/pfam.parquet, see pfam_te_counts.py)
te = pd.read_csv(f"{D}/pfam_te_counts.tsv", sep="\t")
X = X.merge(te, on="species_prefix", how="left")
X["te_per_1000_genes"] = 1000 * X.te_domain_proteins / X.genes
X["te_frac_of_pfam_proteins"] = X.te_domain_proteins / X.proteins_with_pfam.replace(0, np.nan)
X.to_csv(f"{D}/genome_table.tsv.gz", sep="\t", index=False)
print(len(X), "genomes;", X.shape[1], "columns; with BUSCO dup", int(X.busco_duplicated_pct.notna().sum()), "; with Pfam/TE", int(X.te_domain_proteins.notna().sum()))
