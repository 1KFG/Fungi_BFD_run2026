#!/usr/bin/env Rscript
# Nested random-effects model on the taxonomy (no tree available): how much of the variation of each paralog metric sits at phylum / class / order / family / genus / genome,
# with and without assembly covariates. Output: tables/varcomp.tsv (variance shares) and tables/varcomp_fixed.tsv (fixed effects of the covariates).
suppressPackageStartupMessages({library(lme4); library(data.table)})
D <- "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore"
X <- fread(file.path(D, "data/genome_table_final.tsv.gz"))
lg <- function(p) { p <- pmin(pmax(p, 1e-3), 1 - 1e-3); log(p / (1 - p)) }
X[, `:=`(y_H = lg(shareAll_H), y_Y = lg(shareAll_Y), y_MID = lg(shareAll_MID), y_BEYOND = lg(shareAll_BEYOND),
         y_nodes = log(nodes_per_1000_genes), y_paralog = lg(paralog_gene_frac), y_2gene = lg(share_2gene_fam), y_ident = log1p(ident_genes_per_1000))]
X[, `:=`(c_contig = scale(log(contig_count)), c_size = scale(log(total_length_bp)), c_busco = scale(sqrt(busco_duplicated_pct)), c_gc = scale(gc_pct),
         c_genes = scale(log(genes)), c_te = scale(log1p(te_per_1000_genes)))]
for (f in c("phylum", "class", "order", "family_tax", "genus")) X[[f]] <- factor(X[[f]])
resp <- c(H = "y_H", Y = "y_Y", MID = "y_MID", BEYOND = "y_BEYOND", nodes_per_1000 = "y_nodes", paralog_frac = "y_paralog", two_gene_fam = "y_2gene", identical = "y_ident")
form_re <- "(1 | phylum) + (1 | phylum:class) + (1 | phylum:class:order) + (1 | phylum:class:order:family_tax) + (1 | phylum:class:order:family_tax:genus)"
covs <- "c_contig + c_size + c_busco + c_gc + c_genes + c_te"
vc <- list(); fx <- list()
for (nm in names(resp)) {
  y <- resp[[nm]]
  d <- X[is.finite(get(y)) & !is.na(c_contig) & !is.na(c_size) & !is.na(c_busco) & !is.na(c_gc) & !is.na(c_te) & (if (nm %in% c("H", "Y", "MID", "BEYOND")) n_nodes_all >= 100 else TRUE)]
  for (adj in c("no covariates", "with covariates")) {
    f <- as.formula(paste(y, "~", if (adj == "with covariates") covs else "1", "+", form_re))
    m <- tryCatch(lmer(f, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))), error = function(e) NULL)
    if (is.null(m)) next
    v <- as.data.frame(VarCorr(m)); v <- v[, c("grp", "vcov")]
    lev <- c("phylum:class:order:family_tax:genus" = "genus", "phylum:class:order:family_tax" = "family", "phylum:class:order" = "order", "phylum:class" = "class", "phylum" = "phylum", "Residual" = "genome (residual)")
    v$level <- lev[v$grp]; v$share <- v$vcov / sum(v$vcov); v$response <- nm; v$adjustment <- adj; v$n_genomes <- nrow(d)
    vc[[length(vc) + 1]] <- v[, c("response", "adjustment", "level", "vcov", "share", "n_genomes")]
    if (adj == "with covariates") { cf <- as.data.frame(summary(m)$coefficients); cf$term <- rownames(cf); cf$response <- nm; fx[[length(fx) + 1]] <- cf }
    cat(nm, adj, "n =", nrow(d), "\n")
  }
}
fwrite(rbindlist(vc), file.path(D, "tables/varcomp.tsv"), sep = "\t"); fwrite(rbindlist(fx), file.path(D, "tables/varcomp_fixed.tsv"), sep = "\t")
cat("done\n")
