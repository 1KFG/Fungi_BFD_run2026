#!/usr/bin/bash
# Submit the whole exploration pipeline with SLURM dependencies. Run from the login node: bash run_all.sh
# Each job script writes its log to ../data/<name>.<jobid>.log. Edit the input table path in build_features2.py first when the final tables are published.
S=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts
cd $S
j() { sbatch --parsable ${2:+--dependency=afterok:$2} "$1"; }
gs=$(j job_gene_stats.sh)                   # per-genome CDS statistics
f2=$(j job_features2.sh)                    # bins, bootstrap CIs, GMM peaks (reads the node table)
hi=$(j job_hist.sh)                         # log-dS histograms
pt=$(j job_pfam_te.sh)                      # TE Pfam counts per genome
cm=$(j job_compare.sh)                      # yn00 vs codeml (wave-1 genomes)
# genome table needs taxonomy + gene stats + features + TE counts
bt=$(sbatch --parsable --dependency=afterok:$gs:$f2:$pt --wrap="/usr/bin/python3.12 $S/taxonomy.py && /usr/bin/python3.12 $S/build_table.py" -J pf_table -p short --mem 16G -c 2 --time 01:00:00)
st=$(sbatch --parsable --dependency=afterok:$bt:$hi:$cm job_stats.sh)    # clade tables, flags, controls, peak catalogue
vc=$(sbatch --parsable --dependency=afterok:$st job_varcomp.sh)         # nested random-effects model (R, lme4)
sh=$(sbatch --parsable --dependency=afterok:$st job_shape.sh)           # shape-based peaks + BUSCO cross-reference (needs genome_table_final)
ge=$(sbatch --parsable --dependency=afterok:$sh job_gene_enrich.sh)     # gene-level Pfam/CAZy enrichment (needs shape tiers for the clean set)
gp=$(sbatch --parsable --dependency=afterok:$ge job_gene_plots.sh)
pl=$(sbatch --parsable --dependency=afterok:$vc:$sh:$gp job_plots.sh)
sbatch --dependency=afterok:$pl job_report.sh
