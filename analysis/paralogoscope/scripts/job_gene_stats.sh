#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16G --time 02:00:00 -p short -J gene_stats
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/gene_stats.%j.log
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot
/usr/bin/python3.12 $P/eval/explore/scripts/gene_stats.py $P/nf_full/samples_all.csv $P/input/cds $P/eval/explore/data/gene_stats.tsv 8
echo "done $(wc -l < $P/eval/explore/data/gene_stats.tsv) lines"
