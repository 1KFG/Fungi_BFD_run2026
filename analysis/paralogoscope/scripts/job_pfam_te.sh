#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16G --time 02:00:00 -p short -J pfam_te
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/pfam_te.%j.log
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot
/usr/bin/python3.12 $P/eval/explore/scripts/pfam_te_counts.py $P/eval/explore/data/pfam_te_counts.tsv
wc -l $P/eval/explore/data/pfam_te_counts.tsv
