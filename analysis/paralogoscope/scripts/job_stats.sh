#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16G --time 02:00:00 -p short -J pf_stats
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/stats.%j.log
set -e
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore
/usr/bin/python3.12 $E/scripts/stats.py
Rscript $E/scripts/varcomp.R
