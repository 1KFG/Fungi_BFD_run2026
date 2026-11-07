#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8G --time 00:30:00 -p short -J pf_geneplots
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/gene_plots.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/gene_plots.py
