#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 64G --time 04:00:00 -p batch -J upo_perm
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/upo4.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/upo_perm.py
