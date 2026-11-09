#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 64G --time 02:00:00 -p short -J upo_strat
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/upo3.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/upo_strat.py
