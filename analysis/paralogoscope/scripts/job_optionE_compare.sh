#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 48G --time 01:00:00 -p short -J optE_compare
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/optionE_compare.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/optionE_compare.py
