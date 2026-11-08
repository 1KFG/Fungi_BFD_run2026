#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 110G --time 06:00:00 -p batch -J dupmode
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/dupmode.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/dup_mode.py
