#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 120G --time 08:00:00 -p batch -J dupmode_w03
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/dupmode_w03.%j.log
export DUP_RECENT=0.3 DUP_NBR=0.5 DUP_TAG=_w03
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/dup_mode.py
