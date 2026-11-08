#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 64G --time 02:00:00 -p short -J dupmode_tests
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/dupmode_tests.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/dup_mode_tests.py
