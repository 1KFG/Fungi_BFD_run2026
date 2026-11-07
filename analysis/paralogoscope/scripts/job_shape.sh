#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16G --time 01:00:00 -p short -J pf_shape
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/shape.%j.log
S=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts
/usr/bin/python3.12 $S/shape_peaks.py && /usr/bin/python3.12 $S/shape_plots.py
