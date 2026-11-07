#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16G --time 01:00:00 -p short -J cmp_yn00_codeml
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/compare.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/compare_codeml_yn00.py
