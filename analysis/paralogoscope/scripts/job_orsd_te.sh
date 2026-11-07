#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 80G --time 02:00:00 -p short -J orsd_te
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/orsd_te.%j.log
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts/domain_te_context.py PF12013 OrsD
