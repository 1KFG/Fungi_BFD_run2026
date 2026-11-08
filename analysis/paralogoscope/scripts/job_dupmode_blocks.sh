#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 64G --time 02:00:00 -p short -J dupmode_blocks
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/dupmode_blocks_run.%j.log
# usage: DUP_TAG=_w03 sbatch --export=ALL job_dupmode_blocks.sh   (default tag: empty = window dS <= 0.1)
S=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/scripts
/usr/bin/python3.12 $S/dup_mode_tests.py && /usr/bin/python3.12 $S/dup_block_tests.py
