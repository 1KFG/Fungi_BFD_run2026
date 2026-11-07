#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 32G --time 02:00:00 -p short -J pf_features
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/features.%j.log
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot
T=$P/nf_full/tables_yn00
/usr/bin/python3.12 $P/eval/explore/scripts/build_features.py $T/wgd.ks.nodes.parquet $T/wgd.ks.summary.parquet $P/eval/explore/data 8
