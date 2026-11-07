#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16G --time 01:00:00 -p short -J ds_hist
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/hist.%j.log
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot
/usr/bin/python3.12 $P/eval/explore/scripts/build_hist.py $P/nf_full/tables_yn00/wgd.ks.nodes.parquet $P/eval/explore/data/ds_hist.tsv.gz 4
