#!/usr/bin/bash -l
#SBATCH -J optE_codeml -p batch -c 16 --mem 32G --time 06:00:00 --array=1-10
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/option_e_codeml/logs_array.%A_%a.log
O=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/option_e_codeml
bash $O/run_bundle.sh $O/$(printf 'bundle_%03d.txt' $SLURM_ARRAY_TASK_ID)
