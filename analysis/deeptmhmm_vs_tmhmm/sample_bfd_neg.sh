#!/bin/bash
#SBATCH -p batch -c 4 --mem 32G --time 3:00:00 -J dtm_neg -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/sample_bfd_neg.%j.log
set -euo pipefail
module load duckdb
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/sample_bfd_neg.py
sbatch -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/bfd_sample_neg/run_tools.%j.log /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/run_tools.sh /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/bfd_sample_neg
