#!/bin/bash -l
#SBATCH -p exfab -A exfab -N 1 -n 1 -c 8 --mem 48G --gres=gpu:ada6000:1 --time 6:00:00 -J dtm_bench
# usage: sbatch -o DIR/run_tools.%j.log run_tools.sh DIR   (DIR holds bench.fa)
# Runs classic TMHMM 2.0, SignalP 6 (euk, fast; same settings as the BFD SIGNALP module)
# and DeepTMHMM 1.0 (private fast image, --topo-batch 1 = same output as predict.py) on
# DIR/bench.fa. Outputs: DIR/tmhmm_short.tsv, DIR/signalp/, DIR/deeptmhmm/.
set -euo pipefail
D=$(readlink -f $1); cd $D
NCPU=${SLURM_CPUS_PER_TASK:-8}
echo "host $(hostname) $(nvidia-smi --query-gpu=name --format=csv,noheader)"
module load tmhmm
t0=$(date +%s); tmhmm --short --noplot < bench.fa > tmhmm_short.tsv; echo "tmhmm $(( $(date +%s)-t0 )) s"
# same setup as the BFD 'signalp' label (conf/profile_BFD.config): the module
# needs conda's shell hook first, or it fails with "Run 'conda init'".
export CONDA_NO_PLUGINS=true
module load miniconda3
eval "$(conda shell.bash hook)"
module load signalp/6-gpu
rm -rf signalp; t0=$(date +%s)
signalp6 -od signalp -org euk --mode fast -format none -fasta bench.fa --write_procs $NCPU -bs 100 > signalp.log 2>&1
echo "signalp6 $(( $(date +%s)-t0 )) s"
module load apptainer
export TMPDIR=${SCRATCH:?}
S=/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0-fast.sif
rm -rf deeptmhmm; t0=$(date +%s)
apptainer exec --nv --bind $D:$D,$SCRATCH:$SCRATCH $S env NCPU=$NCPU bash -c "cd /opt/deeptmhmm && python3 predict_fast.py --fasta $D/bench.fa --output-dir $SCRATCH/dtm --topo-batch 1" > deeptmhmm.log 2>&1
mkdir -p deeptmhmm; cp $SCRATCH/dtm/TMRs.gff3 $SCRATCH/dtm/predicted_topologies.3line deeptmhmm/
echo "deeptmhmm $(( $(date +%s)-t0 )) s"
