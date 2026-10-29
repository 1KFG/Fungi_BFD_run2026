#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16gb --time 1-00:00:00 -p batch
#SBATCH --job-name=nxf_unannot_merge
#SBATCH --output=logs/unannot_test_merge.%j.log
# 5-genome test 2: same genomes, MERGE_* on, rebased code (30f8afa).
# Own launchDir (sibling of BFD_function) so it does not touch the live run's
# work/, results/ or session lock. Nothing here uses $SCRATCH.
set -euo pipefail
module load nextflow
export NXF_SINGULARITY_CACHEDIR=/bigdata/stajichlab/shared/singularity_cache
export SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
export APPTAINER_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
WORKTREE=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
NXF_OPTS="-Xms512m -Xmx8g" nextflow run $WORKTREE/nextflow/main.nf \
  -c $WORKTREE/nextflow/nextflow.config -profile BFD \
  -resume -params-file $PWD/params_unannot_test_merge.yaml
