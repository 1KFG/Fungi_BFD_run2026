#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16gb --time 2-00:00:00 -p batch
#SBATCH --job-name=nxf_bfd_merge
#SBATCH --output=logs/bfd_merge.%j.log

# Merge-only run (params_merge_only.yaml): refresh stale tables/*.parquet + DuckDB.
# Points explicitly at the worktree's nextflow/main.nf + nextflow.config
# (branch bfd-duckdb-datalake) rather than production Fungi_BFD/nextflow,
# since this run must exercise the branch's code changes, not main.
#
# Run from this directory (BLD_2), a dedicated build/test project root
# separate from do_annotation_asco's live in-progress pipeline.

set -euo pipefail

module load nextflow

export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}

mkdir -p logs/nextflow

WORKTREE=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
PARAMS_FILE=$WORKTREE/nextflow/param_files/params_merge_only.yaml

if [ ! -f "$PARAMS_FILE" ]; then
    echo "ERROR: $PARAMS_FILE not found" >&2
    exit 1
fi

NXF_OPTS="-Xms512m -Xmx12g" \
nextflow run "$WORKTREE/nextflow/main.nf" \
    -c "$WORKTREE/nextflow/nextflow.config" \
    -profile BFD \
    -params-file "$PARAMS_FILE" \
    -resume \
    "$@"
