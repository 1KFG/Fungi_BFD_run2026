#!/usr/bin/bash -l
#SBATCH -p short -c 8 -n 1 -N 1 --mem 8gb --time 2:00:00
#SBATCH -J query_unknowns_mash
#SBATCH -o results/query_unknowns_mash/logs/%x.%j.out
#SBATCH -e results/query_unknowns_mash/logs/%x.%j.err

# Screens samples.csv genomes with no GENUS-level taxonomy against the
# existing results/mash_sketch/ reference cache via mash dist, reporting
# candidate nearby lineages to results/query_unknowns_mash/.
#
# Complements query_orphans.sh (skani, per-CLASS/ORDER-grouped, via the
# query_ANI nextflow workflow): this is a flat, unscoped screen of every
# GENUS-blank genome against the whole mash sketch cache in one job.

mkdir -p results/query_unknowns_mash/logs

module load mash

python3 query_unknowns_mash.py \
    --samples samples.csv \
    --sketch-dir results/mash_sketch \
    --genome-dir input_clean_genomes \
    --outdir results/query_unknowns_mash \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
