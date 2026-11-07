#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 4:00:00 -p batch
#SBATCH --job-name=merge_tf
#SBATCH --output=logs/merge_tf.%j.log

# Standalone MERGE_TF: builds tables/tf_inventory.parquet from the Pfam
# domtblout files, without a Nextflow run. Same steps as
# modules/BFD/MERGE_TF/main.nf:
#   1. keep files whose genome is in samples.csv minus suppress.txt
#      (bin/merge_filter.py, the twin of mergeFileAllowed() in utils.nf)
#   2. bin/merge_tf_domains.py, run on 8 chunks in parallel
#   3. CSV -> parquet with duckdb, then an atomic move into tables/
# Submit from BFD_function/:  sbatch run_merge_tf.sh
#
# Interim. MERGE_TF belongs in the BFD pipeline merge run (run_tf: true with
# params_merge_only.yaml). Retire this script once that run covers it
# (Fungi_BFD todo T-039).

set -euo pipefail

RUN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/BFD_function
NF=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow
PFAM=$RUN/results/function/pfam_hmmscan
TABLES=$RUN/tables
NCHUNK=${SLURM_CPUS_PER_TASK:-8}
TMP="${SCRATCH:?SCRATCH not set}/merge_tf"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load duckdb
# merge_filter.py imports pick_representative_strain, which needs the Python
# duckdb package. The biopython module's python lacks it (job 31806143 failed).
PY=/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin/python3
"$PY" -c 'import duckdb' || { echo "ERROR: $PY has no duckdb module" >&2; exit 1; }

mkdir -p "$TMP"
cd "$RUN"

# 1. Filtered input list.
find -H "$PFAM" -type f -name '*.domtblout.gz' > "$TMP/all.txt"
"$PY" - "$TMP/all.txt" "$TMP/kept.txt" <<EOF
import os, sys
sys.path.insert(0, "$NF/bin")
from merge_filter import allowed_keys, file_allowed
keys = allowed_keys("$RUN/samples.csv", "$RUN/suppress.txt")
src, dst = sys.argv[1], sys.argv[2]
paths = [l.strip() for l in open(src) if l.strip()]
kept = [p for p in paths if file_allowed(keys, os.path.basename(p))]
open(dst, "w").write("\n".join(sorted(kept)) + "\n")
print(f"domtblout files: {len(paths)} found, {len(kept)} kept, {len(paths)-len(kept)} dropped", file=sys.stderr)
EOF

# 2. Parallel extraction, one CSV per chunk.
# Each chunk is passed as one argument list (about 3,000 paths, well under
# ARG_MAX). An over-long list fails with "Argument list too long", and the
# per-PID wait below turns any failure into a job failure.
split -n l/"$NCHUNK" -d "$TMP/kept.txt" "$TMP/chunk."
pids=()
for c in "$TMP"/chunk.??; do
    mapfile -t files < "$c"
    "$PY" "$NF/bin/merge_tf_domains.py" \
        --tf-domains "$NF/assets/fungal_tf_pfam_domains.csv" -o "$c.csv" "${files[@]}" &
    pids+=($!)
done
for p in "${pids[@]}"; do
    wait "$p" || { echo "ERROR: a merge_tf_domains.py chunk failed" >&2; exit 1; }
done

# 3. Concatenate (one header) and convert.
{ head -1 "$TMP/chunk.00.csv"; for f in "$TMP"/chunk.??.csv; do tail -n +2 "$f"; done; } > "$TMP/tf_inventory.csv"
duckdb -c "COPY (SELECT * FROM read_csv_auto('$TMP/tf_inventory.csv', sample_size=-1)) TO '$TMP/tf_inventory.parquet' (FORMAT PARQUET);"
duckdb -c "SELECT count(*) AS rows, count(DISTINCT species_prefix) AS genomes, count(DISTINCT tf_family) AS tf_families FROM '$TMP/tf_inventory.parquet';"

cp "$TMP/tf_inventory.parquet" "$TABLES/.tf_inventory.parquet.tmp"
mv "$TABLES/.tf_inventory.parquet.tmp" "$TABLES/tf_inventory.parquet"
ls -la "$TABLES/tf_inventory.parquet"
rm -rf "$TMP"
