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

# ---------------------------------------------------------------------------
# Result of job 31806285 (2026-11-06): 23,410 domtblout files, 10,946,927
# rows, 23,394 genomes, 21 TF families. 16 genomes have zero TF hits
# (Fungi_BFD todo T-038). In all 16 the protein set has almost no BUSCO
# genes (complete 0.0-0.7%), so gene prediction failed for each one.
# BUSCO = fungi_odb12, C/M = complete/missing %. Proteins = predict_results
# protein count. Assembly stats from tables/asm_stats.parquet.
#
# Group A: good assembly, failed annotation (re-predict)
# LOCUSTAG  Species                              Proteins  gBUSCO C/M  pBUSCO C/M  Mb     Contigs  N50kb
# F00081D8  Amanita vidua AV01                   5834      98.1/1.4    0.4/95.0    56.5   110      2303
# F2E2A438  Amanita subjunquillea LH00001        3073      97.4/1.9    0.3/98.7    53.2   119      839
# F636EC73  Calvatia gigantea CGMCC 5.9          785       98.0/1.9    0.7/99.1    40.1   39       1468
# FAB5C0F4  Lactarius quietus S23C               6084      98.0/1.8    0.2/98.6    115.8  2746     90
# FBDE1BFE  Russula nigricans                    2322      99.2/0.5    0.2/98.0    45.0   11       4242
# FC0DB305  Tricholoma bakamatsutake SF-Tf05     4769      99.7/0.2    0.1/99.8    142.2  14       12427
# FD6195C5  Tricholoma bakamatsutake 20200826-1  3744      99.8/0.1    0.1/99.6    140.4  128      9479
#
# Group B: fragmented or partial assembly (gBUSCO 6.7-74.9%, N50 2.8-14.9 kb)
# F1E1EE23  aff. Porphyrellus sp. UT-M0000367    1844      74.9/16.7   0.3/99.6    38.8   4916     11.5
# F77192FB  Boletus austroedulis NY1393612       1505      72.2/14.8   0.2/99.4    40.5   5840     9.3
# FB348870  Amanita bisporigera Ingham01         2931      70.9/17.5   0.1/99.4    71.1   7368     14.9
# FDAED364  Boletaceae sp. KM206711              2611      41.6/43.6   0.1/99.6    39.7   5472     10.5
# FCF521AA  Boletaceae sp. KM206743              2479      34.6/52.8   0.0/99.6    30.9   5002     8.1
# FA58BA6D  Boletaceae sp. KM206737              1036      30.7/55.9   0.1/99.8    25.9   4533     7.3
# FD8664C9  Austroboletus festivus MICH11833     1712      6.7/86.4    0.0/99.2    9.1    3144     2.8
#
# Group C: special cases
# FBFB4B38  Hamiltosporidium tvaerminnensis      852       12.7/84.1   0.0/99.6    21.6   17       1444  (microsporidian; fungi_odb12 lacks many of its markers)
# FCF75579  Massospora cicadina MCPNR19          415       71.7/24.0   0.0/100.0   1488.9 19691    140   (1.49 Gb genome)
