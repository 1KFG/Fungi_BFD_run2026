#!/bin/bash -l
#SBATCH -p epyc -c 4 --mem 24G --time 6:00:00 -J merge_prof
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile/profile.%j.log
# Profile ways to build gene_exons.parquet from per-genome .gene_exons.csv.gz files.
# Input: exons_1000.manifest (every 23rd genome of the 2026-11-06 MERGE_GENE_STATS
# manifest: ~1000 genomes). Same CPU/memory as the 'merge' label (4 CPUs, 24 GB).
# Outputs are written under this folder or $SCRATCH; nothing in BFD_function is touched.
set -euo pipefail
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile
cd $P; module load duckdb
TMPD=${SCRATCH:?}
M=$P/exons_1000.manifest
COLS="{'exon_id':'VARCHAR','transcript_id':'VARCHAR','order':'BIGINT','chrom':'VARCHAR','start':'BIGINT','end':'BIGINT','strand':'BIGINT','GC_content':'DOUBLE'}"
cut -f1 $M | awk '{printf "%s\x27%s\x27", (NR>1?",":""), $0}' > $TMPD/filelist.sql   # 'a','b',...
echo "host $(hostname); $(wc -l < $M) files; $(cut -f3 $M | paste -sd+ | bc) gz bytes"
run() {  # name, command
  local name=$1; shift
  /usr/bin/time -f "%e s wall, %M KB maxRSS" -o $P/time_$name.txt bash -c "$*"
  echo "$name: $(cat $P/time_$name.txt)"
}
concat() {  # current module loop: zcat each file, keep the first header
  first=1
  while IFS=$'\t' read -r f _m _s; do
    if [ "$first" = 1 ]; then zcat "$f"; first=0; else zcat "$f" | tail -n +2; fi
  done < $M > $TMPD/exons.csv
}
export -f concat; export M TMPD COLS
if [ "${SKIP_AB:-0}" != 1 ]; then
# A: current (2026-11-06) module: concat -> plain CSV on SCRATCH -> read_csv_auto(sample_size=-1) -> parquet on GPFS
run A_concat_csv "concat"
run A_duckdb_autodetect "duckdb -c \"COPY (SELECT * FROM read_csv_auto('$TMPD/exons.csv', sample_size=-1)) TO '$P/A.parquet' (FORMAT PARQUET);\""
# B: same CSV, explicit column types (no full-file type scan)
run B_duckdb_types "duckdb -c \"COPY (SELECT * FROM read_csv('$TMPD/exons.csv', header=true, columns=$COLS)) TO '$P/B.parquet' (FORMAT PARQUET);\""
rm -f $TMPD/exons.csv
fi
# C/D: no concatenation: DuckDB reads the .csv.gz list directly with explicit types.
# The file list goes into a SQL file (one shell argument is limited to 128 KB).
sql() {  # $1 = output parquet
  { printf "COPY (SELECT * FROM read_csv(["; cat $TMPD/filelist.sql
    printf "], header=true, columns=%s)) TO '%s' (FORMAT PARQUET);\n" "$COLS" "$1"; } > $TMPD/build.sql
}
sql $P/C.parquet; run C_duckdb_gzlist_gpfs "duckdb < $TMPD/build.sql"
sql $TMPD/D.parquet; run D_duckdb_gzlist_scratch "duckdb < $TMPD/build.sql"
run D_copy_back "cp $TMPD/D.parquet $P/D.parquet"
# same rows in every output? (row count + order-independent checksum)
for x in A B C D; do
  echo "$x: $(duckdb -csv -noheader -c "SELECT count(*), sum(hash(exon_id, transcript_id, \"order\", chrom, start, \"end\", strand, GC_content)) FROM '$P/$x.parquet'") $(stat -c %s $P/$x.parquet) bytes"
done
