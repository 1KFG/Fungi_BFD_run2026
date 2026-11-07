#!/bin/bash -l
#SBATCH -p stajichlab -c 8 --mem 32G --time 6:00:00 -J strip_seq
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/logs/strip_gene_seq_cols.%j.log
# Write gene_transcripts / gene_proteins without the sequence columns
# (CDS_sequence, peptide) as *.noseq.parquet next to the originals. The
# originals are replaced only after seqdb/ (bgzip + faidx copies) is verified.
set -euo pipefail
module load duckdb
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
duckdb -c "SET threads=8; SET memory_limit='24GB'; SET temp_directory='${SCRATCH:?}';
COPY (SELECT * EXCLUDE (CDS_sequence) FROM 'gene_transcripts.parquet') TO 'gene_transcripts.noseq.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * EXCLUDE (peptide) FROM 'gene_proteins.parquet') TO 'gene_proteins.noseq.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);"
for t in gene_transcripts gene_proteins; do
  a=$(duckdb -csv -noheader -c "SELECT sum(num_rows) FROM parquet_file_metadata('$t.parquet')")
  b=$(duckdb -csv -noheader -c "SELECT sum(num_rows) FROM parquet_file_metadata('$t.noseq.parquet')")
  echo "$t rows: original $a, noseq $b, size $(stat -c %s $t.noseq.parquet)"
  [ "$a" = "$b" ] || { echo "ERROR: row count differs for $t"; exit 1; }
done
