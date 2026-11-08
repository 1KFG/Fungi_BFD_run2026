#!/bin/bash -l
#SBATCH -p stajichlab -c 16 --mem 64G --time 6:00:00 -J introns_db
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/logs/strip_introns_rebuild_db.%j.log
# 1) tables/gene_introns.parquet without seq (row count must be equal);
# 2) db/BFD.duckdb rebuilt as VIEWs over tables/ (replaces the 394 GB copy DB),
#    with the build script from the merge-csv-list branch (Fungi_BFD PR #59).
set -euo pipefail
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
NF=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_wt_csvlist/nextflow
module load duckdb
cd $R/tables
duckdb -c "SET threads=16; SET memory_limit='48GB'; SET temp_directory='${SCRATCH:?}';
COPY (SELECT * EXCLUDE (seq) FROM 'gene_introns.parquet') TO 'gene_introns.noseq.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);"
a=$(duckdb -csv -noheader -c "SELECT sum(num_rows) FROM parquet_file_metadata('gene_introns.parquet')")
b=$(duckdb -csv -noheader -c "SELECT sum(num_rows) FROM parquet_file_metadata('gene_introns.noseq.parquet')")
echo "gene_introns rows: original $a, noseq $b; sizes $(stat -c %s gene_introns.parquet) -> $(stat -c %s gene_introns.noseq.parquet)"
[ "$a" = "$b" ] || { echo "ERROR: row count differs"; rm -f gene_introns.noseq.parquet; exit 1; }
mv gene_introns.noseq.parquet gene_introns.parquet
duckdb -c "DESCRIBE SELECT * FROM 'gene_introns.parquet'"

cd $SCRATCH
TABLES=$R/tables DBDIR=$SCRATCH BIN_DIR=$NF/bin SCHEMA=$NF/../sql/table_schema.json bash $NF/bin/build_BFD_duckDB.sh
ls -l $SCRATCH/BFD.duckdb
cp $SCRATCH/BFD.duckdb $R/db/BFD.duckdb.new && mv $R/db/BFD.duckdb.new $R/db/BFD.duckdb
ls -l $R/db/
duckdb -readonly $R/db/BFD.duckdb -c "SELECT count(*) AS introns FROM gene_introns; SELECT count(*) AS genomes FROM v_species_summary;"
