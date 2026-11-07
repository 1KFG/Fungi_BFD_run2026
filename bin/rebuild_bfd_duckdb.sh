#!/bin/bash -l
#SBATCH -p stajichlab -c 16 --mem 128G --time 24:00:00 -J bfd_duckdb
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/logs/rebuild_bfd_duckdb.%j.log
# Rebuild db/BFD.duckdb from tables/*.parquet outside Nextflow, with the same
# script BUILD_DUCKDB runs (Fungi_BFD/nextflow/bin/build_BFD_duckDB.sh).
# Built in $SCRATCH, then copied to db/ (BFD.duckdb.new -> BFD.duckdb).
set -euo pipefail
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
NF=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow
module load duckdb
cd ${SCRATCH:?}
TABLES=$R/tables DBDIR=$SCRATCH BIN_DIR=$NF/bin SCHEMA=$NF/../sql/table_schema.json bash $NF/bin/build_BFD_duckDB.sh
ls -l $SCRATCH/BFD.duckdb
cp $SCRATCH/BFD.duckdb $R/db/BFD.duckdb.new && mv $R/db/BFD.duckdb.new $R/db/BFD.duckdb
duckdb -readonly $R/db/BFD.duckdb -c "SELECT table_name, table_type FROM information_schema.tables ORDER BY 1" | head -60
