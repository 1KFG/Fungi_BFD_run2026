#!/bin/bash -l
#SBATCH -p stajichlab -c 16 --mem 96G --time 4:00:00 -J sp_check
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/logs/check_swissprot_new.%j.log
# Compare the new SwissProt tables (Fungi_BFD PR #59 design) with the old tables/swissprot.parquet.
set -euo pipefail
module load duckdb
T=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
duckdb -c "SET threads=16; SET memory_limit='80GB'; SET temp_directory='${SCRATCH:?}';
SELECT (SELECT sum(num_rows) FROM parquet_file_metadata('$T/swissprot.parquet')) AS old_rows,
       (SELECT sum(num_rows) FROM parquet_file_metadata('$T/swissprot_hits/swissprot_hits.*.parquet')) AS new_rows;
CREATE TEMP TABLE o AS SELECT DISTINCT protein_id FROM read_parquet('$T/swissprot.parquet');
CREATE TEMP TABLE n AS SELECT DISTINCT protein_id FROM read_parquet('$T/swissprot_hits/swissprot_hits.*.parquet');
SELECT (SELECT count(*) FROM o) AS old_proteins, (SELECT count(*) FROM n) AS new_proteins,
       (SELECT count(*) FROM o ANTI JOIN n USING (protein_id)) AS only_old,
       (SELECT count(*) FROM n ANTI JOIN o USING (protein_id)) AS only_new;
SELECT count(*) AS best_rows, count(DISTINCT protein_id) AS best_proteins FROM '$T/swissprot_best.parquet';
SELECT count(*) AS hits_null_sp_id FROM read_parquet('$T/swissprot_hits/swissprot_hits.*.parquet') WHERE sp_id IS NULL;
SELECT count(*) AS subject_rows FROM '$T/swissprot_subject.parquet';
SELECT count(*) AS summary_rows FROM '$T/swissprot_subject_summary.parquet';
"
# genomes only in old (by protein prefix), sample
duckdb -c "SET threads=16; SET memory_limit='80GB'; SET temp_directory='$SCRATCH';
WITH o AS (SELECT DISTINCT split_part(protein_id,'_',1) g FROM read_parquet('$T/swissprot.parquet')),
     n AS (SELECT DISTINCT split_part(protein_id,'_',1) g FROM read_parquet('$T/swissprot_hits/swissprot_hits.*.parquet'))
SELECT (SELECT count(*) FROM o) old_genomes, (SELECT count(*) FROM n) new_genomes,
       (SELECT count(*) FROM o ANTI JOIN n USING (g)) only_old, (SELECT count(*) FROM n ANTI JOIN o USING (g)) only_new;"
duckdb -readonly /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/db/BFD.duckdb -c "SELECT view_name FROM duckdb_views() WHERE view_name LIKE '%swissprot%'; SELECT * FROM swissprot_best LIMIT 2;"
