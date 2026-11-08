#!/bin/bash -l
#SBATCH -p epyc -c 4 --mem 24G --time 6:00:00 -J merge_prof3
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile/profile_duckdb.%j.log
# DuckDB-only merges (bin/duckdb_merge.py + sql/merge/*.sql) on the ~1,018-genome
# subsets, against the Python paths measured in profile_func.31780317.log.
set -euo pipefail
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile
NF=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_wt_csvlist/nextflow
cd $P; module load duckdb
T=${SCRATCH:?}; O=$P/out_func
run() { local name=$1; shift; /usr/bin/time -f "%e s wall, %M KB maxRSS" -o $O/time_$name.txt bash -c "$*"; echo "$name: $(cat $O/time_$name.txt)"; }
ck() { duckdb -csv -noheader -c "SELECT count(*), sum(hash(COLUMNS(*))) FROM '$1'" | md5sum | cut -c1-12; }
M="--threads 4 --memory-limit 18GB --temp-dir $T"
echo "host $(hostname)"

# PFAM: DuckDB parses domtblout directly; compare with the Python CSV -> typed read (pf1)
run PF3_duckdb_only "python3 $NF/bin/duckdb_merge.py --sql $NF/sql/merge/pfam.sql --files pfam_hmmscan_domtbloutgz.manifest --out $O/pf3.parquet $M"
echo "pfam pf1 (python+types): $(ck $O/pf1.parquet)  pf3 (duckdb only): $(ck $O/pf3.parquet)"
ls -l $O/pf1.parquet $O/pf3.parquet | awk '{print $9, $5}'
# PF2 redo: the Python converter as 4 processes on 4 slices, then a typed read
cut -f1 pfam_hmmscan_domtbloutgz.manifest > $T/pf.list; split -n l/4 -d $T/pf.list $T/pfslice.
run PF2_python_x4 "for s in $T/pfslice.0[0-3]; do python3 $NF/bin/pfamtbl_to_long.py --outfile \$s.csv \$(cat \$s) & done; wait"

# SWISSPROT redesign: narrow hits with sp_id, one task over the 1,018 genomes
run SP3_hits_redesign "python3 $NF/bin/duckdb_merge.py --sql $NF/sql/merge/swissprot_hits.sql --files swissprot_blasttabgz.manifest --out $O/sp3.parquet --param subject=\"'/rhome/jstajich/.claude/jobs/bacad32f/tmp/subject_test.parquet'\" $M"
ls -l $O/sp0.parquet $O/sp3.parquet | awk '{print $9, $5}'
duckdb -c "
SELECT (SELECT count(*) FROM '$O/sp0.parquet') AS rows_old, (SELECT count(*) FROM '$O/sp3.parquet') AS rows_new,
       (SELECT count(DISTINCT protein_id) FROM '$O/sp3.parquet') AS proteins,
       (SELECT max(rk) FROM '$O/sp3.parquet') AS max_rank;
-- value check on a sample of old rows: identity, bitscore, alignment, query coverage
WITH s AS (SELECT * FROM '$O/sp0.parquet' USING SAMPLE 200000 ROWS),
n AS (SELECT h.*, sub.accession FROM '$O/sp3.parquet' h JOIN '/rhome/jstajich/.claude/jobs/bacad32f/tmp/subject_test.parquet' sub USING (sp_id))
SELECT count(*) AS matched,
       sum((round(s.pident*10) = n.pident_x10)::INT) AS pident_same,
       sum((round(s.qcovhsp*10) = n.qcov_x10)::INT) AS qcov_same,
       sum((abs(s.bitscore - n.bitscore) < 0.05)::INT) AS bitscore_same,
       sum((s.aln_length = n.aln_length)::INT) AS aln_same,
       sum((s.evalue = 0 AND n.neg_log10_evalue = 400 OR abs(-log10(nullif(s.evalue,0)) - n.neg_log10_evalue) < 1e-3)::INT) AS evalue_same
FROM s JOIN n ON n.protein_id = s.protein_id AND n.accession = s.swissprot_acc AND n.q_start = s.q_start AND n.s_start = s.s_start;"
