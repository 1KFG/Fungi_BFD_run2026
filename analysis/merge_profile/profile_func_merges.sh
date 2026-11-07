#!/bin/bash -l
#SBATCH -p epyc -c 4 --mem 24G --time 8:00:00 -J merge_prof2
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile/profile_func.%j.log
# Profile MERGE_IDP, MERGE_SWISSPROT and MERGE_PFAM on ~1,018-genome subsets
# (every 23rd file; *.manifest in this folder). 4 CPUs / 24 GB = the 'merge' label.
# Each step is timed alone; outputs are compared by row count and a per-column checksum.
# Nothing in BFD_function is written.
set -euo pipefail
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile
BIN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_wt_csvlist/nextflow/bin
cd $P; module load duckdb
T=${SCRATCH:?}; O=$P/out_func; mkdir -p $O
echo "host $(hostname)"
run() { local name=$1; shift; /usr/bin/time -f "%e s wall, %M KB maxRSS" -o $O/time_$name.txt bash -c "$*"; echo "$name: $(cat $O/time_$name.txt)"; }
ck() { duckdb -csv -noheader -c "SELECT count(*), sum(hash(COLUMNS(*))) FROM '$1'"; }
list_sql() { cut -f1 $1 | awk '{gsub("\x27","\x27\x27"); printf "%s\x27%s\x27", (NR>1?",":""), $0}'; }

# ── IDP ──────────────────────────────────────────────────────────────
concat() { first=1; while IFS=$'\t' read -r f _m _s; do if [ $first = 1 ]; then zcat "$f"; first=0; else zcat "$f" | tail -n +2; fi; done < $1; }
export -f concat
run IDP0_concat "concat aiupred_idpcsvgz.manifest | gzip > $T/idp.csv.gz"
run IDP0_duckdb "duckdb -c \"COPY (SELECT * FROM read_csv_auto('$T/idp.csv.gz', sample_size=-1)) TO '$O/idp0.parquet' (FORMAT PARQUET)\""
run IDP1_direct "python3 $BIN/manifest_to_parquet.py --manifest aiupred_idpcsvgz.manifest --out $O/idp1.parquet --threads 4 --memory-limit 18GB --temp-dir $T"
echo "idp current: $(ck $O/idp0.parquet)"; echo "idp direct:  $(ck $O/idp1.parquet)"
rm -f $T/idp.csv.gz

# ── SWISSPROT ────────────────────────────────────────────────────────
run SP0_python "python3 $BIN/merge_swissprot.py -o $O/swissprot.csv \$(cut -f1 swissprot_blasttabgz.manifest)"
ls -l $O/swissprot.csv | awk '{print "swissprot.csv bytes", $5}'
run SP0_duckdb_autodetect "duckdb -c \"COPY (SELECT * FROM read_csv_auto('$O/swissprot.csv', sample_size=-1)) TO '$O/sp0.parquet' (FORMAT PARQUET)\""
SPCOLS="{'species_prefix':'VARCHAR','protein_id':'VARCHAR','swissprot_acc':'VARCHAR','swissprot_entry':'VARCHAR','swissprot_name':'VARCHAR','pident':'DOUBLE','positive':'BIGINT','nident':'BIGINT','aln_length':'BIGINT','q_len':'BIGINT','s_len':'BIGINT','qcovhsp':'DOUBLE','query_cov':'DOUBLE','hit_cov':'DOUBLE','q_start':'BIGINT','q_end':'BIGINT','s_start':'BIGINT','s_end':'BIGINT','gapopen':'BIGINT','mismatch':'BIGINT','evalue':'DOUBLE','bitscore':'DOUBLE','func_transfer_80_80':'BIGINT'}"
run SP1_duckdb_types "duckdb -c \"COPY (SELECT * FROM read_csv('$O/swissprot.csv', header=true, columns=$SPCOLS)) TO '$O/sp1.parquet' (FORMAT PARQUET)\""
rm -f $O/swissprot.csv
# SP2: no Python. DuckDB reads the 18-column DIAMOND/blastp tables and derives the
# same columns as merge_swissprot.py.
{ echo "SET threads=4; SET memory_limit='18GB'; SET temp_directory='$T';"
  echo "COPY (WITH r AS (SELECT * FROM read_csv([$(list_sql swissprot_blasttabgz.manifest)], delim='\t', header=false, quote='', escape='', columns={'qseqid':'VARCHAR','sseqid':'VARCHAR','pident':'DOUBLE','positive':'BIGINT','nident':'BIGINT','length':'BIGINT','mismatch':'BIGINT','gapopen':'BIGINT','qstart':'BIGINT','qend':'BIGINT','sstart':'BIGINT','send':'BIGINT','evalue':'DOUBLE','bitscore':'DOUBLE','qcovhsp':'DOUBLE','qlen':'BIGINT','slen':'BIGINT','stitle':'VARCHAR'})),"
  echo " d AS (SELECT *, split_part(stitle, ' OS=', 1) AS descr, CASE WHEN qlen <> 0 THEN length::DOUBLE / qlen ELSE 0 END AS qc, CASE WHEN slen <> 0 THEN length::DOUBLE / slen ELSE 0 END AS hc FROM r)"
  echo " SELECT split_part(qseqid, '_', 1) AS species_prefix, qseqid AS protein_id,"
  echo "  CASE WHEN contains(sseqid, '|') THEN split_part(sseqid, '|', 2) ELSE sseqid END AS swissprot_acc,"
  echo "  split_part(sseqid, '|', 3) AS swissprot_entry,"
  echo "  CASE WHEN position(' ' IN descr) > 0 THEN substr(descr, position(' ' IN descr) + 1) ELSE '' END AS swissprot_name,"
  echo "  pident, positive, nident, length AS aln_length, qlen AS q_len, slen AS s_len, qcovhsp,"
  echo "  round(qc, 4) AS query_cov, round(hc, 4) AS hit_cov, qstart AS q_start, qend AS q_end, sstart AS s_start, send AS s_end,"
  echo "  gapopen, mismatch, evalue, bitscore,"
  echo "  CASE WHEN pident >= 80 AND qc >= 0.8 AND hc >= 0.8 THEN 1 ELSE 0 END AS func_transfer_80_80"
  echo " FROM d) TO '$O/sp2.parquet' (FORMAT PARQUET);"; } > $T/sp2.sql
run SP2_duckdb_only "duckdb < $T/sp2.sql"
for x in sp0 sp1 sp2; do echo "swissprot $x: $(ck $O/$x.parquet)"; done
# rows that differ between the Python path and DuckDB-only (full-row compare)
echo "sp0 vs sp2 rows only in one: $(duckdb -csv -noheader -c "SELECT (SELECT count(*) FROM (SELECT * FROM '$O/sp0.parquet' EXCEPT ALL SELECT * FROM '$O/sp2.parquet')), (SELECT count(*) FROM (SELECT * FROM '$O/sp2.parquet' EXCEPT ALL SELECT * FROM '$O/sp0.parquet'))")"

# ── PFAM ─────────────────────────────────────────────────────────────
run PF0_python "python3 $BIN/pfamtbl_to_long.py --outfile $O/pfam.csv \$(cut -f1 pfam_hmmscan_domtbloutgz.manifest)"
ls -l $O/pfam.csv | awk '{print "pfam.csv bytes", $5}'
run PF0_duckdb_autodetect "duckdb -c \"COPY (SELECT * FROM read_csv_auto('$O/pfam.csv', sample_size=-1)) TO '$O/pf0.parquet' (FORMAT PARQUET)\""
PFCOLS="{'protein_id':'VARCHAR','pfam_id':'VARCHAR','pfam_acc':'VARCHAR','pfam_len':'BIGINT','full_seq_e_value':'DOUBLE','full_seq_score':'DOUBLE','full_seq_bias':'DOUBLE','domain_num':'BIGINT','domain_num_of':'BIGINT','domain_c_evalue':'DOUBLE','domain_i_evalue':'DOUBLE','domain_score':'DOUBLE','domain_bias':'DOUBLE','hmm_from':'BIGINT','hmm_to':'BIGINT','ali_from':'BIGINT','ali_to':'BIGINT','env_from':'BIGINT','env_to':'BIGINT'}"
run PF1_duckdb_types "duckdb -c \"COPY (SELECT * FROM read_csv('$O/pfam.csv', header=true, columns=$PFCOLS)) TO '$O/pf1.parquet' (FORMAT PARQUET)\""
rm -f $O/pfam.csv
# PF2: the same Python converter, 4 processes on 4 slices, then DuckDB reads the 4 CSVs with types
split -n l/4 -d <(cut -f1 pfam_hmmscan_domtbloutgz.manifest) $T/pfslice.
run PF2_python_x4 "for s in $T/pfslice.0*; do python3 $BIN/pfamtbl_to_long.py --outfile \$s.csv \$(cat \$s) & done; wait"
run PF2_duckdb_types "duckdb -c \"SET threads=4; COPY (SELECT * FROM read_csv('$T/pfslice.0*.csv', header=true, columns=$PFCOLS)) TO '$O/pf2.parquet' (FORMAT PARQUET)\""
for x in pf0 pf1 pf2; do echo "pfam $x: $(ck $O/$x.parquet)"; done
