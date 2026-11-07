#!/bin/bash -l
#SBATCH -p epyc -c 4 --mem 24G --time 2:00:00 -J m2p_test
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile/m2p_test.%j.log
# Compare bin/manifest_to_parquet.py with the zcat-concatenate path on ~1000-genome subsets.
set -euo pipefail
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile
BIN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_wt_csvlist/nextflow/bin
cd $P; module load duckdb
T=${SCRATCH:?}
ck() { duckdb -csv -noheader -c "SELECT count(*), sum(hash(COLUMNS(*))) FROM '$1'" 2>/dev/null || duckdb -csv -noheader -c "SELECT count(*) FROM '$1'"; }
concat() { first=1; while IFS=$'\t' read -r f _m _s; do if [ $first = 1 ]; then zcat "$f"; first=0; else zcat "$f" | tail -n +2; fi; done < $1; }
/usr/bin/time -f "%e s, %M KB" python3 $BIN/manifest_to_parquet.py --manifest exons_1000.manifest --suffix .gene_exons.csv.gz --out $T/exons.parquet --threads 4 --memory-limit 18GB --temp-dir $T
echo "exons new: $(ck $T/exons.parquet)   reference A: $(ck $P/A.parquet)"
for x in aa ig; do
  /usr/bin/time -f "%e s, %M KB" python3 $BIN/manifest_to_parquet.py --manifest ${x}_1000.manifest --out $T/$x.parquet --threads 4 --memory-limit 18GB --temp-dir $T
  /usr/bin/time -f "concat path %e s, %M KB" bash -c "$(declare -f concat); concat ${x}_1000.manifest > $T/$x.csv; duckdb -c \"COPY (SELECT * FROM read_csv_auto('$T/$x.csv', sample_size=-1)) TO '$T/${x}_ref.parquet' (FORMAT PARQUET)\""
  echo "$x new: $(ck $T/$x.parquet)   concat ref: $(ck $T/${x}_ref.parquet)"
  duckdb -csv -noheader -c "SELECT string_agg(column_name||':'||column_type, ',') FROM (DESCRIBE SELECT * FROM '$T/$x.parquet')" > $T/a; duckdb -csv -noheader -c "SELECT string_agg(column_name||':'||column_type, ',') FROM (DESCRIBE SELECT * FROM '$T/${x}_ref.parquet')" > $T/b; cmp -s $T/a $T/b && echo "$x types identical" || { echo "$x types differ:"; cat $T/a $T/b; }
done
