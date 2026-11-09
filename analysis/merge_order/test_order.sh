#!/bin/bash -l
#SBATCH -p stajichlab -c 4 --mem 24G --time 2:00:00 -J merge_order
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_order/test_order.%j.log
# Row order of merged Parquet: input files in full-path order (sha1 bucket =
# hash order, old) vs file-name order (LOCUSTAG, new). 1,018 genomes of gene_exons.
set -euo pipefail
module load duckdb
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/merge_profile
NEW=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_wt_force/nextflow/bin
OLD=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin
T=${SCRATCH:?}
for v in OLD NEW; do
  B=${!v}
  python3 $B/manifest_to_parquet.py --manifest $P/exons_1000.manifest --suffix .gene_exons.csv.gz --out $T/$v.parquet --threads 4 --memory-limit 18GB --temp-dir $T | cut -c1-120
  echo "$v overlap: $(duckdb -csv -noheader -c "SELECT count(*), sum((mn < lmx)::INT) FROM (SELECT stats_min_value mn, lag(stats_max_value) OVER (ORDER BY row_group_id) lmx FROM parquet_metadata('$T/$v.parquet') WHERE path_in_schema='exon_id')")"
  id=$(duckdb -csv -noheader -c "SELECT exon_id FROM '$T/$v.parquet' USING SAMPLE 1 ROWS (reservoir, 7)")
  for i in 1 2; do s=$(date +%s.%N); duckdb -csv -noheader -c "SELECT count(*) FROM '$T/$v.parquet' WHERE exon_id='$id'" >/dev/null; echo "$v lookup $id: $(echo "$(date +%s.%N)-$s"|bc|cut -c1-5) s"; done
  pre=${id%%_*}
  s=$(date +%s.%N); n=$(duckdb -csv -noheader -c "SELECT count(*) FROM '$T/$v.parquet' WHERE exon_id LIKE '${pre}_%'"); echo "$v genome $pre: $n rows, $(echo "$(date +%s.%N)-$s"|bc|cut -c1-5) s"
done
duckdb -csv -c "SELECT (SELECT count(*) FROM '$T/OLD.parquet') old_rows, (SELECT count(*) FROM '$T/NEW.parquet') new_rows,
 (SELECT count(*) FROM (SELECT * FROM '$T/OLD.parquet' EXCEPT ALL SELECT * FROM '$T/NEW.parquet')) only_old"
