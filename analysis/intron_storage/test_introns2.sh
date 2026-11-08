#!/bin/bash -l
#SBATCH -p stajichlab -c 16 --mem 64G --time 6:00:00 -J introns2
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/intron_storage/introns2.%j.log
# (1) build a VIEW-only BFD.duckdb into $SCRATCH from tables/ and time it;
# (2) 40 genomes: intron FASTA (bgzip+faidx) size vs the Parquet seq column,
#     and whether genome coordinates reproduce the stored intron seq;
# (3) total size of genome_annotation/*/predict_results/*.scaffolds.fa.
set -euo pipefail
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
NF=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_wt_csvlist/nextflow
O=$R/analysis/intron_storage
module load duckdb samtools htslib
T=${SCRATCH:?}; cd $T
echo "== (2) intron sample"
awk -F'\t' 'NR%585==1{print $1"\t"$2}' $R/seqdb/genomes.tsv | head -40 > sample.tsv
wc -l sample.tsv
cut -f1 sample.tsv | sed "s/.*/'&'/" | paste -sd, > lt.txt
duckdb -c "SET threads=16; COPY (SELECT * FROM '$R/tables/gene_introns.parquet' WHERE split_part(transcript_id,'_',1) IN ($(cat lt.txt))) TO 'introns_sample.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * EXCLUDE (seq) FROM 'introns_sample.parquet') TO 'introns_noseq.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);"
echo "parquet with seq $(stat -c %s introns_sample.parquet); without seq $(stat -c %s introns_noseq.parquet)"
duckdb -csv -noheader -c "SELECT count(*), sum(length(seq)) FROM 'introns_sample.parquet'" | sed 's/^/rows,bases: /'
fa_tot=0; fai_tot=0; gzi_tot=0; scaf_tot=0; scafgz_tot=0; ok=0; bad=0
while IFS=$'\t' read -u 3 lt out; do
  duckdb -list -noheader -separator $'\t' -c "SELECT intron_id, seq FROM 'introns_sample.parquet' WHERE split_part(transcript_id,'_',1)='$lt' ORDER BY intron_id" | awk -F'\t' '{print ">"$1; print $2}' > $lt.introns.fa
  bgzip -@4 -f $lt.introns.fa; samtools faidx $lt.introns.fa.gz
  fa_tot=$((fa_tot + $(stat -c %s $lt.introns.fa.gz))); fai_tot=$((fai_tot + $(stat -c %s $lt.introns.fa.gz.fai))); gzi_tot=$((gzi_tot + $(stat -c %s $lt.introns.fa.gz.gzi)))
  sc=$R/genome_annotation/$out/predict_results/$out.scaffolds.fa
  scaf_tot=$((scaf_tot + $(stat -c %s $sc)))
  bgzip -@8 -c $sc > g.fa.gz; samtools faidx g.fa.gz
  scafgz_tot=$((scafgz_tot + $(stat -c %s g.fa.gz) + $(stat -c %s g.fa.gz.fai) + $(stat -c %s g.fa.gz.gzi)))
  # coordinate check on 200 introns: faidx region (reverse complement for strand -1) vs stored seq
  duckdb -csv -noheader -c "SELECT chrom||':'||start||'-'||\"end\", strand, upper(seq) FROM 'introns_sample.parquet' WHERE split_part(transcript_id,'_',1)='$lt' USING SAMPLE 200 ROWS" > chk.csv
  while IFS=, read -u 4 reg st sq; do
    if [ "$st" = "-1" ]; then got=$(samtools faidx -i g.fa.gz "$reg" | grep -v '>' | tr -d '\n' | tr a-z A-Z); else got=$(samtools faidx g.fa.gz "$reg" | grep -v '>' | tr -d '\n' | tr a-z A-Z); fi
    [ "$got" = "$sq" ] && ok=$((ok+1)) || bad=$((bad+1))
  done 4< chk.csv
  rm -f g.fa.gz*
done 3< sample.tsv
echo "intron FASTA bgzip $fa_tot + fai $fai_tot + gzi $gzi_tot = $((fa_tot+fai_tot+gzi_tot)) bytes"
echo "genome scaffolds.fa plain $scaf_tot; bgzip+fai+gzi $scafgz_tot"
echo "coordinate check: $ok match, $bad differ"

echo "genomes processed: $(ls *.introns.fa.gz | wc -l)"
