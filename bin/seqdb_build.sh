#!/bin/bash -l
#SBATCH -p epyc -c 4 --mem 8G --time 2:00:00 -J seqdb
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb/logs/seqdb.%A_%a.log
# seqdb_build.sh -- bgzip + samtools faidx copies of the predicted protein and CDS
# FASTA files, for random access by sequence ID (the BFD gene tables do not
# store sequences). The originals in genome_annotation/ are not touched.
#
# Input: seqdb/genomes.tsv (bin/seqdb_list.py), split into chunks by
#   SLURM_ARRAY_TASK_ID: lines (task*CHUNK+1) .. (task+1)*CHUNK.
# Output: seqdb/{proteins,cds}/<sha1(LOCUSTAG)[:2]>/<LOCUSTAG>.{proteins,cds}.fa.gz
#   + .fai + .gzi. A copy is kept only if its decompressed content has the same
#   md5 as the original; otherwise it is removed and the genome is reported FAIL.
# Re-running skips copies whose .fa.gz, .fai and .gzi all exist.
# usage: sbatch --array=0-N%M seqdb_build.sh [CHUNK]
set -uo pipefail
S=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb
CHUNK=${1:-500}
module load samtools/1.22.1 htslib/1.22.1 pigz 2>/dev/null || module load samtools htslib
T=${SLURM_ARRAY_TASK_ID:-0}
first=$((T * CHUNK + 1)); last=$(((T + 1) * CHUNK))
copy() {  # src dest
  local src=$1 dst=$2
  [ -s "$dst" ] && [ -s "$dst.fai" ] && [ -s "$dst.gzi" ] && { echo "SKIP"; return 0; }
  mkdir -p "$(dirname "$dst")"; rm -f "$dst" "$dst.fai" "$dst.gzi"
  pigz -dc "$src" | bgzip -@ 2 -c > "$dst" || { rm -f "$dst"; echo "FAIL bgzip"; return 1; }
  local a b
  a=$(pigz -dc "$src" | md5sum | cut -c1-32); b=$(bgzip -dc "$dst" | md5sum | cut -c1-32)
  [ "$a" = "$b" ] || { rm -f "$dst"; echo "FAIL md5 $a $b"; return 1; }
  samtools faidx "$dst" 2> "$dst.err" || { echo "FAIL faidx $(head -c 200 "$dst.err")"; rm -f "$dst" "$dst.fai" "$dst.gzi" "$dst.err"; return 1; }
  rm -f "$dst.err"; echo "OK $(wc -l < "$dst.fai")"
}
t0=$(date +%s); n=0
sed -n "${first},${last}p" $S/genomes.tsv | while IFS=$'\t' read -r tag out pep cds; do
  b=$(printf '%s' "$tag" | sha1sum | cut -c1-2)
  r1=$(copy "$pep" "$S/proteins/$b/$tag.proteins.fa.gz")
  r2="none"; [ -n "$cds" ] && r2=$(copy "$cds" "$S/cds/$b/$tag.cds.fa.gz")
  printf '%s\tproteins:%s\tcds:%s\n' "$tag" "$r1" "$r2"
done
echo "task $T lines $first-$last done in $(( $(date +%s) - t0 )) s"
