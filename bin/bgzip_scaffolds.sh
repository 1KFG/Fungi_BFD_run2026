#!/bin/bash -l
#SBATCH -p epyc -c 8 --mem 8G --time 4:00:00 -J bgzip_scaf
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/logs/bgzip_scaffolds/bgzip_scaffolds.%A_%a.log
# bgzip_scaffolds.sh -- backstop: bgzip + samtools faidx every plain
# genome_annotation/*/predict_results/*.scaffolds.fa. New predictions are
# compressed by FUNANNOTATE_PREDICT itself (Fungi_BFD PR #58); this handles the
# older folders.
#
# Steps:
#   1. sbatch bgzip_scaffolds.sh list
#        -> logs/bgzip_scaffolds/files.txt (one path per line)
#   2. DRYRUN=1 sbatch --array=0-N%M bgzip_scaffolds.sh run [CHUNK]
#        -> reports what would be done; changes nothing
#   3. sbatch --array=0-N%M bgzip_scaffolds.sh run [CHUNK]
#   N = ceil(lines / CHUNK) - 1, CHUNK default 500.
#
# Per file (GNU parallel, 4 at a time, bgzip -@ 2 each):
#   - skip when any file in that predict_results/ changed in the last 60 min
#     (a predict task may be writing there) or an *.inflate.* temp file is
#     present (funannotate update/annotate is running on it);
#   - bgzip to <f>.gz.tmp; keep it only when its decompressed size equals the
#     plain file size; give it the plain file's mtime (BFD stale checks compare
#     input mtimes, so a new mtime would rerun every downstream step);
#   - rename to <f>.gz, samtools faidx (.fai + .gzi), then remove the plain file.
#   A failure leaves the plain file in place and removes the partial .gz.
# Output: one TSV line per file in the task log:
#   status <TAB> path <TAB> plain_bytes <TAB> gz_bytes(+fai+gzi)
# input/dna links are switched to the .gz name by SETUP_SYMLINKS at the next BFD run.
set -euo pipefail
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
L=$R/logs/bgzip_scaffolds
mkdir -p $L
MODE=${1:?usage: bgzip_scaffolds.sh list|run [CHUNK]}

if [ "$MODE" = list ]; then
    find $R/genome_annotation -mindepth 3 -maxdepth 3 -path '*/predict_results/*.scaffolds.fa' -type f \
        | LC_ALL=C sort > $L/files.txt.tmp
    mv $L/files.txt.tmp $L/files.txt
    n=$(wc -l < $L/files.txt)
    echo "files: $n; plain bytes: $(tr '\n' '\0' < $L/files.txt | du -cb --files0-from=- | tail -1 | cut -f1)"
    echo "array for CHUNK=500: --array=0-$(( (n + 499) / 500 - 1 ))"
    exit 0
fi

[ "$MODE" = run ] || { echo "unknown mode $MODE" >&2; exit 1; }
CHUNK=${2:-500}
T=${SLURM_ARRAY_TASK_ID:-0}
DRYRUN=${DRYRUN:-0}
module load samtools/1.22.1 htslib/1.22.1
module load parallel
export DRYRUN

one() {
    local f=$1 d n m
    d=$(dirname "$f")
    [ -f "$f" ] || { printf 'GONE\t%s\t0\t0\n' "$f"; return 0; }
    n=$(stat -c %s "$f")
    if [ -n "$(find "$d" -mindepth 1 -maxdepth 1 -newermt '-60 minutes' -print -quit)" ]; then
        printf 'SKIP_RECENT\t%s\t%s\t0\n' "$f" "$n"; return 0
    fi
    if compgen -G "$d/*.inflate.*" >/dev/null; then
        printf 'SKIP_INFLATE\t%s\t%s\t0\n' "$f" "$n"; return 0
    fi
    if [ "$DRYRUN" = 1 ]; then printf 'WOULD\t%s\t%s\t0\n' "$f" "$n"; return 0; fi
    rm -f "$f.gz.tmp"
    if ! bgzip -@ 2 -c "$f" > "$f.gz.tmp"; then
        rm -f "$f.gz.tmp"; printf 'FAIL_BGZIP\t%s\t%s\t0\n' "$f" "$n"; return 0
    fi
    m=$(bgzip -dc "$f.gz.tmp" | wc -c)
    if [ "$m" != "$n" ]; then
        rm -f "$f.gz.tmp"; printf 'FAIL_SIZE\t%s\t%s\t%s\n' "$f" "$n" "$m"; return 0
    fi
    touch -r "$f" "$f.gz.tmp"
    rm -f "$f.gz.fai" "$f.gz.gzi"
    mv "$f.gz.tmp" "$f.gz"
    if ! samtools faidx "$f.gz" 2>/dev/null; then
        rm -f "$f.gz" "$f.gz.fai" "$f.gz.gzi"; printf 'FAIL_FAIDX\t%s\t%s\t0\n' "$f" "$n"; return 0
    fi
    rm -f "$f"
    printf 'OK\t%s\t%s\t%s\n' "$f" "$n" \
        "$(( $(stat -c %s "$f.gz") + $(stat -c %s "$f.gz.fai") + $(stat -c %s "$f.gz.gzi") ))"
}
export -f one

first=$((T * CHUNK + 1)); last=$(((T + 1) * CHUNK))
echo "task $T: lines $first-$last of $L/files.txt; DRYRUN=$DRYRUN; host $(hostname)"
sed -n "${first},${last}p" $L/files.txt \
    | parallel -j $(( ${SLURM_CPUS_PER_TASK:-8} / 2 )) --line-buffer one {} \
    | tee $L/result.task$T.tsv
echo "summary task $T:"
awk -F'\t' '{c[$1]++; p[$1]+=$3; g[$1]+=$4} END {for (s in c) printf "  %s %d files, plain %.2f GB, gz %.2f GB\n", s, c[s], p[s]/1e9, g[s]/1e9}' $L/result.task$T.tsv
