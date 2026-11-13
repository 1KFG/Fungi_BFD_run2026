#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 3-00:00:00 -p epyc
#SBATCH --job-name=rebuild_telomeres
#SBATCH --output=logs/rebuild_telomeres.%j.log
#
# Rebuild the telomere tables for ALL genomes with the batched BATCH_TELOMERES step
# (Fungi_BFD PR #70 and #71, merged 2026-11-12). It replaces the per-genome FIND_TELOMERES
# that had the motif and orientation fix of commit 7426dbc, and adds a de novo finder.
# Modelled on rebuild_gene_introns.sh.
#
# The batched step writes to a NEW store, results/genome_stats/telomeres_v2/, so no stage
# is needed: the old store (telomeres/) is not read by it. Outputs of the old per-genome
# step cannot be reused. The code must come from a Fungi_BFD checkout that has BATCH_TELOMERES
# and the 'set -euo pipefail' fix (PR #73): without it a per-genome error (for example a corrupt
# genome) exits 0 and storeDir keeps a header-only 'no telomeres' file. Set WORKTREE to such a
# checkout. The script refuses to run on one without the fix.
#
# Usage (run from BFD_function/):
#   bash rebuild_telomeres.sh status               # counts only, changes nothing (slow: scans the store)
#   bash rebuild_telomeres.sh stage                # not needed with telomeres_v2 (prints why)
#   N_TEST=20 FORCE=1 sbatch --time 04:00:00 rebuild_telomeres.sh run   # test: first 20 genomes (one batch)
#   bash rebuild_telomeres.sh timing               # BATCH_TELOMERES task times from the latest trace
#   FORCE=1 sbatch rebuild_telomeres.sh run        # all genomes
#   FORCE=1 sbatch --time 1-00:00:00 rebuild_telomeres.sh merge  # parquet tables + pipeline DuckDB
#
# Queues: the head job and every Nextflow task go to one partition, QUEUE (default epyc).
# The pipeline profile sends telomere, setup, merge and duckdb tasks to 'short' and some
# to 'batch'; queue_override_<QUEUE>.config (written by this script) replaces that.
# To change it: QUEUE=stajichlab FORCE=1 sbatch -p stajichlab rebuild_telomeres.sh run
# The #SBATCH header sets the head job's partition (-p epyc); override it on the command line.
#   bash rebuild_telomeres.sh verify               # check files and compare with the old tables
#   bash rebuild_telomeres.sh clean --go           # delete the backup after verify passes
#
# BATCH_TELOMERES uses storeDir (telomeres_v2). A genome with an output there is skipped on a
# rerun. A test run (N_TEST) writes real outputs, and the full run keeps them.
# `clean` deletes the backup of the OLD store (telomeres.pre_both_ends_fix_20261112), which the
# new step does not use.
#
# This script does not touch the dashboard database. After merge, rebuild it with
#   sbatch /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/workflows/build_bfd_db.sh

set -euo pipefail

BFD=/bigdata/stajichlab/shared/projects/BFD
RUNS=$BFD/Fungi_BFD_runs
WORKTREE=${WORKTREE:-$BFD/Fungi_BFD}   # needs BATCH_TELOMERES and the PR #73 fix
RUNDIR=$RUNS/BFD_function
STORE=$RUNS/results/genome_stats/telomeres_v2
OLD_STORE=$RUNS/results/genome_stats/telomeres
BACKUP=$RUNS/results/genome_stats/telomeres.pre_both_ends_fix_20261112
TABLES=$RUNS/BFD_function/tables
FIX_COMMIT=7426dbc
PATTERN='*.telomeres.tsv.gz'
MODE=${1:-status}
GO=${2:-}

# Tables as they were before this rerun (Nov 6 build, read from BFD_dashboard db on 2026-11-12).
OLD_GENOMES=23582
OLD_WITH_TELOMERE=8369
OLD_TRACTS=73242

cd "$RUNDIR"
mkdir -p logs/nextflow

die() { echo "ERROR: $*" >&2; exit 1; }

check_no_live_run() {
    # Do not change files under a live Nextflow run.
    if pgrep -u "$USER" -af "nextflow.*main.nf" | grep -v "pgrep" >/dev/null; then
        pgrep -u "$USER" -af "nextflow.*main.nf" | cut -c1-200 >&2
        [ "${FORCE:-0}" = 1 ] || die "a Nextflow run is live. Wait for it, or set FORCE=1 if it does not use $STORE."
    fi
}

check_fix_in_worktree() {
    git -C "$WORKTREE" merge-base --is-ancestor "$FIX_COMMIT" HEAD \
        || die "$WORKTREE HEAD does not contain $FIX_COMMIT. Update the checkout first."
    grep -q "drop_contained_hits" "$WORKTREE/nextflow/bin/find_telomeres.py" \
        || die "fixed code not found in find_telomeres.py"
    local bt="$WORKTREE/nextflow/modules/BFD/BATCH_TELOMERES/main.nf"
    [ -e "$bt" ] || die "$WORKTREE has no BATCH_TELOMERES (Fungi_BFD PR #70). Update the checkout or set WORKTREE."
    grep -q "pipefail" "$bt" \
        || die "BATCH_TELOMERES has no 'set -euo pipefail' (Fungi_BFD PR #73). A failed genome would be stored as 'no telomeres'. Set WORKTREE to a checkout with the fix."
    grep -Eq "telomere_both_ends[[:space:]]*=[[:space:]]*true" "$WORKTREE/nextflow/conf/profile_BFD.config" \
        || die "telomere_both_ends is not true in conf/profile_BFD.config"
}

export_nf_env() {
    module load nextflow
    export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
    export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
    export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
}

# count_files DIR [NEWER_THAN_ISO]: count telomere files, one find per bucket, 8 in parallel.
count_files() {
    local dir=$1 ts=${2:-}
    [ -d "$dir" ] || { echo 0; return; }
    export PATTERN TS="$ts"
    ls "$dir" | xargs -n 1 -P 8 bash -c '
        d="$0/$1"; [ -d "$d" ] || { echo 0; exit 0; }
        if [ -n "$TS" ]; then find "$d" -maxdepth 1 -name "$PATTERN" -newermt "$TS" | wc -l
        else find "$d" -maxdepth 1 -name "$PATTERN" | wc -l; fi
    ' "$dir" | awk '{n+=$1} END{print n+0}'
}

fix_ts() { git -C "$WORKTREE" show -s --format=%cI "$FIX_COMMIT"; }

QUEUE=${QUEUE:-epyc}
# Write a Nextflow config that puts every process this script can start on one partition.
queue_config() {
    local f="$RUNDIR/queue_override_${QUEUE}.config"
    cat > "$f" <<EOF
// Written by rebuild_telomeres.sh. All tasks of this run go to the '$QUEUE' partition.
process {
    withLabel: 'telomeres' { queue = '$QUEUE' }
    withLabel: 'setup'     { queue = '$QUEUE' }
    withLabel: 'merge'     { queue = '$QUEUE' }
    withLabel: 'duckdb'    { queue = '$QUEUE' }
}
EOF
    echo "$f"
}

case "$MODE" in
status)
    echo "store (telomeres_v2):  $STORE"
    echo "old store (v1):        $OLD_STORE   files: $(count_files "$OLD_STORE")"
    echo "backup of old store:   $BACKUP $( [ -d "$BACKUP" ] && echo "exists, files: $(count_files "$BACKUP")" || echo absent)"
    ts=$(fix_ts)
    echo "telomere files in v2 store:             $(count_files "$STORE")"
    echo "  newer than the 7426dbc commit ($ts): $(count_files "$STORE" "$ts")"
    ls -l --time-style=long-iso "$TABLES"/telomere_summary.parquet "$TABLES"/telomere_tracts.parquet "$RUNDIR/db/BFD.duckdb"
    echo "WORKTREE=$WORKTREE"
    echo "  BATCH_TELOMERES: $( [ -e "$WORKTREE/nextflow/modules/BFD/BATCH_TELOMERES/main.nf" ] && echo present || echo MISSING )"
    echo "  pipefail fix (PR #73): $( grep -qs pipefail "$WORKTREE/nextflow/modules/BFD/BATCH_TELOMERES/main.nf" && echo present || echo MISSING )"
    echo "live Nextflow processes: $(pgrep -u "$USER" -f 'nextflow.*main.nf' | wc -l)"
    ;;

stage)
    echo "stage is not needed. BATCH_TELOMERES writes to $STORE (telomeres_v2), which the old"
    echo "per-genome step never used, so no old file can be reused. Nothing was changed."
    echo "The old store was already emptied on 2026-11-12 (backup: $BACKUP)."
    ;;

run)
    check_fix_in_worktree
    check_no_live_run
    export_nf_env
    # n_test goes in through a params file. A command-line --n_test is a string and
    # fails the pipeline's integer parameter validation.
    params="$RUNDIR/params_telomeres_rerun.yaml"
    if [ -n "${N_TEST:-}" ]; then
        [[ "$N_TEST" =~ ^[0-9]+$ ]] || die "N_TEST must be an integer"
        params="$RUNDIR/params_telomeres_test.yaml"
        sed "s/^n_test:.*/n_test: $N_TEST/" "$RUNDIR/params_telomeres_rerun.yaml" > "$params"
        grep -q "^n_test: $N_TEST\$" "$params" || die "could not set n_test in $params"
        echo "TEST RUN: first $N_TEST genomes only ($params)"
    fi
    # No -resume on purpose: a cached task from an older code version could restore old
    # output into the store. storeDir still skips genomes that finished in this run, so a
    # failed run can be started again as is. Each task runs telomere_batch_size genomes.
    NXF_OPTS="-Xms512m -Xmx12g" \
    nextflow run "$WORKTREE/nextflow/main.nf" \
        -c "$WORKTREE/nextflow/nextflow.config" \
        -c "$(queue_config)" \
        -profile BFD \
        -params-file "$params"
    ;;

timing)
    # Summarise BATCH_TELOMERES task times (one task = telomere_batch_size genomes) from the
    # newest trace file. Check that batches finish well inside the 4 h task limit.
    f=$(ls -t logs/nextflow/BFD_trace.20*.txt 2>/dev/null | head -1)
    [ -n "$f" ] || die "no trace file in logs/nextflow"
    echo "trace: $f"
    /usr/bin/python3.12 - "$f" <<'PYEOF'
import re, sys, statistics as st
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
h = rows[0]; ix = {k: i for i, k in enumerate(h)}
def secs(s):
    t = {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}
    return sum(float(n) * t[u] for n, u in re.findall(r"([\d.]+)(ms|s|m|h|d)\b", s))
sel = [r for r in rows[1:] if "BATCH_TELOMERES" in r[ix["name"]]]
rt = sorted(secs(r[ix["realtime"]]) for r in sel)
ok = [r for r in sel if r[ix["status"]] == "COMPLETED"]
if not rt:
    sys.exit("no BATCH_TELOMERES tasks in this trace")
n = len(rt)
print(f"BATCH_TELOMERES tasks: {n} ({len(ok)} COMPLETED, {n-len(ok)} not)")
print(f"realtime seconds: median {st.median(rt):.1f}, mean {st.mean(rt):.1f}, p90 {rt[int(0.9*(n-1))]:.1f}, max {rt[-1]:.1f}")
print(f"tasks over 3 h: {sum(1 for x in rt if x > 10800)} (task limit is 4 h)")
print(f"23,580 genomes at 100 per task = about 236 tasks; at the mean task time: {st.mean(rt)*236/3600:.0f} CPU hours (estimate from {n} tasks)")
PYEOF
    ;;

merge)
    check_fix_in_worktree
    check_no_live_run
    # The merge reads telomeres_v2 only. Refuse to merge a store that is clearly incomplete.
    total=$(count_files "$STORE")
    echo "v2 store files: $total (the old table had $OLD_GENOMES genomes)"
    [ "$total" -gt 0 ] || die "v2 store is empty"
    [ "$total" -ge $((OLD_GENOMES * 98 / 100)) ] || [ "${FORCE_MERGE:-0}" = 1 ] \
        || die "v2 store has $total files, under 98% of $OLD_GENOMES. Finish run first (FORCE_MERGE=1 to override)."
    export_nf_env
    # The merge reruns when a result file is newer than its table. Rewritten telomere
    # files are newer, so MERGE_TELOMERES rebuilds both parquet files. BUILD_DUCKDB runs
    # after a merge.
    NXF_OPTS="-Xms512m -Xmx12g" \
    nextflow run "$WORKTREE/nextflow/main.nf" \
        -c "$WORKTREE/nextflow/nextflow.config" \
        -c "$(queue_config)" \
        -profile BFD \
        -params-file "$RUNDIR/params_merge_after_telomeres.yaml" \
        -resume
    ;;

verify)
    ts=$(fix_ts)
    new=$(count_files "$STORE")
    newer=$(count_files "$STORE" "$ts")
    echo "v2 store telomere files: $new  (newer than 7426dbc: $newer)"
    [ "$newer" = "$new" ] || echo "WARNING: some store files predate the fix"
    [ "$new" -ge $((OLD_GENOMES * 98 / 100)) ] || echo "WARNING: v2 store is under 98% of the old table's $OLD_GENOMES genomes. Check logs for failed tasks."
    echo "tables:"
    ls -l --time-style=long-iso "$TABLES"/telomere_summary.parquet "$TABLES"/telomere_tracts.parquet
    module load duckdb 2>/dev/null || true
    command -v duckdb >/dev/null || die "duckdb CLI not found"
    duckdb -c "
      SELECT count(*) AS genomes,
             count(*) FILTER (telomere_scaffolds > 0) AS with_telomere,
             count(*) FILTER (telomere_scaffolds_both_ends > 0) AS with_both_ends,
             $OLD_GENOMES AS old_genomes, $OLD_WITH_TELOMERE AS old_with_telomere
      FROM read_parquet('$TABLES/telomere_summary.parquet');
      SELECT count(*) AS tracts, $OLD_TRACTS AS old_tracts, count(DISTINCT tract_id) AS distinct_tract_ids
      FROM read_parquet('$TABLES/telomere_tracts.parquet');
      SELECT end_type, strand, count(*) n FROM read_parquet('$TABLES/telomere_tracts.parquet') GROUP BY ALL ORDER BY ALL;
      SELECT monomer, count(*) n FROM read_parquet('$TABLES/telomere_tracts.parquet') GROUP BY 1 ORDER BY n DESC LIMIT 12;
      -- Control genomes named in commit $FIX_COMMIT. These should now show telomere scaffolds.
      SELECT SPECIES, STRAIN, telomere_scaffolds, telomere_tracts, telomere_top_monomer
      FROM read_parquet('$TABLES/telomere_summary.parquet')
      WHERE SPECIES ILIKE 'Aspergillus flavus%' OR SPECIES ILIKE 'Aspergillus oryzae%'
         OR SPECIES ILIKE 'Cryptococcus neoformans%' OR SPECIES ILIKE 'Tremella fuciformis%'
         OR SPECIES ILIKE 'Wallemia mellicola%' OR SPECIES ILIKE 'Encephalitozoon cuniculi%'
      ORDER BY SPECIES, STRAIN LIMIT 30;
    "
    echo "Duplicate-key checks (ASMID, tract_id) run in the DuckDB build, not here."
    ;;

clean)
    [ -d "$BACKUP" ] || die "no backup dir"
    [ "$GO" = "--go" ] || die "add --go to delete $BACKUP"
    rm -rf "$BACKUP"
    echo "deleted $BACKUP"
    ;;

*)  die "unknown mode: $MODE" ;;
esac
