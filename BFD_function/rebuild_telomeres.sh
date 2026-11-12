#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 1-00:00:00 -p short
#SBATCH --job-name=rebuild_telomeres
#SBATCH --output=logs/rebuild_telomeres.%j.log
#
# Rebuild the telomere tables for ALL genomes after the motif and orientation fix in
# Fungi_BFD commit 7426dbc (2026-11-11): telomere_both_ends = true, TTAGG and CCC[C]+TAA
# added, one G-strand label per motif, overlapping hits dropped.
# Modelled on rebuild_gene_introns.sh.
#
# Usage (run from BFD_function/):
#   bash rebuild_telomeres.sh status               # counts only, changes nothing (slow: scans the store)
#   bash rebuild_telomeres.sh stage                # dry run: shows what would move
#   bash rebuild_telomeres.sh stage --go           # move the old store into the backup dir
#   N_TEST=20 sbatch rebuild_telomeres.sh run      # timing test: first 20 genomes only
#   bash rebuild_telomeres.sh timing               # FIND_TELOMERES run times from the latest trace
#   sbatch -p batch --time 3-00:00:00 rebuild_telomeres.sh run   # all genomes
#   sbatch -p batch --time 1-00:00:00 rebuild_telomeres.sh merge # parquet tables + pipeline DuckDB
#   bash rebuild_telomeres.sh verify               # check files and compare with the old tables
#   bash rebuild_telomeres.sh clean --go           # delete the backup after verify passes
#
# Why stage: FIND_TELOMERES uses storeDir. Nextflow skips a genome when its output exists,
# and a code change does not trigger a rerun, so the old files must be moved away.
# Stage MOVES files (same filesystem, reversible). Only `clean` deletes.
#
# The test run writes to the store. Its outputs come from the fixed code, so the full run
# keeps them (storeDir skips them) and does not need another stage.
#
# This script does not touch the dashboard database. After merge, rebuild it with
#   sbatch /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/workflows/build_bfd_db.sh

set -euo pipefail

BFD=/bigdata/stajichlab/shared/projects/BFD
RUNS=$BFD/Fungi_BFD_runs
WORKTREE=${WORKTREE:-$BFD/Fungi_BFD}
RUNDIR=$RUNS/BFD_function
STORE=$RUNS/results/genome_stats/telomeres
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

case "$MODE" in
status)
    echo "store:   $STORE"
    echo "backup:  $BACKUP $( [ -d "$BACKUP" ] && echo exists || echo absent)"
    [ -e "$BACKUP/.stage_complete" ] && echo "stage:   complete" || echo "stage:   not done"
    ts=$(fix_ts)
    echo "telomere files in store:                $(count_files "$STORE")"
    echo "  newer than the fix commit ($ts): $(count_files "$STORE" "$ts")"
    [ -d "$BACKUP" ] && echo "telomere files in backup:               $(count_files "$BACKUP")"
    ls -l --time-style=long-iso "$TABLES"/telomere_summary.parquet "$TABLES"/telomere_tracts.parquet "$RUNDIR/db/BFD.duckdb"
    echo "live Nextflow processes: $(pgrep -u "$USER" -f 'nextflow.*main.nf' | wc -l)"
    ;;

stage)
    check_fix_in_worktree
    check_no_live_run
    # mv is safe to repeat, so an interrupted stage can be resumed. Stop only after a
    # completed stage, because new files in the store then belong to a later run.
    [ -e "$BACKUP/.stage_complete" ] && die "stage already completed ($BACKUP/.stage_complete). Do not stage again."
    n=$(count_files "$STORE")
    echo "telomere files to stage: $n"
    if [ "$GO" != "--go" ]; then
        echo "Dry run. Re-run with --go to move all telomere files into $BACKUP"
        exit 0
    fi
    mkdir -p "$BACKUP"
    # One mv per bucket directory, 8 in parallel. Layout is preserved.
    export STORE BACKUP
    ls "$STORE" | xargs -n 1 -P 8 bash -c '
        b=$1; src="$STORE/$b"; dst="$BACKUP/$b"
        [ -d "$src" ] || exit 0
        mkdir -p "$dst"
        find "$src" -maxdepth 1 -name "*.telomeres.tsv.gz" -exec mv -t "$dst" {} +
    ' _
    touch "$BACKUP/.stage_complete"
    echo "staged. telomere files left in store: $(count_files "$STORE")"
    echo "telomere files in backup:             $(count_files "$BACKUP")"
    ;;

run)
    check_fix_in_worktree
    check_no_live_run
    [ -e "$BACKUP/.stage_complete" ] || die "run stage --go first"
    export_nf_env
    extra=()
    [ -n "${N_TEST:-}" ] && { extra+=(--n_test "$N_TEST"); echo "TEST RUN: first $N_TEST genomes only"; }
    # No -resume on purpose: FIND_TELOMERES has no code-hash guard, so a cached task from
    # before the fix could restore old output into the store. storeDir still skips
    # genomes that finished in this run, so a failed run can be started again as is.
    NXF_OPTS="-Xms512m -Xmx12g" \
    nextflow run "$WORKTREE/nextflow/main.nf" \
        -c "$WORKTREE/nextflow/nextflow.config" \
        -profile BFD \
        -params-file "$RUNDIR/params_telomeres_rerun.yaml" \
        "${extra[@]}"
    ;;

timing)
    # Summarise FIND_TELOMERES task times from the newest trace file, so you can size the
    # full run (about 1 to 1.5 h of work per job) before launching all genomes.
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
sel = [r for r in rows[1:] if "FIND_TELOMERES" in r[ix["name"]]]
rt = sorted(secs(r[ix["realtime"]]) for r in sel)
ok = [r for r in sel if r[ix["status"]] == "COMPLETED"]
if not rt:
    sys.exit("no FIND_TELOMERES tasks in this trace")
n = len(rt)
print(f"FIND_TELOMERES tasks: {n} ({len(ok)} COMPLETED)")
print(f"realtime seconds: median {st.median(rt):.1f}, mean {st.mean(rt):.1f}, p90 {rt[int(0.9*(n-1))]:.1f}, max {rt[-1]:.1f}")
med = max(st.median(rt), 1.0)
print(f"genomes per 1 h job at the median: about {int(3600/med)}  (this is a guess from {n} tasks)")
print(f"all 23,580 genomes at the mean: about {st.mean(rt)*23580/3600:.1f} CPU hours")
PYEOF
    ;;

merge)
    check_fix_in_worktree
    check_no_live_run
    # The merge must read regenerated files. Refuse if stage or run did not finish.
    [ -e "$BACKUP/.stage_complete" ] || die "stage has not completed. Run stage --go first."
    ts=$(fix_ts)
    total=$(count_files "$STORE")
    newer=$(count_files "$STORE" "$ts")
    old_n=$((total - newer))
    echo "store files: $total, older than the fix: $old_n"
    [ "$old_n" = 0 ] || [ "${FORCE_MERGE:-0}" = 1 ] \
        || die "$old_n telomere files predate the fix. Finish run first (FORCE_MERGE=1 to override)."
    [ "$total" -gt 0 ] || die "store is empty"
    export_nf_env
    # The merge reruns when a result file is newer than its table. Rewritten telomere
    # files are newer, so MERGE_TELOMERES rebuilds both parquet files. BUILD_DUCKDB runs
    # after a merge.
    NXF_OPTS="-Xms512m -Xmx12g" \
    nextflow run "$WORKTREE/nextflow/main.nf" \
        -c "$WORKTREE/nextflow/nextflow.config" \
        -profile BFD \
        -params-file "$RUNDIR/params_merge_after_telomeres.yaml" \
        -resume
    ;;

verify)
    ts=$(fix_ts)
    old=$(count_files "$BACKUP")
    new=$(count_files "$STORE")
    newer=$(count_files "$STORE" "$ts")
    echo "backup telomere files: $old"
    echo "store telomere files:  $new  (newer than fix: $newer)"
    [ "$newer" = "$new" ] || echo "WARNING: some store files predate the fix"
    [ "$new" -ge "$old" ] || echo "WARNING: fewer genomes than before. Check logs for failed tasks."
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
