#!/usr/bin/env python3
"""Repair genome_annotation/<SPECIES>/training symlinks.

Background (see .living/learnings.md 2026-09-10, ../Fungi_BFD/.living):
genome_annotation_training/<SPECIES>/training is the single canonical
training-output directory for a species/strain. do_annotation*/genome_annotation_training
are themselves relative symlinks (../genome_annotation_training) back to that
one canonical directory -- confirmed via matching inodes. genome_annotation/<SPECIES>/training
is supposed to be an absolute symlink to that canonical path.

Two failure modes found on disk (2026-09-10 census):
  - broken symlink: 776 cases, all pointing at .../BFD/Fungi_BFD/do_annotation*/...
    (the sibling project root) instead of .../BFD/Fungi_BFD_runs/... -- recoverable
    by relinking to the canonical absolute path directly.
  - missing entirely: ~10,064 cases -- most have no training data at all (never had
    RNA-seq), but some may have valid canonical training data that was never linked.

This script does NOT touch:
  - entries where training/ is already a real (non-symlink) directory (776 cases)
  - entries where training/ is already a valid symlink resolving to a real directory
    (whether via the canonical top-level path or an indirect do_annotation_X path --
    both resolve to the same canonical inode, so both are left alone)

A species is only (re)linked when the canonical target exists AND passes a
completeness check mirroring FUNANNOTATE_TRAIN's own success markers
(nextflow/modules/funannotate/utils.nf, predict/FUNANNOTATE_TRAIN/main.nf):
  - training/funannotate_train.pasa.gff3 exists and is non-empty (FUNANNOTATE_TRAIN
    `touch`es this file EMPTY on failure/incomplete-Trinity paths, so size>0 is the
    real success signal, not mere existence)
  - training/transcript.alignments.bam exists and is non-empty
  - training/.pasa_train_failed marker absent

Default is a dry run -- prints a report and writes an action-log TSV, makes no
filesystem changes. Pass --apply to actually create/replace symlinks.
"""
import argparse
import csv
import os
import sys
from pathlib import Path


def is_complete_training(training_dir: Path) -> tuple[bool, str]:
    if not training_dir.is_dir():
        return False, "canonical_training_dir_missing"
    failed_marker = training_dir / ".pasa_train_failed"
    if failed_marker.exists():
        return False, "pasa_train_failed_marker_present"
    gff3 = training_dir / "funannotate_train.pasa.gff3"
    if not gff3.exists() or gff3.stat().st_size == 0:
        return False, "pasa_gff3_missing_or_empty"
    bam = training_dir / "transcript.alignments.bam"
    if not bam.exists() or bam.stat().st_size == 0:
        return False, "transcript_bam_missing_or_empty"
    return True, "complete"


def classify_link(link_path: Path) -> str:
    if not link_path.exists() and not link_path.is_symlink():
        return "missing"
    if link_path.is_symlink():
        return "valid_symlink" if link_path.exists() else "broken_symlink"
    if link_path.is_dir():
        return "real_directory"
    return "unexpected_type"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", default=os.getcwd(),
                     help="Fungi_BFD_runs repo root (default: cwd)")
    ap.add_argument("--apply", action="store_true",
                     help="Actually create/replace symlinks. Default is dry-run/report-only.")
    ap.add_argument("--report", default="training_symlink_repair_report.tsv",
                     help="TSV report path (written in both dry-run and apply modes)")
    args = ap.parse_args()

    repo_root = Path(args.repo_root).resolve()
    genome_annotation = repo_root / "genome_annotation"
    genome_annotation_training = repo_root / "genome_annotation_training"

    if not genome_annotation.is_dir() or not genome_annotation_training.is_dir():
        sys.exit(f"expected {genome_annotation} and {genome_annotation_training} to exist under {repo_root}")

    rows = []
    counts = {}

    species_dirs = sorted(p for p in genome_annotation.iterdir() if p.is_dir())
    for sp_dir in species_dirs:
        species = sp_dir.name
        link_path = sp_dir / "training"
        canonical_target = genome_annotation_training / species / "training"

        state = classify_link(link_path)

        if state == "real_directory":
            counts[state] = counts.get(state, 0) + 1
            rows.append([species, state, "", "skipped_real_directory"])
            continue

        if state == "valid_symlink":
            counts[state] = counts.get(state, 0) + 1
            rows.append([species, state, str(link_path.resolve()), "skipped_already_valid"])
            continue

        # state is "missing" or "broken_symlink" -- candidate for (re)linking
        complete, reason = is_complete_training(canonical_target)
        if not complete:
            key = f"{state}_no_valid_source_({reason})"
            counts[key] = counts.get(key, 0) + 1
            rows.append([species, state, str(canonical_target), f"not_linked_{reason}"])
            continue

        counts[f"{state}_relinked"] = counts.get(f"{state}_relinked", 0) + 1
        action = "would_relink" if not args.apply else "relinked"
        rows.append([species, state, str(canonical_target), action])

        if args.apply:
            if link_path.is_symlink() or link_path.exists():
                link_path.unlink()
            link_path.symlink_to(canonical_target)

    with open(args.report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "prior_state", "canonical_target", "action"])
        w.writerows(rows)

    mode = "APPLY" if args.apply else "DRY RUN"
    print(f"[{mode}] {len(species_dirs)} species scanned under {genome_annotation}")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")
    print(f"report written to {args.report}")
    if not args.apply:
        print("Dry run only -- re-run with --apply to create/replace symlinks.")


if __name__ == "__main__":
    main()
