#!/usr/bin/env python3
"""Reconcile the 776 species that have a real (non-symlink) training/ directory
directly under genome_annotation/<SPECIES>/, duplicating the canonical
genome_annotation_training/<SPECIES>/training.

Policy (confirmed with Jason, 2026-09-10): genome_annotation_training/<SPECIES>/training
is canonical. For a given species, if the canonical copy is:
  - complete (funannotate_train.pasa.gff3 + transcript.alignments.bam present and
    non-empty, no .pasa_train_failed marker -- same check as fix_training_symlinks.py),
  - free of broken internal symlinks, and
  - the same age or newer than the genome_annotation/ copy (by pasa.gff3 mtime),
then delete the genome_annotation/<SPECIES>/training real directory and replace it
with a symlink to the canonical path.

Anything that doesn't clear all three checks is left untouched and reported as
ambiguous, for manual investigation -- in particular: canonical older than the
genome_annotation/ copy (canonical.pasa.gff3 mtime < real-dir's pasa.gff3 mtime),
which would mean the genome_annotation/ copy might carry newer/better data that
canonical does not yet have.

Default is a dry run (report only). Pass --apply to actually delete+relink the
resolved (non-ambiguous) cases.
"""
import argparse
import csv
import os
import shutil
import sys
from pathlib import Path


def is_complete_training(training_dir: Path) -> tuple[bool, str]:
    if not training_dir.is_dir():
        return False, "canonical_training_dir_missing"
    if (training_dir / ".pasa_train_failed").exists():
        return False, "pasa_train_failed_marker_present"
    gff3 = training_dir / "funannotate_train.pasa.gff3"
    if not gff3.exists() or gff3.stat().st_size == 0:
        return False, "pasa_gff3_missing_or_empty"
    bam = training_dir / "transcript.alignments.bam"
    if not bam.exists() or bam.stat().st_size == 0:
        return False, "transcript_bam_missing_or_empty"
    return True, "complete"


def broken_internal_symlinks(training_dir: Path) -> list[str]:
    broken = []
    for child in training_dir.iterdir():
        if child.is_symlink() and not child.exists():
            broken.append(child.name)
    return broken


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", default=os.getcwd())
    ap.add_argument("--apply", action="store_true",
                     help="Actually delete genome_annotation/<SPECIES>/training and symlink to canonical for resolved cases.")
    ap.add_argument("--report", default="training_duplicate_reconcile_report.tsv")
    args = ap.parse_args()

    repo_root = Path(args.repo_root).resolve()
    genome_annotation = repo_root / "genome_annotation"
    genome_annotation_training = repo_root / "genome_annotation_training"

    rows = []
    counts = {}

    species_dirs = sorted(p for p in genome_annotation.iterdir() if p.is_dir())
    for sp_dir in species_dirs:
        species = sp_dir.name
        ga_training = sp_dir / "training"
        if ga_training.is_symlink() or not ga_training.is_dir():
            continue  # not one of the "real directory" cases

        canon_training = genome_annotation_training / species / "training"
        complete, reason = is_complete_training(canon_training)
        if not complete:
            counts["ambiguous_canonical_incomplete"] = counts.get("ambiguous_canonical_incomplete", 0) + 1
            rows.append([species, "ambiguous", f"canonical_incomplete:{reason}", ""])
            continue

        broken = broken_internal_symlinks(canon_training)
        if broken:
            counts["ambiguous_canonical_has_broken_internal_symlinks"] = counts.get("ambiguous_canonical_has_broken_internal_symlinks", 0) + 1
            rows.append([species, "ambiguous", "canonical_broken_internal_symlinks:" + ",".join(broken), ""])
            continue

        ga_gff3 = ga_training / "funannotate_train.pasa.gff3"
        canon_gff3 = canon_training / "funannotate_train.pasa.gff3"
        if not ga_gff3.exists() or ga_gff3.stat().st_size == 0:
            counts["ambiguous_real_dir_gff3_missing_or_empty"] = counts.get("ambiguous_real_dir_gff3_missing_or_empty", 0) + 1
            rows.append([species, "ambiguous", "real_dir_gff3_missing_or_empty", ""])
            continue

        ga_mtime = ga_gff3.stat().st_mtime
        canon_mtime = canon_gff3.stat().st_mtime

        if canon_mtime < ga_mtime:
            counts["ambiguous_canonical_older_than_real_dir"] = counts.get("ambiguous_canonical_older_than_real_dir", 0) + 1
            rows.append([species, "ambiguous", f"canonical_older(canon_mtime={canon_mtime:.0f}<real_mtime={ga_mtime:.0f})", ""])
            continue

        # Resolved: canonical is complete, has no broken internal symlinks, and is
        # the same age or newer than the genome_annotation/ copy.
        counts["resolved_canonical_wins"] = counts.get("resolved_canonical_wins", 0) + 1
        action = "would_delete_and_relink" if not args.apply else "deleted_and_relinked"
        rows.append([species, "resolved", str(canon_training), action])

        if args.apply:
            shutil.rmtree(ga_training)
            ga_training.symlink_to(canon_training)

    with open(args.report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "classification", "detail", "action"])
        w.writerows(rows)

    mode = "APPLY" if args.apply else "DRY RUN"
    print(f"[{mode}] {len(rows)} real-directory species evaluated")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")
    print(f"report written to {args.report}")
    if not args.apply:
        print("Dry run only -- re-run with --apply to delete+relink resolved cases.")


if __name__ == "__main__":
    main()
