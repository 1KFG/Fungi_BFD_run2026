#!/usr/bin/env python3
"""Sum SLURM CPU-hours for Nextflow (nf-*) tasks of the rc.6 runs.

Reads `sacct` for the current user. Counts jobs whose WorkDir is under one of
the --workdir prefixes and whose Start is at or after --since.
CPU-hours = CPUTimeRAW / 3600 (allocated CPUs x elapsed time).
Writes a per-day, per-process, per-state TSV and prints a summary.
"""
import argparse
import collections
import getpass
import re
import subprocess
import sys

ROOT = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
FIELDS = "JobID,JobName%60,Start,CPUTimeRAW,State,WorkDir%200"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--since", default="2026-10-03T21:47:00",
                    help="first rc.6 task start (default: first rc.6 launch)")
    ap.add_argument("--workdir", action="append",
                    help="WorkDir substring to include (repeatable); "
                         "default: do_annotation_wave1/work")
    ap.add_argument("--user", default=getpass.getuser())
    ap.add_argument("--sacct-file", help="read saved sacct -P -n output instead of running sacct")
    ap.add_argument("--out", default=f"{ROOT}/docs/rc6_compute_daily.tsv")
    a = ap.parse_args()
    wds = a.workdir or ["do_annotation_wave1/work"]

    if a.sacct_file:
        text = open(a.sacct_file).read()
    else:
        text = subprocess.run(
            ["sacct", "-u", a.user, "-S", a.since[:10], "-E", "now", "-X", "-P", "-n",
             "--format", FIELDS],
            check=True, capture_output=True, text=True).stdout

    agg = collections.defaultdict(lambda: [0, 0.0])
    for line in text.splitlines():
        f = line.split("|")
        if len(f) != 6:
            continue
        _, name, start, cpusec, state, wd = f
        if not name.startswith("nf-") or start in ("None", "Unknown"):
            continue
        if start < a.since or not any(w in wd for w in wds):
            continue
        proc = re.sub(r"_\(.*", "", name)[3:]
        key = (start[:10], proc, state.split()[0])
        agg[key][0] += 1
        agg[key][1] += int(cpusec) / 3600

    with open(a.out, "w") as o:
        o.write("day\tprocess\tstate\tjobs\tcpu_hours\n")
        for k in sorted(agg):
            o.write("%s\t%s\t%s\t%d\t%.1f\n" % (*k, *agg[k]))

    tot = sum(v[1] for v in agg.values())
    print(f"since {a.since}  workdirs {wds}")
    print(f"total jobs {sum(v[0] for v in agg.values())}  CPU-hours {tot:.0f}")
    for title, idx in (("process", 1), ("state", 2), ("day", 0)):
        d = collections.defaultdict(float)
        for k, v in agg.items():
            d[k[idx]] += v[1]
        print(f"-- by {title}")
        for k, v in sorted(d.items(), key=lambda x: (-x[1] if idx else x[0])):
            print(f"  {k:60s} {v:10.0f}")
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
