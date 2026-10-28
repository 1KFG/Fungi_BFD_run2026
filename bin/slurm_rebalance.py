#!/usr/bin/env python3
"""Move pending SLURM jobs to partitions with free capacity.

Why: on UCR HPCC a job waits when its own (account, partition) association is at its group
limit (for example `preempt`: cpu=384, mem=1T), or when its partition has no free node,
while other partitions (`stajichlab`, `epyc`, `highclock`, `exfab`, `short`) have idle CPUs.
A group limit belongs to the association, so a move can change Account and QOS as well as Partition.

What it does, per cycle:
  1. Lists your pending jobs whose Reason is in --reasons and whose name matches --name-regex.
     A job pending on Priority or Resources moves only if its current partition has no node that fits it.
  2. Reads the GrpTRES limit and the current usage of each association with `scontrol show assoc_mgr`.
     Both the (account, partition) row and the parent account row must have headroom.
  3. Reads free CPUs and unallocated memory (Memory - AllocMem) per node with sinfo.
     A node that is in several partitions is counted once.
  4. Places each job on a node that fits.  The first partition in --targets is tried first.
     Jobs of --short-max or less try --short next.  The other partitions follow, most idle fraction first.
  5. Runs `scontrol update JobId=N Partition=P Account=A QOS=Q [Features=F]` (only with --apply).

Each partition has an account and a QOS (--part-map).  For example `preempt` needs Account=preempt,
and `short` needs QOS=short.

Safety rules:
  * Dry run unless --apply is given.
  * Nodes c01-c03 (abu_dhabi) and Ivy Bridge nodes lack BMI2.  Trinity/funannotate die there with SIGILL.
    A job with no --constraint that moves into a partition with such nodes gets
    Features=<safe CPU families>, so it runs only on safe nodes.
  * A job with a --constraint goes only to nodes that match it.
  * A relaxed partition (--relaxed) has verified CPUs but no CPU-family feature.  Features are cleared there.
  * A job moves at most --max-per-job times, and not again within --cooldown minutes.
  * Every move is appended to --log (TSV) so --revert can undo it.

Usage:
  slurm_rebalance.py                       # dry run, one cycle
  slurm_rebalance.py --apply               # move jobs, one cycle
  slurm_rebalance.py --apply --loop 300    # repeat every 5 minutes
  slurm_rebalance.py --short ''            # do not use the short partition
  slurm_rebalance.py --revert moves.tsv --apply   # put moved jobs back (still-pending only)
"""
import argparse
import getpass
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime

UNITS = {"K": 1 / 1024.0, "M": 1.0, "G": 1024.0, "T": 1024.0 * 1024.0}  # to MB
# CPU families that run Trinity/funannotate.  abu_dhabi and Ivy Bridge nodes lack BMI2 (SIGILL).
SAFE_FEATURES = {"milan", "genoa", "rome", "sapphire", "cascade", "broadwell"}
SAFE_CONSTRAINT = "|".join(sorted(SAFE_FEATURES))
BAD_NODE_STATES = ("drain", "down", "maint", "resv", "reserv", "fail", "inval", "unk", "boot",
                   "no_respond", "futur", "power")
DEFAULT_PART_MAP = ("stajichlab=stajichlab:normal,epyc=stajichlab:normal,highclock=stajichlab:normal,"
                    "short=stajichlab:short,preempt=preempt:normal,exfab=exfab:normal")
DEFAULT_RELAXED = "highclock:maxcpu=4:cap=32,exfab:cap=32:unconstrained"
NO_FEATURES = ("", "(null)")
EMPTY = "-"  # log marker for an originally empty value


def run(cmd, timeout=60):
    """Run a command and return stdout.  Return '' and warn on failure."""
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        print("WARN: %s failed: %s" % (" ".join(cmd), exc), file=sys.stderr)
        return ""
    if res.returncode != 0:
        print("WARN: %s exit %d: %s" % (" ".join(cmd), res.returncode, res.stderr.strip()[:200]),
              file=sys.stderr)
    return res.stdout


def mem_to_mb(text, cpus=1):
    """Parse a SLURM memory string such as 12G, 4000M, 4Gc (per CPU), 2T into MB."""
    m = re.match(r"^([\d.]+)\s*([KMGT]?)([cn]?)$", text.strip(), re.I)
    if not m:
        return 0.0
    val = float(m.group(1)) * UNITS.get(m.group(2).upper() or "M", 1.0)
    return val * cpus if m.group(3).lower() == "c" else val


def parse_timelimit(text):
    """Parse D-HH:MM:SS, HH:MM:SS, MM:SS or MM into seconds.  UNLIMITED and NOT_SET give inf."""
    t = text.strip()
    if not t or t.upper() in ("UNLIMITED", "NOT_SET", "INVALID", "PARTITION_LIMIT"):
        return float("inf")
    try:
        if "-" in t:
            days, t = t.split("-", 1)
            h, m, s = ([int(x) for x in t.split(":")] + [0, 0])[:3]
            return int(days) * 86400 + h * 3600 + m * 60 + s
        parts = [int(x) for x in t.split(":")]
    except ValueError:
        return float("inf")
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    return parts[0] * 60


def parse_part_map(text):
    """Parse 'name=account:qos,...' into {partition: (account, qos)}."""
    out = {}
    for item in filter(None, text.split(",")):
        name, _, rest = item.partition("=")
        account, _, qos = rest.partition(":")
        out[name] = (account, qos or "normal")
    return out


def parse_grptres(text):
    """Parse 'cpu=384(251),mem=1048576(226304),node=N(13)' into {'cpu': [limit, used], 'mem': [...]}.
    A limit of None means no limit.  Memory is in MB."""
    out = {}
    for item in text.split(","):
        m = re.match(r"^(cpu|mem)=(N|\d+)\((\d+)\)$", item)
        if m:
            out[m.group(1)] = [None if m.group(2) == "N" else float(m.group(2)), float(m.group(3))]
    return out


def get_assoc(user, accounts):
    """Return {(account, partition): {'cpu': [limit, used], 'mem': [limit, used]}}.
    Your own association rows use the partition name ('' is the row for partitions without their own row).
    The parent account row uses partition None."""
    out = run(["scontrol", "show", "assoc_mgr", "flags=assoc", "users=" + user,
               "accounts=" + ",".join(accounts)])
    assoc = {}
    cur = None
    for line in out.splitlines():
        m = re.match(r"^ClusterName=\S+ Account=(\S+) UserName=(\S*) Partition=(\S*) ", line)
        if m:
            acct, uname, part = m.groups()
            if uname == "":
                cur = (acct, None) if part == "" else None
            elif uname.startswith(user + "("):
                cur = (acct, part)
            else:
                cur = None
            if cur is not None:
                assoc[cur] = {"cpu": [None, 0.0], "mem": [None, 0.0]}
            continue
        m = re.match(r"^\s+GrpTRES=(\S+)", line)
        if m and cur is not None:
            assoc[cur].update(parse_grptres(m.group(1)))
    return assoc


def assoc_key(assoc, account, partition):
    """Pick the partition-specific association if it exists, else the account-wide one."""
    if (account, partition) in assoc:
        return (account, partition)
    return (account, "")


def get_jobs(user, state):
    fmt = "%i\t%j\t%a\t%P\t%r\t%C\t%m\t%f\t%E\t%T\t%l\t%q\t%V"  # tab: Features itself contains '|'
    out = run(["squeue", "-u", user, "-h", "-t", state, "-o", fmt])
    now = datetime.now()
    jobs = []
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 13:
            continue
        cpus = int(p[5])
        try:
            age = (now - datetime.strptime(p[12], "%Y-%m-%dT%H:%M:%S")).total_seconds() / 60.0
        except ValueError:
            age = float("inf")
        jobs.append({"id": p[0], "name": p[1], "account": p[2], "partition": p[3], "reason": p[4],
                     "cpus": cpus, "mem": mem_to_mb(p[6], cpus), "features": p[7],
                     "dep": p[8], "state": p[9], "limit_s": parse_timelimit(p[10]),
                     "qos": p[11], "age_min": age})
    return jobs


def get_nodes(partitions):
    """Return (nodes, part_nodes, part_class).
    nodes is keyed by node name, so a node in several partitions is counted once.
    Each node carries idle CPUs and unallocated memory (Memory - AllocMem, MB).
    part_class[p] is 'safe' when every node has a SAFE_FEATURES CPU family, 'mixed' when some do,
    'none' when none do."""
    out = run(["sinfo", "-h", "-N", "-p", ",".join(partitions), "-O",
               "NodeList:40,PartitionName:30,CPUsState:24,Memory:14,AllocMem:14,StateLong:16,Features:120"])
    nodes = {}
    part_nodes = defaultdict(list)
    for line in out.splitlines():
        p = line.split()
        if len(p) == 6:
            p.append("")
        if len(p) < 7:
            continue
        name, part, cpus, mem, alloc, state, feats = p[:7]
        part = part.rstrip("*")
        feat = set(feats.split(","))
        if name not in nodes:
            usable = not any(state.lower().startswith(bad) for bad in BAD_NODE_STATES)
            idle, total = int(cpus.split("/")[1]), int(cpus.split("/")[3])
            try:
                free_mem = float(mem) - float(alloc)
            except ValueError:
                free_mem = 0.0
            nodes[name] = {"name": name, "cpu": idle if usable else 0, "total": total,
                           "mem": free_mem if usable else 0.0, "features": feat}
        part_nodes[part].append(nodes[name])
    part_class = {}
    for part, plist in part_nodes.items():
        safe = sum(1 for n in plist if n["features"] & SAFE_FEATURES)
        part_class[part] = "safe" if safe == len(plist) else "mixed" if safe else "none"
    return nodes, part_nodes, part_class


def feature_ok(job_features, node_features):
    """Job --constraint 'a|b|c' (OR only).  Anything with & or () is treated as no match."""
    if re.search(r"[&()\[\]*]", job_features):
        return False
    return bool(set(job_features.split("|")) & node_features)


def fits_keys(job, keys, assoc, pad):
    """True when every association key has headroom for the job."""
    for key in keys:
        for res, field in (("cpu", "cpus"), ("mem", "mem")):
            limit, used = assoc[key][res] if key in assoc else (None, 0.0)
            if limit is not None and used + job[field] > limit * (1 - pad):
                return False
    return True


def charge_keys(job, keys, assoc):
    for key in keys:
        if key in assoc:
            assoc[key]["cpu"][1] += job["cpus"]
            assoc[key]["mem"][1] += job["mem"]


def parse_relaxed(text):
    """Parse 'highclock:maxcpu=4:cap=32,exfab:cap=32:unconstrained' into {part: options}.

    A relaxed partition has CPUs verified to run Trinity/funannotate (BMI2 present), but its
    nodes carry no milan/genoa/... feature.  A job moved there with a constraint gets its Features cleared.
      maxcpu=N       only jobs asking for <= N CPUs
      cap=N          cap your total CPUs in that partition (also mem = 8G per CPU)
      unconstrained  only jobs that already have no --constraint
    """
    out = {}
    for item in filter(None, text.split(",")):
        name, *opts = item.split(":")
        spec = {"maxcpu": None, "cap": None, "unconstrained": False}
        for opt in opts:
            if opt == "unconstrained":
                spec["unconstrained"] = True
            else:
                key, _, val = opt.partition("=")
                spec[key] = int(val)
        out[name] = spec
    return out


def read_move_history(path):
    """Return {job id: (move count, time of last ok move)} from the move log."""
    hist = {}
    try:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 7 or f[6] != "ok":
                    continue
                try:
                    when = datetime.strptime(f[0], "%Y-%m-%d %H:%M:%S")
                except ValueError:
                    continue
                count = hist.get(f[1], (0, when))[0]
                hist[f[1]] = (count + 1, when)
    except FileNotFoundError:
        pass
    return hist


def job_fits_here(job, pnodes, keep):
    """True when a node in pnodes could start the job now.  Used to leave well-placed jobs alone."""
    unconstrained = job["features"] in NO_FEATURES
    return any((unconstrained or feature_ok(job["features"], n["features"]))
               and n["cpu"] >= job["cpus"] + keep and n["mem"] >= job["mem"] for n in pnodes)


def pick_node(job, cls, spec, pnodes, keep):
    """Return (node, new_features) for the job in one partition, or None.
    new_features is None to keep the job's Features, '' to clear them, or a constraint string."""
    unconstrained = job["features"] in NO_FEATURES
    new_feat = None
    if spec is not None:  # relaxed partition: CPUs verified, nodes have no family feature
        if spec["maxcpu"] and job["cpus"] > spec["maxcpu"]:
            return None
        if spec["unconstrained"] and not unconstrained:
            return None
        node_ok = lambda n: True
        if not unconstrained:
            new_feat = ""
    elif not unconstrained:
        node_ok = lambda n: feature_ok(job["features"], n["features"])
    elif cls == "safe":
        node_ok = lambda n: True
    elif cls == "mixed":
        node_ok = lambda n: bool(n["features"] & SAFE_FEATURES)
        new_feat = SAFE_CONSTRAINT
    else:
        return None
    fit = [n for n in pnodes if node_ok(n) and n["cpu"] >= job["cpus"] + keep and n["mem"] >= job["mem"]]
    if not fit:
        return None
    return max(fit, key=lambda n: (n["cpu"] - job["cpus"], n["mem"] - job["mem"])), new_feat


def idle_fraction(pnodes, job):
    """Fraction of the partition's CPUs that this job size could still use, limited by free memory."""
    per_cpu = job["mem"] / max(job["cpus"], 1)
    usable = sum(min(n["cpu"], n["mem"] / per_cpu if per_cpu else n["cpu"]) for n in pnodes)
    total = sum(n["total"] for n in pnodes)
    return usable / total if total else 0.0


def plan_moves(args, user):
    part_map = args.part_map
    assoc = get_assoc(user, sorted({a for a, _ in part_map.values()}))
    if not assoc:
        sys.exit("ERROR: no association data from `scontrol show assoc_mgr`; refusing to plan")
    running = get_jobs(user, "R")
    pending = get_jobs(user, "PD")
    for part, spec in args.relaxed.items():
        acct = part_map[part][0]
        if spec["cap"] and (acct, part) not in assoc:  # limit set by this script; count running jobs
            use = [j for j in running if j["account"] == acct and j["partition"] == part]
            assoc[(acct, part)] = {"cpu": [float(spec["cap"]), float(sum(j["cpus"] for j in use))],
                                   "mem": [spec["cap"] * 8 * 1024.0, sum(j["mem"] for j in use)]}
    reasons = set(args.reasons.split(","))
    name_re = re.compile(args.name_regex)
    hist = read_move_history(args.log)
    now = datetime.now()

    skipped = Counter()
    cands = []
    for j in pending:
        if j["reason"] not in reasons or not name_re.search(j["name"]):
            continue
        if j["dep"] not in ("", "(null)"):
            continue
        if j["age_min"] < args.min_age:
            skipped["younger than --min-age"] += 1
        elif j["id"] in hist and (hist[j["id"]][0] >= args.max_per_job or
                                  (now - hist[j["id"]][1]).total_seconds() < args.cooldown * 60):
            skipped["moved recently or too often"] += 1
        else:
            cands.append(j)

    targets = args.targets + [p for p in args.relaxed if p not in args.targets]
    if args.short and args.short not in targets:
        targets.append(args.short)
    cur_parts = {p for j in cands for p in j["partition"].split(",")}
    nodes, part_nodes, part_class = get_nodes(sorted(set(targets) | cur_parts))
    print("%s pending=%d candidates=%d" % (now.strftime("%F %T"), len(pending), len(cands)))
    for part in targets:
        plist = part_nodes.get(part, [])
        print("  part %-10s class=%-6s idle_cpu=%5d free_mem=%7.0fG nodes=%d" %
              (part, "relaxed" if part in args.relaxed else part_class.get(part, "none"),
               sum(n["cpu"] for n in plist), sum(n["mem"] for n in plist) / 1024.0, len(plist)))
    for key in sorted(assoc, key=lambda k: (k[0], k[1] or "~")):
        row = assoc[key]
        if row["cpu"][0] is None and row["mem"][0] is None and not row["cpu"][1]:
            continue
        print("  assoc %-24s use cpu=%5d mem=%6.0fG  limit cpu=%s mem=%s" % (
            "%s|%s" % (key[0], key[1] if key[1] is not None else "<account>"),
            row["cpu"][1], row["mem"][1] / 1024.0,
            "%d" % row["cpu"][0] if row["cpu"][0] is not None else "-",
            "%.0fG" % (row["mem"][0] / 1024.0) if row["mem"][0] is not None else "-"))

    # Largest jobs first so big requests are not starved by small ones.
    cands.sort(key=lambda j: (-j["cpus"], -j["mem"], j["id"]))
    moves = []
    for job in cands:
        if len(moves) >= args.max_moves:
            break
        if job["reason"] in ("Priority", "Resources"):
            here = [n for p in job["partition"].split(",") for n in part_nodes.get(p, [])]
            if job_fits_here(job, here, args.keep_cpus):
                skipped["current partition has a node that fits"] += 1
                continue
        order = [args.targets[0]]
        if args.short and job["limit_s"] <= args.short_max and args.short not in order:
            order.append(args.short)
        rest = [p for p in targets if p not in order]
        rest.sort(key=lambda p: -idle_fraction(part_nodes.get(p, []), job))
        placed = False
        for part in order + rest:
            if part not in part_map:
                continue
            account, qos = part_map[part]
            if job["partition"] == part and job["account"] == account:
                continue
            keys = [assoc_key(assoc, account, part), (account, None)]
            if not fits_keys(job, keys, assoc, args.pad):
                continue
            pick = pick_node(job, part_class.get(part, "none"), args.relaxed.get(part),
                             part_nodes.get(part, []), args.keep_cpus)
            if pick is None:
                continue
            node, new_feat = pick
            node["cpu"] -= job["cpus"]
            node["mem"] -= job["mem"]
            charge_keys(job, keys, assoc)
            moves.append({"job": job, "part": part, "account": account, "qos": qos,
                          "node": node["name"], "features": new_feat})
            placed = True
            break
        if not placed:
            skipped["no partition fits"] += 1
    for reason, count in sorted(skipped.items()):
        print("  skipped %d: %s" % (count, reason))
    by_part = Counter(m["part"] for m in moves)
    print("  planned by partition: %s" % (", ".join("%s=%d" % kv for kv in sorted(by_part.items())) or "none"))
    return moves


def apply_moves(args, moves):
    stamp = datetime.now().strftime("%F %T")
    for mv in moves:
        job = mv["job"]
        cmd = ["scontrol", "update", "JobId=" + job["id"], "Partition=" + mv["part"],
               "Account=" + mv["account"], "QOS=" + mv["qos"]]
        if mv["features"] is not None:
            cmd.append("Features=" + mv["features"])
        feat_note = ""
        if mv["features"] is not None:
            feat_note = " [clear constraint]" if mv["features"] == "" else " [Features=%s]" % mv["features"]
        label = "%s %s %dcpu %.0fG %s|%s -> %s|%s qos=%s (fits %s)%s" % (
            job["id"], job["name"][:50], job["cpus"], job["mem"] / 1024.0,
            job["account"], job["partition"], mv["account"], mv["part"], mv["qos"], mv["node"], feat_note)
        if not args.apply:
            print("  DRY  " + label)
            continue
        res = subprocess.run(cmd, capture_output=True, text=True)
        status = "ok" if res.returncode == 0 else "FAIL " + res.stderr.strip()[:120]
        print("  MOVE " + label + " " + status)
        if mv["features"] is None:
            orig_feat = ""  # Features not touched, so revert leaves them alone
        else:
            orig_feat = EMPTY if job["features"] in NO_FEATURES else job["features"]
        orig_qos = EMPTY if job["qos"] in NO_FEATURES else job["qos"]
        limit = "" if job["limit_s"] == float("inf") else "%d" % job["limit_s"]
        with open(args.log, "a") as fh:
            fh.write("\t".join([stamp, job["id"], job["account"], job["partition"],
                                mv["account"], mv["part"], status, orig_feat, orig_qos, mv["qos"],
                                mv["features"] if mv["features"] is not None else "", limit]) + "\n")


def revert(args, user):
    """Put jobs listed in a move log back to their original account, partition, QOS and Features.
    A job is reverted only if it is still pending and still in the partition and account the log moved it to.
    Old 8-column logs have no QOS: the default QOS of the original partition is used."""
    pending = {j["id"]: j for j in get_jobs(user, "PD")}
    with open(args.revert) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[6] != "ok" or f[1] not in pending:
                continue
            job = pending[f[1]]
            if job["partition"] != f[5] or job["account"] != f[4]:
                print("SKIP   %s no longer in %s|%s" % (f[1], f[4], f[5]))
                continue
            qos = f[8] if len(f) > 8 and f[8] not in ("", EMPTY) else \
                args.part_map.get(f[3].split(",")[0], ("", "normal"))[1]
            cmd = ["scontrol", "update", "JobId=" + f[1], "Partition=" + f[3], "Account=" + f[2],
                   "QOS=" + qos]
            if len(f) > 7 and f[7]:
                cmd.append("Features=" + ("" if f[7] == EMPTY else f[7]))
            if args.apply:
                res = subprocess.run(cmd, capture_output=True, text=True)
                print("REVERT %s -> %s|%s qos=%s %s" % (f[1], f[2], f[3], qos, "ok" if res.returncode == 0
                                                       else res.stderr.strip()))
            else:
                print("DRY    " + " ".join(cmd))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--apply", action="store_true", help="run scontrol update (default is dry run)")
    ap.add_argument("--loop", type=int, default=0, metavar="SEC", help="repeat every SEC seconds")
    ap.add_argument("--targets", default="stajichlab,preempt,epyc,highclock,exfab",
                    help="partitions to fill.  The first is tried first for every job; the others are "
                         "ordered by idle fraction (default %(default)s)")
    ap.add_argument("--short", default="short",
                    help="partition for jobs of --short-max or less, tried after the first target; "
                         "'' turns it off (default %(default)s)")
    ap.add_argument("--short-max", default="2:00:00", help="longest time limit sent to --short (default %(default)s)")
    ap.add_argument("--part-map", default=DEFAULT_PART_MAP,
                    help="account and QOS to set per partition: name=account:qos,... (default %(default)s)")
    ap.add_argument("--relaxed", default=DEFAULT_RELAXED,
                    help="partitions with verified BMI2 CPUs but no CPU-family feature; job "
                         "constraint is cleared on move.  Format: name[:maxcpu=N][:cap=N]"
                         "[:unconstrained],... (default %(default)s)")
    ap.add_argument("--reasons", default="AssocGrpMemLimit,AssocGrpCpuLimit,AssocGrpCPUMinutesLimit,"
                                         "Resources,Priority",
                    help="pending reasons that make a job eligible to move.  Priority and Resources jobs "
                         "move only if their current partition has no node that fits")
    ap.add_argument("--name-regex", default=r"^nf-", help="only move jobs whose name matches")
    ap.add_argument("--max-moves", type=int, default=100, help="cap on moves per cycle")
    ap.add_argument("--min-age", type=float, default=2.0, metavar="MIN",
                    help="leave jobs pending less than MIN minutes alone (default %(default)s)")
    ap.add_argument("--cooldown", type=float, default=30.0, metavar="MIN",
                    help="do not move a job again within MIN minutes of its last move (default %(default)s)")
    ap.add_argument("--max-per-job", type=int, default=3, help="most moves for one job (default %(default)s)")
    ap.add_argument("--pad", type=float, default=0.03,
                    help="fraction of each group limit to leave unused (default %(default)s)")
    ap.add_argument("--keep-cpus", type=int, default=0, help="idle CPUs to leave free on each node")
    ap.add_argument("--log", default="slurm_rebalance_moves.tsv")
    ap.add_argument("--revert", metavar="LOG", help="undo moves recorded in LOG")
    args = ap.parse_args()
    args.targets = args.targets.split(",")
    args.relaxed = parse_relaxed(args.relaxed)
    args.part_map = parse_part_map(args.part_map)
    args.short_max = parse_timelimit(args.short_max)
    for part in args.targets + list(args.relaxed) + ([args.short] if args.short else []):
        if part not in args.part_map:
            ap.error("partition %s has no entry in --part-map" % part)
    user = getpass.getuser()

    if args.revert:
        revert(args, user)
        return
    while True:
        moves = plan_moves(args, user)
        apply_moves(args, moves)
        print("  planned moves: %d%s" % (len(moves), "" if args.apply else " (dry run)"))
        if not args.loop:
            break
        time.sleep(args.loop)


if __name__ == "__main__":
    main()
