# slurm_rebalance.py

Script: `bin/slurm_rebalance.py` (Python 3, standard library only).

## Purpose

A pending job can wait for two reasons.
Its (account, partition) association is at a group limit (`AssocGrpMemLimit`, `AssocGrpCpuLimit`), or its partition has no free node (`Priority`, `Resources`).
Other partitions can have idle CPUs at the same time.

The script moves pending jobs to partitions with free capacity.
It changes `Partition`, `Account` and `QOS` together.
The jobs keep their job IDs, so Nextflow keeps tracking them.

## Partition settings

Set with `--part-map`. Checked 2026-10-28.

| Partition | Account | QOS | Limits that apply |
|---|---|---|---|
| `stajichlab` | `stajichlab` | normal | No user row limit. The parent account cap (cpu=1472, all users) applies. Nodes c01-c03 (abu_dhabi) have no BMI2. |
| `preempt` | `preempt` | normal | User row cpu=384, mem=1T. Running jobs can be preempted (REQUEUE). 29 abu_dhabi nodes are in the partition. |
| `epyc` | `stajichlab` | normal | No own row. Uses the `stajichlab` user row (cpu=384, mem=1T). All nodes are milan. |
| `highclock` | `stajichlab` | normal | User row cpu=32, mem=256G. Nodes have no CPU-family feature (relaxed). |
| `exfab` | `exfab` | normal | User row cpu=32, mem=256G. One GPU node (gpu12). The node has no CPU-family feature (relaxed). |
| `short` | `stajichlab` | `short` | User row cpu=384, mem=1T. QOS `short` has MaxWall 02:00:00 and `DenyOnLimit`. The partition has the same 148 nodes as `preempt`. |

`exfab` jobs under `Account=stajichlab` also run there, but they count against the parent account. The script uses `Account=exfab` so the Slurm-enforced 32-CPU row applies.

## What each cycle does

1. Read limits and usage for the user rows and the parent account rows from `scontrol show assoc_mgr`. A move needs headroom in both.
2. Read free CPUs and unallocated memory (`Memory - AllocMem`) per node from `sinfo`. A node in several partitions is counted once.
3. Select candidate jobs: name matches `--name-regex`, no dependency, pending longer than `--min-age`, and Reason in `--reasons`.
   A job pending on `Priority` or `Resources` is a candidate only if no node in its current partition fits it.
4. Place largest jobs first.
   - The first partition in `--targets` is tried first.
   - A job with a time limit of `--short-max` or less (default 2:00:00) tries `short` next.
   - The other partitions follow, most idle fraction first. The idle fraction is limited by free memory for the job size.
   - In a partition, the node with the most CPUs left after the job is used.
5. With `--apply`, run `scontrol update JobId=N Partition=P Account=A QOS=Q [Features=F]`.

## How to run

Run from the run directory, for example `do_annotation_wave2`. The log file is written to the current directory.

```bash
# 1. Dry run. Prints the planned moves. Changes nothing.
python3 ../bin/slurm_rebalance.py

# 2. Move a small batch first.
python3 ../bin/slurm_rebalance.py --apply --max-moves 20

# 3. Repeat every 5 minutes. Run under nohup, tmux or a small SLURM job.
python3 ../bin/slurm_rebalance.py --apply --loop 300

# 4. Undo moves for jobs that are still pending.
python3 ../bin/slurm_rebalance.py --revert slurm_rebalance_moves.tsv --apply
```

Do not run it from a SLURM job step with `srun`. Use `sbatch` for any new job. `srun` inside a job creates a step in that job.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--apply` | off | Run `scontrol update`. Without it, the script is a dry run. |
| `--loop SEC` | off | Repeat every SEC seconds. |
| `--targets` | `stajichlab,preempt,epyc,highclock,exfab` | Partitions to fill. The first is tried first. The others are ordered by idle fraction. |
| `--short` | `short` | Partition for jobs of `--short-max` or less. `''` turns it off. |
| `--short-max` | `2:00:00` | Longest time limit sent to `--short`. |
| `--part-map` | see table above | `name=account:qos,...` set on moved jobs. |
| `--relaxed` | `highclock:maxcpu=4:cap=32,exfab:cap=32:unconstrained` | Partitions with verified CPUs but no CPU-family feature. See below. |
| `--reasons` | `AssocGrpMemLimit,AssocGrpCpuLimit,AssocGrpCPUMinutesLimit,Resources,Priority` | Pending reasons that make a job eligible. |
| `--name-regex` | `^nf-` | Only jobs whose name matches are moved. |
| `--max-moves` | 100 | Maximum moves per cycle. |
| `--min-age` | 2 | Minutes. Jobs pending less than this are left alone. |
| `--cooldown` | 30 | Minutes. A job is not moved again within this time of its last move. |
| `--max-per-job` | 3 | Most moves for one job. |
| `--pad` | 0.03 | Fraction of each group limit left unused. |
| `--keep-cpus` | 0 | Idle CPUs left free on each node. |
| `--log` | `slurm_rebalance_moves.tsv` | Move log. `--revert` reads it. |

### `--relaxed` format

`name[:maxcpu=N][:cap=N][:unconstrained],...`

- `maxcpu=N`: only jobs that request N CPUs or fewer.
- `cap=N`: maximum total CPUs of yours in that partition, used only when Slurm has no row for it. Memory cap is 8G per CPU.
- `unconstrained`: only jobs that already have no `--constraint`.
- A job with a constraint that moves to a relaxed partition has its `Features` cleared.

## Safety rules

- Nodes c01-c03 (abu_dhabi) and Ivy Bridge nodes have no BMI2. Trinity and funannotate die there with SIGILL.
- A job with a `--constraint` goes only to nodes with a matching feature.
- A job with no constraint can move to any partition with at least one safe node (milan, genoa, rome, sapphire, cascade, broadwell).
  If the partition also has unsafe nodes, the job gets `Features=broadwell|cascade|genoa|milan|rome|sapphire`.
  This keeps it off abu_dhabi nodes. `stajichlab`, `preempt` and `short` are in this class.
- Add a partition to `--relaxed` only after you check its CPUs for BMI2. Test with `sbatch` (not `srun`) and read `/proc/cpuinfo`.
- Write job output to `/bigdata`. Node-local scratch does not exist on the other nodes, and the job fails with exit 53.
- Do not use `scontrol update ExcNodeList=`. It replaces the whole list.
- A job moves into `short` only if its time limit is `--short-max` or less. `DenyOnLimit` rejects longer jobs.
- The script always sets QOS on a move, so a move out of `short` sets `QOS=normal`.

## Move log and revert

Tab-separated columns:
`time, job, orig_account, orig_partition, new_account, new_partition, status, orig_features, orig_qos, new_qos, new_features, time_limit_s`

- `orig_features` is empty when the move did not change Features. `-` means the job had no constraint. `--revert` then sends `Features=`.
- `orig_qos` of `-` means squeue showed no QOS (pending jobs can show `(null)`). `--revert` then uses the default QOS of the original partition (`short` for `short`, `normal` otherwise).
- `--revert` acts only on jobs that are still pending and still in the logged new account and partition.
- `--revert` reads old 8-column logs.

## Known limits

- The script plans from `sinfo`, `squeue` and `assoc_mgr` data. Slurm can still place a job differently.
- Usage in `assoc_mgr` counts all jobs in the association. Moving jobs between partitions of one account does not add CPU headroom under the parent account cap. It helps only where nodes have free memory, or where another account (`preempt`, `exfab`) has headroom.
- Running jobs cannot be moved. Only pending jobs can.
- Jobs moved to `preempt` can be preempted and requeued. 37 jobs ended `PREEMPTED` in 12 hours on 2026-10-28.
- Tested on 2026-10-28 in dry run and revert dry run only. A move with `QOS=` and `Features=` in one `scontrol update` has not been run with `--apply` yet.
- Not verified: that a move to `short` without `QOS=short` fails. The script always sets it.
