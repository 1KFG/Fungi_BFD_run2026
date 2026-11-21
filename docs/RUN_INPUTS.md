# Run inputs

## `samples.csv` is a symlink

`samples.csv` here is a symlink to `../Fungi_BFD/samples.csv`. **Fungi_BFD owns the
file**; this repo only points at it.

It used to be a second tracked copy, and the two drifted: by 2026-11-20 the copy here
and the copy in Fungi_BFD disagreed by 684 genomes, and a third, newer version existed
only as an uncommitted file in the shared checkout. Nothing said which was right.

The file is generated, not edited. To change it, change the curation inputs in
Fungi_BFD and regenerate:

```bash
cd ../Fungi_BFD
python3 scripts/create_samples_file.py        # reads data/curation/*, writes samples.csv
```

`data/curation/README.md` in that repo says which file to use for which kind of fix
(drop a genome, prefer one assembly, correct a field, separate a haplotype pair).
Never edit `samples.csv` by hand: the next regeneration silently reverts it. If edits
were already made, `scripts/extract_hand_edits.py` recovers them as curation rows.

## This repo expects Fungi_BFD as a sibling

The symlink is relative, so a checkout needs the two repos side by side:

```
<parent>/
  Fungi_BFD/           <- owns samples.csv
  Fungi_BFD_runs/      <- samples.csv -> ../Fungi_BFD/samples.csv
```

which is how they sit in `/bigdata/stajichlab/shared/projects/BFD/`. Clone this repo
on its own and the symlink dangles; Nextflow then fails its input validation at
startup rather than running on a stale or missing file.

The Kubernetes path is unaffected: `nextflow/k8/bin/nf-job.sh` copies `samples.csv`
from the Fungi_BFD clone at `/workspace/repo`, which is the owning repo, and `cp`
follows the link.

## `suppress.txt` is a symlink too

`suppress.txt` here is a symlink to `../Fungi_BFD/data/curation/suppress.txt`,
for the same reason and under the same rule as `samples.csv`: **Fungi_BFD owns it**.

It had been a hand-maintained second copy. The two happened to be byte-identical on
2026-11-20, but only because someone had just made the same three edits twice; the
tracked copy here was three decisions behind its own working file. Add a suppression
in `Fungi_BFD/data/curation/suppress.txt` and it applies here at once.

The per-run launch directories under this repo already pointed at the top-level file
(`do_annotation_wave3/suppress.txt -> ../suppress.txt` and about two dozen more), so
they now resolve through it in two hops to the Fungi_BFD copy. A run that needs its
own list keeps a real file instead and names it in its params, as
`do_annotation_wave3` does with `suppress_wave3.txt`.

The Kubernetes path is unaffected: `profile_interproscan6_k8s.config` sets an absolute
`/workspace/suppress.txt`, and `interproscan6_k8s.nf` guards it with `.exists()`.
