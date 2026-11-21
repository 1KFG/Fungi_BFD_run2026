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

## `suppress.txt`

Still a hand-maintained second copy of `Fungi_BFD/data/curation/suppress.txt`
(byte-identical as of 2026-11-20). Same problem, not yet resolved.
