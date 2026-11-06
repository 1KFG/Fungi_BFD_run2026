#!/bin/bash -l
#SBATCH -p batch -c 8 --mem 48G --time 2:00:00 -J dtm_gpi -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/sample_bfd_gpi.%j.log
# GPI stratum (redo): proteins with a PredGPI feature = 'GPI-anchor' (the first draw
# used every protein id in predgpi.parquet, which also lists a 'Chain' row for all
# proteins). 600 drawn by hash order; held-out filtering later keeps 300.
set -euo pipefail
module load duckdb
T=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/bfd_sample_gpi
duckdb -c "SET threads=8; COPY (SELECT 'gpi' AS stratum, g.protein_id, g.species_prefix AS LOCUSTAG, s.PHYLUM, NULL AS len, NULL AS PredHel
  FROM (SELECT DISTINCT protein_id, species_prefix FROM '$T/predgpi.parquet' WHERE feature = 'GPI-anchor') g
  JOIN '$T/species.parquet' s ON s.LOCUSTAG = g.species_prefix
  ORDER BY hash(g.protein_id || '20261106') LIMIT 600) TO 'sample.tsv' (HEADER, DELIMITER '\t');"
/usr/bin/python3.12 - <<'PY'
import csv, gzip, sys
sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin")
from pick_representative_strain import make_sample_tag
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
out = {r["LOCUSTAG"].strip(): make_sample_tag(r["SPECIES"], r["STRAIN"]) for r in csv.DictReader(open(f"{R}/samples.csv"))}
want = {}
for r in csv.DictReader(open("sample.tsv"), delimiter="\t"):
    want.setdefault(r["LOCUSTAG"], set()).add(r["protein_id"])
seqs = {}
for tag, ids in want.items():
    o = out[tag]; cur = None; buf = []
    for line in list(gzip.open(f"{R}/genome_annotation/{o}/predict_results/{o}.proteins.fa.gz", "rt")) + [">END"]:
        if line.startswith(">"):
            if cur in ids: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
with open("bench.fa", "w") as fh:
    for p in sorted(seqs): fh.write(f">{p}\n{seqs[p].rstrip('*')}\n")
print(len(seqs), "sequences from", len(want), "genomes")
PY
module load mmseqs2/17-b804f
mmseqs easy-search bench.fa /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/deeptmhmm_trainset/trainset.fa vs_trainset.m8 ${SCRATCH:?}/mm --threads 8 -e 1e-3 --format-output "query,target,pident,qcov,tcov,evalue" -v 1 >/dev/null
sbatch -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/bfd_sample_gpi/run_tools.%j.log /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/run_tools.sh /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm/bfd_sample_gpi
