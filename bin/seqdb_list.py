#!/usr/bin/env python3
"""seqdb_list.py -- genome list for the indexed sequence copies (seqdb/).

Writes one TSV line per genome in samples.csv that is not in suppress.txt and
has a predicted proteome: LOCUSTAG, out, protein FASTA, CDS FASTA.
usage: seqdb_list.py RUNS_DIR > seqdb/genomes.tsv
"""
import csv, os, sys
R = sys.argv[1] if len(sys.argv) > 1 else "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin")
from pick_representative_strain import make_sample_tag
sup = {l.split(",")[0].strip() for l in open(f"{R}/suppress.txt") if l.strip() and not l.startswith("#")}
seen = set()
for r in csv.DictReader(open(f"{R}/samples.csv")):
    tag = (r.get("LOCUSTAG") or "").strip()
    if not tag or r["ASMID"].strip() in sup or tag in seen:
        continue
    out = make_sample_tag(r["SPECIES"], r["STRAIN"])
    d = f"{R}/genome_annotation/{out}/predict_results"
    pep, cds = f"{d}/{out}.proteins.fa.gz", f"{d}/{out}.cds-transcripts.fa.gz"
    if os.path.exists(pep):
        seen.add(tag)
        print(f"{tag}\t{out}\t{pep}\t{cds if os.path.exists(cds) else ''}")
