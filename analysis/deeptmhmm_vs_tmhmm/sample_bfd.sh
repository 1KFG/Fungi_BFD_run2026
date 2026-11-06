#!/bin/bash
# sample_bfd.sh -- draw the BFD protein sample for the DeepTMHMM vs TMHMM benchmark.
# Reads tables/{tmhmm,signalp.signal_peptide,predgpi,species}.parquet; writes
# bfd_sample/sample.tsv (protein_id, LOCUSTAG, stratum, PHYLUM, len, PredHel, ...)
# and bfd_sample/bench.fa (sequences from genome_annotation/<out>/predict_results).
#
# Strata (deterministic: ordered by hash(protein_id || seed)):
#   tm_pos        2000  TMHMM PredHel >= 1, spread evenly over phyla
#   tm_neg        2000  TMHMM PredHel = 0, spread evenly over phyla
#   nterm_single   600  PredHel = 1, helix starts within the first 60 aa
#   sp_overlap     600  PredHel = 1 within the first 60 aa AND a SignalP 6 signal peptide
#   cterm_single   600  PredHel = 1, helix ends within 40 aa of the C-terminus
#   gpi            600  PredGPI GPI-anchor prediction
#   multi7         600  PredHel >= 7
#   sp_only        600  SignalP 6 signal peptide, PredHel = 0
# A protein can fall in several strata; it is listed once per stratum.
# Draws are 2x the target size: after removing proteins with >=30% identity to
# the DeepTMHMM training set, the first 1000 (or 300) held-out proteins per
# stratum, in hash order, form the held-out sample.
set -euo pipefail
A=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/deeptmhmm_vs_tmhmm
T=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
SEED=20261106
mkdir -p $A/bfd_sample; cd $A/bfd_sample
module load duckdb 2>/dev/null || true
duckdb -c "
SET threads = ${SLURM_CPUS_PER_TASK:-8};
CREATE TEMP TABLE sp AS SELECT DISTINCT protein_id FROM '$T/signalp.signal_peptide.parquet';
CREATE TEMP TABLE gpi AS SELECT DISTINCT protein_id FROM '$T/predgpi.parquet' WHERE feature = 'GPI-anchor';  -- every protein also has a 'Chain' row
CREATE TEMP TABLE tm AS
  SELECT t.*, s.PHYLUM,
         TRY_CAST(regexp_extract(t.Topology, '^[io](\d+)-', 1) AS INT)   AS first_start,
         TRY_CAST(regexp_extract(t.Topology, '-(\d+)[io]\$', 1) AS INT)  AS last_end,
         hash(t.protein_id || '$SEED') AS h
  FROM '$T/tmhmm.parquet' t JOIN '$T/species.parquet' s ON s.LOCUSTAG = t.species_prefix;
CREATE TEMP TABLE pick AS
  WITH pos AS (SELECT *, row_number() OVER (PARTITION BY PHYLUM ORDER BY h) rk FROM tm WHERE PredHel >= 1),
       neg AS (SELECT *, row_number() OVER (PARTITION BY PHYLUM ORDER BY h) rk FROM tm WHERE PredHel = 0)
  (SELECT 'tm_pos' stratum, * EXCLUDE (rk) FROM pos ORDER BY rk, h LIMIT 2000)
  UNION ALL (SELECT 'tm_neg', * EXCLUDE (rk) FROM neg ORDER BY rk, h LIMIT 2000)
  UNION ALL (SELECT 'nterm_single', * FROM tm WHERE PredHel = 1 AND first_start <= 60 ORDER BY h LIMIT 600)
  UNION ALL (SELECT 'sp_overlap', * FROM tm WHERE PredHel = 1 AND first_start <= 60 AND protein_id IN (SELECT protein_id FROM sp) ORDER BY h LIMIT 600)
  UNION ALL (SELECT 'cterm_single', * FROM tm WHERE PredHel = 1 AND last_end >= len - 40 ORDER BY h LIMIT 600)
  UNION ALL (SELECT 'gpi', * FROM tm WHERE protein_id IN (SELECT protein_id FROM gpi) ORDER BY h LIMIT 600)
  UNION ALL (SELECT 'multi7', * FROM tm WHERE PredHel >= 7 ORDER BY h LIMIT 600)
  UNION ALL (SELECT 'sp_only', * FROM tm WHERE PredHel = 0 AND protein_id IN (SELECT protein_id FROM sp) ORDER BY h LIMIT 600);
COPY (SELECT stratum, protein_id, species_prefix AS LOCUSTAG, PHYLUM, len, PredHel, Topology,
             first_start, last_end, protein_id IN (SELECT protein_id FROM sp) AS signalp_sp,
             protein_id IN (SELECT protein_id FROM gpi) AS predgpi
      FROM pick ORDER BY stratum, protein_id) TO 'sample.tsv' (HEADER, DELIMITER '\t');
SELECT stratum, count(*) n, count(DISTINCT PHYLUM) phyla, count(DISTINCT LOCUSTAG) genomes FROM read_csv('sample.tsv') GROUP BY 1 ORDER BY 1;
"
# sequences: one pass per genome over its proteome
/usr/bin/python3.12 - <<'EOF'
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
    for line in gzip.open(f"{R}/genome_annotation/{o}/predict_results/{o}.proteins.fa.gz", "rt"):
        if line.startswith(">"):
            if cur in ids: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if cur in ids: seqs[cur] = "".join(buf)
missing = sum(len(v) for v in want.values()) - len(seqs)
with open("bench.fa", "w") as fh:
    for pid in sorted(seqs):
        fh.write(f">{pid}\n{seqs[pid].rstrip('*')}\n")
print(f"bench.fa: {len(seqs)} unique proteins from {len(want)} genomes; {missing} not found")
EOF
