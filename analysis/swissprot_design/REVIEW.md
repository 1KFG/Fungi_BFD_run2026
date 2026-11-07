# Review: design of `tables/swissprot.parquet`

Date: 2026-11-06. Review only. No files under `tables/` or the pipeline were changed.

## 1. Summary of the recommendation

1. Normalize the subject. Add `sp_id` (INTEGER) to `swissprot_annot` and make it the subject dimension table. The hit table references `sp_id` only.
2. Keep all DIAMOND hits (up to 20 per query) in the hit table, but with 15 narrow columns, ZSTD, sorted by `protein_id, rk`. Estimated size: 36-41 GB (from 146 GB).
3. Add a best-hit table (one row per protein, estimated 3.5 GB) and SQL views for the 80/80 set and the model-fungus hits. Do not store `func_transfer_80_80`.
4. Build the table per genome bucket (256 hive partitions, matching the `storeDir` layout) with DuckDB reading the `.blasttab.gz` files directly. No intermediate CSV.
5. Fix a definition bug while rebuilding: `query_cov` and `hit_cov` are `aln_length / len`, which counts gap columns. They exceed 1.0 in 9.7% and 16.6% of rows. Use `qcovhsp` from DIAMOND and `(s_end - s_start + 1) / s_len`.
6. Load the table into `BFD.duckdb` as a VIEW over the Parquet files. It is not loaded at all today.

If 40 GB is too much, a reduced hit table (top 5 by bitscore, plus all 80/80 hits, plus the best hit per model fungus) is estimated at 15 GB. It keeps 35% of the rows. Section 4 lists what it loses.

## 2. What was measured

Sources:

- Full-table Parquet metadata (`parquet_metadata()` on the 146 GB file, 16 s). No full scan was run.
- A sample of 40 genomes: 10 random buckets x 4 random `.blasttab.gz` files from `BFD_function/results/function/swissprot/`. File list: `sample_files_40genomes.txt` in this folder. The sample was loaded into DuckDB in the exact 23-column shape that `merge_swissprot.py` produces.
- Per-genome raw files: `du` over the whole tree.
- Point lookups on the full file.
- DuckDB 1.1.3 (the only module version). PyArrow 21 for one encoding test.

| Measurement | Value |
|---|---|
| Full table rows | 2,192,440,076 |
| Full table size | 146,373,511,272 bytes = 66.8 bytes/row |
| Row groups | 17,733; 13,100-124,912 rows each (DuckDB default, ~122,880) |
| Codec / encoding | SNAPPY, PLAIN for every numeric column; dictionary only on `protein_id`, `swissprot_name`, `species_prefix` |
| Sample rows / proteins / genomes | 3,974,621 / 285,227 / 40 |
| Sample written in the current form (SNAPPY) | 265.96 MB = 66.9 bytes/row. This matches the full table (66.8), so sizes below scale by 551.6x |
| Hits per protein with >= 1 hit | 13.93 |
| Proteins at the k=20 cap | 153,968 / 285,227 = 54.0% |
| Capped proteins whose 20th hit has bitscore >= 50% of the best | 74,632 (48.5% of capped) |
| Capped proteins whose 20th hit has bitscore >= 90% of the best | 5,045 (3.3% of capped) |
| Duplicate (protein, accession) pairs | 0 |
| Distinct accessions in sample | 138,312 (of 575,503 in SwissProt 2026_02) |
| Accessions with >= 100 rows in sample | 5,695 (4.1%); they hold 26.8% of rows. Top entries are MFS transporters (`YFZ4_SCHPO` 2,162 rows = 54 per genome) |
| Distinct subject organism codes (entry suffix) | 4,897; mean 11.0 distinct organisms among a protein's hits |
| Rows with `func_transfer_80_80 = 1` | 59,911 (1.5%); 22,441 proteins (7.9%) |
| Proteins with >= 1 hit from YEAST, SCHPO, NEUCR or EMENI | 225,862 (79.2%); best hit is from one of these in 128,379 (45.0%) |
| Per organism: any YEAST / SCHPO / NEUCR / EMENI hit | 56.5% / 57.6% / 16.2% / 19.7% of proteins |
| `evalue = 0` (underflow) | 100,371 rows (2.5%) |
| `pident`, `qcovhsp` precision | 1 decimal (831 and 988 distinct values) |
| `bitscore` precision | 1 decimal below 100 (987,623 rows), integer above |
| `query_cov > 1` / `hit_cov > 1` | 385,849 (9.7%) / 660,626 (16.6%) |
| `swissprot_name` != `swissprot_annot.protein_name` | 92,722 rows (2.3%); annot concatenates "Includes" names, stitle gives only the first |
| `swissprot_acc` -> `swissprot_annot.accession` join | 3,974,621 / 3,974,621 (100%) |
| Per-genome `.blasttab.gz` files | 23,410 files, 125 GB (apparent size) |
| Point lookup, full file, by `protein_id` | 3.1 s wall (1.6 s query, rest is reading 17,733 row-group footers) |
| Point lookup, full file, by `species_prefix` | 3.1 s wall, 126,362 rows returned |
| `tables/gene_proteins.parquet` | Absent from `tables/`. Only a 3,786-row copy in `_superseded_20261105/` |
| `swissprot` in `build_BFD_duckDB.sh` / `sql/table_schema.json` | Not referenced. The table is not in `BFD.duckdb` |

Estimates are derived from the sample with the 551.6x factor. They are marked "est." below. The 66.9 vs 66.8 bytes/row match gives confidence for whole-table size, but the sample covers 40 of 23,410 genomes, so per-genome variation is not captured.

## 3. Question 1: normalization and derivable columns

### Where the bytes go now

Per-column compressed size of the full table (measured):

| Column | GB | Note |
|---|---|---|
| evalue | 17.2 | DOUBLE, PLAIN; 2 significant digits of information |
| swissprot_entry | 14.5 | string repeated per hit |
| swissprot_acc | 14.0 | string repeated per hit |
| swissprot_name | 10.4 | string repeated per hit |
| hit_cov, query_cov | 9.3 + 8.6 | DOUBLE; derivable and wrongly defined |
| pident, qcovhsp, bitscore | 6.5 + 6.4 + 6.2 | DOUBLE for 1-decimal values |
| 11 INT64 columns | 51.4 | 8 bytes each for values < 65,536 |
| q_len | 1.6 | |
| protein_id | 1.25 | dictionary, sorted |
| func_transfer_80_80 | 0.9 | INT64 for a flag |
| species_prefix | 0.0 | dictionary, run-length |

Three subject strings cost 38.9 GB (27% of the file). All numeric columns are PLAIN INT64/DOUBLE under SNAPPY, which compresses them 2-3x at best.

### Recommendation

Yes, normalize. Make the subject a dimension table and reference it with an integer.

**`swissprot_subject`** (replaces `swissprot_annot`; one row per SwissProt entry; 575,503 rows; measured 23.8 MB as ZSTD Parquet with all annotation text):

| Column | Type | Source |
|---|---|---|
| sp_id | INTEGER | dense id assigned at build: `row_number() OVER (ORDER BY accession)`; stable for one SwissProt release; rebuild hits when the release changes |
| accession | VARCHAR | `uniprot_sprot.dat` |
| entry_name | VARCHAR | |
| org_code | VARCHAR | suffix of entry_name (`YEAST`, `SCHPO`, ...); cheap organism filter |
| organism | VARCHAR | |
| s_len | INTEGER | sequence length from the `.dat` SQ line or the FASTA; `parse_uniprot_sprot.py` does not emit it yet |
| protein_name, gene_name, gene_synonyms, go_terms, ec_numbers, interpro, pfam | VARCHAR | as today |

**`swissprot_hits`** (fact table):

| Column | Type | Note |
|---|---|---|
| protein_id | VARCHAR | sort key 1 |
| species_prefix | VARCHAR | derivable, but costs ~0 bytes (dictionary + RLE) and matches the other BFD tables; keep |
| sp_id | INTEGER | FK to swissprot_subject |
| rk | UTINYINT | rank within protein by bitscore desc, evalue asc, accession; sort key 2 |
| pident_x10 | USMALLINT | pident * 10, lossless |
| qcov_x10 | USMALLINT | DIAMOND qcovhsp * 10, lossless |
| scov_x10 | USMALLINT | round((s_end - s_start + 1) * 1000 / s_len) |
| aln_length | INTEGER | |
| q_start, q_end, s_start, s_end | INTEGER | |
| q_len | INTEGER | keep until `gene_proteins` exists in `tables/`; est. 0.6 GB |
| bitscore | REAL | 1-decimal values are exact in float32 |
| neg_log10_evalue | REAL | `-log10(evalue)`; store 400.0 when evalue = 0 (underflow; 2.5% of rows) so that `neg_log10_evalue >= 50` works without NULL handling |

Dropped, and why:

| Column | Replacement | Est. saving |
|---|---|---|
| swissprot_acc, swissprot_entry, swissprot_name | join `swissprot_subject` on sp_id | 38.9 GB measured -> 5.8 GB est. for sp_id |
| s_len | swissprot_subject.s_len | 5.7 GB measured -> 0 (3.1 GB est. if kept as INT32 ZSTD) |
| query_cov, hit_cov | qcov_x10, scov_x10 (correct definition) | 17.9 GB measured -> ~6 GB est. for the two USMALLINT columns |
| func_transfer_80_80 | view: `pident_x10 >= 800 AND qcov_x10 >= 800 AND scov_x10 >= 800` | 0.9 GB |
| evalue (DOUBLE) | neg_log10_evalue (REAL) | 17.2 GB measured -> 6.9 GB est. |
| nident, positive, mismatch, gapopen | per-genome raw files | 19.8 GB measured -> 10.6 GB est. if kept; 0 if dropped. They are not used by any of the listed uses |

Note on `float32` for evalue: raw e-values reach 1e-300, which float32 cannot hold. `-log10` as REAL is the right form. Scaled integers (`USMALLINT`) were measured at the same size as REAL for pident/coverage; they were chosen because they are exact and index well.

Note on `query_cov` / `hit_cov`: `merge_swissprot.py` computes `aln_length / q_len`. `aln_length` is the alignment length including gap columns, so coverage goes above 1.0 (9.7% and 16.6% of rows). DIAMOND's `qcovhsp` is `(q_end - q_start + 1) / q_len`. In the sample, 3,111,331 of 3,974,621 rows (78%) differ by more than 0.2 percentage points between the two. The 80/80 flag changes for only 68 of 59,911 flagged rows, so the flag is mostly fine, but the stored coverage values are wrong and should not be kept.

## 4. Question 2: how many hits to keep, and for which use

Row counts from the sample under each rule (all rules keep every protein that has at least one hit):

| Rule | Rows | % of all | Est. full rows | Est. size, 13-col ZSTD |
|---|---|---|---|---|
| all hits (k=20) | 3,974,621 | 100 | 2.19 B | 36-41 GB |
| best 1 | 285,227 | 7.2 | 157 M | 3.5 GB |
| top 3 | 791,822 | 19.9 | 437 M | |
| top 5 | 1,253,769 | 31.5 | 692 M | |
| top 10 | 2,284,143 | 57.5 | 1.26 B | |
| bitscore >= 95% of best | 595,034 | 15.0 | | |
| bitscore >= 90% of best | 883,837 | 22.2 | | |
| bitscore >= 80% of best | 1,430,360 | 36.0 | | |
| bitscore >= 50% of best | 2,794,870 | 70.3 | | |
| 80/80 only | 59,911 | 1.5 | 33 M | |
| best per subject organism | 3,138,210 | 79.0 | | not a reduction |
| top 5 OR 80/80 OR best per model fungus (YEAST, SCHPO, NEUCR, EMENI) | 1,402,826 | 35.3 | 774 M | 15.3 GB |

### Is k=20 too many?

For naming: yes. The best 1-5 hits decide a name. For surveys: k=20 is already too few. 54% of proteins hit the cap, and for 48.5% of those the 20th hit is still within 50% of the best bitscore. So the full table is already a truncated view of "all SwissProt homologs" for half the proteins (tubulins, MFS transporters, kinases, ...). No value of N from 1 to 20 makes the table a complete family survey. Family surveys should come from Pfam or orthology, not from this table.

### Uses and which representation serves each

| Use | Needs | Served by |
|---|---|---|
| Functional annotation transfer / naming | best hit, plus the top few to detect ties and chimeras | best-hit table; `rk <= 5` from the hit table |
| GO / EC transfer with a stringency filter | best hit or 80/80 hits, joined to subject GO/EC | best-hit table and `v_swissprot_80_80` joined to `swissprot_subject` |
| Characterized homologs in model fungi | best hit per protein per org_code | `v_swissprot_model_fungi` over the hit table (all hits needed: the model-fungus hit is not among the top 5 for many proteins; 79% of proteins have one somewhere in the 20) |
| Orthology hints (reciprocal-best style) | best hit per protein and per subject | best-hit table, plus a scan by sp_id |
| Gene-family surveys by SwissProt subject | all hits to one accession across genomes | hit table scanned by sp_id; incomplete beyond rank 20 for 54% of proteins (see above) |
| Raw evidence, re-scoring, alignment coordinates with mismatch/gap counts | everything DIAMOND wrote | per-genome `.blasttab.gz` (125 GB, kept) |

### Recommendation

Keep all hits in the database hit table, in the narrow form (est. 36-41 GB). The reduced table saves an estimated 25 GB more but drops 65% of rows, and the dropped rows are exactly the ones the model-fungus and survey uses need. The raw `.blasttab.gz` files are not queryable without a re-merge, so "full hits in per-genome files" is an archive, not a substitute.

Add the best-hit table as a separate, materialized table (est. 3.5 GB; 157 M rows) with the subject columns joined in (accession, entry_name, protein_name, org_code), because it is the table most ad-hoc queries will touch. The views below cost nothing:

```sql
CREATE VIEW v_swissprot_80_80 AS
SELECT * FROM swissprot_hits
WHERE pident_x10 >= 800 AND qcov_x10 >= 800 AND scov_x10 >= 800;

CREATE VIEW v_swissprot_model_fungi AS
SELECT h.*, s.org_code, s.accession, s.entry_name, s.protein_name
FROM swissprot_hits h JOIN swissprot_subject s USING (sp_id)
WHERE s.org_code IN ('YEAST','SCHPO','NEUCR','EMENI')
QUALIFY row_number() OVER (PARTITION BY h.protein_id, s.org_code ORDER BY h.rk) = 1;
```

If storage must be cut further, use the "top 5 OR 80/80 OR best per model fungus" rule (est. 15 GB). What is lost: hits of rank 6-20 that are not 80/80 and not from the four model fungi. That removes the ability to answer "which BFD proteins hit accession X" for subjects that are typically a weaker hit, and removes homologs in non-model but characterized fungi (ASPFU, CANAL, MAGO7, FUSO4, USTMA, ...). The 125 GB raw archive still holds them.

Option not pursued: a second DIAMOND run against only the four model proteomes (about 15,000 sequences) would give a model-fungus best hit for proteins whose model homolog is outside the top 20. Not measured.

## 5. Question 3: physical layout

### Codec and encoding (measured on the sample, 3,974,621 rows)

| Variant | Sample MB | Est. full GB |
|---|---|---|
| Current 23 columns, SNAPPY (as built) | 265.96 | 146.7 (actual 146.4) |
| Current 23 columns, ZSTD | 152.64 | 84.2 |
| Lean 21 columns (both evalue forms, all counts), ZSTD, sorted by protein_id, rk | 119.91 | 66.1 |
| Same, sorted by sp_id, protein_id | 131.13 | 72.3 |
| Same, ZSTD, ROW_GROUP_SIZE 1,000,000 | 117.07 | 64.6 |
| Same, SNAPPY | 186.56 | 102.9 |
| Minimal 13 columns (protein_id, sp_id, rk, pident_x10, qcov_x10, scov_x10, aln_length, q/s start/end, bitscore, neg_log10_evalue), ZSTD, DuckDB 1.1.3 | 72.39 | 39.9 |
| Same 13 columns, PyArrow 21, ZSTD level 3, default dictionary encoding | 65.30 | 36.0 |
| Same, PyArrow, DELTA_BINARY_PACKED ints + BYTE_STREAM_SPLIT floats, ZSTD 3 | 69.68 | 38.4 |
| Same, ZSTD level 9 | 67.95 | 37.5 |
| Reduced rule (top 5 OR 80/80 OR model fungus), 13 columns, ZSTD | 27.68 | 15.3 |
| Best hit only, 13 columns, ZSTD | 6.31 | 3.5 |
| `swissprot_subject` with all annotation text, ZSTD | 23.82 | 0.024 (full, measured) |

Findings:

- ZSTD alone, with no schema change, cuts the file to 84 GB (est.). The gain is 42%.
- Narrow types and normalization matter more than the codec: 146 -> 40 GB.
- Parquet encodings (delta, byte-stream-split, ZSTD level 9) change the size by under 10% in either direction. Not worth custom writer code. DuckDB's default writer is fine; the 15-column table with `species_prefix` and `q_len` added to the 13-column variant adds an estimated 0.6 GB.
- Sort order: by `protein_id` gives 14.7 MB -> 1.2 MB on `protein_id` and the per-genome row-group pruning measured above (3 s lookups). Sorting by subject makes the file 9% larger and loses the protein/genome pruning. Keep `protein_id, rk`.
- Row group size: 1,000,000 rows is 2.4% smaller than the default and cuts footer entries from ~17,700 to ~2,200 (faster `parquet_metadata`, faster planning). Pruning granularity drops from ~9 to ~75 proteins per row group, which does not matter for point lookups. Recommend 500,000-1,000,000.

### Partitioning

Recommend hive partitioning by the existing 256-way genome bucket (`bucket=XX/`), the same bucket `storeDir` already uses for the raw files, written as `tables/swissprot_hits/bucket=XX/data.parquet`, each file sorted by `protein_id, rk`.

Why:

- The merge becomes 256 independent tasks (est. 8.6 M rows and ~150 MB each) with no intermediate CSV. One bucket can be rebuilt when genomes in it change.
- A genome query touches one file; a protein query touches one file and one or two row groups.
- DuckDB reads the directory as one table (`read_parquet('.../swissprot_hits/*/*.parquet', hive_partitioning=true)`), and `generate_bfd_catalog_sql.py` already supports a glob source.

Partition by phylum was considered and is not recommended. It requires the species join at merge time, changes when taxonomy is corrected, and only helps queries that filter on phylum, which can instead join `species` and prune on `species_prefix` zone maps (the 3 s lookups show this works). If a single file is preferred for distribution, concatenate the bucket files in one final `COPY`; the sort inside each file is preserved and the file-level zone maps still prune.

### Effect on query patterns

| Pattern | Now (146 GB file) | Proposed | Basis |
|---|---|---|---|
| By protein_id | 3.1 s (measured) | similar; fewer footers to read | row-group zone maps on sorted key |
| By genome (species_prefix) | 3.1 s (measured) | one bucket file | sorted; hive partition |
| By accession / subject | scan 14 GB `swissprot_acc` column | scan est. 5.8 GB `sp_id` column over all files; no zone-map help because sp_id is unordered | estimate, not timed |
| By phylum | join species, prune on species_prefix | same | |
| Best hit for every protein | 2.19 B-row window function | read the 3.5 GB best-hit table | |

For subject-first queries, two low-cost additions help: a `swissprot_subject_summary` table (sp_id, n_hits, n_proteins, n_genomes, n_best_hits, max_bitscore; under 600k rows) built in the same merge for survey-style questions, and the best-hit table for reciprocal lookups. A second full copy sorted by `sp_id` (est. 40 GB) was measured for size but is not recommended; the summary table covers most of the need. Newer DuckDB versions write Parquet bloom filters, which would speed `sp_id = ?` lookups; the installed 1.1.3 does not, and this was not tested.

### DuckDB load

`build_BFD_duckDB.sh` does not load `swissprot.parquet` today. Load the new tables as VIEWs over Parquet (the pattern used for species / asm_stats / busco_genome in `sql/table_schema.json`), not as `CREATE TABLE AS` copies. A materialized 2.19 B-row copy with ART indexes inside `BFD.duckdb` was not measured; it would at least duplicate the Parquet volume inside the DB file. Materialize only the best-hit table and `swissprot_subject` if indexed lookups are wanted.

## 6. How the merge would build it

1. `BUILD_SWISSPROT_ANNOT` (per SwissProt release): emit `sp_id`, `org_code`, `s_len` in addition to the current columns. Write `swissprot_subject.parquet` (ZSTD).
2. `MERGE_SWISSPROT` per bucket (256 tasks, label `merge`, one bucket's `*.blasttab.gz` as input; DuckDB reads gzip natively):

```sql
COPY (
  WITH raw AS (
    SELECT * FROM read_csv('<bucket>/*.blasttab.gz', delim='\t', header=false,
      columns={'qseqid':'VARCHAR','sseqid':'VARCHAR','pident':'DOUBLE','positive':'INTEGER',
               'nident':'INTEGER','length':'INTEGER','mismatch':'INTEGER','gapopen':'INTEGER',
               'qstart':'INTEGER','qend':'INTEGER','sstart':'INTEGER','send':'INTEGER',
               'evalue':'DOUBLE','bitscore':'DOUBLE','qcovhsp':'DOUBLE','qlen':'INTEGER',
               'slen':'INTEGER','stitle':'VARCHAR'})
  ),
  j AS (
    SELECT r.qseqid AS protein_id, split_part(r.qseqid,'_',1) AS species_prefix,
           s.sp_id, r.*
    FROM raw r JOIN read_parquet('swissprot_subject.parquet') s
      ON s.accession = split_part(r.sseqid,'|',2)
  )
  SELECT protein_id, species_prefix, sp_id,
         row_number() OVER (PARTITION BY protein_id ORDER BY bitscore DESC, evalue, sp_id)::UTINYINT AS rk,
         round(pident*10)::USMALLINT AS pident_x10,
         round(qcovhsp*10)::USMALLINT AS qcov_x10,
         round((send-sstart+1)*1000.0/slen)::USMALLINT AS scov_x10,
         length AS aln_length, qstart AS q_start, qend AS q_end, sstart AS s_start, send AS s_end,
         qlen AS q_len, bitscore::REAL AS bitscore,
         CASE WHEN evalue = 0 THEN 400.0 ELSE -log10(evalue) END::REAL AS neg_log10_evalue
  FROM j ORDER BY protein_id, rk
) TO 'bucket=<XX>/data.parquet' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000);
```

   Fail the task if any `sseqid` does not join (today the join is 100%, measured on the sample). The coverage check (`coverageCheck`) runs per bucket.
3. Gather step: `swissprot_best.parquet` = `WHERE rk = 1` over all buckets, joined to `swissprot_subject` for accession, entry_name, protein_name, org_code; `swissprot_subject_summary.parquet` = `GROUP BY sp_id`. Both are small.
4. `build_BFD_duckDB.sh`: VIEWs `swissprot_hits` (glob), `swissprot_subject`, `swissprot_best`, `swissprot_subject_summary`, plus `v_swissprot_80_80` and `v_swissprot_model_fungi` above.

I did not time this build. The 9 h 40 min of the current merge is a single Python CSV writer over 2.19 B rows followed by one DuckDB CSV read; the per-bucket plan replaces both and runs 256-wide, but its runtime is an expectation, not a measurement.

## 7. What is lost, in one list

- `nident`, `positive`, `mismatch`, `gapopen`: only in the raw files. Est. 10.6 GB to keep them as INT32 ZSTD if wanted.
- `evalue` below 1e-300 collapses to the 400.0 sentinel (2.5% of rows); `bitscore` still ranks them.
- Subject strings need a join. The best-hit table carries them, so most queries do not need it.
- Reduced-table option only: ranks 6-20 outside 80/80 and outside the four model fungi.
- The `swissprot_name` text from stitle (first name only) is replaced by `protein_name` from the flatfile (all "Includes" names, 2.3% differ). This is a gain, not a loss, but the strings change.

## 8. Files

- This review: `analysis/swissprot_design/REVIEW.md`
- Sample genomes: `analysis/swissprot_design/sample_files_40genomes.txt`
- PyArrow encoding test: `analysis/swissprot_design/pyarrow_encoding_test.py`
- Scratch (temporary, not kept): `/rhome/jstajich/.claude/jobs/bacad32f/tmp/` (sample DuckDB and the `v*.parquet` variants)
