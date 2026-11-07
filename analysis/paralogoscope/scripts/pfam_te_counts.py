#!/usr/bin/python3.12
"""Per genome (locus tag prefix of protein_id): proteins with any Pfam hit, and proteins with a transposable-element Pfam domain.
Reads tables/pfam.parquet with DuckDB (the table is 20 GB; protein_id = <LOCUSTAG>_<number>-T1)."""
import duckdb, sys
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
TE = ["PF00078", "PF07727", "PF13456", "PF00665", "PF00872", "PF20700", "PF01609", "PF13359", "PF13843", "PF03184", "PF13358", "PF02992", "PF04754", "PF10551", "PF14372"]
con = duckdb.connect(); con.execute("SET memory_limit='12GB'; SET threads=4;")
acc = ",".join(f"'{a}'" for a in TE)
q = f"""
COPY (
  SELECT split_part(protein_id,'_',1) AS species_prefix,
         count(DISTINCT protein_id) AS proteins_with_pfam,
         count(DISTINCT CASE WHEN split_part(pfam_acc,'.',1) IN ({acc}) THEN protein_id END) AS te_domain_proteins
  FROM read_parquet('{R}/tables/pfam.parquet')
  GROUP BY 1
) TO '{sys.argv[1]}' (FORMAT CSV, DELIMITER '\t', HEADER)
"""
con.execute(q)
print("done")
