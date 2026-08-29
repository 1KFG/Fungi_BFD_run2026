#!/usr/bin/env python3
"""Screen samples.csv genomes with no GENUS-level taxonomy against the
existing results/mash_sketch/ reference cache via mash, to surface
candidate nearby lineages for manual review.

See run_query_unknowns_mash.sh for the SLURM entry point.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

TAXON_FIELDS = ["PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"]
SKETCH_SUFFIX = ".fa.msh"
GENOME_SUFFIX = ".fa.gz"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", default="samples.csv")
    p.add_argument("--sketch-dir", default="results/mash_sketch")
    p.add_argument("--genome-dir", default="input_clean_genomes")
    p.add_argument("--outdir", default="results/query_unknowns_mash")
    p.add_argument("--kmer", type=int, default=21)
    p.add_argument("--sketch-size", type=int, default=10000)
    p.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    p.add_argument("--top-n", type=int, default=10)
    p.add_argument("--pval", type=float, default=0.05)
    p.add_argument("--mash-bin", default="mash")
    p.add_argument("--scratch-dir", default=os.environ.get("SCRATCH"),
                    help="Node-local scratch dir for intermediate mash files "
                         "(defaults to $SCRATCH; falls back to a system temp dir)")
    return p.parse_args()


def load_samples(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def sketch_for(sketch_dir, asmid):
    return sketch_dir / f"{asmid}{SKETCH_SUFFIX}"


def build_sketch(mash_bin, kmer, size, genome_dir, sketch_dir, asmid):
    """Build a sketch for asmid, matching the naming/internal-ID convention
    of the existing results/mash_sketch cache: sketched from the bare
    filename (cwd=genome_dir) so the internal ID stays '<asmid>.fa.gz' with
    no directory prefix, output written as '<asmid>.fa.msh'."""
    genome_path = genome_dir / f"{asmid}{GENOME_SUFFIX}"
    if not genome_path.exists():
        return False, f"genome fasta not found: {genome_path}"
    out_prefix = (sketch_dir / f"{asmid}.fa").resolve()
    cmd = [mash_bin, "sketch", "-k", str(kmer), "-s", str(size),
           "-o", str(out_prefix), genome_path.name]
    result = subprocess.run(cmd, cwd=genome_dir, capture_output=True, text=True)
    if result.returncode != 0:
        return False, result.stderr.strip()
    return True, None


def strip_suffix(name, suffix):
    return name[: -len(suffix)] if name.endswith(suffix) else name


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    sketch_dir = Path(args.sketch_dir)
    genome_dir = Path(args.genome_dir)

    rows = load_samples(args.samples)
    by_asmid = {r["ASMID"]: r for r in rows}
    queries = [r for r in rows if not r["GENUS"].strip()]
    print(f"[query_unknowns_mash] {len(queries)} genomes missing GENUS out of {len(rows)} total",
          file=sys.stderr)

    query_status = {}
    for r in queries:
        asmid = r["ASMID"]
        if sketch_for(sketch_dir, asmid).exists():
            query_status[asmid] = "existing"
            continue
        ok, err = build_sketch(args.mash_bin, args.kmer, args.sketch_size,
                                genome_dir, sketch_dir, asmid)
        if ok:
            query_status[asmid] = "built"
        else:
            query_status[asmid] = "MISSING_GENOME"
            print(f"[query_unknowns_mash] WARNING: could not sketch {asmid}: {err}",
                  file=sys.stderr)

    usable_queries = [r for r in queries if query_status[r["ASMID"]] in ("existing", "built")]
    print(f"[query_unknowns_mash] {len(usable_queries)} query genomes have usable sketches "
          f"({sum(1 for v in query_status.values() if v == 'built')} newly built)",
          file=sys.stderr)

    scratch_root = args.scratch_dir or tempfile.gettempdir()
    with tempfile.TemporaryDirectory(prefix="query_unknowns_mash.", dir=scratch_root) as tmp:
        tmp = Path(tmp)

        query_list_path = tmp / "query_sketches.txt"
        with open(query_list_path, "w") as fh:
            for r in usable_queries:
                fh.write(str(sketch_for(sketch_dir, r["ASMID"]).resolve()) + "\n")

        query_combined = tmp / "query_combined"
        subprocess.run([args.mash_bin, "paste", str(query_combined), "-l", str(query_list_path)],
                        check=True)
        query_combined_msh = tmp / "query_combined.msh"

        # All sketches currently in the cache are eligible references, including
        # any newly built above and any other GENUS-blank genomes already there
        # (an unfiltered reference set was the explicit design choice here).
        ref_list_path = tmp / "ref_sketches.txt"
        ref_paths = sorted(sketch_dir.glob(f"*{SKETCH_SUFFIX}"))
        with open(ref_list_path, "w") as fh:
            for p in ref_paths:
                fh.write(str(p.resolve()) + "\n")
        print(f"[query_unknowns_mash] comparing against {len(ref_paths)} reference sketches",
              file=sys.stderr)

        raw_dist_path = tmp / "raw_dist.tsv"
        with open(raw_dist_path, "w") as out_fh:
            subprocess.run(
                [args.mash_bin, "dist", "-p", str(args.threads), str(query_combined_msh),
                 "-l", str(ref_list_path)],
                stdout=out_fh, check=True,
            )

        hits_by_query = {}
        with open(raw_dist_path) as fh:
            for line in fh:
                query_id, ref_id, dist, pval, shared = line.rstrip("\n").split("\t")
                query_asmid = strip_suffix(query_id, GENOME_SUFFIX)
                ref_asmid = strip_suffix(ref_id, GENOME_SUFFIX)
                if query_asmid == ref_asmid:
                    continue  # exclude self-hit (query's own sketch is also in the ref list)
                hits_by_query.setdefault(query_asmid, []).append(
                    (float(dist), ref_asmid, float(pval), shared)
                )

    # --- write candidate_lineages.tsv (top-N per query) ---
    cand_path = outdir / "candidate_lineages.tsv"
    cand_header = [
        "query_asmid", "query_species_in", "rank", "ref_asmid",
        "ref_species", "ref_genus", "ref_family", "ref_order",
        "ref_class", "ref_subclass", "ref_phylum", "ref_subphylum",
        "ani_pct", "mash_pval", "shared_hashes", "significant",
    ]
    query_summary = []
    with open(cand_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cand_header)
        for r in usable_queries:
            asmid = r["ASMID"]
            hits = sorted(hits_by_query.get(asmid, []), key=lambda h: h[0])[: args.top_n]
            n_sig = 0
            best_ani = None
            best_ref = None
            best_genus = None
            for rank, (dist, ref_asmid, pval, shared) in enumerate(hits, start=1):
                ref_row = by_asmid.get(ref_asmid, {})
                ani_pct = (1 - dist) * 100
                significant = pval < args.pval
                if significant:
                    n_sig += 1
                if best_ani is None:
                    best_ani, best_ref = ani_pct, ref_asmid
                    best_genus = ref_row.get("GENUS", "")
                w.writerow([
                    asmid, r.get("SPECIES_IN", ""), rank, ref_asmid,
                    ref_row.get("SPECIES", ""), ref_row.get("GENUS", ""),
                    ref_row.get("FAMILY", ""), ref_row.get("ORDER", ""),
                    ref_row.get("CLASS", ""), ref_row.get("SUBCLASS", ""),
                    ref_row.get("PHYLUM", ""), ref_row.get("SUBPHYLUM", ""),
                    f"{ani_pct:.4f}", f"{pval:.6g}", shared, significant,
                ])
            query_summary.append((r, len(hits), n_sig, best_ani, best_ref, best_genus))

    # --- write query_list.tsv (per-query status/summary) ---
    ql_path = outdir / "query_list.tsv"
    with open(ql_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "query_asmid", "query_species_in", "sketch_status", "n_hits_reported",
            "n_significant_hits", "best_ani_pct", "best_ref_asmid", "best_ref_genus",
            "needs_review",
        ])
        for r in queries:
            asmid = r["ASMID"]
            status = query_status[asmid]
            if status == "MISSING_GENOME":
                w.writerow([asmid, r.get("SPECIES_IN", ""), status, 0, 0, "", "", "", True])
                continue
        for r, n_hits, n_sig, best_ani, best_ref, best_genus in query_summary:
            asmid = r["ASMID"]
            w.writerow([
                asmid, r.get("SPECIES_IN", ""), query_status[asmid], n_hits, n_sig,
                f"{best_ani:.4f}" if best_ani is not None else "",
                best_ref or "", best_genus or "",
                n_sig == 0,
            ])

    print(f"[query_unknowns_mash] wrote {cand_path} and {ql_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
