#!/usr/bin/env python3
"""analyze.py -- compare DeepTMHMM 1.0, TMHMM 2.0 and SignalP 6 on the
SwissProt (experimental labels) and BFD (stratified sample) benchmark sets.

Inputs per set folder: bench.fa, tmhmm_short.tsv, signalp/prediction_results.txt,
deeptmhmm/predicted_topologies.3line, vs_trainset.m8 (MMseqs2 vs the DeepTMHMM
training set). SwissProt: labels.tsv and <set>.tsv (UniProt features).
BFD: bfd_sample/sample.tsv and bfd_sample_neg/sample.tsv.

Held out = no DeepTMHMM training protein with >= 30% identity over >= 50% of the
query. BFD strata: the first 1000 (tm_pos, tm_neg) or 300 held-out proteins per
stratum in sha1(protein_id + seed) order.

Writes results/*.tsv and prints the summary tables.
"""
import csv, hashlib, os, re
from collections import Counter, defaultdict

A = os.path.dirname(os.path.abspath(__file__))
OUT = f"{A}/results" + (f"_merge{os.environ['MERGE_GAP']}" if os.environ.get("MERGE_GAP", "0") != "0" else "")
os.makedirs(OUT, exist_ok=True)
SEED = "20261106"
MERGE_GAP = int(os.environ.get("MERGE_GAP", "0"))
OVERLAP = 5   # residues: a predicted helix matches a true helix when they share >= 5


def tmhmm(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        top = f[5].split("=", 1)[1]
        hel = [(int(a), int(b)) for a, b in re.findall(r"(\d+)-(\d+)", top)]
        d[f[0]] = {"n": int(f[4].split("=")[1]), "hel": hel, "nterm": top[0] if hel else None}
    return d


def deeptmhmm(path):
    L = open(path).read().split("\n")
    d = {}
    for i in range(0, len(L) - 2, 3):
        if not L[i].startswith(">"):
            continue
        pid, typ = L[i][1:].split(" | ")
        lab = L[i + 2]
        hel = [(m.start() + 1, m.end()) for m in re.finditer(r"M+", lab)]
        # DeepTMHMM often splits one helix into two with a 1-3 residue "loop"
        # (expert review 2026-11-06). MERGE_GAP > 0 joins such pairs before counting.
        if MERGE_GAP and hel:
            merged = [hel[0]]
            for s_, e_ in hel[1:]:
                if s_ - merged[-1][1] - 1 <= MERGE_GAP:
                    merged[-1] = (merged[-1][0], e_)
                else:
                    merged.append((s_, e_))
            hel = merged
        beta = len(re.findall(r"B+", lab))
        body = lab.lstrip("S")
        nterm = {"I": "i", "O": "o", "P": "o"}.get(body[:1]) if hel else None
        d[pid.strip()] = {"type": typ.strip(), "n": len(hel), "hel": hel, "beta": beta,
                          "sp": lab.startswith("S"), "nterm": nterm}
    return d


def signalp(path):
    d = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        d[f[0].split()[0]] = f[1]
    return d


def heldout(path, ids):
    near = set()
    for line in open(path):
        q, t, pid, qc = line.split("\t")[:4]
        if float(qc) >= 0.5 and float(pid) >= 30:
            near.add(q)
    return {i: i not in near for i in ids}


def pct(n, d):
    return f"{100 * n / d:5.1f}% ({n}/{d})" if d else "   n/a"


def matched(pred, true):
    used, ok = set(), 0
    for s, e in true:
        for j, (ps, pe) in enumerate(pred):
            if j not in used and min(e, pe) - max(s, ps) + 1 >= OVERLAP:
                used.add(j); ok += 1; break
    return ok


# ───────────────────────── SwissProt ─────────────────────────
S = f"{A}/swissprot"
tm, dt, sp = tmhmm(f"{S}/tmhmm_short.tsv"), deeptmhmm(f"{S}/deeptmhmm/predicted_topologies.3line"), signalp(f"{S}/signalp/prediction_results.txt")
lab = {r["accession"]: r["sets"].split(",") for r in csv.DictReader(open(f"{S}/labels.tsv"), delimiter="\t")}
ho = heldout(f"{S}/vs_trainset.m8", lab)

expect = {  # set -> (description, test on (tmhmm, deeptmhmm, signalp) -> correct?)
    "multipass":  ("TM, >= 2 helices",   lambda t, d, s: t["n"] >= 2, lambda d: d["n"] >= 2),
    "singlepass": ("TM, >= 1 helix",     lambda t, d, s: t["n"] >= 1, lambda d: d["n"] >= 1),
    "peripheral": ("no TM helix",        lambda t, d, s: t["n"] == 0, lambda d: d["n"] == 0 and d["beta"] == 0),
    "soluble":    ("no TM helix",        lambda t, d, s: t["n"] == 0, lambda d: d["n"] == 0 and d["beta"] == 0),
    "signal":     ("no TM helix (SP)",   lambda t, d, s: t["n"] == 0, lambda d: d["n"] == 0),
    "secreted":   ("no TM helix (SP)",   lambda t, d, s: t["n"] == 0, lambda d: d["n"] == 0),
}
rows = []
print("\n## SwissProt fungal proteins with experimental labels")
print(f"{'set':11s} {'expected':18s} {'subset':8s} {'TMHMM correct':>22s} {'DeepTMHMM correct':>22s}  extra")
for st, (desc, ft, fd) in expect.items():
    for sub in ("all", "heldout"):
        ids = [a for a, s in lab.items() if st in s and (sub == "all" or ho[a])]
        ct = sum(ft(tm[a], dt[a], sp.get(a)) for a in ids)
        cd = sum(fd(dt[a]) for a in ids)
        extra = ""
        if st in ("signal", "secreted"):
            extra = (f"DeepTMHMM SP call {pct(sum(dt[a]['sp'] for a in ids), len(ids))}; "
                     f"SignalP SP {pct(sum(sp.get(a) == 'SP' for a in ids), len(ids))}")
        if st == "multipass":
            extra = f"DeepTMHMM beta-barrel {sum(dt[a]['beta'] > 0 for a in ids)}"
        print(f"{st:11s} {desc:18s} {sub:8s} {pct(ct, len(ids)):>22s} {pct(cd, len(ids)):>22s}  {extra}")
        rows.append([st, desc, sub, len(ids), ct, cd, extra])

# experimental helix positions and N-terminal side
ft_rows = {r["Entry"]: r for r in csv.DictReader(open(f"{S}/featexp.tsv"), delimiter="\t")}
INSIDE = re.compile(r"Cytoplasmic|Mitochondrial matrix|Nuclear|Chloroplast stroma", re.I)
def true_feats(r):
    hel = [(int(a), int(b)) for a, b in re.findall(r"TRANSMEM (\d+)\.\.(\d+)", r["Transmembrane"])]
    nterm = None
    if hel:
        for a, b, note in re.findall(r'TOPO_DOM (\d+)\.\.(\d+); /note="([^"]+)"', r["Topological domain"]):
            if int(b) < hel[0][0]:
                nterm = "i" if INSIDE.search(note) else "o"
    return hel, nterm
print(f"\n## Experimental helix positions (featexp; helix match = >= {OVERLAP} shared residues)")
for sub in ("all", "heldout"):
    ids = [a for a in ft_rows if a in lab and (sub == "all" or ho[a])]
    ids = [a for a in ids if true_feats(ft_rows[a])[0]]
    res = {}
    for name, pred in (("TMHMM", tm), ("DeepTMHMM", dt)):
        cnt = all_ok = topo_ok = topo_n = 0
        for a in ids:
            th, tn = true_feats(ft_rows[a])
            p = pred[a]
            cnt += p["n"] == len(th)
            all_ok += p["n"] == len(th) and matched(p["hel"], th) == len(th)
            if tn and p["nterm"]:
                topo_n += 1; topo_ok += p["nterm"] == tn
        res[name] = (cnt, all_ok, topo_ok, topo_n)
    print(f"{sub:8s} n={len(ids):3d}  " + "  ".join(
        f"{k}: count {pct(v[0], len(ids))}, all helices placed {pct(v[1], len(ids))}, N-term side {pct(v[2], v[3])}"
        for k, v in res.items()))

with open(f"{OUT}/swissprot_summary.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["set", "expected", "subset", "n", "tmhmm_correct", "deeptmhmm_correct", "extra"]); w.writerows(rows)
with open(f"{OUT}/swissprot_per_protein.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "sets", "heldout", "tmhmm_n", "deeptmhmm_type", "deeptmhmm_n", "deeptmhmm_beta", "signalp"])
    for a in sorted(lab):
        w.writerow([a, ",".join(lab[a]), int(ho[a]), tm[a]["n"], dt[a]["type"], dt[a]["n"], dt[a]["beta"], sp.get(a)])

# ───────────────────────── BFD ─────────────────────────
print("\n## BFD sample (no ground truth: agreement between tools), held-out proteins")
target = {"tm_pos": 1000, "tm_neg": 1000}
per = []
for folder in ("bfd_sample", "bfd_sample_neg", "bfd_sample_gpi"):
    D = f"{A}/{folder}"
    tmb, dtb, spb = tmhmm(f"{D}/tmhmm_short.tsv"), deeptmhmm(f"{D}/deeptmhmm/predicted_topologies.3line"), signalp(f"{D}/signalp/prediction_results.txt")
    smp = list(csv.DictReader(open(f"{D}/sample.tsv"), delimiter="\t"))
    hob = heldout(f"{D}/vs_trainset.m8", {r["protein_id"] for r in smp})
    by = defaultdict(list)
    for r in smp:
        # bfd_sample's first 'gpi' draw used every predgpi.parquet protein (that table
        # also has a 'Chain' row per protein); the redraw is in bfd_sample_gpi.
        if folder == "bfd_sample" and r["stratum"] == "gpi":
            continue
        by[r["stratum"]].append(r)
    for st, lst in sorted(by.items()):
        lst.sort(key=lambda r: hashlib.sha1((r["protein_id"] + SEED).encode()).hexdigest())
        keep = [r for r in lst if hob[r["protein_id"]]][: target.get(st, 300)]
        n = len(keep)
        ids = [r["protein_id"] for r in keep]
        typ = Counter(dtb[i]["type"] for i in ids)
        dn = Counter(min(dtb[i]["n"], 7) for i in ids)
        same_n = sum(dtb[i]["n"] == tmb[i]["n"] for i in ids)
        print(f"\n{st}: {n} held-out of {len(lst)} drawn ({len(lst) - sum(hob[r['protein_id']] for r in lst)} near training set); "
              f"{len({r['LOCUSTAG'] for r in keep})} genomes, {len({r['PHYLUM'] for r in keep})} phyla")
        print("  DeepTMHMM type: " + ", ".join(f"{k} {pct(v, n)}" for k, v in typ.most_common()))
        print(f"  same helix count as TMHMM: {pct(same_n, n)};  DeepTMHMM helices (7 = 7+): "
              + ", ".join(f"{k}:{v}" for k, v in sorted(dn.items())))
        print(f"  SignalP 6 SP: {pct(sum(spb.get(i) == 'SP' for i in ids), n)}")
        for i in ids:
            per.append([st, i, tmb[i]["n"], dtb[i]["type"], dtb[i]["n"], dtb[i]["beta"], spb.get(i)])
with open(f"{OUT}/bfd_per_protein.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["stratum", "protein_id", "tmhmm_n", "deeptmhmm_type", "deeptmhmm_n", "deeptmhmm_beta", "signalp"])
    w.writerows(per)
