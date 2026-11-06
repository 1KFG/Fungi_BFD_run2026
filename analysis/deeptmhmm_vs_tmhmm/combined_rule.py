#!/usr/bin/env python3
"""combined_rule.py -- score "TMHMM + SignalP" against TMHMM and DeepTMHMM on the
SwissProt sets.

Rule: take TMHMM's helices; when SignalP 6 calls SP, drop TMHMM's first helix if
it starts at or before the predicted cleavage site + 5 (the helix is the signal
peptide's hydrophobic core). Everything else is TMHMM unchanged.

Also scores type II (signal-anchor) single-pass proteins separately, the case
where an N-terminal helix must be kept.
Run after analyze.py; writes results/combined_rule.txt.
"""
import csv, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
A = os.path.dirname(os.path.abspath(__file__))
S = f"{A}/swissprot"
import importlib.util
spec = importlib.util.spec_from_file_location("an", f"{A}/analyze.py")

# re-use the parsers without re-running analyze.py's report
src = open(f"{A}/analyze.py").read().split("# ───────────────────────── SwissProt")[0]
ns = {"__file__": f"{A}/analyze.py"}
exec(src, ns)
tmhmm, deeptmhmm, heldout, pct, matched, OVERLAP = (ns[k] for k in ("tmhmm", "deeptmhmm", "heldout", "pct", "matched", "OVERLAP"))

tm = tmhmm(f"{S}/tmhmm_short.tsv")
dt = deeptmhmm(f"{S}/deeptmhmm/predicted_topologies.3line")
cs = {}
for line in open(f"{S}/signalp/prediction_results.txt"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[1] == "SP":
        m = re.search(r"CS pos: (\d+)-", f[4])
        cs[f[0]] = int(m.group(1)) if m else 30
lab = {r["accession"]: r["sets"].split(",") for r in csv.DictReader(open(f"{S}/labels.tsv"), delimiter="\t")}
ho = heldout(f"{S}/vs_trainset.m8", lab)


def combined(a):
    hel = list(tm[a]["hel"])
    if a in cs and hel and hel[0][0] <= cs[a] + 5:
        hel = hel[1:]
    return {"n": len(hel), "hel": hel}


comb = {a: combined(a) for a in lab}
sub = {}
for st in ("singlepass",):
    for r in csv.DictReader(open(f"{S}/{st}.tsv"), delimiter="\t"):
        m = re.findall(r"Single-pass (type [IV]+ |)membrane protein", r["Subcellular location [CC]"])
        sub[r["Entry"]] = (m[0].strip() or "unspecified") if m else "other"

tests = {
    "multipass":  ("TM, >= 2 helices", lambda p: p["n"] >= 2),
    "singlepass": ("TM, >= 1 helix",   lambda p: p["n"] >= 1),
    "peripheral": ("no TM helix",      lambda p: p["n"] == 0),
    "soluble":    ("no TM helix",      lambda p: p["n"] == 0),
    "signal":     ("no TM helix (SP)", lambda p: p["n"] == 0),
    "secreted":   ("no TM helix (SP)", lambda p: p["n"] == 0),
}
out = []
def say(s=""):
    print(s); out.append(s)
say("## TMHMM + SignalP rule vs TMHMM and DeepTMHMM (SwissProt, experimental labels)")
say(f"{'set':11s} {'expected':17s} {'subset':8s} {'TMHMM':>20s} {'TMHMM+SignalP':>20s} {'DeepTMHMM':>20s}")
for st, (desc, ok) in tests.items():
    for s in ("all", "heldout"):
        ids = [a for a, v in lab.items() if st in v and (s == "all" or ho[a])]
        say(f"{st:11s} {desc:17s} {s:8s} {pct(sum(ok(tm[a]) for a in ids), len(ids)):>20s} "
            f"{pct(sum(ok(comb[a]) for a in ids), len(ids)):>20s} {pct(sum(ok(dt[a]) for a in ids), len(ids)):>20s}")
say("\n## Single-pass proteins by type (held out): called >= 1 helix")
for k in ("type I", "type II", "type III", "type IV", "unspecified"):
    ids = [a for a, v in lab.items() if "singlepass" in v and ho[a] and sub.get(a) == k]
    if ids:
        say(f"{k:12s} n={len(ids):3d}  TMHMM {pct(sum(tm[a]['n'] >= 1 for a in ids), len(ids))}  "
            f"TMHMM+SignalP {pct(sum(comb[a]['n'] >= 1 for a in ids), len(ids))}  "
            f"DeepTMHMM {pct(sum(dt[a]['n'] >= 1 for a in ids), len(ids))}")
ft_rows = {r["Entry"]: r for r in csv.DictReader(open(f"{S}/featexp.tsv"), delimiter="\t")}
say(f"\n## Experimental helix positions (featexp, helix match >= {OVERLAP} shared residues)")
for s in ("all", "heldout"):
    ids = [a for a in ft_rows if a in lab and (s == "all" or ho[a])]
    tru = {a: [(int(x), int(y)) for x, y in re.findall(r"TRANSMEM (\d+)\.\.(\d+)", ft_rows[a]["Transmembrane"])] for a in ids}
    ids = [a for a in ids if tru[a]]
    line = f"{s:8s} n={len(ids):3d} "
    for name, pred in (("TMHMM", tm), ("TMHMM+SignalP", comb), ("DeepTMHMM", dt)):
        okc = sum(pred[a]["n"] == len(tru[a]) and matched(pred[a]["hel"], tru[a]) == len(tru[a]) for a in ids)
        line += f" {name}: all helices placed {pct(okc, len(ids))};"
    say(line)
open(f"{A}/results/combined_rule.txt", "w").write("\n".join(out) + "\n")
