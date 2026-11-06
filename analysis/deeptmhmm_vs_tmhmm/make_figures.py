#!/usr/bin/env python3
"""make_figures.py -- figures for REPORT.md (PNG + PDF in figures/).

Numbers come from results/summary.txt, results/combined_rule.txt (written by
analyze.py and combined_rule.py), results/bfd_per_protein.tsv and
review/expert_review.tsv. Palette: the dataviz skill's reference categorical
slots, checked with its validator (light surface #fcfcfb):
TMHMM #2a78d6, TMHMM+SignalP #eb6834, DeepTMHMM #1baf7a -> all checks pass
(contrast WARN for aqua is covered by direct value labels).
"""
import csv, os, re
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = os.path.dirname(os.path.abspath(__file__))
F = f"{A}/figures"; os.makedirs(F, exist_ok=True)
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
TOOLS = {"TMHMM": "#2a78d6", "TMHMM+SignalP": "#eb6834", "DeepTMHMM": "#1baf7a"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{F}/{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ── Figure 1: SwissProt held-out, % correct per set and tool (dot plot) ──
summ = open(f"{A}/results/summary.txt").read()
comb = open(f"{A}/results/combined_rule.txt").read()
def pct_line(text, setname, col):
    for line in text.splitlines():
        if line.startswith(setname) and " heldout " in line:
            vals = re.findall(r"(\d+\.\d)% \((\d+)/(\d+)\)", line)
            return float(vals[col][0]), int(vals[col][2])
    raise KeyError(setname)
def type_line(kind):
    for line in comb.splitlines():
        if line.startswith(kind + " ") or line.startswith(kind.ljust(12)):
            v = re.findall(r"(\d+\.\d)% \((\d+)/(\d+)\)", line)
            return [float(x[0]) for x in v], int(v[0][2])
rows = [("Multi-pass membrane", "multipass", ">= 2 helices"), ("Single-pass membrane", "singlepass", ">= 1 helix")]
data = []
for label, key, exp in rows:
    t, n = pct_line(comb, key, 0); c, _ = pct_line(comb, key, 1); d, _ = pct_line(comb, key, 2)
    data.append((f"{label} (n={n})", exp, t, c, d))
v, n = type_line("type II")
data.append((f"  type II signal anchor (n={n})", ">= 1 helix", *v))
for label, key, exp in [("Peripheral membrane", "peripheral", "no helix"), ("Soluble", "soluble", "no helix"),
                        ("Experimental signal peptide", "signal", "no helix"), ("Secreted", "secreted", "no helix")]:
    t, n = pct_line(comb, key, 0); c, _ = pct_line(comb, key, 1); d, _ = pct_line(comb, key, 2)
    data.append((f"{label} (n={n})", exp, t, c, d))
fig, ax = plt.subplots(figsize=(7.2, 3.9))
ys = list(range(len(data)))[::-1]
for y, (label, exp, t, c, d) in zip(ys, data):
    vals = {"TMHMM": t, "TMHMM+SignalP": c, "DeepTMHMM": d}
    lo, hi = min(vals.values()), max(vals.values())
    ax.plot([lo, hi], [y, y], color=GRID, lw=2, zorder=1)
    for i, (tool, x) in enumerate(vals.items()):
        ax.scatter(x, y + (i - 1) * 0.24, s=46, color=TOOLS[tool], edgecolor=SURF, linewidth=1.5, zorder=3,
                   label=tool if y == ys[0] else None)
        # values near 100% sit in a column right of the axis so they do not collide
        tx = 100.8 if x >= 96 else x + 0.6
        ax.text(tx, y + (i - 1) * 0.24, f"{x:.1f}", va="center", ha="left", fontsize=7,
                color=TOOLS[tool] if x >= 96 else INK2, fontweight="bold" if x >= 96 else "normal")
ax.set_yticks(ys); ax.set_yticklabels([f"{l}\n{e}" for l, e, *_ in data], fontsize=8, color=INK)
ax.set_xlim(60, 103.5); ax.set_xticks(range(60, 101, 5)); ax.set_xlabel("Correct call (%), held-out proteins")
ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.legend(loc="lower left", frameon=False, ncol=3, bbox_to_anchor=(0, 1.0), fontsize=8)
save(fig, "fig1_swissprot_accuracy")

# ── Figure 2: BFD strata, DeepTMHMM call types (100% stacked bars) ──
per = list(csv.DictReader(open(f"{A}/results/bfd_per_protein.tsv"), delimiter="\t"))
order = [("tm_neg", "TMHMM 0 helices"), ("sp_only", "SignalP SP, TMHMM 0 helices"),
         ("gpi", "PredGPI GPI-anchor"), ("sp_overlap", "1 N-term helix + SignalP SP"),
         ("nterm_single", "1 N-terminal helix"), ("cterm_single", "1 C-terminal helix"),
         ("tm_pos", "TMHMM >= 1 helix"), ("multi7", "TMHMM >= 7 helices")]
CALL = [("TM", "#2a78d6"), ("SP+TM", "#1baf7a"), ("SP", "#eb6834"), ("GLOB", "#c9c7bf")]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for y, (st, lab) in zip(range(len(order))[::-1], order):
    c = Counter(r["deeptmhmm_type"] for r in per if r["stratum"] == st)
    n = sum(c.values()); left = 0
    for k, col in CALL:
        w = 100 * c.get(k, 0) / n
        if w:
            ax.barh(y, w, left=left, color=col, height=0.62, edgecolor=SURF, linewidth=2)
            if w >= 7:
                ax.text(left + w / 2, y, f"{w:.0f}", ha="center", va="center", fontsize=7,
                        color="#ffffff" if k in ("TM",) else INK)
        left += w
    ax.text(101, y, f"n={n}", va="center", fontsize=7, color=INK2)
ax.set_yticks(range(len(order))[::-1]); ax.set_yticklabels([l for _, l in order], fontsize=8, color=INK)
ax.set_xlim(0, 108); ax.set_xlabel("DeepTMHMM call (% of proteins in stratum)")
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col) for _, col in CALL], labels=[k for k, _ in CALL],
          loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, frameon=False, fontsize=8)
save(fig, "fig2_bfd_strata")

# ── Figure 3: expert review verdicts per group (stacked counts) ──
rev = list(csv.DictReader(open(f"{A}/review/expert_review.tsv"), delimiter="\t"))
groups = []
for r in rev:
    g = r["case_group"].split(":")[0]
    if g not in groups: groups.append(g)
glabel = {"cterm_single": "1. C-terminal single helix", "nterm_single": "2. N-term helix, DeepTMHMM SP",
          "tm_pos": "3. DeepTMHMM globular", "gpi": "4. PredGPI positive", "multi7": "5. TMHMM >= 7 helices"}
VER = [("TMHMM", TOOLS["TMHMM"]), ("DeepTMHMM", TOOLS["DeepTMHMM"]), ("both wrong", "#8a8880"), ("unclear", "#d9d7d0")]
fig, ax = plt.subplots(figsize=(7.2, 2.9))
for y, g in zip(range(len(groups))[::-1], groups):
    c = Counter(r["verdict"].strip() for r in rev if r["case_group"].split(":")[0] == g)
    left = 0
    for k, col in VER:
        w = c.get(k, 0)
        if w:
            ax.barh(y, w, left=left, color=col, height=0.62, edgecolor=SURF, linewidth=2)
            ax.text(left + w / 2, y, str(w), ha="center", va="center", fontsize=7,
                    color="#ffffff" if k in ("TMHMM", "both wrong") else INK)
        left += w
ax.set_yticks(range(len(groups))[::-1]); ax.set_yticklabels([glabel.get(g, g) for g in groups], fontsize=8, color=INK)
ax.set_xlabel("Proteins (16 per group); verdict = which tool the reviewer judged right")
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col) for _, col in VER], labels=[k for k, _ in VER],
          loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, frameon=False, fontsize=8)
save(fig, "fig3_expert_review")
print("figures written to", F)
