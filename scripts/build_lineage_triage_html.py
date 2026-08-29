#!/usr/bin/env python3
"""Build lineage_triage.html: a filterable/sortable dashboard over
candidate_lineages.tsv + query_list.tsv (see query_unknowns_mash.py).

Usage:
    python3 build_lineage_triage_html.py \
        --candidates results/query_unknowns_mash/candidate_lineages.tsv \
        --query-list results/query_unknowns_mash/query_list.tsv \
        --out        results/query_unknowns_mash/lineage_triage.html
"""
import argparse
import csv
import json
from pathlib import Path

TEMPLATE = r"""<title>Lineage Triage</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Petrona:ital,wght@0,500;0,600;1,500;1,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{
  --bg:#edf1ec; --surface:#ffffff; --surface-2:#f6f8f4; --border:#d7ddd1; --border-strong:#c1cabb;
  --text:#1b211d; --text-muted:#5c665c; --text-faint:#8b968a;
  --accent:#2f6f62; --accent-strong:#1f4e44; --accent-wash:#dcebe6;
  --tier-veryhigh:#3f8f5c; --tier-high:#3e7ca6; --tier-moderate:#a8791f; --tier-low:#b95a26; --tier-verylow:#a83f57;
  --tier-veryhigh-wash:#e0f0e5; --tier-high-wash:#e1edf5; --tier-moderate-wash:#f4ecda; --tier-low-wash:#f6e4d8; --tier-verylow-wash:#f5dfe3;
  --flag-orphan:#6f4fa0; --flag-orphan-wash:#eae3f5;
  --shadow: 0 1px 2px rgba(20,30,20,.06), 0 4px 14px rgba(20,30,20,.05);
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#12160f; --surface:#181e17; --surface-2:#1f261d; --border:#2b332a; --border-strong:#3a4537;
    --text:#e7ece3; --text-muted:#9aab99; --text-faint:#6d7d6b;
    --accent:#5cb8a0; --accent-strong:#7fd4bc; --accent-wash:#1c332d;
    --tier-veryhigh:#6fcb8e; --tier-high:#6fa8d8; --tier-moderate:#dcae52; --tier-low:#df8752; --tier-verylow:#e0728a;
    --tier-veryhigh-wash:#183621; --tier-high-wash:#182b38; --tier-moderate-wash:#332c17; --tier-low-wash:#33231a; --tier-verylow-wash:#331d24;
    --flag-orphan:#b79be0; --flag-orphan-wash:#2a2038;
    --shadow: 0 1px 2px rgba(0,0,0,.4), 0 4px 18px rgba(0,0,0,.35);
  }
}
:root[data-theme="dark"]{
  --bg:#12160f; --surface:#181e17; --surface-2:#1f261d; --border:#2b332a; --border-strong:#3a4537;
  --text:#e7ece3; --text-muted:#9aab99; --text-faint:#6d7d6b;
  --accent:#5cb8a0; --accent-strong:#7fd4bc; --accent-wash:#1c332d;
  --tier-veryhigh:#6fcb8e; --tier-high:#6fa8d8; --tier-moderate:#dcae52; --tier-low:#df8752; --tier-verylow:#e0728a;
  --tier-veryhigh-wash:#183621; --tier-high-wash:#182b38; --tier-moderate-wash:#332c17; --tier-low-wash:#33231a; --tier-verylow-wash:#331d24;
  --flag-orphan:#b79be0; --flag-orphan-wash:#2a2038;
  --shadow: 0 1px 2px rgba(0,0,0,.4), 0 4px 18px rgba(0,0,0,.35);
}

*{box-sizing:border-box;}
body{
  background:var(--bg); color:var(--text);
  font-family:"IBM Plex Sans",system-ui,-apple-system,sans-serif;
  font-size:14px; line-height:1.5;
}
.sci{font-family:"Petrona",Georgia,serif; font-style:italic;}
.mono{font-family:"IBM Plex Mono",ui-monospace,Menlo,monospace; font-variant-numeric:tabular-nums;}

/* ---------- header ---------- */
header{
  position:sticky; top:0; z-index:20;
  background:var(--bg); border-bottom:1px solid var(--border);
  padding:20px 28px 16px;
}
.title-row{display:flex; align-items:baseline; gap:14px; flex-wrap:wrap; margin-bottom:4px;}
h1{
  font-family:"Petrona",Georgia,serif; font-weight:600; font-style:italic;
  font-size:30px; margin:0; text-wrap:balance; color:var(--text);
}
.subtitle{color:var(--text-muted); font-size:13px; max-width:65ch;}

.stats-strip{display:flex; gap:8px; flex-wrap:wrap; margin-top:14px;}
.stat-chip{
  display:flex; align-items:center; gap:7px; padding:6px 12px; border-radius:8px;
  border:1px solid var(--border); background:var(--surface); cursor:pointer;
  font-size:12.5px; color:var(--text-muted); transition:border-color .12s, background .12s;
  user-select:none;
}
.stat-chip:hover{border-color:var(--border-strong);}
.stat-chip.active{background:var(--accent-wash); border-color:var(--accent); color:var(--text); font-weight:500;}
.stat-chip .n{font-family:"IBM Plex Mono",monospace; font-weight:600; color:var(--text);}
.stat-chip .dot{width:9px; height:9px; border-radius:50%; flex:none;}

.controls-row{display:flex; gap:10px; align-items:center; margin-top:14px; flex-wrap:wrap;}
.search-wrap{position:relative; flex:1; min-width:220px; max-width:420px;}
.search-wrap svg{position:absolute; left:10px; top:50%; transform:translateY(-50%); opacity:.5; pointer-events:none;}
input[type="search"]{
  width:100%; padding:8px 12px 8px 32px; border-radius:8px; border:1px solid var(--border);
  background:var(--surface); color:var(--text); font-size:13px; font-family:inherit;
}
input[type="search"]:focus{outline:2px solid var(--accent); outline-offset:-1px;}
input[type="search"]::placeholder{color:var(--text-faint);}

select, .btn{
  padding:8px 12px; border-radius:8px; border:1px solid var(--border); background:var(--surface);
  color:var(--text); font-size:13px; font-family:inherit; cursor:pointer;
}
.btn:hover, select:hover{border-color:var(--border-strong);}
.btn:focus-visible, select:focus-visible, .stat-chip:focus-visible, .row:focus-visible{outline:2px solid var(--accent); outline-offset:1px;}

.result-count{font-size:12.5px; color:var(--text-faint); margin-left:auto; white-space:nowrap;}

/* ---------- list ---------- */
main{padding:16px 28px 60px; max-width:1180px; margin:0 auto;}
.legend{display:flex; gap:16px; flex-wrap:wrap; font-size:11.5px; color:var(--text-faint); margin:2px 0 14px; align-items:center;}
.legend .sw{display:inline-flex; align-items:center; gap:5px;}
.legend .dot{width:8px; height:8px; border-radius:50%;}

.list{display:flex; flex-direction:column; gap:7px;}
.row{
  background:var(--surface); border:1px solid var(--border); border-radius:10px;
  box-shadow:var(--shadow); overflow:hidden;
}
.row-head{
  display:grid; grid-template-columns:1.9fr 1.6fr 1fr 1.5fr auto; gap:14px; align-items:center;
  padding:12px 16px; cursor:pointer;
}
.row-head:hover{background:var(--surface-2);}
.q-id{display:flex; flex-direction:column; gap:2px; min-width:0;}
.q-acc{font-size:12px; color:var(--text-faint);}
.q-name{font-size:14px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;}
.q-status{font-size:10.5px; text-transform:uppercase; letter-spacing:.05em; color:var(--text-faint);}

.best-hit{display:flex; flex-direction:column; gap:2px; min-width:0;}
.best-hit .name{font-size:13.5px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;}
.best-hit .rank-path{font-size:11px; color:var(--text-faint); overflow:hidden; text-overflow:ellipsis; white-space:nowrap;}

.ani-badge{
  display:inline-flex; align-items:baseline; gap:4px; padding:4px 9px; border-radius:999px;
  font-family:"IBM Plex Mono",monospace; font-size:13px; font-weight:600; width:fit-content;
}
.ani-badge .pct-sign{font-size:10px; font-weight:500; opacity:.75;}

.phylum-tag{font-size:11.5px; color:var(--text-muted);}
.orphan-flag{
  display:inline-flex; align-items:center; gap:5px; font-size:10.5px; padding:3px 8px; border-radius:6px;
  background:var(--flag-orphan-wash); color:var(--flag-orphan); font-weight:500; width:fit-content; margin-top:3px;
}

.chevron{transition:transform .15s; color:var(--text-faint); flex:none;}
.row.open .chevron{transform:rotate(90deg);}

.row-detail{display:none; border-top:1px solid var(--border); background:var(--surface-2);}
.row.open .row-detail{display:block;}
.detail-inner{padding:6px 16px 14px; overflow-x:auto;}
table.hits{width:100%; border-collapse:collapse; min-width:720px;}
table.hits th{
  text-align:left; font-size:10.5px; text-transform:uppercase; letter-spacing:.04em; color:var(--text-faint);
  font-weight:500; padding:8px 10px 6px; border-bottom:1px solid var(--border);
}
table.hits td{padding:7px 10px; font-size:12.5px; border-bottom:1px solid var(--border); vertical-align:baseline;}
table.hits tr:last-child td{border-bottom:none;}
table.hits td.rank{color:var(--text-faint); font-family:"IBM Plex Mono",monospace;}
table.hits td.num{font-family:"IBM Plex Mono",monospace; text-align:right;}
.sig-yes{color:var(--tier-veryhigh);}
.sig-no{color:var(--text-faint);}

.empty-state{text-align:center; padding:60px 20px; color:var(--text-faint);}
.empty-state .sci{font-size:16px; display:block; margin-bottom:6px;}

footer{text-align:center; font-size:11.5px; color:var(--text-faint); padding:20px 28px 40px;}

@media (max-width:720px){
  .row-head{grid-template-columns:1fr auto;}
  .best-hit,.phylum-tag{display:none;}
}
</style>

<header>
  <div class="title-row">
    <h1>Lineage Triage</h1>
  </div>
  <p class="subtitle">Mash screen of 341 genomes in samples.csv missing a GENUS assignment against the results/mash_sketch/ reference cache (22,574 sketches, k=21 s=10000). Each card is one unplaced genome; expand it for its top-10 nearest reference sketches by ANI%.</p>

  <div class="stats-strip" id="tierChips"></div>

  <div class="controls-row">
    <div class="search-wrap">
      <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="7"/><line x1="21" y1="21" x2="16.65" y2="16.65"/></svg>
      <input type="search" id="searchInput" placeholder="Search accession, species, or candidate genus…">
    </div>
    <select id="phylumSelect"><option value="">All phyla (best hit)</option></select>
    <button class="btn" id="sortBtn" data-dir="asc">Sort: weakest match first ↑</button>
    <span class="result-count" id="resultCount"></span>
  </div>
</header>

<main>
  <div class="legend" id="legend"></div>
  <div class="list" id="list"></div>
  <div class="empty-state" id="emptyState" style="display:none;">
    <span class="sci">No genomes match this filter</span>
    Try clearing the search or the phylum filter.
  </div>
</main>

<footer>candidate_lineages.tsv → query_unknowns_mash.py · mash dist, k=21 s=10000 · generated from results/query_unknowns_mash/</footer>

<script>
const DATA = __DATA_JSON__;

const TIERS = [
  {key:"veryhigh", label:"≥95% very high", min:95, max:101, color:"var(--tier-veryhigh)", wash:"var(--tier-veryhigh-wash)"},
  {key:"high",     label:"90–95% high",     min:90, max:95,  color:"var(--tier-high)",     wash:"var(--tier-high-wash)"},
  {key:"moderate", label:"80–90% moderate",  min:80, max:90,  color:"var(--tier-moderate)", wash:"var(--tier-moderate-wash)"},
  {key:"low",      label:"70–80% low",       min:70, max:80,  color:"var(--tier-low)",      wash:"var(--tier-low-wash)"},
  {key:"verylow",  label:"<70% very low",    min:-1, max:70,  color:"var(--tier-verylow)",  wash:"var(--tier-verylow-wash)"},
];
function tierOf(ani){
  if (ani==null) return TIERS[4];
  for (const t of TIERS) if (ani>=t.min && ani<t.max) return t;
  return TIERS[4];
}

const state = { search:"", tiers:new Set(), phylum:"", sortDir:"asc", open:new Set() };

// ---- build phylum options from best-hit phylum of each query ----
const phylumCounts = new Map();
DATA.forEach(q=>{
  const p = (q.hits[0]||{}).p || "(none)";
  phylumCounts.set(p, (phylumCounts.get(p)||0)+1);
});
const phylumSelect = document.getElementById("phylumSelect");
[...phylumCounts.entries()].sort((a,b)=>b[1]-a[1]).forEach(([p,n])=>{
  const o = document.createElement("option");
  o.value = p; o.textContent = `${p} (${n})`;
  phylumSelect.appendChild(o);
});
phylumSelect.addEventListener("change", ()=>{ state.phylum = phylumSelect.value; render(); });

// ---- tier stat chips (clickable multi-filter) ----
const tierChipsEl = document.getElementById("tierChips");
function buildChips(){
  const totalChip = document.createElement("div");
  totalChip.className = "stat-chip" + (state.tiers.size===0 ? " active":"");
  totalChip.innerHTML = `<span class="n">${DATA.length}</span> genomes queried`;
  totalChip.onclick = ()=>{ state.tiers.clear(); render(); };
  tierChipsEl.appendChild(totalChip);

  TIERS.forEach(t=>{
    const n = DATA.filter(q=>tierOf(q.best).key===t.key).length;
    const chip = document.createElement("div");
    chip.className = "stat-chip" + (state.tiers.has(t.key) ? " active":"");
    chip.innerHTML = `<span class="dot" style="background:${t.color}"></span><span class="n">${n}</span> ${t.label}`;
    chip.onclick = ()=>{
      if (state.tiers.has(t.key)) state.tiers.delete(t.key); else state.tiers.add(t.key);
      render();
    };
    tierChipsEl.appendChild(chip);
  });

  const orphanN = DATA.filter(q=>!(q.hits[0]||{}).g).length;
  const orphanChip = document.createElement("div");
  orphanChip.className = "stat-chip" + (state.tiers.has("orphan") ? " active":"");
  orphanChip.innerHTML = `<span class="dot" style="background:var(--flag-orphan)"></span><span class="n">${orphanN}</span> best hit also unplaced`;
  orphanChip.onclick = ()=>{
    if (state.tiers.has("orphan")) state.tiers.delete("orphan"); else state.tiers.add("orphan");
    render();
  };
  tierChipsEl.appendChild(orphanChip);
}

// ---- legend ----
document.getElementById("legend").innerHTML = TIERS.map(t=>
  `<span class="sw"><span class="dot" style="background:${t.color}"></span>${t.label}</span>`
).join("") + `<span class="sw"><span class="dot" style="background:var(--flag-orphan)"></span>best hit is itself unplaced</span>`;

// ---- search ----
document.getElementById("searchInput").addEventListener("input", e=>{
  state.search = e.target.value.trim().toLowerCase();
  render();
});

document.getElementById("sortBtn").addEventListener("click", ()=>{
  state.sortDir = state.sortDir==="asc" ? "desc" : "asc";
  render();
});

function fmtPct(v){ return v==null ? "—" : v.toFixed(1); }
function esc(s){ return (s||"").replace(/[&<>"']/g, c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c])); }

function matches(q){
  if (state.tiers.size){
    const isOrphan = !(q.hits[0]||{}).g;
    const tierKey = tierOf(q.best).key;
    const okTier = state.tiers.has(tierKey) || (state.tiers.has("orphan") && isOrphan);
    // require match against ANY selected chip (union), but if only "orphan" style filters chosen, honor as OR
    const anyTierSelected = [...state.tiers].some(k=>k!=="orphan");
    const anyOrphanSelected = state.tiers.has("orphan");
    let pass = false;
    if (anyTierSelected && state.tiers.has(tierKey)) pass = true;
    if (anyOrphanSelected && isOrphan) pass = true;
    if (!pass) return false;
  }
  if (state.phylum){
    const p = (q.hits[0]||{}).p || "(none)";
    if (p !== state.phylum) return false;
  }
  if (state.search){
    const hay = [q.q, q.sp, ...(q.hits[0] ? [q.hits[0].g, q.hits[0].sp, q.hits[0].f] : [])]
      .join(" ").toLowerCase();
    if (!hay.includes(state.search)) return false;
  }
  return true;
}

function rowTemplate(q){
  const best = q.hits[0];
  const tier = tierOf(q.best);
  const isOrphan = best && !best.g;
  const open = state.open.has(q.q);
  return `
  <div class="row ${open ? "open":""}" data-q="${esc(q.q)}">
    <div class="row-head" tabindex="0" role="button" aria-expanded="${open}">
      <div class="q-id">
        <span class="q-acc mono">${esc(q.q)}</span>
        <span class="q-name sci">${esc(q.sp)}</span>
        <span class="q-status">sketch: ${q.st}</span>
      </div>
      <div class="best-hit">
        ${best ? `
          <span class="name sci">${esc(best.g || best.sp)}</span>
          <span class="rank-path">${[best.p,best.c,best.o,best.f].filter(Boolean).join(" › ") || "—"}</span>
          ${isOrphan ? `<span class="orphan-flag">best hit also unplaced</span>` : ""}
        ` : `<span class="rank-path">no reference hits</span>`}
      </div>
      <div>
        <span class="ani-badge" style="background:${tier.wash}; color:${tier.color}">${fmtPct(q.best)}<span class="pct-sign">% ANI</span></span>
      </div>
      <div class="phylum-tag">${best ? esc(best.p||"—") : "—"}</div>
      <svg class="chevron" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="9 18 15 12 9 6"/></svg>
    </div>
    <div class="row-detail">
      <div class="detail-inner">
        <table class="hits">
          <thead><tr>
            <th>#</th><th>Candidate genus / species</th><th>Family</th><th>Order</th><th>Class</th><th>Phylum</th>
            <th style="text-align:right">ANI %</th><th style="text-align:right">p-value</th><th style="text-align:right">shared k-mers</th><th>sig.</th>
          </tr></thead>
          <tbody>
            ${q.hits.map(h=>`
              <tr>
                <td class="rank">${h.r}</td>
                <td class="sci">${esc(h.g || h.sp || "—")}${h.g && h.sp ? ` <span class="mono" style="font-style:normal;color:var(--text-faint)">(${esc(h.sp)})</span>`:""}</td>
                <td>${esc(h.f)||"—"}</td>
                <td>${esc(h.o)||"—"}</td>
                <td>${esc(h.c)||"—"}</td>
                <td>${esc(h.p)||"—"}</td>
                <td class="num">${h.ani.toFixed(2)}</td>
                <td class="num">${Number(h.pv).toExponential(2)}</td>
                <td class="num">${esc(h.sh)}</td>
                <td class="${h.sig?"sig-yes":"sig-no"}">${h.sig?"yes":"no"}</td>
              </tr>`).join("")}
          </tbody>
        </table>
      </div>
    </div>
  </div>`;
}

function render(){
  tierChipsEl.innerHTML = "";
  buildChips();

  let rows = DATA.filter(matches);
  rows.sort((a,b)=>{
    const av = a.best==null?-1:a.best, bv = b.best==null?-1:b.best;
    return state.sortDir==="asc" ? av-bv : bv-av;
  });

  document.getElementById("sortBtn").textContent =
    state.sortDir==="asc" ? "Sort: weakest match first ↑" : "Sort: strongest match first ↓";
  document.getElementById("resultCount").textContent = `${rows.length} of ${DATA.length}`;

  const list = document.getElementById("list");
  document.getElementById("emptyState").style.display = rows.length ? "none" : "block";
  list.innerHTML = rows.map(rowTemplate).join("");

  list.querySelectorAll(".row-head").forEach(el=>{
    const toggle = ()=>{
      const rowEl = el.closest(".row");
      const qid = rowEl.dataset.q;
      if (state.open.has(qid)) state.open.delete(qid); else state.open.add(qid);
      rowEl.classList.toggle("open");
      el.setAttribute("aria-expanded", rowEl.classList.contains("open"));
    };
    el.addEventListener("click", toggle);
    el.addEventListener("keydown", e=>{ if (e.key==="Enter"||e.key===" "){ e.preventDefault(); toggle(); } });
  });
}

render();
</script>
"""


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", default="results/query_unknowns_mash/candidate_lineages.tsv")
    p.add_argument("--query-list", default="results/query_unknowns_mash/query_list.tsv")
    p.add_argument("--out", default="results/query_unknowns_mash/lineage_triage.html")
    return p.parse_args()


def build_data(candidates_path, query_list_path):
    hits_by_query = {}
    with open(candidates_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            hits_by_query.setdefault(row["query_asmid"], []).append({
                "r": int(row["rank"]),
                "a": row["ref_asmid"],
                "sp": row["ref_species"],
                "g": row["ref_genus"],
                "f": row["ref_family"],
                "o": row["ref_order"],
                "c": row["ref_class"],
                "p": row["ref_phylum"],
                "ani": round(float(row["ani_pct"]), 2),
                "pv": row["mash_pval"],
                "sh": row["shared_hashes"],
                "sig": row["significant"] == "True",
            })

    queries = []
    with open(query_list_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asmid = row["query_asmid"]
            queries.append({
                "q": asmid,
                "sp": row["query_species_in"],
                "st": row["sketch_status"],
                "best": round(float(row["best_ani_pct"]), 2) if row["best_ani_pct"] else None,
                "hits": hits_by_query.get(asmid, []),
            })
    return queries


def main():
    args = parse_args()
    data = build_data(args.candidates, args.query_list)
    data_json = json.dumps(data, separators=(",", ":"))
    html = TEMPLATE.replace("__DATA_JSON__", data_json)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    print(f"[build_lineage_triage_html] wrote {out_path} "
          f"({len(data)} queries, {sum(len(q['hits']) for q in data)} hit rows, {len(html):,} bytes)")


if __name__ == "__main__":
    main()
