"""Plan 9 summary: python breadth/plan9_summary.py <dir with selection.json, S*_info.json, b_S*_k*.json>
-> breadth/out/plan9_summary.json, paper2/plan9nums.tex, paper2/table_breadth.tex (rules exactly as breadth/PLAN9.md 9.4)."""
import os, sys, json, glob, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = sys.argv[1]; P2 = os.path.join(HERE, "..", "paper2"); os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
sel = json.load(open(os.path.join(D, "selection.json"))); M = ("PB", "KX", "XDSM", "XDSMh"); rows = []
# 9.1 excludes every collection already used in this project. The keyword match missed the PsychAD collection (its name lacks the word);
# it is excluded here by collection id and reported separately (not counted as a new cohort).
USED_COLL = {"84ce6837-548d-4a1f-919f-0bc0d9a3952f": "PsychAD (already used; keyword match missed it)"}
for s in sel["settings"]:
    i = s["i"]; fi = os.path.join(D, f"S{i}_info.json")
    fs = sorted(glob.glob(os.path.join(D, f"b_S{i}_k*.json")), key=lambda f: int(re.search(r"_k(\d+)\.json", f).group(1)))
    if not os.path.exists(fi) or not fs: rows.append(dict(i=i, cell_type=s["cell_type"], collection=s["collection_name"], status="not run")); continue
    info = json.load(open(fi)); R = [json.load(open(f)) for f in fs]; r = dict(i=i, info=info, n_runs=len(R), status="complete" if len(R) == 20 else f"{len(R)} runs")
    for m in M:
        c = np.array([x[m]["calls"] for x in R]); t = np.array([x[m]["t05"] for x in R])
        r[m] = dict(runs_with_call=int((c > 0).sum()), mean_calls=float(c.mean()), mean_t05=float(t.mean()), holds=bool((c > 0).sum() <= 2 and t.mean() <= 0.055 and len(R) == 20))
    r["xdsm_consistency"] = int(sum(x["XDSM"]["calls"] == x["XDSMcheck"]["calls"] for x in R)); rows.append(r)
for r, s_ in zip(rows, sel["settings"]): r["used_cohort"] = USED_COLL.get(s_["collection_id"])
allc = [r for r in rows if r.get("status") == "complete"]; done = [r for r in allc if not r["used_cohort"]]
out = dict(plan="breadth/PLAN9.md", plan_sha256=open(os.path.join(HERE, "PLAN9.sha256")).read().split()[0], n_selected=len(sel["settings"]), n_complete=len(done),
           excluded_collections=sel.get("excluded_collections"), n_eligible_settings=sel.get("n_eligible_settings"), settings=rows,
           held={m: int(sum(r[m]["holds"] for r in done)) for m in M}, pb_failed=int(sum(not r["PB"]["holds"] for r in done)))
json.dump(out, open(os.path.join(HERE, "out", "plan9_summary.json"), "w"), indent=1)
fmt = lambda x: f"{x:,.0f}".replace(",", "{,}") if x >= 100 else (f"{x:.1f}" if x >= 1 else f"{x:.2f}")
def esc(t): return str(t).replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")
lines = []
for r in allc:
    I = r["info"]; ct = esc(I["cell_type"]); coll = esc(I["collection_name"])[:58] + ("\\ldots" if len(I["collection_name"]) > 58 else "")
    cells = " & ".join(f"{fmt(r[m]['mean_calls'])} ({r[m]['runs_with_call']})" for m in M)
    coll = ("PsychAD$^{*}$" if r["used_cohort"] else coll)
    lines.append(f"{coll}; {ct} & {I['donors']} & {I['median_nd']:.0f} & {cells}\\\\")
tab = ("\\begin{table}[h]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n\\caption{Plan 9 (post hoc, hashed before any data were read): donor-centred null in the " + str(len(done)) +
       " CELLxGENE Census collections other than those already used that have a cell type with $\\ge200$ donors (one cell type per collection, the one with most donors; at most 200 cells per donor; 400-gene panels; 20 reassignments). Entries: mean BH false calls per run (runs of 20 with $\\ge1$ call). Donors: donors with $\\ge4$ cells; $n_d$: median cells per donor (capped at 200). $^{*}$Same cohort as our PsychAD data (MSSM astrocytes); selected because the exclusion keyword missed it, reported but not counted.}\\label{tab:breadth}\n"
       "\\begin{tabular}{@{}p{0.42\\linewidth}rrrrrr@{}}\\toprule\nCohort; cell type & Donors & $n_d$ & PB & KX & XDS-M & XDS-Mh\\\\\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\\end{table}\n")
open(os.path.join(P2, "table_breadth.tex"), "w").write(tab)
mac = dict(BrN=str(len(done)), BrPBfail=str(out["pb_failed"]), BrKXheld=str(out["held"]["KX"]), BrXheld=str(out["held"]["XDSM"]), BrHheld=str(out["held"]["XDSMh"]),
           BrDonors=(f"{min(r['info']['donors'] for r in done)}--{max(r['info']['donors'] for r in done):,}".replace(",", "{,}") if done else "--"),
           BrNd=(f"{min(r['info']['median_nd'] for r in done):.0f}--{max(r['info']['median_nd'] for r in done):.0f}" if done else "--"),
           BrPBcalls=(f"{fmt(min(r['PB']['mean_calls'] for r in done))}--{fmt(max(r['PB']['mean_calls'] for r in done))}" if done else "--"),
           BrXtI=(f"{min(r['XDSM']['mean_t05'] for r in done):.3f}--{max(r['XDSM']['mean_t05'] for r in done):.3f}" if done else "--"),
           BrHtI=(f"{min(r['XDSMh']['mean_t05'] for r in done):.3f}--{max(r['XDSMh']['mean_t05'] for r in done):.3f}" if done else "--"),
           BrSha=out["plan_sha256"][:16])
open(os.path.join(P2, "plan9nums.tex"), "w").write("% generated by breadth/plan9_summary.py; do not edit\n" + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in mac.items()))
print(json.dumps({k: v for k, v in out.items() if k != "settings"}, indent=1)); print(tab); print(mac)
for r in rows: print(r.get("i"), r.get("status"), (r.get("info") or {}).get("cell_type"), {m: (r[m]["runs_with_call"], round(r[m]["mean_calls"], 1), round(r[m]["mean_t05"], 4), r[m]["holds"]) for m in M if m in r}, r.get("xdsm_consistency"))
