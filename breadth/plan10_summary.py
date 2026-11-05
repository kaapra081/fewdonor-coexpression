"""Plan 10 summary: python breadth/plan10_summary.py <outdir A> [<outdir B> ...]  (results from Modal and Colab; rules as breadth/PLAN10.md)
-> breadth/out/plan10_summary.json, paper2/plan10nums.tex, paper2/table_breadth10.tex"""
import os, sys, json, glob, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); P2 = os.path.join(HERE, "..", "paper2"); DIRS = sys.argv[1:]; M = ("PB", "KX", "XDSM", "XDSMh")
sels = [json.load(open(os.path.join(d, "selection.json"))) for d in DIRS]
assert all(json.dumps(s["settings"], sort_keys=True) == json.dumps(sels[0]["settings"], sort_keys=True) for s in sels), "selections differ between machines"
sel = sels[0]; rows = []; overlap = []
for s in sel["settings"]:
    i = s["i"]; per = {}
    for d in DIRS:
        fs = {int(re.search(r"_k(\d+)\.json", f).group(1)): json.load(open(f)) for f in glob.glob(os.path.join(d, f"b_S{i}_k*.json"))}
        fi = os.path.join(d, f"S{i}_info.json")
        if len(fs) == 20 and os.path.exists(fi): per[d] = (fs, json.load(open(fi)))
    if not per: rows.append(dict(i=i, cell_type=s["cell_type"], collection=s["collection_name"], status="not run")); continue
    ds = list(per)
    if len(ds) > 1:
        a, b = per[ds[0]][0], per[ds[1]][0]
        overlap.append(dict(i=i, identical_calls=int(sum(all(a[k][m]["calls"] == b[k][m]["calls"] for m in M) for k in range(20))), n=20))
    R = [per[ds[0]][0][k] for k in range(20)]; info = per[ds[0]][1]; r = dict(i=i, info=info, status="complete", machine=os.path.basename(os.path.dirname(ds[0].rstrip("/"))) or ds[0])
    for m in M:
        c = np.array([x[m]["calls"] for x in R]); t = np.array([x[m]["t05"] for x in R])
        r[m] = dict(runs_with_call=int((c > 0).sum()), mean_calls=float(c.mean()), mean_t05=float(t.mean()), holds=bool((c > 0).sum() <= 2 and t.mean() <= 0.055))
    r["xdsm_consistency"] = int(sum(x["XDSM"]["calls"] == x["XDSMcheck"]["calls"] for x in R)); rows.append(r)
done = [r for r in rows if r["status"] == "complete"]; colls = sorted(set(r["info"]["collection_name"] for r in done))
out = dict(plan="breadth/PLAN10.md", plan_sha256=open(os.path.join(HERE, "PLAN10.sha256")).read().split()[0], n_selected=len(sel["settings"]), n_complete=len(done), n_collections=len(colls),
           n_eligible_settings=sel.get("n_eligible_settings"), excluded_collections=sel.get("excluded_collections"), overlap=overlap, settings=rows,
           held={m: int(sum(r[m]["holds"] for r in done)) for m in M}, pb_failed=int(sum(not r["PB"]["holds"] for r in done)),
           pb_runs_all20=int(sum(r["PB"]["runs_with_call"] == 20 for r in done)))
os.makedirs(os.path.join(HERE, "out"), exist_ok=True); json.dump(out, open(os.path.join(HERE, "out", "plan10_summary.json"), "w"), indent=1)
fmt = lambda x: f"{x:,.0f}".replace(",", "{,}") if x >= 100 else (f"{x:.1f}" if x >= 1 else f"{x:.2f}")
def esc(t): return str(t).replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")
SHORT = {"Asian Immune Diversity Atlas (AIDA)": "AIDA blood", "The integrated Human Lung Cell Atlas": "Lung atlas (HLCA)", "Single cell atlas of large B-cell lymphoma": "B-cell lymphoma",
         "Harmonized single-cell landscape, intercellular crosstalk and tumor architecture of glioblastoma": "Glioblastoma (GBmap)"}
lines = []
for r in done:
    I = r["info"]; cn = SHORT.get(I["collection_name"], I["collection_name"]); cn = esc(cn[:34] + ("\\ldots" if len(cn) > 34 else ""))
    lines.append(f"{cn}; {esc(I['cell_type'][:34])} & {I['donors']} & {I['median_nd']:.0f} & " + " & ".join(f"{fmt(r[m]['mean_calls'])} ({r[m]['runs_with_call']})" for m in M) + "\\\\")
tab = ("\\begin{table}[h]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n\\caption{Plan 10 (post hoc, hashed before any of its data were read): donor-centred null in the largest CELLxGENE Census cohorts not otherwise used with a cell type of $\\ge100$ donors (9 of 18 qualifying; up to three cell types per cohort, the three with most donors; at most 24 settings, ranked by donors; all cells; 400-gene panels; 20 reassignments). Entries: mean BH false calls per run (runs of 20 with $\\ge1$ call). Donors: donors with $\\ge4$ cells; $n_d$: median cells per donor."
       + (f" Settings computed on both machines gave identical call counts in {sum(o['identical_calls'] for o in overlap)} of {sum(o['n'] for o in overlap)} runs." if overlap else "") + "}\\label{tab:breadth10}\n"
       "\\begin{tabular}{@{}p{0.44\\linewidth}rrrrrr@{}}\\toprule\nCohort; cell type & Donors & $n_d$ & PB & KX & XDS-M & XDS-Mh\\\\\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\\end{table}\n")
open(os.path.join(P2, "table_breadth10.tex"), "w").write(tab)
rng = lambda v, f: f"{f(min(v))}--{f(max(v))}"
mac = dict(BtN=str(len(done)), BtSel=str(len(sel["settings"])), BtColl=str(len(colls)), BtPBfail=str(out["pb_failed"]), BtPBall=str(out["pb_runs_all20"]), BtKXheld=str(out["held"]["KX"]),
           BtXheld=str(out["held"]["XDSM"]), BtHheld=str(out["held"]["XDSMh"]), BtDonors=rng([r["info"]["donors"] for r in done], str), BtNd=rng([r["info"]["median_nd"] for r in done], lambda x: f"{x:.0f}"),
           BtCells=rng([r["info"]["cells"] for r in done], lambda x: f"{x:,}".replace(",", "{,}")), BtPBcalls=rng([r["PB"]["mean_calls"] for r in done], fmt),
           BtXtI=rng([r["XDSM"]["mean_t05"] for r in done], lambda x: f"{x:.3f}"), BtHtI=rng([r["XDSMh"]["mean_t05"] for r in done], lambda x: f"{x:.3f}"), BtSha=out["plan_sha256"][:16],
           BtXfail=", ".join(esc(SHORT.get(r["info"]["collection_name"], r["info"]["collection_name"])[:30]) + " " + esc(r["info"]["cell_type"]) for r in done if not r["XDSM"]["holds"]) or "none",
           BtHfail=", ".join(esc(SHORT.get(r["info"]["collection_name"], r["info"]["collection_name"])[:30]) + " " + esc(r["info"]["cell_type"]) for r in done if not r["XDSMh"]["holds"]) or "none")
from scipy.stats import beta as _be, binom as _bi
for m, tg in (("XDSM", "X"), ("XDSMh", "H")):
    k = sum(r[m]["runs_with_call"] for r in done); n = 20 * len(done); lo = _be.ppf(0.025, k, n - k + 1) if k else 0.0; hi = _be.ppf(0.975, k + 1, n - k)
    mac[f"Bt{tg}k"], mac[f"Bt{tg}n"], mac[f"Bt{tg}pct"], mac[f"Bt{tg}ci"] = str(k), str(n), f"{100*k/n:.1f}\\%", f"{100*lo:.1f}--{100*hi:.1f}\\%"
mac["BtExpFail"] = f"{len(done) * (1 - _bi.cdf(2, 20, 0.05)):.1f}"; mac["BtOverlap"] = f"{sum(o['identical_calls'] for o in overlap)} of {sum(o['n'] for o in overlap)}"
open(os.path.join(P2, "plan10nums.tex"), "w").write("% generated by breadth/plan10_summary.py; do not edit\n" + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in mac.items()))
print(json.dumps({k: v for k, v in out.items() if k != "settings"}, indent=1)); print(mac)
for r in rows: print(r["i"], r["status"], (r.get("info") or {}).get("cell_type", r.get("cell_type")), {m: (r[m]["runs_with_call"], round(r[m]["mean_calls"], 1), round(r[m]["mean_t05"], 4), r[m]["holds"]) for m in M if m in r}, r.get("xdsm_consistency"))
