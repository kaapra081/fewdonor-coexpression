"""Plan 8 summary: python rev5/plan8_summary.py <dir with f_*.json>  -> rev5/out/plan8_summary.json (rules exactly as rev5/PLAN8.md).
Reproduction checks: XDS-M against Plan 7 (lupus/results/lupus_results_plan7.tgz), Plan 5 R2 (rev4/out/R2_*.json), Plan 6d S (rev5/out/d_S_*.json)."""
import os, sys, json, glob, tarfile, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, ".."); D = sys.argv[1]
J = lambda f: json.load(open(f))
out = {"plan": "rev5/PLAN8.md", "plan_sha256": open(os.path.join(HERE, "PLAN8.sha256")).read().split()[0]}
tf = tarfile.open(os.path.join(R, "lupus", "results", "lupus_results_plan7.tgz"))
p7 = {os.path.basename(m.name): json.load(tf.extractfile(m)) for m in tf.getmembers() if m.name.endswith(".json") and not os.path.basename(m.name).startswith("._")}


def runs_summary(rows, methods, max_runs, key="all"):
    o = {}
    for m in methods:
        calls = np.array([r[key][m]["calls"] for r in rows]); t05 = np.array([r[key][m]["t05"] for r in rows])
        o[m] = dict(n_runs=len(rows), runs_with_call=int((calls > 0).sum()), mean_calls=float(calls.mean()), max_calls=int(calls.max()), mean_t05=float(t05.mean()),
                    holds=bool((calls > 0).sum() <= max_runs and t05.mean() <= 0.055))
    return o


# 8A lupus nulls (rule: <= 2 of 20 runs with a call and mean per-pair type-I <= 0.055)
out["8A"] = {}
for mode, p7m in (("LA", "A"), ("LD", "D")):
    for sub in ("LM", "FC"):
        fs = sorted(glob.glob(os.path.join(D, f"f_{mode}_{sub}_k*.json")), key=lambda f: int(f.split("_k")[-1][:-5])); rows = [J(f) for f in fs]
        if not rows: continue
        s = runs_summary(rows, ("XDSM", "XDSMh"), 2)
        rep = [(r["all"]["XDSM"]["calls"], p7[f"{p7m}_{sub}_k{int(f.split('_k')[-1][:-5])}.json"]["XDSM"]["calls"]) for r, f in zip(rows, fs)]
        s["repro_XDSM_vs_plan7"] = dict(identical_runs=int(sum(a == b for a, b in rep)), n=len(rep))
        out["8A"][f"{mode}_{sub}"] = s
for sub in ("LM", "FC"):
    f = os.path.join(D, f"f_LR_{sub}.json")
    if os.path.exists(f): out["8A"][f"LR_{sub}"] = J(f)
out["8A"]["held_out_donor_centred_both"] = bool(all(out["8A"].get(f"LA_{s}", {}).get("XDSMh", {}).get("holds", False) for s in ("LM", "FC")))

# 8B brain-bank-stratified null (rule: <= 4 of 40 runs with a call and mean per-pair type-I <= 0.055)
out["8B"] = {}
for sub in ("P1", "P2", "P3"):
    fs = sorted(glob.glob(os.path.join(D, f"f_B_{sub}_k*.json")), key=lambda f: int(f.split("_k")[-1][:-5])); rows = [J(f) for f in fs]
    if not rows: continue
    s = runs_summary(rows, ("XDSM", "XDSMh", "XDSMhb"), 4)
    rep = [(r["all"]["XDSM"]["calls"], J(os.path.join(R, "rev4", "out", f"R2_{sub}_k{int(f.split('_k')[-1][:-5])}.json"))["all"]["XDSM"]["V"]) for r, f in zip(rows, fs)]
    s["repro_XDSM_vs_plan5_R2"] = dict(identical_runs=int(sum(a == b for a, b in rep)), n=len(rep)); out["8B"][sub] = s
for m in ("XDSMh", "XDSMhb"):
    out["8B"][f"{m}_holds_all3"] = bool(all(out["8B"].get(s, {}).get(m, {}).get("holds", False) for s in ("P1", "P2", "P3")))

# 8C matched-FDR power (descriptive)
out["8C"] = {}
for sub in ("P1", "P2", "P3"):
    fs = sorted(glob.glob(os.path.join(D, f"f_S_{sub}_k*.json")), key=lambda f: int(f.split("_k")[-1][:-5])); rows = [J(f) for f in fs]
    if not rows: continue
    o = {}
    for m in ("PB", "KX", "XDSM", "XDSMh"):
        g = lambda k_: np.array([r[m][k_] for r in rows])
        o[m] = {k_: dict(mean=float(g(k_).mean()), sd=float(g(k_).std(ddof=1))) for k_ in ("FDP", "power", "power_at_fdp0.05", "power_at_fdp0.10")}
    rep = []
    for r, f in zip(rows, fs):
        ref = os.path.join(R, "rev5", "out", f"d_S_{sub}_k{int(f.split('_k')[-1][:-5])}.json")
        if os.path.exists(ref): rep.append(all(r[m]["calls"] == J(ref)[m]["calls"] for m in ("XDSM", "XDSMh")))
    o["n_runs"] = len(rows); o["repro_vs_plan6d_S"] = dict(identical_runs=int(sum(rep)), n=len(rep)); out["8C"][sub] = o

json.dump(out, open(os.path.join(HERE, "out", "plan8_summary.json"), "w"), indent=1)
print(json.dumps(out, indent=1)[:6000])
