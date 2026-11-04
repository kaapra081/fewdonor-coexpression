"""PLAN7: evaluate every pre-specified result against its rule. python lupus_summary.py <dir with result JSONs> -> lupus_summary.json"""
import os, sys, json, glob, numpy as np
OD = sys.argv[1] if len(sys.argv) > 1 else "out"; SUBS = ("LM", "FC"); NRUN = 20
holds = lambda calls, t05: bool(sum(c > 0 for c in calls) <= 2 and np.mean(t05) <= 0.055)


def load(mode, sub):
    fs = sorted(glob.glob(os.path.join(OD, f"{mode}_{sub}_k*.json")), key=lambda f: int(f.rsplit("_k", 1)[1][:-5])); return [json.load(open(f)) for f in fs]


def block(R, methods, nested=None):
    out = {}
    for m in methods:
        rr = [(r[nested][m] if nested else r[m]) for r in R]; c = [x["calls"] for x in rr]; t = [x["t05"] for x in rr]
        out[m] = dict(n_runs=len(R), mean_calls=float(np.mean(c)), sd_calls=float(np.std(c, ddof=1)) if len(c) > 1 else None, min_calls=int(min(c)), max_calls=int(max(c)), runs_ge1=int(sum(x > 0 for x in c)),
                      mean_t05=float(np.mean(t)), holds=holds(c, t), complete=bool(len(R) >= NRUN))
    return out


S = dict(plan_sha256=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "PLAN7.sha256")).read().split()[0] if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "PLAN7.sha256")) else None)
for sub in SUBS:
    d = dict()
    RA = load("A", sub); RA2 = load("A2", sub); RB = load("B", sub); RD = load("D", sub); rf = os.path.join(OD, f"R_{sub}.json"); RR = json.load(open(rf)) if os.path.exists(rf) else None
    if RA: d["A"] = dict(D=RA[0]["D"], nd_median=RA[0]["nd_median"], methods=block(RA, ("PB", "PBw", "KX", "XDSM")))
    if RA2:
        d["A2"] = dict(D=RA2[0]["D"], methods=block(RA2, ("PB", "PBw", "KX", "XDSM")))
        if RA: d["A2"]["ratio_PB_full_over_100"] = d["A"]["methods"]["PB"]["mean_calls"] / max(d["A2"]["methods"]["PB"]["mean_calls"], 1e-9)
    if RB:
        d["B"] = dict(n1=RB[0]["n1"], n0=RB[0]["n0"], methods=block(RB, ("PBDC", "XDSMgu")))
        if RR:
            for m, nm in (("PBDC", "rho_PB"), ("XDSMgu", "rho_XDS")):
                real = RR["C2"]["calls"][m]; d["B"][nm] = (d["B"]["methods"][m]["mean_calls"] / real) if real > 0 else None; d["B"]["real_calls_" + m] = real
    if RD:
        d["D"] = dict(D=RD[0]["D"], D_lib_mean=float(np.mean([r["D_lib"] for r in RD])), D_lib_min=int(min(r["D_lib"] for r in RD)), methods=block(RD, ("PB", "PBw", "KX", "XDSM", "XDSMlib")), DC=block(RD, ("PBDC", "XDSMgu"), nested="DC"))
    if RR: d["R"] = RR
    S[sub] = d
r = {}
if all("A" in S[s] for s in SUBS):
    xh = {s: S[s]["A"]["methods"]["XDSM"]["holds"] for s in SUBS}; pf = {s: not S[s]["A"]["methods"]["PB"]["holds"] for s in SUBS}
    r["A-a"] = dict(XDSM_holds=xh, outcome=bool(all(xh.values()))); r["A-b"] = dict(PB_fails=pf)
    r["A-c"] = dict(replicates=bool(all(xh.values()) and any(pf.values())), fully=bool(all(xh.values()) and all(pf.values())))
    r["A-d"] = dict(t05_PB={s: S[s]["A"]["methods"]["PB"]["mean_t05"] for s in SUBS}, FC_gt_LM=bool(S["FC"]["A"]["methods"]["PB"]["mean_t05"] > S["LM"]["A"]["methods"]["PB"]["mean_t05"]))
    r["A-other"] = {s: {m: S[s]["A"]["methods"][m]["holds"] for m in ("PBw", "KX")} for s in SUBS}
if all("B" in S[s] for s in SUBS):
    xg = {s: S[s]["B"]["methods"]["XDSMgu"]["holds"] for s in SUBS}
    r["B-a"] = dict(XDSMgu_holds=xg, outcome=bool(all(xg.values()))); r["B-b"] = dict(PBDC_fails={s: not S[s]["B"]["methods"]["PBDC"]["holds"] for s in SUBS})
    r["B-c"] = dict(rho_PB={s: S[s]["B"].get("rho_PB") for s in SUBS}, not_interpretable={s: (S[s]["B"].get("rho_PB") is not None and S[s]["B"]["rho_PB"] >= 0.5) for s in SUBS})
if all("D" in S[s] for s in SUBS):
    xm = {s: S[s]["D"]["methods"]["XDSM"]["holds"] for s in SUBS}; xl = {s: S[s]["D"]["methods"]["XDSMlib"]["holds"] for s in SUBS}
    r["D-a"] = dict(XDSM_holds=xm, outcome=bool(all(xm.values()))); r["D-b"] = dict(XDSMlib_holds=xl, outcome=bool(all(xl.values())))
    r["D-c"] = dict(applies=bool(not all(xm.values())), XDSMlib_recommended=bool((not all(xm.values())) and all(xl.values())))
    r["D-other"] = {s: {m: S[s]["D"]["methods"][m]["holds"] for m in ("PB", "PBw", "KX")} for s in SUBS}
S["rules"] = r; json.dump(S, open("lupus_summary.json", "w"), indent=1); print(json.dumps(r, indent=1))
