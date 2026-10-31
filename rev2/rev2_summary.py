"""Summarize REV2 outputs into one compact JSON (runs on the compute host). python rev2_summary.py"""
import os, sys, glob, json, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
O = os.path.join(HERE, "out"); S = {"A": {}, "C": {}, "D1": {}, "D3": {}, "D2": None, "E": {}, "F": {}}
def runlevel(files, getter):
    agg = {}
    for f in files:
        for key, m, V, t in getter(json.load(open(f))):
            a = agg.setdefault(key, {}).setdefault(m, dict(n=0, any=0, V=[], t05=[])); a["n"] += 1; a["any"] += int(V >= 1); a["V"].append(V); a["t05"].append(t)
    for key in agg:
        for m, a in agg[key].items():
            a["meanV"] = float(np.mean(a["V"])); a["medV"] = float(np.median(a["V"])); a["maxV"] = int(np.max(a["V"]))
            a["mean_t05"] = float(np.mean(a["t05"])); a["sd_t05"] = float(np.std(a["t05"], ddof=1)) if len(a["t05"]) > 1 else None; del a["V"], a["t05"]
    return agg
for s in ("P1", "P2", "P3"):
    fs = sorted(glob.glob(os.path.join(O, f"A_{s}_k*.json")))
    if fs: S["A"][s] = runlevel(fs, lambda r: [(D, m, r[D][m]["V"], r[D][m]["t05"]) for D in r for m in r[D]])
    fc = os.path.join(O, f"C_{s}.json")
    if os.path.exists(fc):
        c = json.load(open(fc)); sx = c.pop("sexpairs", [])
        xy = [p for p in sx if p["xist"]]; yy = [p for p in sx if not p["xist"]]
        def cnt(L, sgn, key_p, key_s):
            return dict(n=len(L), sig_right_sign=int(sum((p[key_p] < 0.05) and ((p[key_s] < 0) if sgn < 0 else (p[key_s] > 0)) for p in L)), max_p=float(max([p[key_p] for p in L], default=np.nan)))
        c["sex_XDSM"] = dict(xist_y=cnt(xy, -1, "p_XDSM", "s_XDSM"), y_y=cnt(yy, 1, "p_XDSM", "s_XDSM"))
        c["sex_PB"] = dict(xist_y=cnt(xy, -1, "p_PB", "r_PB"), y_y=cnt(yy, 1, "p_PB", "r_PB")); S["C"][s] = c
    for f in sorted(glob.glob(os.path.join(O, f"D1_{s}_*.json"))): S["D1"].setdefault(s, {})[os.path.basename(f)[len(f"D1_{s}_"):-5]] = json.load(open(f))
    f3 = os.path.join(O, f"D3_{s}.json")
    if os.path.exists(f3): S["D3"][s] = json.load(open(f3))
if os.path.exists(os.path.join(O, "D2_summary.json")): S["D2"] = json.load(open(os.path.join(O, "D2_summary.json")))
EO = os.path.join(HERE, "..", "onek", "out")
for t in ("B1", "B2", "B3"):
    fs = sorted(glob.glob(os.path.join(EO, f"E_{t}_k*.json")))
    if fs: S["E"][t] = runlevel(fs, lambda r: [("all", m, r[m]["V"], r[m]["t05"]) for m in r])["all"]
for d in ("L200", "L500"):
    od = os.path.join(O, "F_" + d); reps = sorted(glob.glob(os.path.join(od, "rep*.npz")))
    if not reps: continue
    tr = np.load(os.path.join(od, "truth.npz")); alt = tr["altB"]; wn = tr["cls"].astype(str) == "QQ_cross"; res = {}
    Z = [dict(np.load(f)) for f in reps]
    for m in ("PB", "LCP", "LCPd", "XDS", "XDSM", "KX"):
        f_, p_, w_ = [], [], []
        for z in Z:
            p = z[m].astype(float); rj = bh(p); n = rj.sum(); f_.append((rj & ~alt).sum() / max(n, 1)); p_.append((rj & alt).sum() / alt.sum()); w_.append((p[wn] < 0.05).mean())
        R = len(Z); se = lambda v: float(np.std(v, ddof=1) / np.sqrt(R)) if R > 1 else None
        res[m] = dict(FDR=float(np.mean(f_)), FDR_mcse=se(f_), power=float(np.mean(p_)), power_mcse=se(p_), WN=float(np.mean(w_)), WN_mcse=se(w_), R=R)
    S["F"][d] = res
json.dump(S, open(os.path.join(O, "REV2_SUMMARY.json"), "w"), indent=1); print(len(json.dumps(S)))
