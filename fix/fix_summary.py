import os, sys, glob, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
out = {"null": {}, "sim": {}}
for st in ("S1", "S2"):
    Z = [dict(np.load(f)) for f in sorted(glob.glob(os.path.join(HERE, f"null_{st}_k*.npz")))]
    r = {}
    for m in ("PB", "LCP", "LCP2"):
        calls = [int(bh(z[f"p_{m}"].astype(float)).sum()) for z in Z]
        t1 = np.mean([(z[f"p_{m}"] < 0.05).mean() for z in Z])
        top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
        t1top = np.concatenate([z[f"p_{m}"][tp] < 0.05 for z, tp in zip(Z, top)])
        r[m] = dict(mean_calls=float(np.mean(calls)), calls=calls, frac_any=float(np.mean([c > 0 for c in calls])),
                    typeI=float(t1), typeI_top1=float(t1top.mean()), typeI_top1_se=float(np.sqrt(0.05 * 0.95 / len(t1top))))
    r["decision_leak_real"] = "YES" if (r["PB"]["mean_calls"] >= 50 and r["PB"]["mean_calls"] >= 10 * max(r["LCP"]["mean_calls"], 1)) else "NO"
    for m in ("LCP", "LCP2"):
        r[f"{m}_holds"] = bool(r[m]["typeI"] <= 0.055 and r[m]["typeI_top1"] <= 0.05 + 3 * r[m]["typeI_top1_se"])
    r["K"] = len(Z); out["null"][st] = r
for od in sorted(glob.glob(os.path.join(HERE, "runs", "F_*"))):
    tr = np.load(os.path.join(od, "truth.npz")); alt = tr["altB"]; cls = tr["cls"]; res = {}
    reps = [np.load(f) for f in sorted(glob.glob(os.path.join(od, "rep*.npz")))]
    for m in ("LCP", "LCP2", "PB"):
        f_, p_, s_ = [], [], []
        for z in reps:
            p = z[m].astype(float); rj = bh(p); n = rj.sum()
            f_.append((rj & ~alt).sum() / max(n, 1)); p_.append((rj & alt).sum() / alt.sum()); s_.append((p[cls == "SS"] < 0.05).mean())
        R = len(reps); se = lambda v: float(np.std(v, ddof=1) / np.sqrt(R))
        res[m] = dict(FDR=float(np.mean(f_)), FDR_mcse=se(f_), power=float(np.mean(p_)), power_mcse=se(p_), SS=float(np.mean(s_)), R=R)
    out["sim"][os.path.basename(od)] = res
S = out["sim"]; ok = lambda d, m: S[d][m]["FDR"] <= 0.05 + 2 * S[d][m]["FDR_mcse"]
a = all(ok(d, "LCP2") for d in ("F_D24", "F_D40", "F_S2", "F_D16_n100")); b = ok("F_D16", "LCP2") or ok("F_S1", "LCP2")
c = all(S[d]["LCP2"]["power"] >= 0.9 * S[d]["LCP"]["power"] for d in ("F_D24", "F_D40"))
out["LCP2_decision"] = dict(a=a, b=b, c=c, adopt=bool(a and b and c))
json.dump(out, open(os.path.join(HERE, "fix_results.json"), "w"), indent=1)
print(json.dumps(out["null"], indent=0)[:2500])
for d, r in S.items(): print(d, {m: (round(v["FDR"], 3), round(v["FDR_mcse"], 3), round(v["power"], 3), round(v["SS"], 3)) for m, v in r.items()})
print(out["LCP2_decision"])
