import os, sys, glob, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
out = {"sim": {}, "null": {}}
for od in sorted(glob.glob(os.path.join(HERE, "runs", "F_*"))):
    dn = os.path.basename(od); fx = os.path.join(HERE, "..", "fix", "runs", dn)
    tr = np.load(os.path.join(fx, "truth.npz")); alt = tr["altB"]; cls = tr["cls"].astype(str); wn = cls == "QQ_cross"
    files = sorted(glob.glob(os.path.join(od, "rep*.npz"))); res = {}
    reps = [(dict(np.load(f)), dict(np.load(os.path.join(fx, os.path.basename(f))))) for f in files]
    for m in ("LCPX", "LCPXS", "LCP", "LCP2", "PB"):
        f_, p_, w_ = [], [], []
        for x, y in reps:
            p = (x if m in x else y)[m].astype(float); rj = bh(p); n = rj.sum()
            f_.append((rj & ~alt).sum() / max(n, 1)); p_.append((rj & alt).sum() / alt.sum()); w_.append((p[wn] < 0.05).mean())
        R = len(reps); se = lambda v: float(np.std(v, ddof=1) / np.sqrt(R))
        res[m] = dict(FDR=float(np.mean(f_)), FDR_mcse=se(f_), power=float(np.mean(p_)), power_mcse=se(p_), WN=float(np.mean(w_)), WN_mcse=se(w_), R=R)
    out["sim"][dn] = res
for st in ("S1", "S2"):
    Z = [dict(np.load(f)) for f in sorted(glob.glob(os.path.join(HERE, f"xnull_{st}_k*.npz")))]
    if not Z: continue
    r = {}
    for m in ("LCPX", "LCPXS"):
        calls = [int(bh(z[f"p_{m}"].astype(float)).sum()) for z in Z]
        t1 = np.mean([(z[f"p_{m}"] < 0.05).mean() for z in Z])
        top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
        t1top = np.concatenate([z[f"p_{m}"][tp] < 0.05 for z, tp in zip(Z, top)])
        se = float(np.sqrt(0.05 * 0.95 / len(t1top)))
        r[m] = dict(mean_calls=float(np.mean(calls)), calls=calls, typeI=float(t1), typeI_top1=float(t1top.mean()), typeI_top1_se=se,
                    holds=bool(t1 <= 0.055 and t1top.mean() <= 0.05 + 3 * se))
    r["K"] = len(Z); out["null"][st] = r
S = out["sim"]
if len(S) == 6 and len(out["null"]) == 2:
    ok = lambda d, m: S[d][m]["FDR"] <= 0.05 + 2 * S[d][m]["FDR_mcse"]
    A = all(ok(d, "LCPX") for d in S); Bn = all(out["null"][s]["LCPX"]["holds"] for s in ("S1", "S2"))
    C = all(S[d]["LCPX"]["power"] >= 0.8 * S[d]["LCP"]["power"] for d in ("F_D24", "F_D40"))
    W = all(S[d]["LCPXS"]["WN"] <= 0.05 + 2 * S[d]["LCPXS"]["WN_mcse"] and ok(d, "LCPXS") for d in ("F_D24", "F_D40", "F_S1", "F_S2"))
    out["decision"] = dict(a=A, b=Bn, c=C, adopt=bool(A and Bn and C), weak_null_claim=bool(W))
json.dump(out, open(os.path.join(HERE, "xs_results.json"), "w"), indent=1)
for d, r in S.items(): print(d, r["LCPX"]["R"], {m: (round(v["FDR"], 3), round(v["FDR_mcse"], 3), round(v["power"], 3), round(v["WN"], 3)) for m, v in r.items()})
print(json.dumps(out["null"])); print(out.get("decision"))
