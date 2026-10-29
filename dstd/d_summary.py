import os, sys, glob, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
out = {"sim": {}, "null": {}}
for od in sorted(glob.glob(os.path.join(HERE, "runs", "F_*"))):
    dn = os.path.basename(od); fx = os.path.join(HERE, "..", "fix", "runs", dn); xs = os.path.join(HERE, "..", "xsplit", "runs", dn)
    tr = np.load(os.path.join(fx, "truth.npz")); alt = tr["altB"]; cls = tr["cls"].astype(str); wn = cls == "QQ_cross"
    files = sorted(glob.glob(os.path.join(od, "rep*.npz"))); res = {}
    reps = []
    for f in files:
        z = dict(np.load(f)); z.update(dict(np.load(os.path.join(fx, os.path.basename(f))))); z.update(dict(np.load(os.path.join(xs, os.path.basename(f))))); reps.append(z)
    for m in ("LCPd", "LCPXd", "PBlin", "LCP", "LCP2", "LCPX", "PB"):
        f_, p_, w_ = [], [], []
        for z in reps:
            p = z[m].astype(float); rj = bh(p); n = rj.sum()
            f_.append((rj & ~alt).sum() / max(n, 1)); p_.append((rj & alt).sum() / alt.sum()); w_.append((p[wn] < 0.05).mean())
        R = len(reps); se = lambda v: float(np.std(v, ddof=1) / np.sqrt(R))
        res[m] = dict(FDR=float(np.mean(f_)), FDR_mcse=se(f_), power=float(np.mean(p_)), power_mcse=se(p_), WN=float(np.mean(w_)), R=R)
    out["sim"][dn] = res
for st in ("S1", "S2"):
    Z = [dict(np.load(f)) for f in sorted(glob.glob(os.path.join(HERE, f"dnull_{st}_k*.npz")))]
    if not Z: continue
    r = {}
    for m in ("LCPd", "LCPXd", "PBlin"):
        calls = [int(bh(z[f"p_{m}"].astype(float)).sum()) for z in Z]
        t1 = np.mean([(z[f"p_{m}"] < 0.05).mean() for z in Z])
        top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
        t1top = np.concatenate([z[f"p_{m}"][tp] < 0.05 for z, tp in zip(Z, top)]); se = float(np.sqrt(0.05 * 0.95 / len(t1top)))
        r[m] = dict(mean_calls=float(np.mean(calls)), calls=calls, typeI=float(t1), typeI_top1=float(t1top.mean()), typeI_top1_se=se,
                    holds=bool(t1 <= 0.055 and t1top.mean() <= 0.05 + 3 * se))
    r["K"] = len(Z); out["null"][st] = r
S = out["sim"]; ok = lambda d, m: S[d][m]["FDR"] <= 0.05 + 2 * S[d][m]["FDR_mcse"]
for m in ("LCPd", "LCPXd"):
    a = all(ok(d, m) for d in ("F_D24", "F_D40", "F_S2", "F_D16_n100"))
    c = all(S[d][m]["power"] >= 0.9 * S[d]["LCP"]["power"] for d in ("F_D24", "F_D40"))
    b = all(out["null"].get(s, {}).get(m, {}).get("holds", None) for s in ("S1", "S2")) if len(out["null"]) == 2 else None
    out[f"decision_{m}"] = dict(a=a, b=b, c=c, all_six_FDR_ok=all(ok(d, m) for d in S))
json.dump(out, open(os.path.join(HERE, "d_results.json"), "w"), indent=1)
for d, r in S.items(): print(d, {m: (round(v["FDR"], 3), round(v["FDR_mcse"], 3), round(v["power"], 3)) for m, v in r.items()})
print(json.dumps(out["null"])); print({k: v for k, v in out.items() if k.startswith("decision")})
