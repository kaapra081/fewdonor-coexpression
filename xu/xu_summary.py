import os, sys, glob, json, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
out = {"E3": {}, "E4": {}, "E1": {}, "E2": {}}
def simsum(od, methods, extra=None):
    tr = np.load(os.path.join(od, "truth.npz")); alt = tr["altB"]; cls = tr["cls"].astype(str); wn = cls == "QQ_cross"
    reps = []
    for f in sorted(glob.glob(os.path.join(od, "rep*.npz"))):
        z = dict(np.load(f))
        if extra: z.update({k: v for k, v in extra(os.path.basename(f)).items()})
        reps.append(z)
    res = {}
    for m in methods:
        if m not in reps[0]: continue
        f_, p_, w_ = [], [], []
        for z in reps:
            p = z[m].astype(float); rj = bh(p); n = rj.sum()
            f_.append((rj & ~alt).sum() / max(n, 1)); p_.append((rj & alt).sum() / alt.sum()); w_.append((p[wn] < 0.05).mean())
        R = len(reps); se = lambda v: float(np.std(v, ddof=1) / np.sqrt(R))
        res[m] = dict(FDR=float(np.mean(f_)), FDR_mcse=se(f_), power=float(np.mean(p_)), power_mcse=se(p_), WN=float(np.mean(w_)), WN_mcse=se(w_), R=R)
    return res
for d in ("L200", "L500", "L1000"):
    od = os.path.join(HERE, "runs", d)
    if os.path.isdir(od) and glob.glob(os.path.join(od, "rep*.npz")): out["E3"][d] = simsum(od, ("PB", "LCP", "LCPd", "XDS", "XU"))
R2 = os.path.join(HERE, "..")
for d in ("F_D16", "F_D24", "F_D40", "F_S1", "F_S2", "F_D16_n100"):
    od = os.path.join(HERE, "runs", d)
    if not (os.path.isdir(od) and glob.glob(os.path.join(od, "rep*.npz"))): continue
    def ex(fn, d=d):
        z = dict(np.load(os.path.join(R2, "fix", "runs", d, fn))); z.update({"XDS": np.load(os.path.join(R2, "dstd", "runs", d, fn))["LCPXd"]}); return z
    out["E4"][d] = simsum(od, ("PB", "LCP", "LCP2", "XDS", "XU"), ex)
for f in sorted(glob.glob(os.path.join(HERE, "out", "XU_*_D*.npz"))):
    sub, k, D = re.search(r"XU_(P\d)_k(\d+)_D(\w+)\.npz", f).groups(); out["E1"].setdefault(sub, {}).setdefault(D, []).append(dict(np.load(f)))
for sub in list(out["E1"]):
    for D, Z in list(out["E1"][sub].items()):
        r = {"K": len(Z)}
        for m in ("PB", "LCP", "XDS", "XU"):
            if m not in Z[0]: continue
            calls = [int(bh(z[m].astype(float)).sum()) for z in Z]; t1 = float(np.mean([(z[m] < 0.05).mean() for z in Z]))
            top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
            tt = np.concatenate([z[m][tp] < 0.05 for z, tp in zip(Z, top)]); se = float(np.sqrt(0.05 * 0.95 / len(tt)))
            r[m] = dict(mean_calls=float(np.mean(calls)), calls=calls, typeI=t1, typeI_top1=float(tt.mean()), se=se,
                        holds=bool(np.mean(calls) <= 1 and t1 <= 0.055 and tt.mean() <= 0.05 + 3 * se))
        out["E1"][sub][D] = r
for st in ("S1", "S2"):
    Z = [dict(np.load(f)) for f in sorted(glob.glob(os.path.join(HERE, "out", f"xunull_{st}_k*.npz")))]
    if not Z: continue
    calls = [int(bh(z["p_XU"].astype(float)).sum()) for z in Z]; t1 = float(np.mean([(z["p_XU"] < 0.05).mean() for z in Z]))
    top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
    tt = np.concatenate([z["p_XU"][tp] < 0.05 for z, tp in zip(Z, top)]); se = float(np.sqrt(0.05 * 0.95 / len(tt)))
    out["E2"][st] = dict(K=len(Z), mean_calls=float(np.mean(calls)), calls=calls, typeI=t1, typeI_top1=float(tt.mean()), se=se,
                         holds=bool(t1 <= 0.055 and tt.mean() <= 0.05 + 3 * se))
E1 = out["E1"]; E3 = out["E3"]
if len(E1) == 3 and len(E3) == 3:
    h = all(E1[s][D]["XU"]["holds"] for s in E1 for D in E1[s] if D == "all" or int(D) >= 200)
    ok = lambda d, m: E3[d][m]["FDR"] <= 0.05 + 2 * E3[d][m]["FDR_mcse"]
    rec = h and all(ok(d, "XU") for d in E3) and all(E3[d]["XU"]["power"] >= 1.5 * E3[d]["XDS"]["power"] for d in ("L500", "L1000"))
    wn = all(E3[d]["XU"]["WN"] <= 0.05 + 2 * E3[d]["XU"]["WN_mcse"] for d in E3)
    out["decision"] = dict(holds_pop=h, recommended=bool(rec), weak_null=bool(wn))
json.dump(out, open(os.path.join(HERE, "xu_results.json"), "w"), indent=1)
for k in ("E3", "E4"):
    for d, r in out[k].items(): print(k, d, r[list(r)[0]]["R"], {m: (round(v["FDR"], 3), round(v["FDR_mcse"], 3), round(v["power"], 3), round(v["WN"], 3)) for m, v in r.items()})
for s, r in E1.items():
    for D, v in r.items(): print("E1", s, D, v["K"], {m: (round(v[m]["mean_calls"], 1), round(v[m]["typeI"], 3), round(v[m]["typeI_top1"], 3)) for m in v if m != "K"})
print("E2", json.dumps(out["E2"])); print(out.get("decision"))
