import os, sys, glob, json, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
MS = ("PB", "PBlin", "LCP", "LCPX", "LCPXS", "LCPd", "LCPXd"); out = {"A": {}, "rep_null": {}, "B": {}}
for sub in ("P1", "P2", "P3"):
    fs = sorted(glob.glob(os.path.join(HERE, "out", f"AD_{sub}_k*_D*.npz")))
    byD = {}
    for f in fs:
        k, D = re.search(r"_k(\d+)_D(\w+)\.npz", f).groups(); byD.setdefault(D, []).append(dict(np.load(f)))
    out["A"][sub] = {}
    for D, Z in byD.items():
        r = {"K": len(Z)}
        for m in MS:
            calls = [int(bh(z[m].astype(float)).sum()) for z in Z]
            top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
            t1 = float(np.mean([(z[m] < 0.05).mean() for z in Z])); tt = np.concatenate([z[m][tp] < 0.05 for z, tp in zip(Z, top)])
            se = float(np.sqrt(0.05 * 0.95 / len(tt)))
            r[m] = dict(mean_calls=float(np.mean(calls)), calls=calls, typeI=t1, typeI_top1=float(tt.mean()), se_top1=se,
                        holds=bool(np.mean(calls) <= 1 and t1 <= 0.055 and tt.mean() <= 0.05 + 3 * se))
        out["A"][sub][D] = r
    # replication of false calls across pseudo-donor halves
    rr = {}
    for f0 in sorted(glob.glob(os.path.join(HERE, "out", f"ADhalf_{sub}_k*_h0.npz"))):
        h0 = dict(np.load(f0)); h1 = dict(np.load(f0.replace("_h0", "_h1"))); k = re.search(r"_k(\d+)_", f0).group(1)
        for m, sg in (("PB", "sign_PB"), ("LCP", "sign_LCP"), ("LCPXd", None), ("LCPd", "sign_LCPd")):
            c0 = bh(h0[m].astype(float)); same = (h1[m] < 0.05) if sg is None else ((h1[m] < 0.05) & (h1[sg] == h0[sg]))
            rr.setdefault(m, []).append(int((c0 & same).sum()))
    out["rep_null"][sub] = {m: dict(replicated_false=v, mean=float(np.mean(v))) for m, v in rr.items()}
    ft = os.path.join(HERE, "out", f"B_{sub}_full.npz")
    if os.path.exists(ft):
        z = dict(np.load(ft)); b = {m: int(bh(z[m].astype(float)).sum()) for m in MS}
        cX = bh(z["LCPXd"].astype(float)); cP = bh(z["PB"].astype(float)); rw = np.abs(z["r_w"])
        b["PB_and_LCPXd"] = int((cX & cP).sum()); b["PB_only"] = int((cP & ~cX).sum()); b["LCPXd_only"] = int((cX & ~cP).sum())
        b["median_rw_PBonly"] = float(np.nanmedian(rw[cP & ~cX])) if (cP & ~cX).any() else None
        b["median_rw_shared"] = float(np.nanmedian(rw[cP & cX])) if (cP & cX).any() else None
        h0f = ft.replace("_full", "_h0")
        if os.path.exists(h0f) and os.path.exists(ft.replace("_full", "_h1")):
            h0 = dict(np.load(h0f)); h1 = dict(np.load(ft.replace("_full", "_h1")))
            for m, sg in (("PB", "sign_PB"), ("LCPXd", None), ("LCP", "sign_LCP")):
                c0 = bh(h0[m].astype(float)); same = (h1[m] < 0.05) if sg is None else ((h1[m] < 0.05) & (h1[sg] == h0[sg]))
                b[f"rep_{m}"] = dict(calls_h0=int(c0.sum()), replicated=int((c0 & same).sum()))
        out["B"][sub] = b
json.dump(out, open(os.path.join(HERE, "pop_results.json"), "w"), indent=1)
for sub, r in out["A"].items():
    for D, v in sorted(r.items(), key=lambda x: (x[0] != "all", -int(x[0]) if x[0] != "all" else 0)):
        print(sub, D, "K", v["K"], {m: (round(v[m]["mean_calls"]), round(v[m]["typeI"], 3), round(v[m]["typeI_top1"], 3)) for m in MS})
print(json.dumps(out["rep_null"])); print(json.dumps(out["B"], indent=0))
