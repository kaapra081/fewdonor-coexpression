"""REV3_PREREG A2/A2x summary: per data set, null, D, method, stratum: pooled observed/expected tail ratios, run-level BH, frozen decision."""
import os, sys, glob, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); OD = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "out")
THR = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]; MS = ("PB", "XDS", "XDSM1", "XDSM5", "XDSM10")
bound = lambda n: 0.05 + 2 * math.sqrt(0.05 * 0.95 / n)
out = {}
for ds in ("B1", "B2", "B1x", "P3"):
    for nl in (0, 1):
        fs = sorted(glob.glob(os.path.join(OD, f"A2_{ds}_n{nl}_k*.json")))
        R = [json.load(open(f)) for f in fs]; R = [r for r in R if r.get("valid", True)]
        if not R: continue
        key = f"{ds}_n{nl}"; out[key] = {"runs": len(R), "invalid": len(fs) - len(R), "D_all": R[0]["D_all"], "maxcheck": max(max(r["check"].values()) for r in R)}
        for D in ("200", "all"):
            for m in MS:
                for si in map(str, range(5)):
                    st = [r[D][m]["strata"][si] for r in R]; npair = st[0]["n"]
                    if npair == 0: continue
                    exc = [sum(s["exc"][i] for s in st) for i in range(5)]; pos = [sum(s["pos"][i] for s in st) for i in range(5)]
                    expc = [npair * len(R) * t for t in THR]
                    ratio = [e / x for e, x in zip(exc, expc)]
                    anyV = sum(s["V"] >= 1 for s in st); runok = anyV / len(R) <= bound(len(R))
                    ass4, ass5 = expc[1] >= 20, expc[2] >= 20
                    tail_ok = (ratio[1] <= 1.5 if ass4 else None, ratio[2] <= 1.5 if ass5 else None)
                    dec = None if (tail_ok[0] is None) else (tail_ok[0] and (tail_ok[1] in (True, None)) and runok)
                    out[key].setdefault(D, {}).setdefault(m, {})[si] = dict(pairs=npair, exc=exc, pos=pos, exp=[round(x, 1) for x in expc], ratio=[round(x, 3) for x in ratio],
                        anyV=anyV, runs=len(R), run_ok=runok, calibrated=dec, assessable_1e5=ass5, meanVall=sum(s["Vall"] for s in st) / len(R))
                out[key][D][m]["Vtotal_mean"] = sum(r[D][m]["V"] for r in R) / len(R); out[key][D][m]["anyV_total"] = sum(r[D][m]["V"] >= 1 for r in R)
json.dump(out, open(os.path.join(OD, "A2_SUMMARY.json"), "w"), indent=1)
for k, v in out.items():
    print(k, "runs", v["runs"], "D", v["D_all"], "maxcheck %.1e" % v["maxcheck"])
    for D in ("200", "all"):
        for m in ("PB", "XDS", "XDSM1", "XDSM5", "XDSM10"):
            x = v[D][m]; print(" ", D, m, "V=%.1f any=%d" % (x["Vtotal_mean"], x["anyV_total"]), " | ".join(f"s{si}:r4={x[si]['ratio'][1]:.2f},r5={x[si]['ratio'][2]:.2f},any={x[si]['anyV']},{'OK' if x[si]['calibrated'] else ('NA' if x[si]['calibrated'] is None else 'FAIL')}" for si in map(str, range(5)) if si in x))
