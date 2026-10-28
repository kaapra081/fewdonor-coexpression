"""POST HOC ablation summaries with MCSE. python analyze_ablate.py --runs runs_ablate --out ablate_results"""
import os, sys, glob, argparse
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from methods_v2 import bh  # noqa: E402
NULLS = ["WW", "SS", "HH", "NN", "QQ_cross"]; ALTS = ["BB", "MM", "OO_same", "OO_opp", "QQ_same"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", default="runs_ablate"); ap.add_argument("--out", default="ablate_results")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rows = []
    for dd in sorted(glob.glob(os.path.join(a.runs, "*"))):
        tr = np.load(os.path.join(dd, "truth.npz")); cls, alt = tr["cls"], tr["altB"]; null = ~alt
        for fn in sorted(glob.glob(os.path.join(dd, "rep*.npz"))):
            z = np.load(fn); rep = int(os.path.basename(fn)[3:7])
            for m in z.files:
                p = z[m].astype(float); rej = bh(p); nd = rej.sum()
                r = dict(design=os.path.basename(dd), rep=rep, method=m, fdp=(rej & null).sum() / max(nd, 1), power=(rej & alt).sum() / alt.sum())
                for c in NULLS + ALTS:
                    r["rej05_" + c] = float((p[cls == c] < 0.05).mean())
                rows.append(r)
    df = pd.DataFrame(rows); df.to_csv(os.path.join(a.out, "per_rep.csv"), index=False)
    out = []
    for (d, m), g in df.groupby(["design", "method"]):
        R = len(g); r = dict(design=d, method=m, R=R)
        for k in ["fdp", "power"] + [c for c in g.columns if c.startswith("rej05_")]:
            r[k] = g[k].mean(); r[k + "_mcse"] = g[k].std(ddof=1) / np.sqrt(R)
        out.append(r)
    s = pd.DataFrame(out).rename(columns={"fdp": "FDR", "fdp_mcse": "FDR_mcse"}); s.to_csv(os.path.join(a.out, "summary.csv"), index=False)
    print(s[["design", "method", "R", "FDR", "FDR_mcse", "power", "power_mcse", "rej05_SS", "rej05_NN"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
