"""Per-replicate metrics and design summaries with MCSE (DEV_PREREG.md section 5).
python analyze_dev.py --runs runs_dev --out dev_results   ->  per_rep.csv, summary.csv, coverage.csv"""
import os, sys, glob, argparse
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from methods_v2 import bh  # noqa: E402

TESTS = ["LCP", "LCP_stud", "CMB_mom", "CMB_naive", "PB_logcpm", "E5", "ORACLE", "BOUND"]
NULLS = ["WW", "SS", "HH", "NN", "QQ_cross"]
ALTS = ["BB", "MM", "OO_same", "OO_opp", "QQ_same"]
INTERVALS = {"EQ": ("EQ_lo", "EQ_hi"), "CMB_jack": ("CMB_lo", "CMB_hi")}


def per_rep(design_dir):
    tr = np.load(os.path.join(design_dir, "truth.npz"))
    cls, alt, sb = tr["cls"], tr["altB"], tr["sb"]
    null = ~alt
    rows, crow = [], []
    for fn in sorted(glob.glob(os.path.join(design_dir, "rep*.npz"))):
        z = np.load(fn); rep = int(os.path.basename(fn)[3:7])
        for m in TESTS:
            p = z[m].astype(float); rej = bh(p)
            nd = rej.sum(); fd = (rej & null).sum()
            r = dict(rep=rep, method=m, n_disc=int(nd), fdp=fd / max(nd, 1), power=(rej & alt).sum() / alt.sum())
            for c in NULLS + ALTS:
                r["rej05_" + c] = float((p[cls == c] < 0.05).mean())
            rows.append(r)
        for nm, (lo, hi) in list(INTERVALS.items()) + ([("BOOT_T", ("BT_lo", "BT_hi")), ("ORACLE_BOOT_T", ("OBT_lo", "OBT_hi"))] if "BT_lo" in z.files else []):
            cov = (z[lo] <= sb) & (sb <= z[hi])
            r = dict(rep=rep, interval=nm)
            for c in NULLS + ALTS:
                r["cov_" + c] = float(cov[cls == c].mean())
            crow.append(r)
        rows[-1]["runtime_total_s"] = float(z["runtime_total_s"])
        rows[-1]["runtime_LCP_s"] = float(z["runtime_LCP_s"])
    return pd.DataFrame(rows), pd.DataFrame(crow)


def summarize(df, cov):
    out = []
    for m, g in df.groupby("method"):
        R = len(g)
        r = dict(method=m, R=R)
        for k in ["fdp", "power"] + [c for c in g.columns if c.startswith("rej05_")]:
            r[k] = g[k].mean(); r[k + "_mcse"] = g[k].std(ddof=1) / np.sqrt(R)
        out.append(r)
    s = pd.DataFrame(out).rename(columns={"fdp": "FDR", "fdp_mcse": "FDR_mcse"})
    cs = []
    for nm, g in cov.groupby("interval"):
        R = len(g); r = dict(interval=nm, R=R)
        for k in [c for c in g.columns if c.startswith("cov_")]:
            r[k] = g[k].mean(); r[k + "_mcse"] = g[k].std(ddof=1) / np.sqrt(R)
        cs.append(r)
    return s, pd.DataFrame(cs)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", default="runs_dev"); ap.add_argument("--out", default="dev_results")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    P, C, Sm, Cs = [], [], [], []
    for dd in sorted(glob.glob(os.path.join(a.runs, "*"))):
        if not os.path.isdir(dd) or not glob.glob(os.path.join(dd, "rep*.npz")):
            continue
        name = os.path.basename(dd)
        df, cov = per_rep(dd); df.insert(0, "design", name); cov.insert(0, "design", name)
        s, cs = summarize(df.drop(columns="design"), cov.drop(columns="design"))
        s.insert(0, "design", name); cs.insert(0, "design", name)
        P.append(df); C.append(cov); Sm.append(s); Cs.append(cs)
    pd.concat(P).to_csv(os.path.join(a.out, "per_rep.csv"), index=False)
    pd.concat(C).to_csv(os.path.join(a.out, "per_rep_coverage.csv"), index=False)
    pd.concat(Sm).to_csv(os.path.join(a.out, "summary.csv"), index=False)
    pd.concat(Cs).to_csv(os.path.join(a.out, "coverage.csv"), index=False)
    print(pd.concat(Sm)[["design", "method", "R", "FDR", "FDR_mcse", "power", "power_mcse"]].to_string(index=False))


if __name__ == "__main__":
    main()
