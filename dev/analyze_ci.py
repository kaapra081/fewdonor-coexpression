"""Coverage of BOOT_T and ORACLE_BOOT_T (change 1), per design and pair class, with MCSE."""
import os, glob, argparse
import numpy as np, pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("--runs", default="runs_ci"); ap.add_argument("--out", default="dev_results")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rows = []
for dd in sorted(glob.glob(os.path.join(a.runs, "*"))):
    tr = np.load(os.path.join(dd, "truth.npz")); cls, sb = tr["cls"], tr["sb"]
    for fn in sorted(glob.glob(os.path.join(dd, "rep*.npz"))):
        z = np.load(fn)
        for nm, lo, hi in (("BOOT_T", "BT_lo", "BT_hi"), ("ORACLE_BOOT_T", "OBT_lo", "OBT_hi")):
            cov = (z[lo] <= sb) & (sb <= z[hi]); w = z[hi] - z[lo]
            for c in np.unique(cls):
                rows.append(dict(design=os.path.basename(dd), interval=nm, cls=c, rep=fn[-8:-4],
                                 cov=cov[cls == c].mean(), width=np.median(w[cls == c])))
df = pd.DataFrame(rows)
s = df.groupby(["design", "interval", "cls"]).agg(R=("cov", "size"), coverage=("cov", "mean"),
                                                   mcse=("cov", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
                                                   median_width=("width", "median")).reset_index()
s.to_csv(os.path.join(a.out, "coverage_boot_t.csv"), index=False)
w = s.pivot_table(index=["design", "cls"], columns="interval", values="coverage").reset_index()
w["diff"] = w["BOOT_T"] - w["ORACLE_BOOT_T"]
w.to_csv(os.path.join(a.out, "coverage_boot_t_wide.csv"), index=False)
print(w.round(3).to_string(index=False))
