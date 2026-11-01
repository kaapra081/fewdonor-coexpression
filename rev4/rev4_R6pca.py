"""REV4 R6: PCA coordinates (PC1, PC2; PC1 sign = large for few cells) appended to donor tables. Run from repo/."""
import numpy as np, pandas as pd, os
Z = np.load("rev4/out/R6_P1.npz", allow_pickle=True)
for tag in ("null", "real"):
    df = pd.read_csv(f"rev4/out/R6/donors_{tag}.csv")
    for prof in ("L", "R"):
        X = Z[f"{tag}_{prof}"].astype(float); Xc = X - X.mean(0); U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
        sc = U[:, :2] * s[:2]
        for i in range(2):
            v = sc[:, i]; v = v if np.corrcoef(v, np.log(df.cells))[0, 1] <= 0 or i else -v   # PC1 sign: large for few cells
            df[f"pc{i+1}_{prof}"] = np.round(v, 5)
    df.to_csv(f"rev4/out/R6/donors_{tag}.csv", index=False)
print("ok")
