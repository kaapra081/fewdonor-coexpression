"""REV4 R6: PC1 loading signs and gene-mean versus log cell number. Run from repo/."""
import numpy as np, json
Z = np.load("rev4/out/R6_P1.npz", allow_pickle=True); out = {}
for tag in ("null", "real"):
    for prof in ("L", "R"):
        X = Z[f"{tag}_{prof}"].astype(float); Xc = X - X.mean(0); U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
        v = Vt[0]; sc = U[:, 0] * s[0]; sg = np.sign(np.corrcoef(sc, np.log(Z[f"{tag}_nd"]))[0, 1]); v = v * sg
        out[f"{tag}_{prof}"] = dict(frac_same_sign_PC1=float(max((v > 0).mean(), (v < 0).mean())), corr_PC1_rowmean=float(abs(np.corrcoef(sc, Xc.mean(1))[0, 1])), corr_rowmean_logn=float(np.corrcoef(Xc.mean(1), np.log(Z[f"{tag}_nd"]))[0, 1]))
json.dump(out, open("rev4/out/R6/R6_pc1.json", "w"), indent=1); print(json.dumps(out))
