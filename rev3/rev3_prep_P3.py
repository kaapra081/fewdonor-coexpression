"""REV3_PREREG A2-P3: ADAM_GPNMB raw counts for genes ranked 1-2,000 by mean CP10k on a 5,000-nucleus random subsample (seed 46,000,000),
donors with >= 4 nuclei (as POP_PREREG). Writes rev3/data/a2_P3_{X,meta}.npz and a2_P3_info.json."""
import os, sys, json, numpy as np, anndata as ad, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); OD = os.path.join(HERE, "data"); os.makedirs(OD, exist_ok=True)
DATA = os.environ.get("FEWDONOR_DATA", "/Users/kaavinprasanna/Documents/recomb27/data/")
a = ad.read_h5ad(os.path.join(DATA, "pfc_myeloid_1469.h5ad"), backed="r"); o = a.obs
sym = np.asarray(a.raw.var["feature_name"] if "feature_name" in a.raw.var else a.var["feature_name"]).astype(str)
rows = np.where(o["subtype"].astype(str).values == "ADAM_GPNMB")[0]; don = o["donor_id"].astype(str).values[rows]
u, c = np.unique(don, return_counts=True); keep = np.isin(don, u[c >= 4]); rows, don = rows[keep], don[keep]
rng = np.random.default_rng(46_000_000); sub = np.sort(rng.choice(len(rows), size=min(5000, len(rows)), replace=False))
blocks, lib = [], np.zeros(len(rows))
for i in range(0, len(rows), 3000):
    X = sp.csr_matrix(a.raw.X[rows[i:i + 3000]]); lib[i:i + 3000] = np.asarray(X.sum(1)).ravel(); blocks.append(X)
Xall = sp.vstack(blocks).tocsr(); del blocks
Xs = Xall[sub]; cp = np.asarray(Xs.multiply(1e4 / lib[sub][:, None]).mean(0)).ravel()
top = np.sort(np.argsort(-cp)[:2000])
Xp = Xall[:, top].tocsr().astype(np.float32); umi = np.asarray(Xp.mean(0)).ravel()
sp.save_npz(os.path.join(OD, "a2_P3_X.npz"), Xp)
np.savez_compressed(os.path.join(OD, "a2_P3_meta.npz"), lib=lib, donor=don, sym=sym[top], umi=umi)
info = dict(subtype="ADAM_GPNMB", cells=int(len(rows)), donors=int(len(u[c >= 4])), G=int(len(top)),
            strata=[int(((umi > lo) & (umi <= hi)).sum()) for lo, hi in ((1, 1e9), (0.3, 1), (0.1, 0.3), (0.03, 0.1), (-1, 0.03))])
json.dump(info, open(os.path.join(OD, "a2_P3_info.json"), "w"), indent=1); print(info)
