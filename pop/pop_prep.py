"""POP_PREREG.md: extract raw counts of one subtype for the 417-gene panel. python pop_prep.py P1"""
import os, sys, json, numpy as np, anndata as ad, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "real"))
from real_prep import FORCE  # noqa: E402
DATA = "/Users/kaavinprasanna/Documents/recomb27/data/"
SUB = {"P1": "Homeo_CECR2", "P2": "PVM_CD163", "P3": "ADAM_GPNMB"}[sys.argv[1]]
a = ad.read_h5ad(DATA + "pfc_myeloid_1469.h5ad", backed="r")
o = a.obs; sym = np.asarray(a.raw.var["feature_name"] if "feature_name" in a.raw.var else a.var["feature_name"]).astype(str)
rows = np.where(o["subtype"].astype(str).values == SUB)[0]
don = o["donor_id"].astype(str).values[rows]
u, c = np.unique(don, return_counts=True); keep = np.isin(don, u[c >= 4]); rows, don = rows[keep], don[keep]
lib = np.zeros(len(rows)); cps = np.zeros(a.raw.shape[1]); blocks = []
for i in range(0, len(rows), 4000):
    X = sp.csr_matrix(a.raw.X[rows[i:i + 4000]]); l = np.asarray(X.sum(1)).ravel(); lib[i:i + 4000] = l
    cps += np.asarray(X.multiply(1e4 / l[:, None]).sum(0)).ravel(); blocks.append(X)
cps /= len(rows)
order = np.argsort(-cps); panel = list(order[:400])
for f in FORCE:
    w = np.where(sym == f)[0]
    if len(w) and w[0] not in panel and cps[w[0]] > 0: panel.append(int(w[0]))
panel = np.array(sorted(panel))
Xp = sp.vstack([b[:, panel] for b in blocks]).tocsr().astype(np.float32)
sp.save_npz(os.path.join(HERE, f"pop_{sys.argv[1]}_X.npz"), Xp)
np.savez_compressed(os.path.join(HERE, f"pop_{sys.argv[1]}_meta.npz"), lib=lib, donor=don, panel=panel, sym=sym[panel])
json.dump(dict(subtype=SUB, nuclei=int(len(rows)), donors=int(len(np.unique(don))), median_nd=float(np.median(np.unique(don, return_counts=True)[1])),
               G=int(len(panel)), forced=[f for f in FORCE if f in set(sym[panel])]), open(os.path.join(HERE, f"pop_{sys.argv[1]}_info.json"), "w"), indent=1)
print(open(os.path.join(HERE, f"pop_{sys.argv[1]}_info.json")).read())
