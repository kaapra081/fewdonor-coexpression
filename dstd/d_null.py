"""DSTD_PREREG.md E2: LCPd/LCPXd/PBlin on the FIX_PREREG Item 1: real-cell null. Donor labels permuted across the nuclei of a supertype (n_d kept), so Sigma_B = 0 exactly.
python fix_null.py --set S1 --k 0"""
import os, sys, json, argparse, time
import numpy as np, anndata as ad, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("real", "dev"): sys.path.insert(0, os.path.join(HERE, "..", p))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "fix"))
from real_prep import donor_stats  # noqa: E402
from methods_v2 import paule_mandel, bh  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
import torch
ap = argparse.ArgumentParser(); ap.add_argument("--set", required=True); ap.add_argument("--k", type=int, required=True); ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args(); torch.set_num_threads(a.threads)
DATA = os.environ.get("FEWDONOR_DATA", "/Users/kaavinprasanna/Documents/recomb27/data/")
SUP = {"S1": "Micro-PVM_1", "S2": "Micro-PVM_3-SEAAD"}[a.set]
rs = np.load(os.path.join(HERE, "..", "rare", f"rare_stats_{a.set}.npz"), allow_pickle=True)
gk = rs["gene_idx"].astype(int); ids = rs["gene_ids"].astype(str)
full = ad.read_h5ad(DATA + "seaad_mtg_microglia.h5ad", backed="r").obs
mask = ((full["assay"].astype(str) == "10x 3' v3") & (full["Neurotypical reference"].astype(str) != "True")).values
obs = full[mask]; A = ad.read_h5ad(DATA + "scdemon_in_seaad_v3_84donors.h5ad")
X = sp.csr_matrix(A.X); lib = np.asarray(X.sum(1)).ravel()
X = X[:, np.array([A.var_names.get_loc(g) for g in ids])].tocsr()
sup = obs["Supertype"].astype(str).values; don = obs["donor_id"].astype(str).values
rows = np.where(sup == SUP)[0]; keepd = set(rs["donors"].astype(str)); rows = rows[np.isin(don[rows], list(keepd))]
rng = np.random.default_rng(31_000_000 + a.k + (0 if a.set == "S1" else 1000))
pd_lab = rng.permutation(don[rows])                     # pseudo-donor labels: same n_d, random nuclei
Xs = X[rows]; ls = lib[rows]
d = donor_stats(Xs, ls, pd_lab); d["C"] = d["C"].astype(np.float64)
G = d["m"].shape[1]; iu = np.triu_indices(G, 1)
t0 = time.time()
cells = ((Xs[np.where(pd_lab == u)[0]].toarray().astype(np.float64), ls[np.where(pd_lab == u)[0]].astype(np.float64)) for u in d["donors"])
dd = donor_level_cells(cells, 31_000_000 + a.k + 7)
o = methods_d(dd, iu, B=100_000, seed=31_000_000 + a.k)
ref = np.load(os.path.join(HERE, "..", "fix", f"null_{a.set}_k{a.k}.npz"))
res = dict(p_LCPd=o["LCPd"], p_LCPXd=o["LCPXd"], p_PBlin=o["PBlin"], r_w=ref["r_w"])
np.savez_compressed(os.path.join(HERE, f"dnull_{a.set}_k{a.k}.npz"), **{k: np.asarray(v, np.float32) for k, v in res.items()})
print(a.set, a.k, "D", len(d["donors"]), "calls", {m: int(bh(res[f"p_{m}"].astype(float)).sum()) for m in ("LCPd", "LCPXd", "PBlin")}, round(time.time() - t0), "s", flush=True)
