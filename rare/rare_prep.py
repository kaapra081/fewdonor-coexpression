"""Rare-supertype analysis (RARE_PREREG.md): per-donor count-scale statistics for SEA-AD Micro-PVM_1 (S1), Micro-PVM_3-SEAAD (S2)
and Micro-PVM_2 (REF), on the 417 genes of the main analysis. Donors with < 4 nuclei of a supertype are excluded."""
import os, sys, json
import numpy as np, anndata as ad, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "real"))
from real_prep import donor_stats  # noqa: E402
DATA = os.environ.get("FEWDONOR_DATA", "/Users/kaavinprasanna/Documents/recomb27/data/")
STATS = os.environ.get("FEWDONOR_STATS", os.path.join(HERE, "..", "real"))
SETS = {"S1": "Micro-PVM_1", "S2": "Micro-PVM_3-SEAAD", "REF": "Micro-PVM_2"}

ref = np.load(os.path.join(STATS, "real_stats_seaad.npz"), allow_pickle=True)
ids, sym = ref["gene_ids"].astype(str), ref["gene_sym"].astype(str)
full = ad.read_h5ad(DATA + "seaad_mtg_microglia.h5ad", backed="r").obs
mask = ((full["assay"].astype(str) == "10x 3' v3") & (full["Neurotypical reference"].astype(str) != "True")).values
obs = full[mask]
a = ad.read_h5ad(DATA + "scdemon_in_seaad_v3_84donors.h5ad")
assert a.n_obs == len(obs) and (np.asarray(a.obs["donor_id"]).astype(str) == obs["donor_id"].astype(str).values).all()
X = sp.csr_matrix(a.X); lib = np.asarray(X.sum(1)).ravel()
X = X[:, np.array([a.var_names.get_loc(g) for g in ids])].tocsr()
sup = obs["Supertype"].astype(str).values; don = obs["donor_id"].astype(str).values
sexmap = obs.drop_duplicates("donor_id").set_index("donor_id")["sex"].astype(str).to_dict()
meta = {}
for key, st in SETS.items():
    rows = np.where(sup == st)[0]
    cnt = {d: int((don[rows] == d).sum()) for d in np.unique(don[rows])}
    keep_d = np.array([d for d, c in cnt.items() if c >= 4])
    rows = rows[np.isin(don[rows], keep_d)]
    tot = np.asarray(X[rows].sum(0)).ravel(); gk = np.where(tot >= 20)[0]
    s = donor_stats(X[rows][:, gk].tocsr(), lib[rows], don[rows])
    np.savez_compressed(os.path.join(HERE, f"rare_stats_{key}.npz"), **s, gene_ids=ids[gk], gene_sym=sym[gk], gene_idx=gk,
                        sex=np.array([sexmap.get(d, "NA") for d in s["donors"]]))
    meta[key] = dict(supertype=st, D=int(len(s["donors"])), cells=int(len(rows)), n_median=float(np.median(s["n"])),
                     n_q10_q90=[float(np.quantile(s["n"], 0.1)), float(np.quantile(s["n"], 0.9))], G=int(len(gk)), n_list=s["n"].astype(int).tolist())
    print(key, {k: v for k, v in meta[key].items() if k != "n_list"}, flush=True)
json.dump(meta, open(os.path.join(HERE, "rare_prep_meta.json"), "w"), indent=1)
