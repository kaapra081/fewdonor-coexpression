"""Per-donor count-scale statistics for the real-data analyses (SEA-AD MTG microglia v3; Lee 2026 microglia 80x1600).
Gene panel: the 2,007-gene Lee/SEA-AD universe, restricted to the G most expressed genes in BOTH cohorts plus forced
control genes (heat-shock, PVM markers, sex-linked). Output: real_stats_<cohort>.npz with m (D x G), psi, kappa, C (D x G x G),
S, n, donor ids, sex, gene ids/symbols."""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "3")
import numpy as np, anndata as ad, scipy.sparse as sp
DATA = os.environ.get("FEWDONOR_DATA", "data/")
G_TOP = 400
FORCE = ["HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1", "HSPD1", "HSPE1", "DNAJB1", "DNAJA1", "HSPA8", "HSPB1",
         "CD163", "MRC1", "LYVE1", "F13A1", "XIST", "UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14", "KDM5C", "KDM6A"]


def donor_stats(X, lib, donor, cap=None, seed=0):
    ud = np.unique(donor); out = dict(m=[], psi=[], C=[], kappa=[], S=[], n=[], donors=ud)
    rng = np.random.default_rng(seed)
    for d in ud:
        r = np.where(donor == d)[0]
        if cap is not None and len(r) > cap:
            r = np.sort(rng.choice(r, size=cap, replace=False))
        s = lib[r].astype(float); Y = X[r].toarray().astype(float)
        S = s.sum(); x = Y / s[:, None]; m = Y.sum(0) / S; kap = (s ** 2).sum() / S ** 2
        xc = x - x.mean(0); n = len(r)
        sig2 = (xc ** 2).sum(0) / (n - 1) - (x / s[:, None]).mean(0)
        out["m"].append(m); out["psi"].append(np.maximum(kap * np.maximum(sig2, 0) + m / S, 1e-30))
        out["C"].append((xc.T @ xc / (n - 1)).astype(np.float32)); out["kappa"].append(kap); out["S"].append(S); out["n"].append(n)
    return {k: np.array(v) for k, v in out.items()}


def main():
    a = ad.read_h5ad(DATA + "scdemon_in_seaad_v3_84donors.h5ad")
    b = ad.read_h5ad(DATA + "P1_lee2026_microglia_80x1600.h5ad")
    ids = np.array(b.var_names); sym = np.array(b.var["feature_name"]).astype(str)
    XA = sp.csr_matrix(a.X); libA = np.asarray(XA.sum(1)).ravel()
    colA = np.array([a.var_names.get_loc(g) for g in ids])
    XA = XA[:, colA]
    XB = sp.csr_matrix(b.X); libB = np.asarray(b.obs["lib_size"]).astype(float)
    cpA = np.asarray(XA.multiply(1e4 / libA[:, None]).mean(0)).ravel()
    cpB = np.asarray(XB.multiply(1e4 / libB[:, None]).mean(0)).ravel()
    score = np.minimum(cpA / np.median(cpA), cpB / np.median(cpB))
    order = np.argsort(-score)
    keep = list(order[:G_TOP])
    for f in FORCE:
        w = np.where(sym == f)[0]
        if len(w) and w[0] not in keep and cpA[w[0]] > 0 and cpB[w[0]] > 0:
            keep.append(int(w[0]))
    keep = np.array(sorted(keep))
    meta = dict(G=len(keep), forced_present=[f for f in FORCE if f in set(sym[keep])],
                forced_missing=[f for f in FORCE if f not in set(sym[keep])])
    sexA_src = ad.read_h5ad(DATA + "seaad_mtg_microglia.h5ad", backed="r").obs[["donor_id", "sex"]]
    sexA = sexA_src.drop_duplicates("donor_id").set_index("donor_id")["sex"].astype(str).to_dict()
    sexB = b.obs.drop_duplicates("donor_id").set_index("donor_id")["sex"].astype(str).to_dict()
    for name, X, lib, don, sexmap in (("seaad", XA, libA, np.asarray(a.obs["donor_id"]).astype(str), sexA),
                                      ("lee", XB, libB, np.asarray(b.obs["donor_id"]).astype(str), sexB)):
        st = donor_stats(X[:, keep].tocsr(), lib, don)
        np.savez_compressed(f"real_stats_{name}.npz", **st, gene_ids=ids[keep], gene_sym=sym[keep],
                            sex=np.array([sexmap.get(d, "NA") for d in st["donors"]]))
        meta[name] = dict(D=int(len(st["donors"])), n_median=float(np.median(st["n"])), cells=int(st["n"].sum()))
        for rep in range(3):                                    # rare-cell analogue: <= 40 cells per donor
            stc = donor_stats(X[:, keep].tocsr(), lib, don, cap=40, seed=1000 + rep)
            np.savez_compressed(f"real_stats_{name}_cap40_r{rep}.npz", **stc, gene_ids=ids[keep], gene_sym=sym[keep],
                                sex=np.array([sexmap.get(d, "NA") for d in stc["donors"]]))
    json.dump(meta, open("real_prep_meta.json", "w"), indent=1); print(json.dumps(meta))


if __name__ == "__main__":
    main()
