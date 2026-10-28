"""Leak experiment without selection bias: split each donor's cells at random into halves A and B.
Half A (all its cells) is used ONLY to select within-dominant pairs; half B, capped at 40 cells per donor (all of half B
if it has fewer), is used ONLY to test. Outputs real_stats_<coh>_halfA.npz and real_stats_<coh>_halfB_cap40_r{0,1,2}.npz."""
import os, numpy as np, anndata as ad, scipy.sparse as sp
from real_prep import donor_stats, DATA
keep_ids = np.load("real_stats_seaad.npz", allow_pickle=True)["gene_ids"]


def run(name, X, lib, don, ids):
    col = np.array([np.where(ids == g)[0][0] for g in keep_ids])
    Xs = X[:, col].tocsr()
    rng = np.random.default_rng(777)
    half = np.zeros(len(don), int)
    for d in np.unique(don):
        r = np.where(don == d)[0]; perm = rng.permutation(r); half[perm[len(r) // 2:]] = 1
    A = np.where(half == 0)[0]; B = np.where(half == 1)[0]
    st = donor_stats(Xs[A], lib[A], don[A])
    np.savez_compressed(f"real_stats_{name}_halfA.npz", **st, gene_ids=keep_ids, gene_sym=np.load("real_stats_seaad.npz", allow_pickle=True)["gene_sym"], sex=np.array(["NA"] * len(st["donors"])))
    for rep in range(3):
        stc = donor_stats(Xs[B], lib[B], don[B], cap=40, seed=2000 + rep)
        np.savez_compressed(f"real_stats_{name}_halfB_cap40_r{rep}.npz", **stc, gene_ids=keep_ids,
                            gene_sym=np.load("real_stats_seaad.npz", allow_pickle=True)["gene_sym"], sex=np.array(["NA"] * len(stc["donors"])))
    print(name, "A cells", len(A), "B cells", len(B), "min B cells/donor", int(np.bincount(np.unique(don[B], return_inverse=True)[1]).min()))


a = ad.read_h5ad(DATA + "scdemon_in_seaad_v3_84donors.h5ad")
XA = sp.csr_matrix(a.X); run("seaad", XA, np.asarray(XA.sum(1)).ravel(), np.asarray(a.obs["donor_id"]).astype(str), np.array(a.var_names))
b = ad.read_h5ad(DATA + "P1_lee2026_microglia_80x1600.h5ad")
run("lee", sp.csr_matrix(b.X), np.asarray(b.obs["lib_size"]).astype(float), np.asarray(b.obs["donor_id"]).astype(str), np.array(b.var_names))
