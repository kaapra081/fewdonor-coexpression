"""Build a real-cell null for between-donor co-expression from an h5ad file (one command).

  python examples/make_null.py --h5ad cells.h5ad --donor-col donor_id \
      [--type-col subtype --type Homeo_CECR2] [--mode donor-centred|pooled|depth-matched] [--seed 0] --out null.h5ad

pooled         donor labels of the cells are permuted; pseudo-donors keep the real cell numbers. Under this null no gene pair
               co-varies between pseudo-donors, but real between-donor differences become within-pseudo-donor structure.
donor-centred  each donor's shrunken (Paule-Mandel) mean effect is divided out of its cells first (factor clipped to
               [0.1, 10]), then labels are permuted, as in the paper (rev4/rev4_C.py, centre()).
depth-matched  as donor-centred, but labels are permuted within 10 equal-count bins of log library size, so each pseudo-donor
               keeps its real donor's depth distribution (Plan 6 D).
Counts are read from .raw.X if present, else .X; library size is the per-cell total over all genes. The output h5ad holds the
(centred) counts in X, obs['pseudo_donor'] and obs['library_size']; any between-donor test run on it should make no calls.
"""
import argparse, os, sys, numpy as np, scipy.sparse as sp, anndata as ad
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dev"))
from methods_v2 import paule_mandel  # noqa: E402  (same estimator as every analysis in the paper)


def donor_stats(X, lib, rows):
    """Ratio mean m and design sampling variance psi of one donor (as dstd._ms)."""
    Y = X[rows].toarray().astype(np.float64); s = lib[rows]; S = s.sum(); x = Y / s[:, None]; m = Y.sum(0) / S
    kap = (s ** 2).sum() / S ** 2; n = len(s); xc = x - x.mean(0)
    sig2 = (xc ** 2).sum(0) / max(n - 1, 1) - (x / s[:, None]).mean(0)
    return m, np.maximum(kap * np.maximum(sig2, 0) + m / S, 1e-30)


def centre(X, lib, don, lo=0.1, hi=10.0):
    ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]; st = [donor_stats(X, lib, r) for r in rows]
    m = np.array([a for a, _ in st]); ps = np.array([b for _, b in st]); tau2, mu, _ = paule_mandel(m, ps)
    th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), lo, hi)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for i, r in enumerate(rows): idx[r] = i
    return sp.csr_matrix(X.multiply(f[idx]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True); ap.add_argument("--donor-col", required=True); ap.add_argument("--type-col"); ap.add_argument("--type")
    ap.add_argument("--mode", default="donor-centred", choices=("pooled", "donor-centred", "depth-matched")); ap.add_argument("--min-cells", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--out", required=True); a = ap.parse_args()
    A = ad.read_h5ad(a.h5ad); keep = np.ones(A.n_obs, bool)
    if a.type_col: keep &= (A.obs[a.type_col].astype(str).values == a.type)
    X = (A.raw.X if A.raw is not None else A.X); X = sp.csr_matrix(X)[np.where(keep)[0]]; var = (A.raw.var if A.raw is not None else A.var)
    lib = np.asarray(X.sum(1)).ravel().astype(np.float64); don = A.obs[a.donor_col].astype(str).values[keep]
    u, c = np.unique(don, return_counts=True); ok = np.isin(don, u[c >= a.min_cells]) & (lib > 0); X, lib, don = X[np.where(ok)[0]], lib[ok], don[ok]
    rng = np.random.default_rng(a.seed); lab = don.copy()
    if a.mode == "pooled": Xn = X.astype(np.float64); lab = rng.permutation(don)
    else:
        Xn = centre(X, lib, don)
        if a.mode == "donor-centred": lab = rng.permutation(don)
        else:
            b = np.digitize(np.log(lib), np.quantile(np.log(lib), np.linspace(0, 1, 11)[1:-1]))
            for k in range(10):
                ix = np.where(b == k)[0]; lab[ix] = rng.permutation(don[ix])
    out = ad.AnnData(X=Xn, var=var.copy()); out.obs["pseudo_donor"] = lab; out.obs["library_size"] = lib; out.obs["source_donor"] = don
    out.uns["null"] = dict(mode=a.mode, seed=a.seed, source=a.h5ad, cells=int(len(lab)), donors=int(len(np.unique(lab))))
    out.write_h5ad(a.out); print(f"{a.mode} null: {len(lab)} cells, {len(np.unique(lab))} pseudo-donors -> {a.out}")


if __name__ == "__main__":
    main()
