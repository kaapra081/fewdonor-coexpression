"""Post hoc real-data diagnostics (not in any analysis plan).
(1) Depth-state coupling: pooled within-donor Pearson correlation between log library size and x_ij = Y_ij / lib_i for the 417 analysed genes
    (same statistic as diag_sim.py; zero in expectation when depth is independent of cell state).
(2) Low-rank donor factor: share of variance of standardized donor means r_dj on the first principal component, and its correlation with
    donor-level covariates (mean log library size, cell number, sex; SEA-AD also mean mitochondrial fraction)."""
import os, sys, json
import numpy as np, anndata as ad, scipy.sparse as sp, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import paule_mandel  # noqa: E402
DATA = os.environ.get("FEWDONOR_DATA", "data/"); RS = os.environ.get("FEWDONOR_STATS", os.path.join(HERE, "..", "real"))
out = {}
for coh, fn in (("seaad", "scdemon_in_seaad_v3_84donors.h5ad"), ("lee", "P1_lee2026_microglia_80x1600.h5ad")):
    z = np.load(os.path.join(RS, f"real_stats_{coh}.npz"), allow_pickle=True)
    ids = z["gene_ids"].astype(str); donors = z["donors"].astype(str)
    a = ad.read_h5ad(DATA + fn)
    X = sp.csr_matrix(a.X)
    lib = np.asarray(X.sum(1)).ravel() if coh == "seaad" else np.asarray(a.obs["lib_size"]).astype(float)
    X = X[:, np.array([a.var_names.get_loc(g) for g in ids])].tocsr()
    don = np.asarray(a.obs["donor_id"]).astype(str)
    lc_all, num, den_x, den_l, sdl = [], np.zeros(len(ids)), np.zeros(len(ids)), 0.0, []
    for d in donors:
        r = np.where(don == d)[0]
        ls = np.log(lib[r]); lc = ls - ls.mean(); sdl.append(ls.std(ddof=1))
        x = X[r].toarray() / lib[r][:, None]; xc = x - x.mean(0)
        num += lc @ xc; den_x += (xc ** 2).sum(0); den_l += (lc ** 2).sum()
    c = num / np.sqrt(den_l * den_x)
    # donor factor
    m, psi = z["m"], z["psi"]; tau2, mu, _ = paule_mandel(m, psi); R = (m - mu) / np.sqrt(tau2[None, :] + psi)
    Rc = R - R.mean(0); U, Sv, Vt = np.linalg.svd(Rc, full_matrices=False); pc1 = U[:, 0] * Sv[0]
    share = float(Sv[0] ** 2 / (Sv ** 2).sum()); share2 = float(Sv[1] ** 2 / (Sv ** 2).sum())
    cov = pd.DataFrame(dict(donor=donors, mloglib=[np.log(lib[don == d]).mean() for d in donors], n=z["n"].astype(float),
                            male=(z["sex"].astype(str) == "male").astype(float)))
    if coh == "seaad":
        ob = ad.read_h5ad(DATA + "seaad_mtg_microglia.h5ad", backed="r").obs
        col = [k for k in ob.columns if "mito" in k.lower()][0]
        mt = ob.groupby("donor_id", observed=True)[col].mean().astype(float)
        cov["mito"] = cov.donor.map(mt).astype(float)
    corr = {k: float(np.corrcoef(pc1, cov[k])[0, 1]) for k in cov.columns if k != "donor" and cov[k].std() > 0}
    out[coh] = dict(coupling_median=float(np.median(c)), coupling_absmedian=float(np.median(np.abs(c))),
                    coupling_q=[float(v) for v in np.quantile(c, [0.05, 0.25, 0.75, 0.95])],
                    frac_abs_gt_0045=float(np.mean(np.abs(c) > 0.045)), frac_abs_gt_01=float(np.mean(np.abs(c) > 0.1)),
                    sd_loglib_median=float(np.median(sdl)), n_cells=int(sum((don == d).sum() for d in donors)),
                    pc1_share=share, pc2_share=share2, pc1_corr=corr)
    pd.DataFrame(dict(gene=z["gene_sym"].astype(str), coupling=c)).to_csv(os.path.join(HERE, f"coupling_{coh}.csv"), index=False)
    print(coh, json.dumps(out[coh]), flush=True)
json.dump(out, open(os.path.join(HERE, "diag_real.json"), "w"), indent=1)
