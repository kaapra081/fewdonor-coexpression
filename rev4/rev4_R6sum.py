"""REV4 Amendment 1, R6 summaries and Mantis map tables from out/R6_P1.npz -> out/R6/ (CSV + R6_summary.json)."""
import os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import cross_val_score, KFold
import umap
HERE = os.path.dirname(os.path.abspath(__file__)); OD = os.path.join(HERE, "out", "R6"); os.makedirs(OD, exist_ok=True)
Z = np.load(os.path.join(HERE, "out", "R6_P1.npz"), allow_pickle=True); sym = Z["sym"].astype(str); G = len(sym); iu = np.triu_indices(G, 1); A_, B_ = iu
src = json.load(open(os.path.join(HERE, "donor_source.json"))); ad = json.load(open(os.path.join(HERE, "..", "rev2", "adlab.json")))
S = {}


def sp(x, y):
    f = np.isfinite(x) & np.isfinite(y); return float(spearmanr(x[f], y[f]).correlation)


for tag in ("null", "real"):
    nd = Z[f"{tag}_nd"]; dp = Z[f"{tag}_depth"]; don = Z[f"{tag}_donor"].astype(str)
    df = pd.DataFrame({"donor": [f"pseudo-donor {i + 1:04d}" for i in range(len(nd))] if tag == "null" else don, "cells": nd.astype(int), "log10_cells": np.log10(nd), "mean_depth_umi": np.round(dp, 1)})
    if tag == "real":
        df["brain_bank"] = [src.get(d_, "NA") for d_ in don]; df["AD"] = [{"1.0": "AD", "0.0": "control"}.get(ad.get(d_, "nan"), "NA") for d_ in don]
    for prof in ("L", "R", "Lz"):   # L (centred log-CPM) and R are the frozen profiles; Lz (gene-standardized log-CPM, as Pearson sees it) is an added, labelled variant
        X = Z[f"{tag}_{prof[0]}"].astype(np.float64); Xc = X - X.mean(0)
        if prof == "Lz": Xc = Xc / np.maximum(Xc.std(0), 1e-12)
        nrm = np.sqrt((Xc ** 2).sum(1)); emb = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=0).fit_transform(Xc)
        knn = cross_val_score(KNeighborsRegressor(n_neighbors=15), Xc, np.log(nd), cv=KFold(5, shuffle=True, random_state=0), scoring="r2").mean()
        U, s, Vt = np.linalg.svd(Xc, full_matrices=False); pcs = U[:, :5] * s[:5]
        r2pc = [float(np.corrcoef(pcs[:, i], np.log(nd))[0, 1] ** 2) for i in range(5)]
        S[f"{tag}_{prof}"] = dict(spearman_norm_ncells=sp(nrm, nd), knn_cv_r2_logn=float(knn), r2_logn_PC1to5=r2pc, varfrac_PC1to5=(s[:5] ** 2 / (s ** 2).sum()).tolist(), spearman_norm_depth=sp(nrm, dp))
        df[f"x_{prof}"] = np.round(emb[:, 0], 4); df[f"y_{prof}"] = np.round(emb[:, 1], 4); df[f"profile_norm_{prof}"] = np.round(nrm, 3)
    df.to_csv(os.path.join(OD, f"donors_{tag}.csv"), index=False)
rw = Z["r_w"]; S["pairs"] = dict(n=int(len(rw)), finite_rw=int(np.isfinite(rw).sum()),
                                  spearman_nullPB_rw=sp(Z["null_rPB"], rw), spearman_realPB_rw=sp(Z["real_rPB"], rw), spearman_realXDSMz_rw=sp(Z["real_zXDSM"], rw),
                                  spearman_nullPB_realPB=sp(Z["null_rPB"], Z["real_rPB"]), spearman_nullPB_realXDSMz=sp(Z["null_rPB"], Z["real_zXDSM"]),
                                  null_calls=int(Z["null_cPB"].sum()), real_PB_calls=int(Z["real_cPB"].sum()), real_XDSM_calls=int(Z["real_cXDSM"].sum()))
cP = Z["real_cPB"].astype(bool); cX = Z["real_cXDSM"].astype(bool); c0 = Z["null_cPB"].astype(bool)
# among real pseudobulk calls: fraction that are also null false calls, by class
S["pairs"].update(PBonly=int((cP & ~cX).sum()), shared=int((cP & cX).sum()), XDSMonly=int((cX & ~cP).sum()),
                  frac_PBonly_in_null=float(c0[cP & ~cX].mean()), frac_shared_in_null=float(c0[cP & cX].mean()), frac_XDSMonly_in_null=float(c0[cX & ~cP].mean()),
                  med_rw_PBonly=float(np.nanmedian(rw[cP & ~cX])), med_rw_shared=float(np.nanmedian(rw[cP & cX])), med_rw_XDSMonly=float(np.nanmedian(rw[cX & ~cP])), med_rw_all=float(np.nanmedian(rw)))
def deg(mask): return np.bincount(np.concatenate([A_[mask], B_[mask]]), minlength=G)
arw = np.abs(rw); mrw = np.array([np.nanmean(np.concatenate([arw[A_ == g], arw[B_ == g]])) for g in range(G)])
gt = pd.DataFrame({"gene": sym, "null_false_degree": deg(c0), "PB_degree": deg(cP), "XDSM_degree": deg(cX), "PB_only_degree": deg(cP & ~cX), "XDSM_only_degree": deg(cX & ~cP),
                   "mean_umi_per_cell": np.round(Z["mean_umi"], 4), "detection_fraction": np.round(Z["detect"], 4), "mean_abs_within_donor_r": np.round(mrw, 4)})
gt.to_csv(os.path.join(OD, "genes_numeric.csv"), index=False)
S["genes"] = dict(spearman_nullfalse_PBonly=sp(gt.null_false_degree.values.astype(float), gt.PB_only_degree.values.astype(float)),
                  spearman_nullfalse_meanabsrw=sp(gt.null_false_degree.values.astype(float), mrw), spearman_nullfalse_meanumi=sp(gt.null_false_degree.values.astype(float), Z["mean_umi"]),
                  spearman_PBonly_meanabsrw=sp(gt.PB_only_degree.values.astype(float), mrw), spearman_XDSM_meanabsrw=sp(gt.XDSM_degree.values.astype(float), mrw))
json.dump(S, open(os.path.join(OD, "R6_summary.json"), "w"), indent=1); print(json.dumps(S, indent=1))
