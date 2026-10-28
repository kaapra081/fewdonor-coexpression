"""Summaries of the real-data runs (small CSVs only).
python real_summarize.py --runs real_out --stats-dir . --out real_summary"""
import os, glob, argparse
import numpy as np, pandas as pd
from scipy import stats
from methods_v2 import bh

HS = ["HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1", "HSPD1", "HSPE1", "DNAJB1", "DNAJA1", "HSPA8", "HSPB1"]
SEX = ["XIST", "UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "KDM6A"]
METH = {"LCP": ("p_LCP", "T_LCP"), "PB": ("p_PB", "r_PB"), "MC": ("p_MC", "est_MC"), "naive": ("p_naive", None)}


def calls(z, m):
    p = z[METH[m][0]].astype(float); rej = bh(p)
    s = np.sign(z[METH[m][1]]) if METH[m][1] else np.sign(z["S_pb"])
    return rej, s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", default="real_out")
    ap.add_argument("--stats-dir", default="."); ap.add_argument("--out", default="real_summary")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    sym = np.load(os.path.join(a.stats_dir, "real_stats_seaad.npz"), allow_pickle=True)["gene_sym"].astype(str)
    G = len(sym); iu = np.triu_indices(G, 1)
    Z = {os.path.basename(f)[:-4]: dict(np.load(f)) for f in sorted(glob.glob(os.path.join(a.runs, "*.npz")))}
    rows = []
    for name, z in Z.items():
        coh, part = name.split("_")[0], name.split("_")[1]
        r = dict(name=name, cohort=coh, part=part, D=len(z["rows"]), B=int(z["B"]))
        share = np.where(np.abs(z["S_pb"]) > 0, z["leak"] / z["S_pb"], np.nan)
        for m in METH:
            rej, _ = calls(z, m); r[f"n_{m}"] = int(rej.sum())
            if m in ("LCP", "PB"):
                r[f"leakshare_med_{m}"] = float(np.nanmedian(share[rej])) if rej.any() else np.nan
                r[f"frac_leak_gt_half_{m}"] = float(np.nanmean(share[rej] > 0.5)) if rej.any() else np.nan
        pl = z["p_LCP"].astype(float); rej = bh(pl)
        r["LCP_bh_threshold"] = float(pl[rej].max()) if rej.any() else np.nan
        r["perm_floor"] = 1.0 / (int(z["B"]) + 1)
        rows.append(r)
    runs = pd.DataFrame(rows); runs.to_csv(os.path.join(a.out, "runs_summary.csv"), index=False)

    # subsampling against the full cohort (reference = full-cohort LCP calls; itself an estimate)
    sub = []
    for name, z in Z.items():
        coh, part = name.split("_")[0], name.split("_")[1]
        if not part.startswith("sub") or f"{coh}_full" not in Z:
            continue
        full = Z[f"{coh}_full"]; ref, rs = calls(full, "LCP")
        for m in ("LCP", "PB", "MC"):
            rej, s = calls(z, m)
            same = rej & ref & (s == rs)
            sub.append(dict(name=name, cohort=coh, D=int(part[3:]), method=m, n=int(rej.sum()),
                            precision_vs_fullLCP=(same.sum() / rej.sum()) if rej.sum() else np.nan,
                            recall_of_fullLCP=same.sum() / max(ref.sum(), 1),
                            spearman_T_vs_full=stats.spearmanr(z["T_LCP"], full["T_LCP"])[0] if m == "LCP" else np.nan))
    pd.DataFrame(sub).to_csv(os.path.join(a.out, "subsample_vs_full.csv"), index=False)

    # cross-cohort replication (full cohorts)
    cc = []
    if "seaad_full" in Z and "lee_full" in Z:
        A, Bz = Z["seaad_full"], Z["lee_full"]
        for m in ("LCP", "PB", "MC"):
            ra, sa = calls(A, m); rb, sb = calls(Bz, m)
            key = METH[m][1]
            cc.append(dict(method=m, n_seaad=int(ra.sum()), n_lee=int(rb.sum()), n_both_same_sign=int((ra & rb & (sa == sb)).sum()),
                           repl_seaad_in_lee=float((ra & rb & (sa == sb)).sum() / max(ra.sum(), 1)),
                           repl_lee_in_seaad=float((ra & rb & (sa == sb)).sum() / max(rb.sum(), 1)),
                           spearman_stat=float(stats.spearmanr(A[key], Bz[key])[0])))
        # replication of LCP-only and PB-only calls
        ra, sa = calls(A, "LCP"); rp, sp = calls(A, "PB"); rbL, sbL = calls(Bz, "LCP")
        pbonly = rp & ~ra
        cc.append(dict(method="PB_only_in_seaad->LCP_lee", n_seaad=int(pbonly.sum()), n_lee=int(rbL.sum()),
                       n_both_same_sign=int((pbonly & rbL & (sp == sbL)).sum()),
                       repl_seaad_in_lee=float((pbonly & rbL & (sp == sbL)).sum() / max(pbonly.sum(), 1))))
    pd.DataFrame(cc).to_csv(os.path.join(a.out, "cross_cohort.csv"), index=False)

    # cell-capping (rare-cell analogue) vs full
    cap = []
    for name, z in Z.items():
        coh, part = name.split("_")[0], name.split("_")[1]
        if part != "cap40" or f"{coh}_full" not in Z:
            continue
        full = Z[f"{coh}_full"]
        rw = full["r_w"]; thr = np.nanquantile(np.abs(rw), 0.99)
        within_only = (np.abs(rw) >= thr) & (full["p_LCP"] > 0.05) & (full["p_PB"] > 0.05)
        for m in ("LCP", "PB", "naive", "MC"):
            rej, _ = calls(z, m); reff, _ = calls(full, "LCP")
            cap.append(dict(name=name, cohort=coh, method=m, n=int(rej.sum()),
                            n_within_only_ref=int(within_only.sum()),
                            callrate_within_only=float(rej[within_only].mean()) if within_only.any() else np.nan,
                            frac_calls_in_fullLCP=float((rej & reff).sum() / max(rej.sum(), 1))))
    pd.DataFrame(cap).to_csv(os.path.join(a.out, "cap40_vs_full.csv"), index=False)

    # split-half leak experiment: select on half A, test on half B capped at 40 cells (no reuse of cells)
    sp_ = []
    for name, z in Z.items():
        coh, part = name.split("_")[0], name.split("_")[1]
        if not part.startswith("halfBcap40") or f"{coh}_halfA" not in Z:
            continue
        A = Z[f"{coh}_halfA"]; rw = A["r_w"]; thr = np.nanquantile(np.abs(rw), 0.99)
        wo = (np.abs(rw) >= thr) & (A["p_LCP"] > 0.05) & (A["p_PB"] > 0.05)
        share = np.where(np.abs(z["S_pb"]) > 0, z["leak"] / z["S_pb"], np.nan)
        for m in ("LCP", "PB", "naive", "MC"):
            rej, _ = calls(z, m)
            sp_.append(dict(name=name, cohort=coh, method=m, n_calls=int(rej.sum()), n_within_only=int(wo.sum()),
                            n_called_within_only=int(rej[wo].sum()), callrate_within_only=float(rej[wo].mean()) if wo.any() else np.nan,
                            leakshare_med_PBcalls=float(np.nanmedian(share[calls(z, "PB")[0]]))))
    pd.DataFrame(sp_).to_csv(os.path.join(a.out, "split_leak.csv"), index=False)
    # cross-method replication (non-circular): SEA-AD calls of method X replicated by method Y in Lee
    cm = []
    if "seaad_full" in Z and "lee_full" in Z:
        for X in ("LCP", "PB"):
            rx, sx = calls(Z["seaad_full"], X)
            for Y in ("LCP", "PB"):
                ry, sy = calls(Z["lee_full"], Y)
                cm.append(dict(seaad_method=X, lee_method=Y, n_seaad=int(rx.sum()), repl=float((rx & ry & (sx == sy)).sum() / max(rx.sum(), 1))))
    pd.DataFrame(cm).to_csv(os.path.join(a.out, "cross_method_replication.csv"), index=False)

    # genes of interest
    goi = []
    idx = {g: i for i, g in enumerate(sym)}
    for name in ("seaad_full", "lee_full"):
        if name not in Z:
            continue
        z = Z[name]
        pos = {(iu[0][k], iu[1][k]): k for k in range(len(iu[0]))}
        for setname, genes in (("heat_shock", HS), ("sex", SEX)):
            gs = [g for g in genes if g in idx]
            for i1 in range(len(gs)):
                for i2 in range(i1 + 1, len(gs)):
                    j, k = sorted((idx[gs[i1]], idx[gs[i2]])); q = pos[(j, k)]
                    goi.append(dict(cohort=name.split("_")[0], set=setname, g1=sym[j], g2=sym[k],
                                    p_LCP=float(z["p_LCP"][q]), T_LCP=float(z["T_LCP"][q]), rhoB=float(z["rhoB"][q]),
                                    p_PB=float(z["p_PB"][q]), r_PB=float(z["r_PB"][q]), r_within=float(z["r_w"][q]),
                                    leak_share=float(z["leak"][q] / z["S_pb"][q]) if z["S_pb"][q] != 0 else np.nan,
                                    LCP_BH=bool(bh(z["p_LCP"].astype(float))[q]), PB_BH=bool(bh(z["p_PB"].astype(float))[q])))
    pd.DataFrame(goi).to_csv(os.path.join(a.out, "genes_of_interest.csv"), index=False)
    print(runs.to_string(index=False))


if __name__ == "__main__":
    main()
