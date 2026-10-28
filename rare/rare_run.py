"""Rare-supertype analysis (RARE_PREREG.md). LCP (B = 1e5, CPU torch, no refinement, seed 28102026) and PB on S1, S2, REF; Q1 and Q2."""
import os, sys, json
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "real"))
from real_gpu import load, lcp_gpu, comparators  # noqa: E402
from methods_v2 import bh  # noqa: E402
B, SEED, G0 = 100_000, 28102026, 417
HS = {"HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1", "HSPD1", "HSPE1", "DNAJB1", "DNAJA1", "HSPA8", "HSPB1"}
sym0 = np.load(os.path.join(HERE, "..", "real", "real_stats_seaad.npz"), allow_pickle=True)["gene_sym"].astype(str)


def analyse(key):
    d = load(os.path.join(HERE, f"rare_stats_{key}.npz")); gi = d["gene_idx"].astype(int)
    o = lcp_gpu(d, B, SEED, "cpu"); c = comparators(d, o["iu"]); a, b = o["iu"]
    P = pd.DataFrame(dict(g1=gi[a], g2=gi[b], p_LCP=o["p"], T_LCP=o["T"], p_PB=c["p_PB"], r_PB=c["r_PB"], r_w=c["r_w"],
                          share=np.where(np.abs(c["S_pb"]) > 0, c["leak"] / c["S_pb"], np.nan)))
    P["LCP"] = bh(P.p_LCP.values); P["PB"] = bh(P.p_PB.values)
    return P.set_index(["g1", "g2"]), d


def top_hubs(edges, k=20):
    deg = np.zeros(G0, int)
    for g1, g2 in edges: deg[g1] += 1; deg[g2] += 1
    order = sorted(range(G0), key=lambda g: (-deg[g], g)); return set(order[:k]), deg


out = {}
REF, dref = analyse("REF"); out["REF"] = dict(D=int(len(dref["donors"])), n_LCP=int(REF.LCP.sum()), n_PB=int(REF.PB.sum()))
for key in ("S1", "S2"):
    P, d = analyse(key)
    J_ = P.join(REF[["r_w", "LCP", "T_LCP"]].rename(columns={"r_w": "r_w_ref", "LCP": "LCP_ref", "T_LCP": "T_ref"}), how="left")
    pbo = J_[J_.PB & ~J_.LCP]; both = J_[J_.PB & J_.LCP]; lco = J_[~J_.PB & J_.LCP]
    sgn = lambda X, col: np.sign(X[col].values)
    def repl(X, col):
        Y = X.dropna(subset=["T_ref"]); return float(((Y.LCP_ref.values) & (np.sign(Y.T_ref.values) == np.sign(Y[col].values))).mean()) if len(Y) else np.nan
    w = lambda X: np.abs(X.r_w_ref.dropna().values)
    q1 = dict(n_PBonly=len(pbo), n_both=len(both), n_LCPonly=len(lco),
              med_absrw_PBonly=float(np.median(w(pbo))) if len(pbo) else np.nan, med_absrw_both=float(np.median(w(both))) if len(both) else np.nan,
              med_absrw_LCPonly=float(np.median(w(lco))) if len(lco) else np.nan, med_absrw_all=float(np.median(np.abs(J_.r_w_ref.dropna()))))
    q1["p_MW_vs_both"] = float(mannwhitneyu(w(pbo), w(both), alternative="greater").pvalue) if len(pbo) and len(both) else np.nan
    q1["p_MW_vs_LCPonly"] = float(mannwhitneyu(w(pbo), w(lco), alternative="greater").pvalue) if len(pbo) and len(lco) >= 50 else np.nan
    q1["repl_ref_PBonly"] = repl(pbo, "r_PB"); q1["repl_ref_both"] = repl(both, "r_PB"); q1["repl_ref_LCPonly"] = repl(lco, "T_LCP")
    thr = np.nanquantile(np.abs(J_.r_w_ref), 0.99)
    q1["frac_top1pct_rw_PBonly"] = float((w(pbo) >= thr).mean()) if len(pbo) else np.nan; q1["frac_top1pct_rw_both"] = float((w(both) >= thr).mean()) if len(both) else np.nan
    q1["leakshare_med_PBonly"] = float(np.nanmedian(pbo.share)) if len(pbo) else np.nan; q1["leakshare_med_both"] = float(np.nanmedian(both.share)) if len(both) else np.nan
    cond_b = (q1["med_absrw_PBonly"] > q1["med_absrw_both"]) and q1["p_MW_vs_both"] < 0.01 and \
             (len(lco) < 50 or ((q1["med_absrw_PBonly"] > q1["med_absrw_LCPonly"]) and q1["p_MW_vs_LCPonly"] < 0.01))
    cond_c = (q1["repl_ref_both"] - q1["repl_ref_PBonly"]) >= 0.10
    q1["decision"] = "not assessable" if len(pbo) < 50 else ("YES" if (cond_b and cond_c) else "NO"); q1["cond_b"] = bool(cond_b); q1["cond_c"] = bool(cond_c)
    # Q2
    E = list(J_[J_.PB].index); Eb = list(both.index)
    H, deg = top_hubs(E); Hb, degb = top_hubs(Eb); Jobs = len(H & Hb) / len(H | Hb)
    rng = np.random.default_rng(SEED); Jr = []
    for _ in range(200):
        keep = rng.choice(len(E), size=len(E) - len(pbo), replace=False); Hr, _ = top_hubs([E[i] for i in keep]); Jr.append(len(H & Hr) / len(H | Hr))
    q2 = dict(J=float(Jobs), J_random_p05=float(np.quantile(Jr, 0.05)), J_random_median=float(np.median(Jr)),
              decision="YES" if (Jobs < 0.5 and Jobs < np.quantile(Jr, 0.05)) else "NO",
              hubs_PB=[sym0[g] for g in sorted(H, key=lambda g: (-deg[g], g))], hubs_noPBonly=[sym0[g] for g in sorted(Hb, key=lambda g: (-degb[g], g))],
              hubs_dropped=[sym0[g] for g in sorted(H - Hb)], hubs_added=[sym0[g] for g in sorted(Hb - H)])
    gl = pd.Series([sym0[g] for p in pbo.index for g in p]).value_counts()
    out[key] = dict(D=int(len(d["donors"])), cells=int(d["n"].sum()), n_median=float(np.median(d["n"])), G=int(len(d["gene_idx"])),
                    n_pairs=int(len(J_)), n_LCP=int(J_.LCP.sum()), n_PB=int(J_.PB.sum()), Q1=q1, Q2=q2,
                    PBonly_top_genes=gl.head(15).to_dict(), PBonly_heatshock_pairs=int(sum(1 for p in pbo.index if sym0[p[0]] in HS and sym0[p[1]] in HS)))
    J_.reset_index().to_csv(os.path.join(HERE, f"rare_pairs_{key}.csv.gz"), index=False)
    print(key, json.dumps(out[key], default=str)[:3000], flush=True)
json.dump(out, open(os.path.join(HERE, "rare_results.json"), "w"), indent=1, default=str)
