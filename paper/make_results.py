"""Generate results.tex (LaTeX macros), appendix_tables.tex and the simulation tables from the result CSVs,
so that every number in the paper is copied programmatically from the saved results."""
import os, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root
V2 = ROOT
P = os.path.join(ROOT, "paper")
s = pd.read_csv(os.path.join(ROOT, "results", "dev", "summary_all.csv"))
ho_path = os.path.join(ROOT, "results", "heldout", "summary_all.csv")
ho = pd.read_csv(ho_path) if os.path.exists(ho_path) else None
R = os.path.join(ROOT, "results", "real")
cal = json.load(open(os.path.join(ROOT, "results", "real_between_variance_calibration.json")))
deff = json.load(open(os.path.join(R, "deff_real.json")))
runs = pd.read_csv(os.path.join(R, "runs_summary.csv")); sub = pd.read_csv(os.path.join(R, "subsample_vs_full.csv"))
cap = pd.read_csv(os.path.join(R, "split_leak.csv")); xm = pd.read_csv(os.path.join(R, "cross_method_replication.csv")); cc = pd.read_csv(os.path.join(R, "cross_cohort.csv"))
goi = pd.read_csv(os.path.join(R, "genes_of_interest.csv"))

M = {}
f3 = lambda x: f"{x:.3f}"
f2 = lambda x: f"{x:.2f}"


def get(df, d, m, col):
    v = df[(df.design == d) & (df.method == m)][col]
    assert len(v) == 1, (d, m, col, len(v)); return float(v.iloc[0])


def fm(df, d, m, col):          # value with MCSE
    return f"{get(df, d, m, col):.3f}\\,({get(df, d, m, col + '_mcse'):.3f})"


for d, nm in (("D8", "Deight"), ("D12", "Dtwelve"), ("D16", "Dsixteen"), ("D24", "Dtwentyfour"), ("D40", "Dforty")):
    M[f"LCPFDR{nm}"] = fm(s, d, "LCP", "FDR")
    M[f"LCPpow{nm}"] = f3(get(s, d, "LCP", "power"))
    M[f"BOUNDpow{nm}"] = f3(get(s, d, "BOUND", "power"))
    M[f"ORACLEpow{nm}"] = f3(get(s, d, "ORACLE", "power"))
    M[f"MCpow{nm}"] = f3(get(s, d, "CMB_mom", "power"))
pbf = [get(s, d, "PB_logcpm", "FDR") for d in ("D24", "D40")]
M["PBFDRrange"] = f"{min(pbf):.2f}--{max(pbf):.2f}"
of = [get(s, d, mm, "FDR") for d in ("D24", "D40") for mm in ("E5", "CMB_naive")]
M["OtherFDRrange"] = f"{min(of):.2f}--{max(of):.2f}"
qq = lambda m: [get(s, d, m, "rej05_QQ_cross") for d in ("D24", "D40")]
M["QQLCP"] = f"{min(qq('LCP')):.3f}--{max(qq('LCP')):.3f}"; M["QQORACLE"] = f"{min(qq('ORACLE')):.3f}--{max(qq('ORACLE')):.3f}"
M["STRESSLCPFDR"] = fm(s, "D16_STRESS_depth", "LCP", "FDR"); M["STRESSPBFDR"] = f3(get(s, "D16_STRESS_depth", "PB_logcpm", "FDR"))
M["PBSSmax"] = f2(s[s.method == "PB_logcpm"].rej05_SS.max())
M["SEAADnsForty"] = f2(cal["seaad"]["noise_signal_n40_q"][2]); M["LEEnsForty"] = f2(cal["lee"]["noise_signal_n40_q"][2])
dj = json.load(open(os.path.join(R, "deff_real.json")))
M["SEAADnsFortyTop"] = f2(dj["real_stats_seaad_cap40_r0.npz"]["median_noise_to_signal"]); M["LEEnsFortyTop"] = f2(dj["real_stats_lee_cap40_r0.npz"]["median_noise_to_signal"])
_dt = [16 * dj[k]["Deff_over_D_median"] for k in ("real_stats_seaad_cap40_r0.npz", "real_stats_lee_cap40_r0.npz")]
M["DeffSixteenTop"] = f"{min(_dt):.1f}--{max(_dt):.1f}"
_du = [16 / (1 + cal[c]["noise_signal_n40_q"][2]) ** 2 for c in ("seaad", "lee")]   # pair of two median universe genes
M["DeffSixteenUni"] = f"{min(_du):.1f}--{max(_du):.1f}"
_ds = json.load(open(os.path.join(ROOT, "results", "diag", "diag_sim.json")))
_are = [v["ARE_median"] for k, v in _ds.items() if "STRESS" not in k]
M["AREmin"] = f"{min(_are):.3f}"; M["AREmax"] = f"{max(_are):.3f}"
_dsim = [_ds[k]["Deff_over_D_median"] for k in ("D16", "D24", "D40")]
M["DeffSimMin"] = f2(min(_dsim)); M["DeffSimMax"] = f2(max(_dsim))
_dc = sorted({f2(dj[k]["Deff_over_D_median"]) for k in ("real_stats_seaad_cap40_r0.npz", "real_stats_lee_cap40_r0.npz")})
M["DeffCapRange"] = "--".join(_dc)
_dr = json.load(open(os.path.join(ROOT, "results", "diag", "diag_real.json"))); _dn = json.load(open(os.path.join(ROOT, "results", "diag", "diag_null.json")))
for c, nm in (("seaad", "SEAAD"), ("lee", "LEE")):
    M[f"cplMed{nm}"] = f3(_dr[c]["coupling_median"]); M[f"cplAbs{nm}"] = f3(_dr[c]["coupling_absmedian"])
    M[f"cplFracS{nm}"] = f2(_dr[c]["frac_abs_gt_0045"]); M[f"cplFracT{nm}"] = f2(_dr[c]["frac_abs_gt_01"])
    M[f"sdLib{nm}"] = f2(_dr[c]["sd_loglib_median"])
    M[f"pcOne{nm}"] = f2(_dr[c]["pc1_share"]); M[f"pcNull{nm}"] = f3(_dn[c]["pc1_null_max"])
    M[f"pcCorrMax{nm}"] = f2(max(abs(v) for v in _dr[c]["pc1_corr"].values()))
M["cplStressW"] = f3(_ds["D16_STRESS_depth"]["coupling_W_median"]); M["sdLibNull"] = f2(_dn["D16"]["sd_loglib_median"]); M["sdLibStress"] = f2(_dn["D16_STRESS_depth"]["sd_loglib_median"])
M["DeffSEAADcap"] = f2(deff["real_stats_seaad_cap40_r0.npz"]["Deff_over_D_median"])
M["DeffLEEcap"] = f2(deff["real_stats_lee_cap40_r0.npz"]["Deff_over_D_median"])
M["DeffSEAADfull"] = f2(deff["real_stats_seaad.npz"]["Deff_over_D_median"])
M["DeffLEEfull"] = f2(deff["real_stats_lee.npz"]["Deff_over_D_median"])

# real data
G = 417; npairs = G * (G - 1) // 2
M["RealG"] = str(G); M["RealPairs"] = f"{npairs:,}".replace(",", "{,}")
for coh, nm in (("seaad", "SEAAD"), ("lee", "LEE")):
    r = runs[runs.name == f"{coh}_full"].iloc[0]
    for m in ("LCP", "PB", "MC"):
        M[f"{nm}full{m}"] = f"{int(r['n_' + m]):,}".replace(",", "{,}")
    M[f"{nm}fullLCPfrac"] = f2(r["n_LCP"] / npairs)
    M[f"{nm}leakFull"] = f3(r["leakshare_med_PB"])
    rc = runs[(runs.cohort == coh) & (runs.part == "cap40")]
    M[f"{nm}leakCap"] = f3(rc["leakshare_med_PB"].mean())
for D in (24, 40):
    for m in ("LCP", "PB", "MC"):
        g = sub[(sub.D == D) & (sub.method == m)]
        M[f"prec{m}D{'twentyfour' if D == 24 else 'forty'}"] = f2(g.precision_vs_fullLCP.mean())
        M[f"rec{m}D{'twentyfour' if D == 24 else 'forty'}"] = f2(g.recall_of_fullLCP.mean())

M["capNwithin"] = f"{int(cap[cap.cohort == 'seaad'].n_within_only.iloc[0])} (SEA-AD) and {int(cap[cap.cohort == 'lee'].n_within_only.iloc[0])} (Lee)"
for coh, nm in (("seaad", "SEAAD"), ("lee", "LEE")):
    for m in ("LCP", "PB", "naive", "MC"):
        g = cap[(cap.cohort == coh) & (cap.method == m)]
        M[f"cap{m}{nm}"] = f3(g.callrate_within_only.mean())
    M[f"{nm}leakSplit"] = f3(cap[(cap.cohort == coh)].leakshare_med_PBcalls.mean())
r_ = lambda X, Y: f2(float(xm[(xm.seaad_method == X) & (xm.lee_method == Y)].repl.iloc[0]))
M["xLCPbyPB"] = r_("LCP", "PB"); M["xPBbyLCP"] = r_("PB", "LCP")
for m in ("LCP", "PB"):
    x = cc[cc.method == m].iloc[0]
    M[f"repl{m}"] = f2(x.repl_seaad_in_lee); M[f"rho{m}"] = f2(x.spearman_stat)
x = cc[cc.method == "PB_only_in_seaad->LCP_lee"].iloc[0]
M["replPBonly"] = f2(x.repl_seaad_in_lee); M["nPBonly"] = f"{int(x.n_seaad):,}".replace(",", "{,}")
sx = goi[(goi.set == "sex") & ((goi.g1 == "XIST") | (goi.g2 == "XIST")) & ~((goi.g1 == "KDM6A") | (goi.g2 == "KDM6A"))]
M["sexRhoRange"] = f"{sx.rhoB.min():.2f} to {sx.rhoB.max():.2f}"; M["sexWithinMax"] = f2(sx.r_within.abs().max())
hs = goi[goi.set == "heat_shock"]
for coh, nm in (("seaad", "SEAAD"), ("lee", "LEE")):
    h = hs[hs.cohort == coh]
    M[f"hsRho{nm}"] = f2(h.rhoB.median()); M[f"hsW{nm}"] = f2(h.r_within.median()); M[f"hsLCP{nm}"] = f2(h.LCP_BH.mean())

# held-out
if ho is not None:
    ALLHO = ("HO_seaad24", "HO_seaad24_cap40", "HO_lee40", "HO_D32", "HO_weak40")
    hod = [d for d in ALLHO if d in set(ho.design) and get(ho, d, "LCP", "R") >= 200]
    fd = [get(ho, d, "LCP", "FDR") for d in hod]
    ratio = [get(ho, d, "LCP", "power") / max(get(ho, d, "BOUND", "power"), 1e-9) for d in hod]
    pbh = [get(ho, d, "PB_logcpm", "FDR") for d in hod]
    # frozen validity claim (HELDOUT_SPEC.md): FDR <= 0.05 + 2 MCSE; type-I <= 0.05 + 3 MCSE in WW, SS, HH, NN
    fails = []
    for d in hod:
        if get(ho, d, "LCP", "FDR") > 0.05 + 2 * get(ho, d, "LCP", "FDR_mcse"):
            fails.append(f"{d} FDR")
        for c in ("WW", "SS", "HH", "NN"):
            if get(ho, d, "LCP", "rej05_" + c) > 0.05 + 3 * get(ho, d, "LCP", "rej05_" + c + "_mcse"):
                fails.append(f"{d} {c}")
    pw = [r >= 0.8 for r in ratio]
    M["HOn"] = str(len(hod))
    M["HOLCPFDRrange"] = f"{min(fd):.3f}--{max(fd):.3f}"
    M["HOratio"] = f"{min(ratio):.2f}--{max(ratio):.2f}"
    M["HOPBFDRrange"] = f"{min(pbh):.2f}--{max(pbh):.2f}"
    M["QQHO"] = f"{ho[(ho.method == 'LCP') & ho.design.isin(hod)].rej05_QQ_cross.min():.3f}--{ho[(ho.method == 'LCP') & ho.design.isin(hod)].rej05_QQ_cross.max():.3f}"
    sent = (f"With a new gene panel and new seeds, the frozen validity criterion "
            + ("held in all " + str(len(hod)) + " held-out designs with at least 24 donors" if not fails else
               "failed in " + "; ".join(x.replace("_", "\\_") for x in fails))
            + f": LCP's FDR was {M['HOLCPFDRrange']} and its power {M['HOratio']} times NMO's"
            + ("" if all(pw) else " (below the frozen 0.8 threshold in " + str(sum(not p for p in pw)) + " design(s))")
            + f", while pseudobulk's FDR was {M['HOPBFDRrange']} (the upper end from the weak-signal design, where pseudobulk makes few discoveries)."
            + (f" In the one design with real SEA-AD cell numbers (18--1{{,}}068 per donor), pseudobulk's FDR is {get(ho, 'HO_seaad24', 'PB_logcpm', 'FDR'):.3f}, and its power exceeds LCP's"
               f" ({get(ho, 'HO_seaad24', 'PB_logcpm', 'power'):.3f} vs {get(ho, 'HO_seaad24', 'LCP', 'power'):.3f})." if "HO_seaad24" in hod else ""))
    if "HO_D16_n100" in set(ho.design):
        M["HOsixteenFDR"] = fm(ho, "HO_D16_n100", "LCP", "FDR")
        sent += (f" With 16 donors but about 100 cells each, LCP's FDR was {M['HOsixteenFDR']}"
                 f" (pseudobulk {get(ho, 'HO_D16_n100', 'PB_logcpm', 'FDR'):.2f}), so the 16-donor failure in development"
                 " depends on cells per donor as well as on the number of donors.")
    if len(hod) < len(ALLHO):
        sent += " Held-out designs still running: " + ", ".join(d.replace("_", "\\_") for d in ALLHO if d not in hod) + "."
    M["HELDOUTsentence"] = sent
else:
    M["HELDOUTsentence"] = "are pending."
with open(os.path.join(P, "results.tex"), "w") as f:
    for k, v in M.items():
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

# simulation table (development)
designs = ["D8", "D12", "D16", "D24", "D40", "D16_strong", "D16_STRESS_depth"]
names = {"D8": "8", "D12": "12", "D16": "16", "D24": "24", "D40": "40", "D16_strong": "16, strong", "D16_STRESS_depth": "16, stress$^\\dagger$"}
meths = [("LCP", "LCP"), ("BOUND", "NMO"), ("ORACLE", "ORACLE"), ("PB_logcpm", "PB"), ("E5", "E5"), ("CMB_naive", "naive"), ("CMB_mom", "MC")]


def table(df, dl, nm, caption, label, meths=meths, size="\\small"):
    lines = ["\\begin{table}[t]\\centering" + size + ("\\setlength{\\tabcolsep}{2.5pt}" if "footnotesize" in size else ""), f"\\caption{{{caption}}}\\label{{{label}}}",
             "\\begin{tabular}{l" + "cc" * len(meths) + "}", "\\toprule",
             "& " + " & ".join(f"\\multicolumn{{2}}{{c}}{{{m[1]}}}" for m in meths) + " \\\\",
             "Design & " + " & ".join(["FDR & Pow"] * len(meths)) + " \\\\", "\\midrule"]
    for d in dl:
        if d not in set(df.design):
            continue
        row = [nm.get(d, d.replace("_", "\\_"))]
        for m, _ in meths:
            row += [f"{get(df, d, m, 'FDR'):.3f}", f"{get(df, d, m, 'power'):.3f}"]
        lines.append(" & ".join(row) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    return "\n".join(lines)


MAIN = [("LCP", "LCP"), ("BOUND", "NMO"), ("ORACLE", "ORACLE"), ("PB_logcpm", "Pseudobulk"), ("CMB_mom", "Moment-corr.")]
open(os.path.join(P, "table_sim_main.tex"), "w").write(table(s, designs, names, "Development simulations ($R=200$ per design, 3{,}160 pairs, ${\\sim}40$ cells per donor). FDR and power at BH $q=0.05$. MCSE of FDR $\\le$0.031 and of power $\\le$0.011. $^\\dagger$Depth--state coupling and doubled within-donor depth spread (stress design). All methods in Table~\\ref{tab:simfull}.", "tab:sim", meths=MAIN, size="") + "\n")
tab = table(s, designs, names, "Development simulations ($R=200$ per design; panel of 3{,}160 pairs). BH at $q=0.05$; MCSE of FDR $\\le$0.031 and of power $\\le$0.011 throughout (full values in the repository).", "tab:simfull", size="\\footnotesize")
open(os.path.join(P, "table_sim.tex"), "w").write(tab + "\n")
if ho is not None:
    hn = {"HO_seaad24": "SEA-AD 24", "HO_seaad24_cap40": "SEA-AD 24, 40 cells", "HO_lee40": "Lee 40", "HO_D32": "32",
          "HO_weak40": "40, weak", "HO_D16_n100": "16, 100 cells"}
    capm = (f"Frozen held-out evaluation ($R=200$ per design; new gene panel and seeds). SEA-AD and Lee designs use real donors' cell numbers or depths; "
            f"MCSE of FDR $\\le${ho.FDR_mcse.max():.3f} and of power $\\le${ho.power_mcse.max():.3f}. All methods in Table~\\ref{{tab:heldout}}.")
    open(os.path.join(P, "table_heldout_main.tex"), "w").write(table(ho, list(hn), hn, capm, "tab:heldoutmain", meths=MAIN, size="") + "\n")
    tabh = table(ho, list(hn), hn, f"Frozen held-out evaluation ($R=200$ per design; new gene panel and seeds; the SEA-AD and Lee designs use real donors' cell numbers or depths). MCSE of FDR $\\le${ho.FDR_mcse.max():.3f} and of power $\\le${ho.power_mcse.max():.3f}. Weak-null type-I error of LCP at $\\alpha=0.05$: {ho[ho.method=='LCP'].rej05_QQ_cross.min():.3f}--{ho[ho.method=='LCP'].rej05_QQ_cross.max():.3f}.", "tab:heldout", size="\\footnotesize")
else:
    tabh = "\\begin{table}[t]\\centering\\caption{Held-out evaluation pending.}\\label{tab:heldout}\\end{table}"
open(os.path.join(P, "table_heldout.tex"), "w").write(tabh + "\n")
print(len(M), "macros")

# ---- interval coverage table (appendix) ----
cv = pd.read_csv(os.path.join(ROOT, "results", "dev", "coverage_boot_t_wide.csv"))
cdes = ["D8", "D12", "D16", "D16_strong", "D24", "D40"]
clab = {"BB": "between only", "MM": "both levels, same sign", "OO_same": "both levels, signs aligned",
        "OO_opp": "both levels, signs opposed", "WW": "within only (moderate)", "SS": "within only (strong)",
        "HH": "sign-heterogeneous within", "QQ_same": "weak null, same block", "QQ_cross": "weak null, cross block",
        "NN": "independent"}
lines = [r"\begin{table}[h]\centering\footnotesize", r"\caption{Coverage of nominal 95\% donor bootstrap-$t$ intervals for $\SB_{jk}$ (deployed BOOT\_T, with the same interval computed from the true donor means in parentheses). $R=200$; MCSE $\le$0.02.}\label{tab:cov}",
         r"\begin{tabular}{l" + "c" * len(cdes) + "}", r"\hline", "Pair class & " + " & ".join(d.replace("_strong", " strong").replace("D", "$D=") + "$" for d in cdes) + r"\\", r"\hline"]
for k, lab in clab.items():
    row = [lab]
    for d in cdes:
        g = cv[(cv.design == d) & (cv.cls == k)]
        row.append(f"{g.BOOT_T.iloc[0]:.2f} ({g.ORACLE_BOOT_T.iloc[0]:.2f})" if len(g) else "--")
    lines.append(" & ".join(row) + r"\\")
lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
open(os.path.join(P, "table_cov.tex"), "w").write("\n".join(lines) + "\n")
nbad = int((cv[cv.design.isin(cdes + ["D16_STRESS_depth"])]["diff"] < -0.02).sum())
pass

# ---- post hoc ablations (appendix) ----
abp = os.path.join(ROOT, "results", "ablate", "summary.csv")
if os.path.exists(abp):
    ab = pd.read_csv(abp)
    AV = [("LCP_gpu", "LCP"), ("LCP_noleak", "no leak term"), ("LCP_truev", "true variances"), ("LCP_unstd", "unstandardized"), ("PB_perm", "pseudobulk + perm.")]
    AD = [d for d in ("D16", "D24", "D40") if d in set(ab.design)]
    L = [r"\begin{table}[h]\centering\footnotesize\setlength{\tabcolsep}{3pt}",
         r"\caption{Post hoc ablations (not in the analysis plan), on the development replicates ($R=%d$ per design, ${\sim}40$ cells per donor). All variants use the same $10^5$ donor permutations without refinement. FDR and power at BH $q=0.05$; SS, type-I error at $\alpha=0.05$ for strong within-only pairs. MCSE of FDR $\le$%.3f, of power $\le$%.3f, of SS $\le$%.3f.}\label{tab:ablate}" % (int(ab.R.min()), ab.FDR_mcse.max(), ab.power_mcse.max(), ab.rej05_SS_mcse.max()),
         r"\begin{tabular}{l" + "ccc" * len(AD) + "}", r"\toprule",
         "& " + " & ".join(r"\multicolumn{3}{c}{$D=%s$}" % d[1:] for d in AD) + r" \\",
         "Variant & " + " & ".join(["FDR & Pow & SS"] * len(AD)) + r" \\", r"\midrule"]
    for k, lab in AV:
        row = [lab]
        for d in AD:
            row += [f"{get(ab, d, k, 'FDR'):.3f}", f"{get(ab, d, k, 'power'):.3f}", f"{get(ab, d, k, 'rej05_SS'):.3f}"]
        L.append(" & ".join(row) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(P, "table_ablate.tex"), "w").write("\n".join(L) + "\n")
    _ss = [get(ab, d, "LCP_gpu", "rej05_SS") for d in AD]
    with open(os.path.join(P, "results.tex"), "a") as f:
        f.write(f"\\newcommand{{\\ABfullSSrange}}{{{min(_ss):.3f}--{max(_ss):.3f}}}\n")
        for d in AD:
            for k, _ in AV:
                for col, tag in (("FDR", "F"), ("power", "P"), ("rej05_SS", "S")):
                    nm = "AB" + k.replace("_", "") + tag + {"D16": "sixteen", "D24": "twentyfour", "D40": "forty"}[d]
                    f.write(f"\\newcommand{{\\{nm}}}{{{get(ab, d, k, col):.3f}}}\n")
else:
    open(os.path.join(P, "table_ablate.tex"), "w").write("% ablations pending\n")

# ---- design map (main text): where LCP is calibrated, by donors x cells per donor ----
_rs = json.load(open(os.path.join(ROOT, "results", "rare", "rare_sim_summary.json"))); _rd = json.load(open(os.path.join(ROOT, "results", "diag", "diag_rare_deff.json")))
_dsj = json.load(open(os.path.join(ROOT, "results", "diag", "diag_sim.json")))
rows = []
def _row(lab, D, cells, src, df, key, deff):
    if src == "rare":
        g = _rs[key]; rows.append((lab, D, cells, deff, g["LCP"]["fdp"][0], g["LCP"]["fdp"][1], g["PB"]["fdp"][0], g["LCP"]["pw"][0]))
    else:
        rows.append((lab, D, cells, deff, get(df, key, "LCP", "FDR"), get(df, key, "LCP", "FDR_mcse"), get(df, key, "PB_logcpm", "FDR"), get(df, key, "LCP", "power")))
_row("SEA-AD Micro-PVM\\_1 counts$^a$", _rs["sim_S1"]["LCP"]["R"] and len(json.load(open(os.path.join(ROOT, "results", "rare", "rare_prep_meta.json")))["S1"]["n_list"]), "15", "rare", None, "sim_S1", _rd["S1"]["Deff_over_D_median"])
_row("Development", 16, "${\\sim}40$", "dev", s, "D16", _dsj["D16"]["Deff_over_D_median"])
_row("Development", 24, "${\\sim}40$", "dev", s, "D24", _dsj["D24"]["Deff_over_D_median"])
if ho is not None:
    _row("Held-out, SEA-AD depths", 24, "${\\sim}40$", "ho", ho, "HO_seaad24_cap40", _dsj["HO_seaad24_cap40"]["Deff_over_D_median"])
    _row("Held-out", 32, "${\\sim}40$", "ho", ho, "HO_D32", _dsj["HO_D32"]["Deff_over_D_median"])
_row("Development", 40, "${\\sim}40$", "dev", s, "D40", _dsj["D40"]["Deff_over_D_median"])
_row("SEA-AD Micro-PVM\\_3 counts$^a$", len(json.load(open(os.path.join(ROOT, "results", "rare", "rare_prep_meta.json")))["S2"]["n_list"]), "43", "rare", None, "sim_S2", _rd["S2"]["Deff_over_D_median"])
if ho is not None:
    _row("Held-out, Lee depths", 40, "${\\sim}60$", "ho", ho, "HO_lee40", _dsj["HO_lee40"]["Deff_over_D_median"])
    _row("Held-out", 16, "${\\sim}100$", "ho", ho, "HO_D16_n100", _dsj["HO_D16_n100"]["Deff_over_D_median"])
    _row("Held-out, SEA-AD counts", 24, "18--1{,}068", "ho", ho, "HO_seaad24", _dsj["HO_seaad24"]["Deff_over_D_median"])
L = [r"\begin{table}[t]\centering\small\setlength{\tabcolsep}{4pt}",
     r"\caption{Where LCP is calibrated. Simulated designs ordered by cells per donor ($R=100$ for the two SEA-AD supertype designs, $R=200$ otherwise). FDR at BH $q=0.05$ with MCSE; \checkmark: FDR $\le0.05+2\,$MCSE. $\Deff/D$: median over between-donor pairs. $^a$Exact per-donor cell numbers of the real supertype (Section~\ref{sec:real}).}\label{tab:map}",
     r"\begin{tabular}{lrcccccc}", r"\toprule", r"Design & $D$ & cells/donor & $\Deff/D$ & LCP FDR & & PB FDR & LCP power \\", r"\midrule"]
for lab, D, cells, deff, fd, fse, pbf, pw in rows:
    ok = r"\checkmark" if fd <= 0.05 + 2 * fse else r"$\times$"
    L.append(f"{lab} & {D} & {cells} & {deff:.2f} & {fd:.3f}\\,({fse:.3f}) & {ok} & {pbf:.3f} & {pw:.3f} \\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
open(os.path.join(P, "table_map.tex"), "w").write("\n".join(L) + "\n")
_areall = [v["ARE_median"] for k, v in _dsj.items() if "STRESS" not in k] + [v["ARE_median"] for v in _rd.values()]
with open(os.path.join(P, "results.tex"), "a") as f:
    f.write(f"\\renewcommand{{\\AREmin}}{{{min(_areall):.3f}}}\\renewcommand{{\\AREmax}}{{{max(_areall):.3f}}}\n")
    f.write(f"\\newcommand{{\\SoneFDR}}{{{_rs['sim_S1']['LCP']['fdp'][0]:.3f}}}\\newcommand{{\\SonePBFDR}}{{{_rs['sim_S1']['PB']['fdp'][0]:.2f}}}\n")
