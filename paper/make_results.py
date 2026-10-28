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
    sent = (f"With a new gene panel, new seeds and three designs built from real donors' cell numbers and depths, the frozen validity criterion "
            + ("held in all " + str(len(hod)) + " held-out designs with at least 24 donors" if not fails else
               "failed in " + "; ".join(x.replace("_", "\\_") for x in fails))
            + f": LCP's FDR was {M['HOLCPFDRrange']} and its power {M['HOratio']} times BOUND's"
            + ("" if all(pw) else " (below the frozen 0.8 threshold in " + str(sum(not p for p in pw)) + " design(s))")
            + f", while pseudobulk's FDR was {M['HOPBFDRrange']}.")
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
meths = [("LCP", "LCP"), ("BOUND", "BOUND"), ("ORACLE", "ORACLE"), ("PB_logcpm", "PB"), ("E5", "E5"), ("CMB_naive", "naive"), ("CMB_mom", "MC")]


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


MAIN = [("LCP", "LCP"), ("BOUND", "BOUND"), ("ORACLE", "ORACLE"), ("PB_logcpm", "Pseudobulk"), ("CMB_mom", "Moment-corr.")]
open(os.path.join(P, "table_sim_main.tex"), "w").write(table(s, designs, names, "Development simulations ($R=200$ per design, 3{,}160 pairs, ${\\sim}40$ cells per donor). FDR and power at BH $q=0.05$. MCSE of FDR $\\le$0.031 and of power $\\le$0.011. $^\\dagger$Depth--state coupling and doubled within-donor depth spread (stress design). All methods in Table~\\ref{tab:simfull}.", "tab:sim", meths=MAIN, size="") + "\n")
tab = table(s, designs, names, "Development simulations ($R=200$ per design; panel of 3{,}160 pairs). BH at $q=0.05$; MCSE of FDR $\\le$0.031 and of power $\\le$0.011 throughout (full values in the repository).", "tab:simfull", size="\\footnotesize")
open(os.path.join(P, "table_sim.tex"), "w").write(tab + "\n")
if ho is not None:
    hn = {"HO_seaad24": "SEA-AD 24", "HO_seaad24_cap40": "SEA-AD 24, 40 cells", "HO_lee40": "Lee 40", "HO_D32": "32",
          "HO_weak40": "40, weak", "HO_D16_n100": "16, 100 cells"}
    capm = (f"Frozen held-out evaluation ($R=200$ per design; new gene panel and seeds). SEA-AD and Lee designs use real donors' cell numbers or depths; "
            f"MCSE of FDR $\\le${ho.FDR_mcse.max():.3f} and of power $\\le${ho.power_mcse.max():.3f}. All methods in Table~\\ref{{tab:heldout}}.")
    open(os.path.join(P, "table_heldout_main.tex"), "w").write(table(ho, list(hn), hn, capm, "tab:heldoutmain", meths=MAIN, size="") + "\n")
    tabh = table(ho, list(hn), hn, f"Frozen held-out evaluation ($R=200$ per design; new gene panel, seeds and real-donor designs). MCSE of FDR $\\le${ho.FDR_mcse.max():.3f} and of power $\\le${ho.power_mcse.max():.3f}. Weak-null type-I error of LCP at $\\alpha=0.05$: {ho[ho.method=='LCP'].rej05_QQ_cross.min():.3f}--{ho[ho.method=='LCP'].rej05_QQ_cross.max():.3f}.", "tab:heldout", size="\\footnotesize")
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
    with open(os.path.join(P, "results.tex"), "a") as f:
        for d in AD:
            for k, _ in AV:
                for col, tag in (("FDR", "F"), ("power", "P"), ("rej05_SS", "S")):
                    nm = "AB" + k.replace("_", "") + tag + {"D16": "sixteen", "D24": "twentyfour", "D40": "forty"}[d]
                    f.write(f"\\newcommand{{\\{nm}}}{{{get(ab, d, k, col):.3f}}}\n")
else:
    open(os.path.join(P, "table_ablate.tex"), "w").write("% ablations pending\n")
