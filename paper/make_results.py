"""Generate results.tex (LaTeX macros), appendix_tables.tex and the simulation tables from the result CSVs,
so that every number in the paper is copied programmatically from the saved results."""
import os, json
import numpy as np, pandas as pd

V2 = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = os.path.join(V2, "paper")
s = pd.read_csv(os.path.join(V2, "results", "dev", "summary_all.csv"))
ho_path = os.path.join(V2, "results", "heldout", "summary_all.csv")
ho = pd.read_csv(ho_path) if os.path.exists(ho_path) else None
R = os.path.join(V2, "results", "real")
cal = json.load(open(os.path.join(V2, "results", "real_between_variance_calibration.json")))
deff = json.load(open(os.path.join(R, "deff_real.json")))
runs = pd.read_csv(os.path.join(R, "runs_summary.csv")); sub = pd.read_csv(os.path.join(R, "subsample_vs_full.csv"))
cap = pd.read_csv(os.path.join(R, "cap40_vs_full.csv")); cc = pd.read_csv(os.path.join(R, "cross_cohort.csv"))
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
for m in ("LCP", "PB", "MC", "naive"):
    g = cap[cap.method == m]
    M[f"cap{m}"] = f3(g.callrate_within_only.mean())
M["capNwithin"] = f"{int(cap[cap.cohort == 'seaad'].n_within_only_ref.iloc[0])} (SEA-AD) and {int(cap[cap.cohort == 'lee'].n_within_only_ref.iloc[0])} (Lee)"
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
    hod = [d for d in ("HO_seaad24", "HO_seaad24_cap40", "HO_lee40", "HO_D32", "HO_weak40") if d in set(ho.design)]
    fd = [get(ho, d, "LCP", "FDR") for d in hod]
    ratio = [get(ho, d, "LCP", "power") / max(get(ho, d, "BOUND", "power"), 1e-9) for d in hod]
    pbh = [get(ho, d, "PB_logcpm", "FDR") for d in hod]
    M["HOn"] = str(len(hod))
    M["HOLCPFDRrange"] = f"{min(fd):.3f}--{max(fd):.3f}"
    M["HOratio"] = f"{min(ratio):.2f}--{max(ratio):.2f}"
    M["HOPBFDRrange"] = f"{min(pbh):.2f}--{max(pbh):.2f}"
    if "HO_D16_n100" in set(ho.design):
        M["HOsixteenFDR"] = fm(ho, "HO_D16_n100", "LCP", "FDR")
    M["HELDOUTsentence"] = ("confirm this on new gene panels, seeds and real-donor designs: across the "
                            f"{len(hod)} held-out designs with at least 24 donors, LCP's FDR is {M['HOLCPFDRrange']} and its power is "
                            f"{M['HOratio']} times BOUND's, while pseudobulk's FDR is {M['HOPBFDRrange']}.")
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
    lines = ["\\begin{table}[t]\\centering" + size, f"\\caption{{{caption}}}\\label{{{label}}}",
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
open(os.path.join(P, "table_sim_main.tex"), "w").write(table(s, designs, names, "Development simulations ($R=200$ per design, 3{,}160 pairs, ${\\sim}40$ cells per donor). FDR and power at BH $q=0.05$. MCSE of FDR $\\le$0.03 and of power $\\le$0.011. $^\\dagger$Depth--state coupling and doubled within-donor depth spread (stress design). All methods in Table~\\ref{tab:simfull}.", "tab:sim", meths=MAIN, size="") + "\n")
tab = table(s, designs, names, "Development simulations ($R=200$ per design; panel of 3{,}160 pairs). BH at $q=0.05$; MCSE of FDR $\\le$0.03 and of power $\\le$0.013 throughout (full values in the repository).", "tab:simfull", size="\\footnotesize")
open(os.path.join(P, "table_sim.tex"), "w").write(tab + "\n")
if ho is not None:
    hn = {"HO_seaad24": "SEA-AD 24", "HO_seaad24_cap40": "SEA-AD 24, 40 cells", "HO_lee40": "Lee 40", "HO_D32": "32",
          "HO_weak40": "40, weak", "HO_D16_n100": "16, 100 cells"}
    tabh = table(ho, list(hn), hn, "Frozen held-out evaluation ($R=200$ per design; new panel, seeds and real-donor designs).", "tab:heldout")
else:
    tabh = "\\begin{table}[t]\\centering\\caption{Held-out evaluation pending.}\\label{tab:heldout}\\end{table}"
open(os.path.join(P, "table_heldout.tex"), "w").write(tabh + "\n")
print(len(M), "macros")
