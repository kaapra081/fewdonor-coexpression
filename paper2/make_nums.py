"""Generates nums.tex, table_sim.tex, table_null.tex, table_ablate.tex and paper2/dc_summary.json from result files."""
import os, sys, json, glob, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..")
sys.path.insert(0, os.path.join(R, "dev")); from methods_v2 import bh  # noqa: E402
J = lambda p: json.load(open(os.path.join(R, p)))
pop = J("pop/pop_results.json"); dres = J("dstd/d_results.json"); xs = J("xsplit/xs_results.json"); fx = J("fix/fix_results.json")
xu = J("xu/xu_results.json"); sp = J("split/split_results.json"); dg = J("diag/diag_real.json")
info = {s: J(f"pop/pop_{s}_info.json") for s in ("P1", "P2", "P3")}
S = ("P1", "P2", "P3"); SUBN = {"P1": "Homeo\\_CECR2", "P2": "PVM\\_CD163", "P3": "ADAM\\_GPNMB"}
c = lambda x: f"{int(round(x)):,}".replace(",", "{,}")
f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"
def rng(vals, fmt=c):
    lo, hi = min(vals), max(vals); return fmt(lo) if fmt(lo) == fmt(hi) else f"{fmt(lo)}--{fmt(hi)}"
A = pop["A"]; M = {}
M["NnucPsy"] = "289{,}279"
M["DpsyRange"] = rng([info[s]["donors"] for s in S]); M["NdMedRange"] = rng([info[s]["median_nd"] for s in S])
for i, s in enumerate(S):
    M["D" + ["Pone", "Ptwo", "Pthree"][i]] = c(info[s]["donors"]); M["Med" + ["Pone", "Ptwo", "Pthree"][i]] = c(info[s]["median_nd"])
full = {s: A[s]["all"] for s in S}
M["PBnullCalls"] = rng([full[s]["PB"]["mean_calls"] for s in S]); M["PBnullTI"] = rng([full[s]["PB"]["typeI"] for s in S], f2)
M["LCPnullCalls"] = rng([full[s]["LCP"]["mean_calls"] for s in S]); M["LCPnullTI"] = rng([full[s]["LCP"]["typeI"] for s in S], f2)
M["XDSnullTI"] = rng([full[s]["LCPXd"]["typeI"] for s in S], f3); M["XDSnullTop"] = rng([full[s]["LCPXd"]["typeI_top1"] for s in S], f3)
M["XDSnullMaxCall"] = c(max(max(full[s]["LCPXd"]["calls"]) for s in S))
M["PBrepFalse"] = rng([v for s in S for v in pop["rep_null"][s]["PB"]["replicated_false"]])
M["XDSrepFalse"] = c(max(v for s in S for v in pop["rep_null"][s]["LCPXd"]["replicated_false"]))
M["PBgrowA"] = rng([A[s]["25"]["PB"]["mean_calls"] for s in S]); M["PBgrowC"] = rng([A[s]["200"]["PB"]["mean_calls"] for s in S])
M["LCPgrowB"] = rng([A[s][D]["LCP"]["mean_calls"] for s in S for D in ("25", "50", "100")])
M["XUnullCalls"] = rng([xu["E1"][s]["all"]["XU"]["mean_calls"] for s in S])
M["XDSfreshCalls"] = rng([xu["E1"][s]["all"]["XDS"]["mean_calls"] for s in S], lambda x: f"{x:.1f}")
E3 = xu["E3"]; E4 = xu["E4"]
for k, d in zip("abc", ("L200", "L500", "L1000")):
    M["LPBFDR" + k] = f2(E3[d]["PB"]["FDR"]); M["LLCPFDR" + k] = f3(E3[d]["LCP"]["FDR"])
M["LXDSFDRrange"] = rng([E3[d]["XDS"]["FDR"] for d in E3], f3); M["LXDSpowRange"] = rng([E3[d]["XDS"]["power"] for d in E3], f3)
M["SXDSFDRrange"] = rng([E4[d]["XDS"]["FDR"] for d in E4], f3)
M["SXDSpowDforty"] = f3(E4["F_D40"]["XDS"]["power"]); M["SLCPpowDforty"] = f3(E4["F_D40"]["LCP"]["power"])
dn = dres["null"]
M["SeaXDSTI"] = rng([dn[s]["LCPXd"]["typeI"] for s in ("S1", "S2")], f3); M["SeaXDSTop"] = rng([dn[s]["LCPXd"]["typeI_top1"] for s in ("S1", "S2")], f3)
M["SeaLCPTop"] = rng([fx["null"][s]["LCP"]["typeI_top1"] for s in ("S1", "S2")], f3)
M["SeaPBcalls"] = rng([fx["null"][s]["PB"]["mean_calls"] for s in ("S1", "S2")])
B = pop["B"]; G = {"P1": 89676, "P2": 87990, "P3": 88831}
M["RealXDSfrac"] = rng([100 * B[s]["LCPXd"] / G[s] for s in S], lambda x: f"{x:.0f}") + "\\%"
M["RealXDSrep"] = rng([100 * B[s]["rep_LCPXd"]["replicated"] / B[s]["rep_LCPXd"]["calls_h0"] for s in S], lambda x: f"{x:.0f}") + "\\%"
M["RealPBrep"] = rng([100 * B[s]["rep_PB"]["replicated"] / B[s]["rep_PB"]["calls_h0"] for s in S], lambda x: f"{x:.0f}") + "\\%"
M["RealPBonly"] = rng([B[s]["PB_only"] for s in S])
# coupling and transformation diagnostics
M["CouplingSEAAD"] = f3(dg["seaad"]["coupling_median"]); M["CouplingFrac"] = f"{100 * dg['seaad']['frac_abs_gt_0045']:.0f}\\%"
M["CouplingFracLee"] = f"{100 * dg['lee']['frac_abs_gt_0045']:.0f}\\%"
M["Rcorrln"] = "0.155"; M["Mcorrln"] = "0.000"   # from the disclosed exploratory diagnostic in DSTD_PREREG.md (P3, k = 0)
M["PBposFrac"] = rng([100 * v for s in S for v in sp[s]["PB_false_frac_positive"]], lambda x: f"{x:.0f}") + "\\%"
cp = J("xu/couple_results.json")
fm = lambda m: ", ".join(f3(cp[d][m]["FDR"]) for d in ("LS200", "LS500", "LS1000"))
for m in ("PB", "LCP", "XDS", "XU"): M["Couple" + m] = fm(m)
# differential co-expression
dc = {}
for s in S:
    nul = [np.load(f) for f in sorted(glob.glob(os.path.join(R, "dc", "out", f"DC_{s}_k*.npz")))]
    real = np.load(os.path.join(R, "dc", "out", f"DC_{s}_real.npz"))
    dc[s] = {m: dict(null_calls=[int(bh(z[m].astype(float)).sum()) for z in nul], null_typeI=[float((z[m] < 0.05).mean()) for z in nul],
                     real_calls=int(bh(real[m].astype(float)).sum())) for m in ("PB", "LCP", "XDS")}
    dc[s]["n1"] = int(real["n1"][0]); dc[s]["n0"] = int(real["n0"][0])
json.dump(dc, open(os.path.join(H, "dc_summary.json"), "w"), indent=1)
pbm = {s: np.mean(dc[s]["PB"]["null_calls"]) for s in S}
dcyes = sum(pbm[s] >= 100 for s in S) >= 2
xdsok = all(np.mean(dc[s]["XDS"]["null_calls"]) <= 1 and max(dc[s]["XDS"]["null_typeI"]) <= 0.055 for s in S)
M["DCtext"] = (
    "Cell numbers differ between AD cases and controls in all three PsychAD subtypes (medians 21 versus 36 nuclei in Homeo\\_CECR2, "
    "23 versus 18 in PVM\\_CD163 and 17 versus 10 in ADAM\\_GPNMB; Mann--Whitney $p\\le8\\times10^{-5}$). By Propositions~\\ref{prop:leak} and~\\ref{prop:transform}, "
    "groups with different cell numbers have different biases, so comparing co-expression between groups should produce false "
    "differences. We tested this in the real-cell null, with each pseudo-donor carrying the AD label of the real donor whose cell number it has "
    "(three reassignments per subtype; frozen plan DC\\_PREREG). Pseudobulk differential co-expression (difference of Fisher-transformed correlations) made "
    f"{rng([v for s in S for v in dc[s]['PB']['null_calls']])} false calls per reassignment (per-pair type-I {rng([v for s in S for v in dc[s]['PB']['null_typeI']], f2)}), "
    f"so the pre-registered decision that pseudobulk differential co-expression is confounded by cell number is {'YES' if dcyes else 'NO'}. "
    f"XDS differential co-expression (Welch test of the per-donor cross-split products between groups) made {rng([v for s in S for v in dc[s]['XDS']['null_calls']])} "
    f"(type-I {rng([v for s in S for v in dc[s]['XDS']['null_typeI']], f3)}; decision {'YES' if xdsok else 'NO'}), and the leak-corrected version "
    f"made {rng([v for s in S for v in dc[s]['LCP']['null_calls']])} (type-I {rng([v for s in S for v in dc[s]['LCP']['null_typeI']], f3)}). "
    f"On real donors pseudobulk reported {c(dc['P1']['PB']['real_calls'])}, {c(dc['P2']['PB']['real_calls'])} and {c(dc['P3']['PB']['real_calls'])} AD-associated "
    f"co-expression changes in Homeo\\_CECR2, PVM\\_CD163 and ADAM\\_GPNMB. With only cell-number differences it makes {c(pbm['P1'])}, {c(pbm['P2'])} and {c(pbm['P3'])} on average, "
    f"that is {100 * pbm['P1'] / dc['P1']['PB']['real_calls']:.0f}\\%, {100 * pbm['P2'] / dc['P2']['PB']['real_calls']:.0f}\\% and {100 * pbm['P3'] / dc['P3']['PB']['real_calls']:.0f}\\% of those numbers, "
    f"so in the first two subtypes pseudobulk differential co-expression cannot be distinguished from the cell-number artifact. XDS reported {c(dc['P1']['XDS']['real_calls'])}, {c(dc['P2']['XDS']['real_calls'])} and {c(dc['P3']['XDS']['real_calls'])}. "
    "The caveat of Section~\\ref{sec:null} applies: the pooled null has somewhat stronger within-donor structure than real donors.")
M["SplitText"] = ("XDS uses one random split per donor. On real PsychAD donors we repeated XDS with three split seeds (post hoc). The numbers of BH calls were "
    + "; ".join(f"{SUBN[s]} {', '.join(c(v) for v in sp[s]['n_calls'])}" for s in S)
    + ", and the Jaccard index between call sets from different seeds was "
    + "; ".join(f"{rng(sp[s]['jaccard'], f2)} ({SUBN[s]})" for s in S)
    + ". Call counts are stable, but a fifth to a third of individual calls depend on the seed; aggregating several splits (for example by a Cauchy combination of $p$-values) is a natural extension that we have not evaluated.")
with open(os.path.join(H, "nums.tex"), "w") as fh:
    for k, v in M.items(): fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
# ---------------- tables
def cell(d, m, key="FDR"):
    v = d[m]; return f"{v[key]:.3f}\\,({v[key + '_mcse']:.3f})"
rows = []
for lab, src, d in (("16 donors, $\\sim$40 cells", E4, "F_D16"), ("24 donors, $\\sim$40 cells", E4, "F_D24"), ("40 donors, $\\sim$40 cells", E4, "F_D40"),
                    ("79 donors, Micro-PVM\\_1 cells (median 15)", E4, "F_S1"), ("83 donors, Micro-PVM\\_3 cells (median 43)", E4, "F_S2"), ("16 donors, $\\sim$100 cells", E4, "F_D16_n100"),
                    ("200 donors, median 20 cells", E3, "L200"), ("500 donors, median 20 cells", E3, "L500"), ("1,000 donors, median 20 cells", E3, "L1000")):
    r = src[d]; rows.append(" & ".join([lab] + [cell(r, m) for m in ("PB", "LCP", "XDS")] + [f"{r[m]['power']:.3f}" for m in ("PB", "LCP", "XDS")]) + " \\\\")
open(os.path.join(H, "table_sim.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{4pt}\\caption{Simulations with exact truth. FDR (MCSE) at BH $q=0.05$ and power. $R=100$ replicates for designs with up to 83 donors, $R=50$ for larger designs. "
    "Pseudobulk and LCP lose FDR control as donors are added; XDS keeps it, and its power deficit disappears from about 500 donors.}\\label{tab:sim}\n"
    "\\begin{tabular}{lcccccc}\\toprule & \\multicolumn{3}{c}{FDR} & \\multicolumn{3}{c}{Power}\\\\\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
    "Design & PB & LCP & XDS & PB & LCP & XDS\\\\\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\\end{tabular}\\end{table}\n")
rows = []
for s in S:
    a = full[s]; rp = pop["rep_null"][s]
    rows.append(f"{SUBN[s]} ({c(info[s]['donors'])}; {c(info[s]['median_nd'])}) & " + " & ".join(
        f"{c(a[m]['mean_calls'])} & {a[m]['typeI']:.3f} & {a[m]['typeI_top1']:.3f}" for m in ("PB", "LCP", "LCPXd")) +
        f" & {c(np.mean(rp['PB']['replicated_false']))} & {c(np.mean(rp['LCPXd']['replicated_false']))} \\\\")
open(os.path.join(H, "table_null.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3.5pt}\\caption{PsychAD real-cell null at full cohort size (pseudo-donors; median nuclei), mean over three reassignments. "
    "Calls: BH false discoveries among $\\sim$89{,}000 pairs. TI: per-pair type-I error at 0.05, over all pairs and over the 1\\% most within-correlated pairs. Rep.: false calls replicated (same sign, $p<0.05$) in an independent half of the pseudo-donors.}\\label{tab:null}\n"
    "\\begin{tabular}{lccccccccccc}\\toprule & \\multicolumn{3}{c}{PB} & \\multicolumn{3}{c}{LCP} & \\multicolumn{3}{c}{XDS} & \\multicolumn{2}{c}{Rep.}\\\\"
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-12}\n Subtype & Calls & TI & top 1\\% & Calls & TI & top 1\\% & Calls & TI & top 1\\% & PB & XDS\\\\\\midrule\n"
    + "\n".join(rows) + "\n\\bottomrule\\end{tabular}\\end{table}\n")
# ablation table: which bias each method leaves, with evidence
ab = [("PB", "log-CPM, own cells", "--", "leak, transform", "PB"),
      ("PBlin", "ratio mean, own cells", "--", "leak", "PBlin"),
      ("LCP", "self-standardized", "subtract estimate", "transform, leak-estimate error", "LCP"),
      ("LCPd", "design-only", "subtract estimate", "leak-estimate error", "LCPd"),
      ("LCPX", "self-standardized", "disjoint halves", "transform (larger)", "LCPX"),
      ("XU", "design-only weights", "exact U-statistic", "depth--state coupling", "XU"),
      ("XDS", "design-only", "disjoint halves", "none detected", "LCPXd")]
rows = []
for nm, stdz, leak, left, key in ab:
    if key == "XU":
        calls = rng([xu["E1"][s]["all"]["XU"]["mean_calls"] for s in S]); ti = rng([xu["E1"][s]["all"]["XU"]["typeI"] for s in S], f3)
    else:
        calls = rng([full[s][key]["mean_calls"] for s in S]); ti = rng([full[s][key]["typeI"] for s in S], f3)
    if key in ("PB", "LCP"): sea = rng([fx["null"][s][key]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key == "LCPX": sea = rng([xs["null"][s]["LCPX"]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key in ("LCPd", "LCPXd", "PBlin"): sea = rng([dres["null"][s][key]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key == "XU": sea = rng([xu["E2"][s]["typeI_top1"] for s in ("S1", "S2")], f3)
    rows.append(f"{nm} & {stdz} & {leak} & {left} & {calls} & {ti} & {sea} \\\\")
open(os.path.join(H, "table_ablate.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3pt}\\caption{Each correction and the bias it leaves. Calls and TI: PsychAD real-cell null at full size (range over subtypes; XU on fresh reassignments). "
    "SEA-AD: type-I for the 1\\% most within-correlated pairs in the SEA-AD real-cell null (79 and 83 pseudo-donors). LCP2 (not shown) had SEA-AD top-1\\% type-I "
    + rng([fx["null"][s]["LCP2"]["typeI_top1"] for s in ("S1", "S2")], f3) + ".}\\label{tab:ablate}\n"
    "\\resizebox{\\linewidth}{!}{\\begin{tabular}{lllllll}\\toprule Method & Standardization & Leak handling & Bias left & Calls & TI & SEA-AD\\\\\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\\end{tabular}}\\end{table}\n")
print(len(M), "macros"); print({k: v for k, v in M.items() if k not in ("DCtext", "SplitText")}); print(M["DCtext"][-700:])
