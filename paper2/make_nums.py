"""Generates nums.tex, table_sim.tex, table_null.tex, table_ablate.tex and paper2/dc_summary.json from result files."""
import os, sys, json, glob, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..")
sys.path.insert(0, os.path.join(R, "dev")); from methods_v2 import bh  # noqa: E402
def J(p):
    """Local layout (recomb_v2/<dir>/<file>) or repository layout (results/<dir>/<file>)."""
    for cand in (os.path.join(R, p), os.path.join(R, "results", p)):
        if os.path.exists(cand): return json.load(open(cand))
    raise FileNotFoundError(p)
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

# ---------------- review-response analyses (REV_PREREG, ONEK_PREREG)
import math
rv = os.path.join(R, "rev")
r1 = {}
for s_ in S:
    Z = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(rv, "out", f"R1_{s_}_k*.json")))]
    if not Z: continue
    r1[s_] = {m: dict(n=len(Z), runs_any=int(sum(z[m]["V"] > 0 for z in Z)), frac=float(np.mean([z[m]["V"] > 0 for z in Z])),
                      se=float(math.sqrt(0.05 * 0.95 / len(Z))), n5=float(np.mean([z[m]["n5"] for z in Z])), n4=float(np.mean([z[m]["n4"] for z in Z])),
                      t05=float(np.mean([z[m]["t05"] for z in Z])), V=float(np.mean([z[m]["V"] for z in Z]))) for m in ("PB", "XDS", "XDSM", "KX")}
    r1[s_]["P"] = Z[0]["P"]
    for m in ("XDS", "XDSM", "KX"): r1[s_][m]["holds"] = bool(r1[s_][m]["frac"] <= 0.05 + 2 * math.sqrt(0.05 * 0.95 / len(Z)))
json.dump(r1, open(os.path.join(H, "r1_summary.json"), "w"), indent=1)
if r1:
    M["RoneXDSruns"] = rng([r1[s_]["XDS"]["runs_any"] for s_ in r1]); M["RoneN"] = rng([r1[s_]["XDS"]["n"] for s_ in r1])
    M["RoneXDSMruns"] = rng([r1[s_]["XDSM"]["runs_any"] for s_ in r1]); M["RoneKXruns"] = rng([r1[s_]["KX"]["runs_any"] for s_ in r1])
    M["RoneXDSnfive"] = rng([r1[s_]["XDS"]["n5"] for s_ in r1], lambda x: f"{x:.2f}"); M["RoneExpnfive"] = rng([r1[s_]["P"] * 1e-5 for s_ in r1], lambda x: f"{x:.2f}")
    M["RoneXDSnfour"] = rng([r1[s_]["XDS"]["n4"] for s_ in r1], lambda x: f"{x:.1f}"); M["RoneExpnfour"] = rng([r1[s_]["P"] * 1e-4 for s_ in r1], lambda x: f"{x:.1f}")
r6 = {}
for s_ in S:
    Z = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(rv, "out", f"R6_{s_}_k*.json")))]
    if Z: r6[s_] = Z
if r6:
    M["RsixPBcalls"] = rng([np.mean([z["PB"]["V"] for z in r6[s_]]) for s_ in r6]); M["RsixPBTI"] = rng([np.mean([z["PB"]["t05"] for z in r6[s_]]) for s_ in r6], f2)
    M["RsixXDScalls"] = rng([z["XDS"]["V"] for s_ in r6 for z in r6[s_]]); M["RsixXDSTI"] = rng([z["XDS"]["t05"] for s_ in r6 for z in r6[s_]], f3)
    M["RsixXDSTop"] = rng([z["XDS"]["top1"] for s_ in r6 for z in r6[s_]], f3); M["RsixK"] = rng([len(r6[s_]) for s_ in r6])
    M["RsixRw"] = rng([np.mean([z["median_abs_rw"] for z in r6[s_]]) for s_ in r6], f3)
pool = J("rev/pooling.json")
M["PoolRwReal"] = rng([pool[s_]["median_abs_real"] for s_ in S], f3); M["PoolRwPseudo"] = rng([pool[s_]["median_abs_pseudo"] for s_ in S], f3)
M["PoolRatio"] = rng([100 * (pool[s_]["ratio_median"] - 1) for s_ in S], lambda x: f"{x:.0f}") + "\\%"
r35 = J("rev/rev_dc_r5.json")
M["FalseRmed"] = rng([r35[s_]["false_r_median"] for s_ in S], f2); M["FalseRgtOne"] = rng([100 * r35[s_]["false_r_gt01"] for s_ in S], lambda x: f"{x:.0f}") + "\\%"
M["FalseRgtTwo"] = rng([100 * r35[s_]["false_r_gt02"] for s_ in S], lambda x: f"{x:.0f}") + "\\%"
dcg = {}
for s_ in S:
    nul = [np.load(f) for f in sorted(glob.glob(os.path.join(rv, "out", f"DCg_{s_}_k*.npz")))]
    real = np.load(os.path.join(rv, "out", f"DCg_{s_}_real.npz"))
    dcg[s_] = dict(null_calls=[int(bh(z["XDSg"].astype(float)).sum()) for z in nul], null_typeI=[float((z["XDSg"] < 0.05).mean()) for z in nul],
                   real_calls=int(bh(real["XDSg"].astype(float)).sum()))
M["DCgNull"] = rng([v for s_ in S for v in dcg[s_]["null_calls"]]); M["DCgTI"] = rng([v for s_ in S for v in dcg[s_]["null_typeI"]], f3)
M["DCgReal"] = ", ".join(c(dcg[s_]["real_calls"]) for s_ in S)
M["DCovP"] = ", ".join(f"{r35[s_]['overlap']:,} versus {r35[s_]['expected']:,.0f}".replace(",", "{,}") for s_ in S)
M["DCgHubs"] = ", ".join("\\emph{" + g + "}" for g, _ in r35["P1"]["XDSg_top"][:6])
sx = {s_: json.load(open(os.path.join(rv, "out", f"real_{s_}.json"))) for s_ in S}
xy = [r for s_ in S for r in sx[s_]["sexpairs"] if "XIST" in r["pair"]]; yy = [r for s_ in S for r in sx[s_]["sexpairs"] if "XIST" not in r["pair"]]
M["PosXYn"] = c(len(xy)); M["PosYYn"] = c(len(yy)); M["PosXYneg"] = c(sum(r["sign_XDS"] < 0 for r in xy)); M["PosYYpos"] = c(sum(r["sign_XDS"] > 0 for r in yy))
def sci(x):
    m_, e_ = f"{x:.0e}".split("e"); return f"{m_}\\times10^{{{int(e_)}}}"
M["PosXYmaxp"] = sci(max(r['p_XDS'] for r in xy)); M["PosYYmaxp"] = sci(max(r['p_XDS'] for r in yy))
rs = J("rev/rev_sim_results.json")
M["XDSMLFDR"] = rng([rs[d]["XDSM"]["FDR"] for d in ("L200", "L500", "L1000", "LS200", "LS500", "LS1000")], f3)
M["XDSMLpow"] = f3(rs["L200"]["XDSM"]["power"]); M["XDSLpowTwo"] = f3(E3["L200"]["XDS"]["power"])
M["XDSMWN"] = rng([rs[d]["XDSM"]["WN"] for d in ("L200", "L500", "L1000")], f3); M["XDSWN"] = rng([E3[d]["XDS"]["WN"] for d in ("L200", "L500", "L1000")], f3)
M["XDSMsmallFDR"] = rng([rs[d]["XDSM"]["FDR"] for d in ("F_D16", "F_D24", "F_D16_n100")], f3)
M["KXLpow"] = rng([rs[d]["KX"]["power"] for d in ("L200", "L500", "L1000")], f3); M["KXsmallFDR"] = rng([rs[d]["KX"]["FDR"] for d in ("F_D16", "F_D24", "F_D40", "F_D16_n100")], f2)



on = {}
for t_ in ("B1", "B2", "B3"):
    fi = os.path.join(R, "onek", f"onek_{t_}_info.json")
    if not os.path.exists(fi): continue
    inf = json.load(open(fi)); res = {"info": inf}
    for D in ("all", "200"):
        Z = [np.load(f) for f in sorted(glob.glob(os.path.join(R, "onek", "out", f"ON_{t_}_k[0-9]_D{D}.npz")))]
        if not Z: continue
        res[D] = {}
        for m in ("PB", "LCP", "LCPXd", "XDSM"):
            if m not in Z[0]: continue
            calls = [int(bh(z[m].astype(float)).sum()) for z in Z]; t1 = float(np.mean([(z[m] < 0.05).mean() for z in Z]))
            top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
            tt = np.concatenate([z[m][tp] < 0.05 for z, tp in zip(Z, top)]); se = float(np.sqrt(0.05 * 0.95 / len(tt)))
            res[D][m] = dict(K=len(Z), mean_calls=float(np.mean(calls)), calls=calls, typeI=t1, top1=float(tt.mean()), holds=bool(np.mean(calls) <= 1 and t1 <= 0.055 and tt.mean() <= 0.05 + 3 * se))
    on[t_] = res
json.dump(on, open(os.path.join(H, "onek_summary.json"), "w"), indent=1)
ok_ = [t_ for t_ in on if "all" in on[t_] and "XDSM" in on[t_]["all"]]
if ok_:
    nm = {"B1": "transitional B cells", "B2": "naive B cells", "B3": "naive CD4 T cells"}
    parts = [f"{nm[t_]} ({c(on[t_]['info']['donors'])} donors, median {c(on[t_]['info']['median_n'])} cells)" for t_ in ok_]
    M["OnekText"] = ("In OneK1K blood \\citep{Yazar2022}, queried from CELLxGENE Census release 2025-11-08, we repeated the null in " + ", ".join(parts) +
        " (frozen plan ONEK\\_PREREG, three reassignments each). At full cohort size pseudobulk made " +
        rng([on[t_]["all"]["PB"]["mean_calls"] for t_ in ok_]) + " false calls per reassignment (type-I " + rng([on[t_]["all"]["PB"]["typeI"] for t_ in ok_], f2) +
        "), LCP " + rng([on[t_]["all"]["LCP"]["mean_calls"] for t_ in ok_]) + ", XDS " + rng([on[t_]["all"]["LCPXd"]["mean_calls"] for t_ in ok_], lambda x: f"{x:.1f}") +
        " and XDS-M " + rng([on[t_]["all"]["XDSM"]["mean_calls"] for t_ in ok_], lambda x: f"{x:.1f}") + " (type-I " + rng([on[t_]["all"]["XDSM"]["typeI"] for t_ in ok_], f3) +
        "). Pseudobulk's false calls at full size were " + ", ".join(f"{on[t_]['all']['PB']['mean_calls'] / on[t_]['200']['PB']['mean_calls']:.1f}" for t_ in ok_)
        + " times those at 200 donors, so the frozen growth criterion (at least 2) was met in " + str(sum(on[t_]['all']['PB']['mean_calls'] >= 2 * on[t_]['200']['PB']['mean_calls'] and on[t_]['all']['PB']['mean_calls'] >= 100 for t_ in ok_))
        + " of " + str(len(ok_)) + " cell types; XDS met its criterion in " + str(sum(on[t_]['all']['LCPXd']['holds'] and on[t_]['200']['LCPXd']['holds'] for t_ in ok_)) + " and XDS-M in "
        + str(sum(on[t_]['all']['XDSM']['holds'] and on[t_]['200']['XDSM']['holds'] for t_ in ok_))
        + "; XDS-M made no false call but exceeded the top-1\\% bound at 200 transitional B cell donors (" + f3(on['B1']['200']['XDSM']['top1']) + "), a reminder that its $t$ reference needs many donors. Even with a median of 244 cells per donor, pseudobulk fails at 981 donors; the leak is small there, but LCP's residual false calls ("
        + c(on['B3']['all']['LCP']['mean_calls']) + ") and pseudobulk's far larger count point to the log-transform bias of Proposition~\\ref{prop:transform}, which does not vanish as fast as the leak.")
    M["OnekAbs"] = ("In three OneK1K blood cell types, pseudobulk made " + rng([on[t_]["all"]["PB"]["mean_calls"] for t_ in ok_]) + " false calls per reassignment, even with a median of 244 cells per donor, and XDS-M none.")
else:
    M["OnekText"] = ""; M["OnekAbs"] = ""

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
    f"so in the first two subtypes pseudobulk differential co-expression cannot be distinguished from the cell-number artifact. XDS with pooled centring reported {c(dc['P1']['XDS']['real_calls'])}, {c(dc['P2']['XDS']['real_calls'])} and {c(dc['P3']['XDS']['real_calls'])}. "
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
dsim = dres["sim"]
def get(d, m):
    if m == "XDSM": return rs[d]["XDSM"]
    if m == "LCPd": return (E3[d]["LCPd"] if d in E3 else dsim[d]["LCPd"])
    return (E3[d] if d in E3 else E4[d])[m]
MS_ = ("PB", "LCP", "LCPd", "XDS", "XDSM")
for lab, d in (("16 donors, $\\sim$40 cells", "F_D16"), ("24 donors, $\\sim$40 cells", "F_D24"), ("40 donors, $\\sim$40 cells", "F_D40"),
               ("16 donors, $\\sim$100 cells", "F_D16_n100"), ("79 donors, median 15 cells", "F_S1"), ("83 donors, median 43 cells", "F_S2"),
               ("200 donors, median 20", "L200"), ("500 donors, median 20", "L500"), ("1,000 donors, median 20", "L1000")):
    rows.append(" & ".join([lab] + [f"{get(d, m)['FDR']:.3f}" for m in MS_] + [f"{get(d, m)['power']:.2f}" for m in MS_]) + " \\\\")
mcse = max(get(d, m)["FDR_mcse"] for d in ("F_D16", "F_D24", "F_D40", "F_D16_n100", "F_S1", "F_S2", "L200", "L500", "L1000") for m in MS_)
open(os.path.join(H, "table_sim.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3pt}\\caption{Simulations with exact truth: FDR at BH $q=0.05$ and power. $R=100$ replicates for designs with up to 83 donors, $R=50$ for larger designs; the largest FDR MCSE is "
    + f"{mcse:.3f}" + " (all MCSEs in the repository). XDS-M uses a normal reference and is intended for large cohorts.}\\label{tab:sim}\n"
    "\\begin{tabular}{lcccccccccc}\\toprule & \\multicolumn{5}{c}{FDR} & \\multicolumn{5}{c}{Power}\\\\\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}\n"
    "Design & PB & LCP & LCPd & XDS & XDS-M & PB & LCP & LCPd & XDS & XDS-M\\\\\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\\end{tabular}\\end{table}\n")
rows = []
for s in S:
    a = full[s]; rp = pop["rep_null"][s]
    rows.append(f"{SUBN[s]} ({c(info[s]['donors'])}; {c(info[s]['median_nd'])}) & " + " & ".join(
        f"{c(a[m]['mean_calls'])} & {a[m]['typeI']:.3f} & {a[m]['typeI_top1']:.3f}" for m in ("PB", "LCP", "LCPXd")) +
        f" & {c(np.mean(rp['PB']['replicated_false']))} & {r1[s]['XDS']['runs_any']}/{r1[s]['XDS']['n']} & {r1[s]['XDSM']['runs_any']}/{r1[s]['XDSM']['n']} \\\\")
open(os.path.join(H, "table_null.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3.5pt}\\caption{PsychAD real-cell null at full cohort size (pseudo-donors; median nuclei), mean over three reassignments. "
    "Calls: BH false discoveries among $\\sim$89{,}000 pairs. TI: per-pair type-I error at 0.05, over all pairs and over the 1\\% most within-correlated pairs. Rep.: PB false calls replicated (same sign, $p<0.05$) in an independent half of the pseudo-donors. Runs: reassignments with at least one false call (equal to the FDR under this complete null) among further reassignments, for XDS and XDS-M.}\\label{tab:null}\n"
    "\\resizebox{\\linewidth}{!}{\\begin{tabular}{lcccccccccccc}\\toprule & \\multicolumn{3}{c}{PB} & \\multicolumn{3}{c}{LCP} & \\multicolumn{3}{c}{XDS} & Rep. & \\multicolumn{2}{c}{Runs $V\\ge1$}\\\\"
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){12-13}\n Subtype & Calls & TI & top 1\\% & Calls & TI & top 1\\% & Calls & TI & top 1\\% & PB & XDS & XDS-M\\\\\\midrule\n"
    + "\n".join(rows) + "\n\\bottomrule\\end{tabular}}\\end{table}\n")
# ablation table: which bias each method leaves, with evidence
ab = [("PB", "log-CPM, own cells", "--", "leak, transform", "PB"),
      ("PBlin", "ratio mean, own cells", "--", "leak", "PBlin"),
      ("LCP", "self-standardized", "subtract estimate", "transform, leak-estimate error", "LCP"),
      ("LCPd", "design-only", "subtract estimate", "leak-estimate error", "LCPd"),
      ("LCPX", "self-standardized", "disjoint halves", "transform (larger)", "LCPX"),
      ("XU", "design-only weights", "exact U-statistic", "depth--state coupling", "XU"),
      ("XDS", "design-only", "disjoint halves", "tail of normal reference", "LCPXd"),
      ("XDS-M", "design-only", "disjoint halves, 5 splits", "none detected (needs large $D$)", "XDSM"),
      ("KX", "design-only", "regress on $\\kappa_d$", "one subtype fails; low power", "KX")]
rows = []
for nm, stdz, leak, left, key in ab:
    if key in ("XDSM", "KX"):
        calls = rng([r1[s][key]["V"] for s in S], lambda x: f"{x:.2f}"); ti = rng([r1[s][key]["t05"] for s in S], f3)
    elif key == "XU":
        calls = rng([xu["E1"][s]["all"]["XU"]["mean_calls"] for s in S]); ti = rng([xu["E1"][s]["all"]["XU"]["typeI"] for s in S], f3)
    else:
        calls = rng([full[s][key]["mean_calls"] for s in S]); ti = rng([full[s][key]["typeI"] for s in S], f3)
    if key in ("PB", "LCP"): sea = rng([fx["null"][s][key]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key == "LCPX": sea = rng([xs["null"][s]["LCPX"]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key in ("LCPd", "LCPXd", "PBlin"): sea = rng([dres["null"][s][key]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key == "XU": sea = rng([xu["E2"][s]["typeI_top1"] for s in ("S1", "S2")], f3)
    elif key in ("XDSM", "KX"):
        vals = []
        for st_ in ("S1", "S2"):
            Zs = [np.load(f) for f in sorted(glob.glob(os.path.join(R, "rev", "out", f"seanull_{st_}_k*.npz")))]
            if Zs:
                top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Zs]
                vals.append(float(np.concatenate([z["p_" + key][tp] < 0.05 for z, tp in zip(Zs, top)]).mean()))
        sea = rng(vals, f3) if vals else "--"
    rows.append(f"{nm} & {stdz} & {leak} & {left} & {calls} & {ti} & {sea} \\\\")
open(os.path.join(H, "table_ablate.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3pt}\\caption{Each correction and the bias it leaves. Calls and TI: PsychAD real-cell null at full size, mean calls per reassignment and per-pair type-I (range over subtypes; XU, XDS-M and KX on further reassignments). "
    "SEA-AD: type-I for the 1\\% most within-correlated pairs in the SEA-AD real-cell null (79 and 83 pseudo-donors). LCP2 (not shown) had SEA-AD top-1\\% type-I "
    + rng([fx["null"][s]["LCP2"]["typeI_top1"] for s in ("S1", "S2")], f3) + ".}\\label{tab:ablate}\n"
    "\\resizebox{\\linewidth}{!}{\\begin{tabular}{lllllll}\\toprule Method & Standardization & Leak handling & Bias left & Calls & TI & SEA-AD\\\\\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\\end{tabular}}\\end{table}\n")
print(len(M), "macros"); print({k: v for k, v in M.items() if k not in ("DCtext", "SplitText")}); print(M["DCtext"][-700:])
