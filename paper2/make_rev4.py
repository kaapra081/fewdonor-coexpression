"""REV4 summaries -> rev4/results/REV4_SUMMARY.json, paper2/rev4nums.tex, table_pbfix.tex, table_calib.tex, table_R1.tex, table_R2.tex.
Run from repo/: python paper2/make_rev4.py   (reads rev4/out/*.json, rev2/out/D1_*.json, rev2/interim/REV2_SUMMARY_latest.json, rev3/results/A3B3_SUMMARY.json)"""
import os, glob, json, math, numpy as np
RR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = lambda *a: os.path.join(RR, *a)
SUBS = ("P1", "P2", "P3"); SN = {"P1": "Homeo\\_CECR2", "P2": "PVM\\_CD163", "P3": "ADAM\\_GPNMB"}


def bound(n): return 0.05 + 2 * math.sqrt(0.05 * 0.95 / n)
def holds(k, n): return k / n <= bound(n) + 1e-12
def fmt(x): return f"{int(round(x)):,}".replace(",", "{,}")
def rng(v): v = [int(round(x)) for x in v]; return fmt(min(v)) if min(v) == max(v) else f"{fmt(min(v))}--{fmt(max(v))}"
def load(pat): return {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(P(pat)))}


S = {}; M = {}
# ---------- R1: small D, B = 2e6 ----------
R1 = load("rev4/out/R1_P*_k*.json"); S["R1"] = {}
for s in SUBS:
    runs = [v for k, v in R1.items() if k.startswith(f"R1_{s}_")]; S["R1"][s] = {"n": len(runs)}
    for D in ("25", "50", "100"):
        S["R1"][s][D] = {m: dict(any=sum(r[D][m]["V"] > 0 for r in runs), n=len(runs), meanV=float(np.mean([r[D][m]["V"] for r in runs])) if runs else None,
                             floor_mean=float(np.mean([r[D][m]["nfloor"] for r in runs])) if runs else None) for m in ("PB", "LCP", "LCPd", "XDSperm", "XDS", "XDSM", "KX")}
# ---------- R2: bank-stratified null ----------
R2 = load("rev4/out/R2_P*_k*.json"); S["R2"] = {}
for s in SUBS:
    runs = [v for k, v in R2.items() if k.startswith(f"R2_{s}_")]
    S["R2"][s] = {m: dict(any=sum(r["all"][m]["V"] > 0 for r in runs), n=len(runs), meanV=float(np.mean([r["all"][m]["V"] for r in runs])) if runs else None,
                          t05=float(np.mean([r["all"][m]["t05"] for r in runs])) if runs else None) for m in ("PB", "LCP", "LCPd", "XDS", "XDSM", "KX")}
# ---------- R3 / R3o: pseudobulk family ----------
R3 = load("rev4/out/R3_P*_k*.json"); R3o = load("rev4/out/R3o_*_k*.json"); S["R3"] = {}; S["R3o"] = {}
for s in SUBS:
    runs = [v for k, v in R3.items() if k.startswith(f"R3_{s}_")]
    S["R3"][s] = {m: dict(any=sum(r[m]["V"] > 0 for r in runs), n=len(runs), meanV=float(np.mean([r[m]["V"] for r in runs])), t05=float(np.mean([r[m]["t05"] for r in runs]))) for m in ("PB", "PBres", "PBw", "PC5", "PC10")}
for t in ("B1", "B2", "B3"):
    runs = [v for k, v in R3o.items() if k.startswith(f"R3o_{t}_")]
    S["R3o"][t] = {m: dict(any=sum(r[m]["V"] > 0 for r in runs), n=len(runs), meanV=float(np.mean([r[m]["V"] for r in runs])), t05=float(np.mean([r[m]["t05"] for r in runs]))) for m in ("PB", "PBw", "PC5", "PC10")}
# ---------- R4 ----------
S["R4"] = {s: json.load(open(P(f"rev4/out/R4_{s}.json"))) for s in SUBS}
# ---------- R5: DC null, 20 runs ----------
D1old = json.load(open(P("rev2/interim/REV2_SUMMARY_latest.json")))["D1"]; S["R5"] = {}
for s in SUBS:
    runs = {f"k{k}": D1old[s][f"k{k}"] for k in (10, 11, 12)}
    for f in glob.glob(P(f"rev2/out/D1_{s}_k*.json")): runs[os.path.basename(f).split("_")[-1][:-5]] = json.load(open(f))
    ks = sorted(runs, key=lambda x: int(x[1:])); assert len(ks) == 20 and ks[0] == "k10" and ks[-1] == "k29", ks
    S["R5"][s] = {m: dict(any=sum(runs[k][m]["V"] > 0 for k in ks), n=20, meanV=float(np.mean([runs[k][m]["V"] for k in ks])), t05=float(np.mean([runs[k][m]["t05"] for k in ks]))) for m in ("PB", "XDSg", "XDSMg", "XDSMgu")}
# ---------- R6 ----------
r6 = P("rev4/out/R6/R6_summary.json")
if os.path.exists(r6): S["R6"] = json.load(open(r6))

# ---------- macros ----------
def mac(name, val): M[name] = val
r3 = S["R3"]; r3o = S["R3o"]
mac("RfPBmean", rng([r3[s]["PB"]["meanV"] for s in SUBS])); mac("RfPBres", rng([r3[s]["PBres"]["meanV"] for s in SUBS]))
mac("RfPBwOne", fmt(r3["P1"]["PBw"]["meanV"])); mac("RfPBwTwo", fmt(r3["P2"]["PBw"]["meanV"])); mac("RfPBwThree", fmt(r3["P3"]["PBw"]["meanV"]))
mac("RfPCadj", rng([r3[s][m]["meanV"] for s in SUBS for m in ("PC5", "PC10")]))
mac("RfPCanyAll", "every" if all(r3[s][m]["any"] == r3[s][m]["n"] for s in SUBS for m in ("PC5", "PC10")) else "not every")
mac("RfOnekPB", rng([r3o[t]["PB"]["meanV"] for t in r3o])); mac("RfOnekPBw", rng([r3o[t]["PBw"]["meanV"] for t in r3o])); mac("RfOnekPC", rng([r3o[t][m]["meanV"] for t in r3o for m in ("PC5", "PC10")]))
mac("RfOnekPBwTI", rng([1000 * r3o[t]["PBw"]["t05"] for t in r3o]).replace("{,}", "")); M["RfOnekPBwTI"] = "--".join(f"{x:.3f}" for x in (min(r3o[t]["PBw"]["t05"] for t in r3o), max(r3o[t]["PBw"]["t05"] for t in r3o)))
r5 = S["R5"]; mac("RfDCPB", rng([r5[s]["PB"]["meanV"] for s in SUBS])); mac("RfDCgu", ", ".join(str(r5[s]["XDSMgu"]["any"]) for s in SUBS)); mac("RfDCg", ", ".join(str(r5[s]["XDSMg"]["any"]) for s in SUBS))
mac("RfDCguHolds", str(sum(holds(r5[s]["XDSMgu"]["any"], 20) for s in SUBS))); mac("RfDCgHolds", str(sum(holds(r5[s]["XDSMg"]["any"], 20) for s in SUBS)))
r4 = S["R4"]
mac("RfRealPB", ", ".join(fmt(r4[s]["calls"]["PB"]) for s in SUBS)); mac("RfRealXDSM", ", ".join(fmt(r4[s]["calls"]["XDSM"]) for s in SUBS))
mac("RfRealPBonly", ", ".join(fmt(r4[s]["classes"]["PBonly"]["n"]) for s in SUBS))
mac("RfRealPBonlyPct", rng([100 * r4[s]["classes"]["PBonly"]["n"] / r4[s]["calls"]["PB"] for s in SUBS]) + "\\%")
mac("RfRealXonly", ", ".join(fmt(r4[s]["classes"]["XDSMonly"]["n"]) for s in SUBS)); mac("RfRealShared", ", ".join(fmt(r4[s]["classes"]["shared"]["n"]) for s in SUBS))
mac("RfHubOverlap", rng([r4[s]["hub_overlap_top20"] for s in SUBS])); mac("RfTopOverlap", rng([100 * r4[s]["top5000_pair_overlap"] / 5000 for s in SUBS]) + "\\%")
if "R6" in S:
    r6 = S["R6"]; f2 = lambda x: f"{x:.2f}".replace("-", "$-$")
    mac("RsixNullCalls", fmt(r6["pairs"]["null_calls"])); mac("RsixPairs", fmt(r6["pairs"]["n"]))
    mac("RsixNormPB", f2(r6["null_L"]["spearman_norm_ncells"])); mac("RsixNormR", f2(r6["null_R"]["spearman_norm_ncells"]))
    mac("RsixPConeVar", f"{100 * r6['null_L']['varfrac_PC1to5'][0]:.1f}\\%"); mac("RsixPConeR", f"{r6['null_L']['r2_logn_PC1to5'][0]:.2f}")
    mac("RsixPConeRR", f"{max(r6['null_R']['r2_logn_PC1to5']):.3f}"); mac("RsixRealPCtwoR", f"{r6['real_L']['r2_logn_PC1to5'][1]:.2f}")
    mac("RsixRealNormPB", f2(r6["real_L"]["spearman_norm_ncells"])); mac("RsixRealNormR", f2(r6["real_R"]["spearman_norm_ncells"]))
    pc1 = json.load(open(P("rev4/out/R6/R6_pc1.json")))
    mac("RsixPConeSign", f"{100 * pc1['null_L']['frac_same_sign_PC1']:.0f}\\%"); mac("RsixRowMeanLogn", f"{pc1['null_L']['corr_rowmean_logn']:.2f}"); mac("RsixRowMeanLognR", f2(pc1['null_R']['corr_rowmean_logn']))
    import pandas as pd; from sklearn.cluster import KMeans
    dr = pd.read_csv(P("rev4/out/R6/donors_real.csv")); km = KMeans(2, n_init=10, random_state=0).fit_predict(dr[["x_L", "y_L"]].values)
    ok = dr.sex.isin(["male", "female"]).values; tab = pd.crosstab(km[ok], dr.sex.values[ok]).values; conc = max(tab[0, 0] + tab[1, 1], tab[0, 1] + tab[1, 0])
    mac("RsixSexConc", f"{fmt(conc)} of {fmt(ok.sum())}"); S["R6"]["sex_cluster_concordance"] = [int(conc), int(ok.sum())]
    mac("RsixPairNull", f2(r6["pairs"]["spearman_nullPB_rw"])); mac("RsixPairRealPB", f2(r6["pairs"]["spearman_realPB_rw"])); mac("RsixPairRealX", f2(r6["pairs"]["spearman_realXDSMz_rw"]))
# R1 / R2 pooled + holds
def pooled(block, D, m):
    sv = [s for s in SUBS if block[s].get("n", 0)]
    k = sum(block[s][D][m]["any"] for s in sv); n = sum(block[s][D][m]["n"] for s in sv); h = sum(holds(block[s][D][m]["any"], block[s][D][m]["n"]) for s in sv); return k, n, h
mac("RfRoneN", rng([S["R1"][s]["n"] for s in SUBS]))
from scipy.stats import beta as _beta
_up = lambda k, n: 100 * _beta.ppf(0.975, k + 1, n - k)
_k400 = sum(RA_[x]["400"]["XDSM"]["any"] for x in SUBS) if (RA_ := json.load(open(P("rev2/interim/REV2_SUMMARY_latest.json")))["A"]) else 0
mac("RfUpperPooled", f"{_up(_k400, 300):.0f}\\%"); mac("RfUpperThree", f"{_up(RA_['P3']['400']['XDSM']['any'], 100):.0f}\\%")
_secs = [v["sec"] for k, v in R1.items() if "sec" in v]; mac("RfRoneMin", f"{np.median(_secs) / 60:.0f}" if _secs else "?")
mc = json.load(open(P("rev4/out/R6/R6_mantis_clusters.json"))); S.setdefault("R6", {})["mantis_clusters"] = mc
fp = lambda x: f"{x:.2g}" if x >= 0.001 else f"{x:.1g}"
mac("McLeaves", str(mc["n_leaves"])); mac("McTotal", str(mc["n_clusters_total"])); mac("McLevels", str(mc["n_levels"])); mac("McTopLabel", mc["top_PBonly_cluster"]["label"])
mac("McTopN", str(mc["top_PBonly_cluster"]["n"])); mac("McTopPBonly", f"{mc['top_PBonly_cluster']['median']:.0f}"); mac("McTopNull", f"{mc['top_null_cluster']['median']:.0f}"); mac("McTopXDSM", f"{mc['top_PBonly_cluster']['xdsm_median']:.0f}")
mac("McMedPBonly", f"{mc['gene_medians']['PB_only_degree']:.0f}"); mac("McMedNull", f"{mc['gene_medians']['null_false_degree']:.0f}"); mac("McMedXDSM", f"{mc['gene_medians']['XDSM_degree']:.0f}")
mac("McKWPBonly", fp(mc["kw"]["PB_only_degree"]["p"])); mac("McKWXDSM", fp(mc["kw"]["XDSM_degree"]["p"])); mac("McKWNull", fp(mc["kw"]["null_false_degree"]["p"]))
mac("McSameTop", "yes" if mc["top_PBonly_cluster"]["label"] == mc["top_null_cluster"]["label"] else "no")
hs = json.load(open(P("rev4/out/R6/R6_hs.json"))); S.setdefault("R6", {})["heatshock_posthoc"] = hs
mac("RsixHSn", str(hs["n_hs_genes"])); mac("RsixHSpairs", str(hs["hs_pairs"])); mac("RsixHScalled", str(hs["hs_pairs_called"])); mac("RsixHSpct", f"{100 * hs['hs_pairs_called'] / hs['hs_pairs']:.0f}\\%")
mac("RsixAllpct", f"{100 * hs['frac_all_called']:.1f}\\%"); mac("RsixHSrw", f"{hs['median_rw_hs_pairs']:.3f}"); mac("RsixAllrw", f"{hs['median_rw_all']:.3f}"); mac("RsixYcalled", f"{hs['y_pairs_called']} of {hs['y_pairs']}")
sea = json.load(open(P("paper2/seanull_summary.json"))); mac("SeaXDSMTop", f"{min(sea['XDSM']):.3f}--{max(sea['XDSM']):.3f}"); mac("SeaKXTop", f"{min(sea['KX']):.3f}--{max(sea['KX']):.3f}")
mac("RfKXpoolThree", str(json.load(open(P("paper2/r1_summary.json")))["P3"]["KX"]["runs_any"]))
if any(S["R1"][s]["n"] for s in SUBS):
    for D, w in (("25", "TwentyFive"), ("50", "Fifty"), ("100", "Hundred")):
        for m, mn in (("LCP", "LCP"), ("LCPd", "LCPd"), ("XDSperm", "XDSperm"), ("XDSM", "XDSM"), ("KX", "KX")):
            k, n, h = pooled(S["R1"], D, m); mac(f"Rf{mn}{w}", f"{k}/{n}"); mac(f"Rf{mn}{w}Holds", str(h))
if all(S["R2"][s]["PB"]["n"] for s in SUBS):
    for m in ("PB", "LCP", "LCPd", "XDS", "XDSM", "KX"):
        mac(f"RfBank{m}", ", ".join(f"{S['R2'][s][m]['any']}/{S['R2'][s][m]['n']}" for s in SUBS)); mac(f"RfBank{m}Holds", str(sum(holds(S['R2'][s][m]['any'], S['R2'][s][m]['n']) for s in SUBS)))
        mac(f"RfBank{m}Mean", rng([S["R2"][s][m]["meanV"] for s in SUBS]))
with open(P("paper2/rev4nums.tex"), "w") as fh:
    fh.write("% generated by make_rev4.py from rev4/results/REV4_SUMMARY.json\n")
    for k, v in M.items(): fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

# ---------- tables ----------
def cell(k, n): return f"{k}/{n}" + ("" if holds(k, n) else "$^\\dagger$")
# Table 1: pseudobulk and its corrections
L = ["\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{4pt}",
     "\\caption{Pseudobulk and its corrections on real-cell nulls at full cohort size: mean BH false calls per reassignment (per-pair type-I error at 0.05 in parentheses). Every reassignment of every variant made at least one false call. PsychAD: donor-centred null, 20 reassignments per subtype ($\\sim$89{,}000 pairs). OneK1K: pooled null, 10 reassignments per cell type (B1 transitional, B2 naive B, B3 naive CD4$^+$ T; $\\sim$88{,}000 pairs). PBres: log-CPM adjusted for log cell number and log depth; PBw: donors weighted by cell number; PC$k$: top $k$ principal components removed.}\\label{tab:pbfix}",
     "\\begin{tabular}{lcccccc}\\toprule & \\multicolumn{3}{c}{PsychAD (donor-centred)} & \\multicolumn{3}{c}{OneK1K (pooled)}\\\\\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
     "Method & Homeo & PVM & ADAM & B1 & B2 & B3\\\\\\midrule"]
for m, lab in (("PB", "Pseudobulk"), ("PBres", "PBres"), ("PBw", "PBw"), ("PC5", "PC5"), ("PC10", "PC10")):
    row = [lab]
    for s in SUBS: row.append(f"{fmt(r3[s][m]['meanV'])} ({r3[s][m]['t05']:.2f})")
    for t in ("B1", "B2", "B3"): row.append(f"{fmt(r3o[t][m]['meanV'])} ({r3o[t][m]['t05']:.2f})" if m in r3o[t] else "--")
    L.append(" & ".join(row) + "\\\\")
L += ["\\bottomrule\\end{tabular}\\end{table}"]
open(P("paper2/table_pbfix.tex"), "w").write("\n".join(L) + "\n")
# Appendix: R1 per subtype, R2 per subtype
if all(S["R1"][s]["n"] for s in SUBS):
    L = ["\\begin{table}[h]\\centering\\small\\setlength{\\tabcolsep}{3pt}",
         f"\\caption{{Small cohorts with $2\\times10^6$ donor permutations (REV4 R1): reassignments with at least one BH call, of {S['R1']['P1']['n']}, {S['R1']['P2']['n']} and {S['R1']['P3']['n']} per subtype, donor-centred null (clip $[0.1,10]$). $^\\dagger$: outside the frozen bound. LCP, LCPd and XDS (perm.) use donor-permutation $p$-values (smallest attainable $p=5\\times10^{{-7}}$, below the BH first-call threshold); XDS (normal), XDS-M and KX use normal references.}}\\label{tab:R1}",
         "\\begin{tabular}{lccccccccc}\\toprule & \\multicolumn{3}{c}{Homeo\\_CECR2} & \\multicolumn{3}{c}{PVM\\_CD163} & \\multicolumn{3}{c}{ADAM\\_GPNMB}\\\\",
         "Method & 25 & 50 & 100 & 25 & 50 & 100 & 25 & 50 & 100\\\\\\midrule"]
    for m, lab in (("PB", "PB"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDSperm", "XDS (perm.)"), ("XDS", "XDS (normal)"), ("XDSM", "XDS-M"), ("KX", "KX")):
        L.append(lab + " & " + " & ".join(cell(S["R1"][s][D][m]["any"], S["R1"][s][D][m]["n"]) for s in SUBS for D in ("25", "50", "100")) + "\\\\")
    L += ["\\bottomrule\\end{tabular}\\end{table}"]; open(P("paper2/table_R1.tex"), "w").write("\n".join(L) + "\n")
if all(S["R2"][s]["PB"]["n"] for s in SUBS):
    L = ["\\begin{table}[h]\\centering\\small\\setlength{\\tabcolsep}{4pt}",
         "\\caption{Brain-bank-stratified null (REV4 R2): donor-centred cells reassigned only among donors of the same brain bank, full cohort size. Reassignments with at least one BH call (mean false calls). $^\\dagger$: outside the frozen bound.}\\label{tab:R2}",
         "\\begin{tabular}{lccc}\\toprule Method & Homeo\\_CECR2 & PVM\\_CD163 & ADAM\\_GPNMB\\\\\\midrule"]
    for m, lab in (("PB", "PB"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDS", "XDS (one split)"), ("XDSM", "XDS-M"), ("KX", "KX")):
        L.append(lab + " & " + " & ".join(f"{cell(S['R2'][s][m]['any'], S['R2'][s][m]['n'])} ({fmt(S['R2'][s][m]['meanV']) if S['R2'][s][m]['meanV'] >= 1 else round(S['R2'][s][m]['meanV'], 2)})" for s in SUBS) + "\\\\")
    L += ["\\bottomrule\\end{tabular}\\end{table}"]; open(P("paper2/table_R2.tex"), "w").write("\n".join(L) + "\n")

# ---------- Table 2: run-level calibration, pooled over the three PsychAD subtypes ----------
import re
A3 = json.load(open(P("rev3/results/A3B3_SUMMARY.json")))["A3"]; RA = json.load(open(P("rev2/interim/REV2_SUMMARY_latest.json")))["A"]
r2n = open(P("paper2/rev2nums.tex")).read()
def macro(name): return re.search(r"\\newcommand\{\\" + name + r"\}\{([^}]*)\}", r2n).group(1)
oneK = {"PB": macro("OnekEanyPB"), "LCP": macro("OnekEanyLCP"), "XDSM": macro("OnekEanyXDSM")}
R2c = load("rev4/out/R2c_P*_k*.json"); S["R2c"] = {}
for s_ in SUBS:
    runs = [v for k, v in R2c.items() if k.startswith(f"R2c_{s_}_")]
    S["R2c"][s_] = {m: dict(any=sum(r[m]["V"] > 0 for r in runs), n=len(runs), meanV=float(np.mean([r[m]["V"] for r in runs])) if runs else None) for m in ("PB", "PBb", "XDSM", "XDSMb")}
    if runs: S["R2c"][s_]["bank_shift_median"] = {b: float(np.median([r["bank_shift_median"][b] for r in runs])) for b in runs[0]["bank_shift_median"]}
def pool(cells):   # cells: list of (any, n) per subtype -> "k/n" with dagger if any cell is outside its bound
    k = sum(a for a, n in cells); n = sum(n for a, n in cells)
    if n == 0: return "--"
    return f"{k}/{n}" + ("" if all(holds(a, n_) for a, n_ in cells if n_) else "$^\dagger$")
rows = [("PB", "Pseudobulk"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDS", "XDS, one split"), ("XDSM", "XDS-M"), ("KX", "KX")]
perm = {"LCP", "LCPd", "XDS"}
L = ["\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3.2pt}",
     "\\caption{Run-level error on real-cell nulls: reassignments with at least one BH call, pooled over the three PsychAD subtypes (under a complete null this fraction is the FDR). $^\\dagger$: outside the frozen bound ($0.05+2$\\,SE) in at least one subtype or cell type. Donor-centred null: 25--100 pseudo-donors from REV3 (100 runs per subtype) for pseudobulk, XDS-M and KX, and from REV4 with $2\\times10^6$ donor permutations (" + str(S["R1"]["P1"]["n"]) + "--" + str(max(S["R1"][x]["n"] for x in SUBS)) + " runs per subtype) for LCP, LCPd and XDS, which use permutation $p$-values at $D\\le100$; 200 to full size from REV2 (100 runs per subtype, normal references). Bank: brain-bank-stratified null at full size (REV4, 40 runs per subtype); the bank-centred row is post hoc (Section~\\ref{sec:calib}). OneK1K: pooled null at full size, three cell types (REV2; 30--33 runs each).}\\label{tab:calib}",
     "\\begin{tabular}{lcccccc|c|c}\\toprule & \\multicolumn{6}{c|}{Donor-centred null, pseudo-donors} & Bank & OneK1K\\\\",
     "Method & 25 & 50 & 100 & 200 & 400 & full & full & full\\\\\\midrule"]
for m, lab in rows:
    row = [lab]
    for D in ("25", "50", "100"):
        if m in perm: mm = "XDSperm" if m == "XDS" else m; row.append(pool([(S["R1"][x][D][mm]["any"], S["R1"][x][D][mm]["n"]) for x in SUBS]))
        else: row.append(pool([(A3[x][D][m]["any"], A3[x][D][m]["n"]) for x in SUBS]))
    for D in ("200", "400", "all"): row.append(pool([(RA[x][D][m]["any"], RA[x][D][m]["n"]) for x in SUBS]))
    row.append(pool([(S["R2"][x][m]["any"], S["R2"][x][m]["n"]) for x in SUBS]) if all(S["R2"][x][m]["n"] for x in SUBS) else "--")
    if m in oneK: row.append(pool([tuple(int(v) for v in c.strip().split("/")) for c in oneK[m].split(",")]))
    else: row.append("--")
    L.append(" & ".join(row) + "\\\\")
if all(S["R2c"][x]["XDSMb"]["n"] for x in SUBS):
    L.append("XDS-M, bank-centred & & & & & & & " + pool([(S["R2c"][x]["XDSMb"]["any"], S["R2c"][x]["XDSMb"]["n"]) for x in SUBS]) + " & \\\\")
L += ["\\bottomrule\\end{tabular}\\end{table}"]
open(P("paper2/table_calib.tex"), "w").write("\n".join(L) + "\n")
if all(S["R2c"][x]["XDSMb"]["n"] for x in SUBS):
    mac("RfBankXDSMb", ", ".join(f"{S['R2c'][x]['XDSMb']['any']}/{S['R2c'][x]['XDSMb']['n']}" for x in SUBS)); mac("RfBankXDSMbHolds", str(sum(holds(S['R2c'][x]['XDSMb']['any'], S['R2c'][x]['XDSMb']['n']) for x in SUBS)))
    mac("RfBankXDSMc", ", ".join(f"{S['R2c'][x]['XDSM']['any']}/{S['R2c'][x]['XDSM']['n']}" for x in SUBS)); mac("RfBankPBb", rng([S["R2c"][x]["PBb"]["meanV"] for x in SUBS]))
    sh = [v for x in SUBS for v in S["R2c"][x]["bank_shift_median"].values()]; mac("RfBankShift", f"{min(sh):.2f}--{max(sh):.2f}")
with open(P("paper2/rev4nums.tex"), "w") as fh:
    fh.write("% generated by make_rev4.py from rev4/results/REV4_SUMMARY.json\n")
    for k, v in M.items(): fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
os.makedirs(P("rev4/results"), exist_ok=True); json.dump(S, open(P("rev4/results/REV4_SUMMARY.json"), "w"), indent=1)
print(len(M), "macros;", {k: S["R1"][k]["n"] for k in SUBS}, {k: S["R2"][k]["PB"]["n"] for k in SUBS})
