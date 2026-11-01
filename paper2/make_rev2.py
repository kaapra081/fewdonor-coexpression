"""Run-level table and macros from REV2_SUMMARY.json (REV2_PREREG). python make_rev2.py <REV2_SUMMARY.json>"""
import sys, json, math, os
S = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "rev2_summary.json")))
HERE = os.path.dirname(os.path.abspath(__file__))
# REV2_PREREG E: XDS-M run-level fraction over all 33 OneK1K runs (3 ONEK_PREREG reassignments + 30 REV2 reassignments)
_on = os.path.join(HERE, "onek_summary.json")
if S.get("E") and os.path.exists(_on):
    _o = json.load(open(_on))
    for _t, _v in S["E"].items():
        if "XDSM" in _v and _t in _o and "XDSM" in _o[_t].get("all", {}):
            _c = _o[_t]["all"]["XDSM"]["calls"]; _v["XDSM"] = dict(_v["XDSM"]); _v["XDSM"]["any"] += sum(c >= 1 for c in _c); _v["XDSM"]["n"] += len(_c)
def holds(a): return a["any"] / a["n"] <= 0.05 + 2 * math.sqrt(0.05 * 0.95 / a["n"])
MS = [("PB", "PB"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDS", "XDS"), ("XDSM", "XDS-M"), ("KX", "KX")]
cols = [(s, D) for s in ("P1", "P2", "P3") for D in ("200", "400", "all") if s in S["A"]]
ecols = [t for t in ("B1", "B2", "B3") if t in S.get("E", {})]
def cell(a):
    if a is None: return "--"
    x = f"{a['any']}"; return x if holds(a) else x + "$^\\dagger$"
L = ["\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{2.6pt}",
     "\\caption{Run-level error on real-cell nulls (frozen plan REV2\\_PREREG): number of reassignments with at least one BH call (under a complete null, the FDR equals the fraction of such reassignments). PsychAD, donor-centred null: P1 Homeo\\_CECR2, P2 PVM\\_CD163, P3 ADAM\\_GPNMB, " + ", ".join(sorted({str(S["A"][s]["all"]["PB"]["n"]) for s in S["A"]})) + " reassignments each at nested subsets of 200, 400 and all pseudo-donors. OneK1K, pooled null at full size: B1 transitional B cells, B2 naive B cells, B3 naive CD4 T cells; " + (", ".join(sorted({str(S["E"][t]["PB"]["n"]) for t in ecols})) if ecols else "--") + " reassignments (XDS-M: " + (", ".join(sorted({str(S["E"][t]["XDSM"]["n"]) for t in ecols})) if ecols else "--") + ", including the three ONEK\\_PREREG reassignments, as the plan specifies). $^\\dagger$ fails the frozen rule (fraction $>0.05+2$\\,SE).}\\label{tab:runlevel}",
     "\\begin{tabular}{l" + "c" * len(cols) + ("|" + "c" * len(ecols) if ecols else "") + "}\\toprule"]
hdr1 = " & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{{s}}}" for s in ("P1", "P2", "P3") if s in S["A"]) + ((" & \\multicolumn{%d}{c}{OneK1K}" % len(ecols)) if ecols else "") + "\\\\"
hdr2 = "Method & " + " & ".join("all" if D == "all" else D for _, D in cols) + ((" & " + " & ".join(ecols)) if ecols else "") + "\\\\\\midrule"
L += [hdr1, hdr2]
for m, lab in MS:
    row = [cell(S["A"][s][D].get(m)) for s, D in cols] + [cell(S["E"][t].get(m)) for t in ecols if True]
    L.append(lab + " & " + " & ".join(row) + "\\\\")
L.append("\\bottomrule\\end{tabular}\\end{table}")
open(os.path.join(HERE, "table_runlevel.tex"), "w").write("\n".join(L) + "\n")
# macros
M = {}
def runs(m, ss=("P1", "P2", "P3"), Ds=("200", "400", "all")):
    v = [S["A"][s][D][m] for s in ss if s in S["A"] for D in Ds]; return v
for m, lab in MS:
    v = runs(m); k = sum(a["any"] for a in v); n = sum(a["n"] for a in v); nh = sum(holds(a) for a in v)
    key = lab.replace("-", "")
    M[f"RtwoAny{key}"] = f"{k}/{n}"; M[f"RtwoHolds{key}"] = f"{nh} of {len(v)}"
    M[f"RtwoRange{key}"] = f"{min(a['any'] for a in v)}--{max(a['any'] for a in v)}"
    M[f"RtwoMeanV{key}"] = f"{min(a['meanV'] for a in v):,.1f}--{max(a['meanV'] for a in v):,.1f}" if max(a['meanV'] for a in v) < 100 else f"{min(a['meanV'] for a in v):,.0f}--{max(a['meanV'] for a in v):,.0f}"
from scipy.stats import beta as _beta
def _cp(k, n):
    lo = 0.0 if k == 0 else _beta.ppf(0.025, k, n - k + 1); hi = 1.0 if k == n else _beta.ppf(0.975, k + 1, n - k); return lo, hi
for m, lab in MS:
    v = runs(m); k = sum(a["any"] for a in v); n = sum(a["n"] for a in v); lo, hi = _cp(k, n); key = lab.replace("-", "")
    M[f"RtwoPooled{key}"] = f"{100 * k / n:.1f}\\% (95\\% CI {100 * lo:.1f}--{100 * hi:.1f}\\%)"
    M[f"RtwoCells{key}"] = str(len(v)); M[f"RtwoRuns{key}"] = str(n)
    va = runs(m, Ds=("all",)); ka = sum(a["any"] for a in va); na = sum(a["n"] for a in va); lo, hi = _cp(ka, na)
    M[f"RtwoPooledAll{key}"] = f"{ka}/{na} ({100 * ka / na:.1f}\\%, 95\\% CI {100 * lo:.1f}--{100 * hi:.1f}\\%)"
    M[f"RtwoPctAll{key}"] = f"{100 * ka / na:.0f}\\%"
for t in ecols:
    for m, lab in MS:
        if m in S["E"][t]: M[f"Onek{t.replace('1','One').replace('2','Two').replace('3','Three')}Any{lab.replace('-', '')}"] = f"{S['E'][t][m]['any']}/{S['E'][t][m]['n']}"

# REV3 B2 (stronger pseudobulk baselines), if summarized
_b2 = os.path.join(HERE, "..", "rev3", "results", "B2_summary.json")
if os.path.exists(_b2):
    B2 = json.load(open(_b2)); f0 = lambda x: f"{x:,.0f}"
    for m, lab in (("PB", "PB"), ("PBres", "PBres"), ("PBw", "PBw")):
        v = [B2[s][m]["meanV"] for s in ("P1", "P2", "P3")]
        M[f"Btwo{lab}"] = f"{f0(min(v))}--{f0(max(v))}"
        for s in ("P1", "P2", "P3"):
            M[f"Btwo{lab}{s.replace('P1','One').replace('P2','Two').replace('P3','Three')}"] = f0(B2[s][m]["meanV"])
            M[f"Btwo{lab}Any{s.replace('P1','One').replace('P2','Two').replace('P3','Three')}"] = f"{B2[s][m]['any']}/{B2[s][m]['n']}"
            M[f"Btwo{lab}TI{s.replace('P1','One').replace('P2','Two').replace('P3','Three')}"] = f"{B2[s][m]['t05']:.3f}"
    M["BtwoPBwRangeLow"] = f"{min(B2[s]['PBw']['minV'] for s in ('P1','P3'))}--{max(B2[s]['PBw']['maxV'] for s in ('P1','P3'))}"

# REV2 C (real donors with XDS-M; A5 diagnostic; positive control)
if S.get("C"):
    C = S["C"]; ss = [s for s in ("P1", "P2", "P3") if s in C]; pct = lambda x: f"{100 * x:.0f}\\%"
    M["RealXDSMfrac"] = ", ".join(pct(C[s]["seeds"]["frac_pairs"]) for s in ss)
    M["RealXDSMcalls"] = ", ".join(f"{C[s]['seeds']['n_calls'][0]:,}" for s in ss)
    J = [x for s in ss for x in C[s]["seeds"]["jaccard"]]; M["RealXDSMjac"] = f"{min(J):.2f}--{max(J):.2f}"
    R = [C[s]["halves"]["replicated"] / max(C[s]["halves"]["calls_h0"], 1) for s in ss]; M["RealXDSMrep"] = f"{100 * min(R):.0f}--{100 * max(R):.0f}\\%"
    r2 = [x for s in ss for x in C[s]["A5"]["mean_R2"]]; M["AfiveRtwo"] = f"{min(r2):.2f}--{max(r2):.2f}"
    M["PosXYnM"] = str(sum(C[s]["sex_XDSM"]["xist_y"]["n"] for s in ss)); M["PosYYnM"] = str(sum(C[s]["sex_XDSM"]["y_y"]["n"] for s in ss))
    M["PosXYokM"] = str(sum(C[s]["sex_XDSM"]["xist_y"]["sig_right_sign"] for s in ss)); M["PosYYokM"] = str(sum(C[s]["sex_XDSM"]["y_y"]["sig_right_sign"] for s in ss))
    def sci(x):
        e = math.floor(math.log10(x)); c = math.ceil(x / 10 ** e - 1e-12)
        c, e = (1, e + 1) if c == 10 else (c, e); return f"{c}\\times10^{{{e}}}"
    M["PosXYmaxpM"] = sci(max(C[s]["sex_XDSM"]["xist_y"]["max_p"] for s in ss)); M["PosYYmaxpM"] = sci(max(C[s]["sex_XDSM"]["y_y"]["max_p"] for s in ss))
    M["PosXYokPB"] = str(sum(C[s]["sex_PB"]["xist_y"]["sig_right_sign"] for s in ss)); M["PosYYokPB"] = str(sum(C[s]["sex_PB"]["y_y"]["sig_right_sign"] for s in ss))
    M["PosYYmaxpPB"] = sci(max(C[s]["sex_PB"]["y_y"]["max_p"] for s in ss)); M["PosXYmaxpPB"] = sci(max(C[s]["sex_PB"]["xist_y"]["max_p"] for s in ss))

# REV2 D (differential co-expression): D1 donor-centred null and real, D3 halves, D2 true-DC simulation
if S.get("D1"):
    D1 = S["D1"]; ss = [s for s in ("P1", "P2", "P3") if s in D1]
    nul = {m: [D1[s][k][m] for s in ss for k in D1[s] if k.startswith("k")] for m in ("PB", "XDSg", "XDSMg", "XDSMgu")}
    rg = lambda v, f=lambda x: f"{x:,}": f"{f(min(v))}--{f(max(v))}" if min(v) != max(v) else f(min(v))
    M["DCdcPB"] = rg([a["V"] for a in nul["PB"]]); M["DCdcPBTI"] = rg([a["t05"] for a in nul["PB"]], lambda x: f"{x:.2f}")
    xs = nul["XDSg"] + nul["XDSMg"] + nul["XDSMgu"]; M["DCdcXDS"] = rg([a["V"] for a in xs]); M["DCdcXDSTI"] = rg([a["t05"] for a in xs], lambda x: f"{x:.3f}")
    M["DCdcRuns"] = str(len(nul["PB"]) // max(len(ss), 1))
    for m, lab in (("PB", "PB"), ("XDSMg", "M"), ("XDSMgu", "U"), ("XDSg", "G")):
        if all("real" in D1[s] for s in ss): M[f"DCreal{lab}"] = ", ".join(f"{D1[s]['real'][m]['V']:,}" for s in ss)
if S.get("D3") and "P1" in S["D3"]:
    d3 = S["D3"]["P1"]
    for m, lab in (("PB", "PB"), ("XDSMgu", "U"), ("XDSMg", "M")):
        M[f"DCthree{lab}"] = f"{d3[m]['replicated']:,} of {d3[m]['calls_h0']:,}"
if S.get("D2"):
    d2 = S["D2"]; M["DCtwoAlt"] = f"{d2['n_alt']:,}"; M["DCtwoPairs"] = f"{d2['n_pairs']:,}"; M["DCtwoR"] = str(d2["PB"]["R"])
    for m, lab in (("PB", "PB"), ("XDSMg", "M"), ("XDSMgu", "U"), ("XDSg", "G")):
        if m in d2: M[f"DCtwoFDR{lab}"] = f"{d2[m]['FDR']:.3f}"; M[f"DCtwoPow{lab}"] = f"{d2[m]['power']:.3f}"; M[f"DCtwoFDRse{lab}"] = f"{d2[m]['FDR_mcse']:.3f}"
# OneK1K E (30 further reassignments)
if S.get("E"):
    E = S["E"]; tt = [t_ for t_ in ("B1", "B2", "B3") if t_ in E]
    M["OnekEruns"] = ", ".join(f"{E[t_]['PB']['n']}" for t_ in tt)
    for m, lab in (("PB", "PB"), ("PBlin", "PBlin"), ("LCP", "LCP"), ("XDSM", "XDSM")):
        M[f"OnekE{lab}"] = ", ".join(f"{E[t_][m]['meanV']:,.0f}" if E[t_][m]['meanV'] >= 10 else f"{E[t_][m]['meanV']:.1f}" for t_ in tt)
        M[f"OnekEany{lab}"] = ", ".join(f"{E[t_][m]['any']}/{E[t_][m]['n']}" for t_ in tt)
    _nm = {"B1": "transitional B cells", "B2": "naive B cells", "B3": "naive CD4 T cells"}
    fl = [_nm[t_] for t_ in tt if not holds(E[t_]["XDSM"])]; ok = [_nm[t_] for t_ in tt if holds(E[t_]["XDSM"])]
    M["OnekEfailXDSM"] = (" and ".join(fl) if fl else "none of the cell types"); M["OnekEholdXDSM"] = ("all three cell types" if len(ok) == 3 else (" and ".join(ok) if ok else "none"))

# REV3 A2 (tail calibration by expression stratum), if summarized
_a2 = os.path.join(HERE, "..", "rev3", "results", "A2_SUMMARY.json")
if os.path.exists(_a2):
    A2 = json.load(open(_a2))
    kk = [k for k in A2 if k.startswith("B1_")]
    if kk:
        npairs = sum(A2[kk[0]]["all"]["PB"][s]["pairs"] for s in A2[kk[0]]["all"]["PB"] if s.isdigit())
        fr = [A2[k]["all"]["PB"]["Vtotal_mean"] / npairs for k in kk]
        M["AtwoPBfrac"] = f"{100 * min(fr):.0f}--{100 * max(fr):.0f}\\%"; M["AtwoPairs"] = f"{npairs / 1e6:.1f} million"
        M["AtwoPBcalls"] = f"{min(A2[k]['all']['PB']['Vtotal_mean'] for k in kk) / 1e6:.2f}--{max(A2[k]['all']['PB']['Vtotal_mean'] for k in kk) / 1e6:.2f} million"
        M["AtwoXDSMcalls"] = ", ".join(f"{A2[k]['all']['XDSM5']['Vtotal_mean']:.1f}" for k in kk)
        M["AtwoRuns"] = ", ".join(str(A2[k]["runs"]) for k in kk)
    def _rr(keys, m, strata, thr_i, Ds=("200", "all")):
        v = [A2[k][D][m][s]["ratio"][thr_i] for k in keys for D in Ds for s in strata if s in A2[k][D][m] and A2[k][D][m][s]["exp"][thr_i] >= 20]
        return (min(v), max(v)) if v else None
    def _dec(keys, m, strata, Ds=("200", "all")):
        c = [A2[k][D][m][s]["calibrated"] for k in keys for D in Ds for s in strata if s in A2[k][D][m]]
        return sum(x is True for x in c), sum(x is not None for x in c)
    hi_keys = [k for k in A2 if k.split("_")[0] in ("B1", "B2")]; lo_keys = [k for k in A2 if k.startswith("B1x")]
    fmt = lambda t: ("n/a" if not t else (f"{t[0]:,.0f}--{t[1]:,.0f}" if t[0] >= 100 else f"{t[0]:.2f}--{t[1]:.2f}"))
    for m, lab in (("XDSM5", "M"), ("XDS", "S"), ("XDSM1", "Mone"), ("XDSM10", "Mten"), ("PB", "PB")):
        M[f"AtwoHiFour{lab}"] = fmt(_rr(hi_keys, m, ("0", "1", "2"), 1)); M[f"AtwoHiFive{lab}"] = fmt(_rr(hi_keys, m, ("0", "1", "2"), 2))
        M[f"AtwoLoFour{lab}"] = fmt(_rr(lo_keys, m, ("3",), 1)); M[f"AtwoLoFive{lab}"] = fmt(_rr(lo_keys, m, ("3",), 2))
        a_, b_ = _dec(hi_keys, m, ("0", "1", "2")); M[f"AtwoHiDec{lab}"] = f"{a_} of {b_}"
        a_, b_ = _dec(lo_keys, m, ("3",)); M[f"AtwoLoDec{lab}"] = f"{a_} of {b_}"
    if "P3_n0" in A2 or "P3_n1" in A2: M["AtwoHasPthree"] = "yes"
    pk = [k for k in A2 if k.startswith("P3_")]
    if pk:
        def _dec_all(keys, m):
            c = [A2[k][D][m][s]["calibrated"] for k in keys for D in ("200", "all") for s in A2[k][D][m] if s.isdigit()]
            return sum(x is True for x in c), sum(x is not None for x in c)
        a_, b_ = _dec_all(pk, "XDSM5"); M["AtwoPthreeDecM"] = f"{a_} of {b_}"
        a_, b_ = _dec_all(pk, "XDS"); M["AtwoPthreeDecS"] = f"{a_} of {b_}"
        a_, b_ = _dec_all(pk, "PB"); M["AtwoPthreeDecPB"] = f"{a_} of {b_}"
        pv = [A2[k][D]["PB"]["Vtotal_mean"] for k in pk for D in ("200", "all")]
        M["AtwoPthreePBcalls"] = f"{min(pv) / 1e3:,.0f},000--{max(pv) / 1e3:,.0f},000"
        M["AtwoPthreeRuns"] = str(min(A2[k]["runs"] for k in pk))
    M["AtwoRunsAll"] = "; ".join(f"{k.replace('_n0', ' pooled').replace('_n1', ' donor-centred')}: {A2[k]['runs']}" for k in sorted(A2))

# REV2 F (A5 overdispersion stress), if summarized
if S.get("F"):
    for d in ("L200", "L500"):
        if d in S["F"]:
            lab = {"L200": "Ftwo", "L500": "Ffive"}[d]
            for m, ml in (("PB", "PB"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDS", "XDS"), ("XDSM", "XDSM"), ("KX", "KX")):
                v = S["F"][d][m]; M[f"{lab}FDR{ml}"] = f"{v['FDR']:.3f}"; M[f"{lab}se{ml}"] = f"{v['FDR_mcse']:.3f}"; M[f"{lab}Pow{ml}"] = f"{v['power']:.2f}"
            M[f"{lab}R"] = str(S["F"][d]["PB"]["R"])

# REV3 A3 (small-D run-level) and B3 (clip sensitivity)
_a3 = os.path.join(HERE, "..", "rev3", "results", "A3B3_SUMMARY.json")
if os.path.exists(_a3):
    AB = json.load(open(_a3)); A3 = {s: v for s, v in AB["A3"].items() if v}; B3 = AB["B3"].get("P3")
    MS3 = [("PB", "PB"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDSperm", "XDS (perm.)"), ("XDSM", "XDS-M"), ("KX", "KX")]
    subs = [s for s in ("P1", "P2", "P3") if s in A3]; Ds3 = ("25", "50", "100")
    L3 = ["\\begin{table}[h]\\centering\\small\\setlength{\\tabcolsep}{3pt}",
          "\\caption{Run-level error at small cohort sizes on the donor-centred null (frozen plan REV3\\_PREREG, analysis A3): reassignments with at least one BH call, out of " + ", ".join(f"{s}: {A3[s]['n']}" for s in subs) + ", at nested subsets of 25, 50 and 100 pseudo-donors. XDS (perm.), LCP and LCPd use donor-permutation $p$-values with $10^5$ permutations, whose smallest attainable value ($10^{-5}$) exceeds the first BH threshold among $\\sim$89{,}000 pairs, so their zeros here are uninformative; Table~\\ref{tab:R1} repeats them with $2\\times10^6$ permutations. In the last block (D $\\ge$ 200) XDS uses the normal-moment reference. Last block: P3 with the donor-centring factors clipped to $[0.5,2]$ instead of $[0.1,10]$ (analysis B3, " + (str(B3["n"]) if B3 else "--") + " reassignments); with the narrower clip single-split XDS fails less often than at $[0.1,10]$, so part of its failure depends on the clip. $^\\dagger$ fails the frozen rule.}\\label{tab:smallD}",
          "\\begin{tabular}{l" + "ccc" * len(subs) + ("|ccc" if B3 else "") + "}\\toprule",
          " & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{{s}}}" for s in subs) + (" & \\multicolumn{3}{c}{P3, clip $[0.5,2]$}" if B3 else "") + "\\\\",
          "Method & " + " & ".join(D for _ in subs for D in Ds3) + (" & 200 & 400 & all" if B3 else "") + "\\\\\\midrule"]
    for m, lab in MS3:
        row = [cell(A3[s][D].get(m)) for s in subs for D in Ds3] + ([cell(B3[D].get(m)) for D in ("200", "400", "all")] if B3 else [])
        L3.append(lab + " & " + " & ".join(row) + "\\\\")
    L3.append("\\bottomrule\\end{tabular}\\end{table}")
    open(os.path.join(HERE, "table_smallD.tex"), "w").write("\n".join(L3) + "\n")
    def rng3(m, Ds):
        v = [A3[s][D][m]["any"] / A3[s][D][m]["n"] for s in subs for D in Ds]; a_, b_ = round(100 * min(v)), round(100 * max(v))
        return f"{a_}\\%" if a_ == b_ else f"{a_}--{b_}\\%"
    M["AthreeSubs"] = ", ".join(f"{s} ({A3[s]['n']} runs)" for s in subs)
    M["AthreeXDSpermAny"] = str(sum(A3[s][D]["XDSperm"]["any"] for s in subs for D in Ds3)); M["AthreeXDSpermN"] = str(sum(A3[s][D]["XDSperm"]["n"] for s in subs for D in Ds3))
    M["AthreeXDSMtwentyfive"] = rng3("XDSM", ("25",)); M["AthreeXDSMfifty"] = rng3("XDSM", ("50",)); M["AthreeXDSMhundred"] = rng3("XDSM", ("100",))
    M["AthreeLCPsmall"] = rng3("LCP", ("25", "50")); M["AthreeLCPhundred"] = rng3("LCP", ("100",)); M["AthreeKXtwentyfive"] = rng3("KX", ("25",)); M["AthreeLCPdtwentyfive"] = rng3("LCPd", ("25",))
    M["AthreeXDSMhundredHolds"] = "all" if all(holds(A3[s]["100"]["XDSM"]) for s in subs) else "not all"
    _w = {0: "none", 1: "one", 2: "two", 3: "three"}
    for m_, ml in (("XDSM", "XDSM"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("KX", "KX"), ("XDSperm", "XDSperm")):
        for D, dl in (("25", "TwentyFive"), ("50", "Fifty"), ("100", "Hundred")):
            M[f"AthreeHolds{ml}{dl}"] = f"{_w.get(sum(holds(A3[s][D][m_]) for s in subs), str(sum(holds(A3[s][D][m_]) for s in subs)))} of {_w.get(len(subs), str(len(subs)))}"
    M["AthreeLCPfailCells"] = ", ".join(f"{s} at {D}" for s in subs for D in Ds3 if not holds(A3[s][D]["LCP"])) or "none"
    if B3:
        M["BthreeXDSM"] = ", ".join(f"{B3[D]['XDSM']['any']}/{B3[D]['XDSM']['n']}" for D in ("200", "400", "all"))
        M["BthreeLCPd"] = ", ".join(f"{B3[D]['LCPd']['any']}/{B3[D]['LCPd']['n']}" for D in ("200", "400", "all"))
        M["BthreeXDS"] = ", ".join(f"{B3[D]['XDS']['any']}/{B3[D]['XDS']['n']}" for D in ("200", "400", "all")); M["BthreeR"] = str(B3["n"])

# Outcome strings for the pre-registration appendix
try:
    hc = {lab: f"{sum(holds(a) for a in runs(m))} of {len(runs(m))}" for m, lab in MS}
    fails = [f"{s} at {D}" for s in ("P1", "P2", "P3") if s in S["A"] for D in ("200", "400", "all") if not holds(S["A"][s][D]["XDSM"])]
    o2 = (f"(A) cells within the run-level bound: XDS-M {hc['XDS-M']}" + (f" (fails {', '.join(fails)})" if fails else "") + f", KX {hc['KX']}, LCPd {hc['LCPd']}, single-split XDS {hc['XDS']}, pseudobulk and LCP {hc['PB']}; LCPd is therefore not co-recommended. ")
    if S.get("C"): o2 += "(B, C) no trend of the standardized mean square with cell number; real-donor XDS-M results in Section~\\ref{sec:real}. "
    if S.get("D1"): o2 += f"(D) pseudobulk differential co-expression confounded on the donor-centred null (\\DCdcPB{{}} false calls); group-centred XDS-M variants \\DCdcXDS{{}}. "
    if S.get("E"): o2 += f"(E) XDS-M within the bound over 33 runs in \\OnekEholdXDSM; PBlin made thousands of false calls, so the transformation is not the main cause. "
    if S.get("F"): o2 += "(F) reported in Section~\\ref{sec:null}."
    M["REVtwoOutcomes"] = o2
except Exception as e:
    print("outcome2 failed", e)
try:
    o3 = []
    if "AtwoHiDecM" in M: o3.append(f"(A2) XDS-M (S = 5) met the calibration rule in {M['AtwoHiDecM']} combinations for genes above 0.1 UMI per cell and in {M['AtwoLoDecM']} for 0.03--0.1 UMI per cell; pseudobulk in none")
    if "AtwoPthreeDecM" in M: o3.append(f"(A2-P3, ADAM\\_GPNMB, 2{{,}}000 genes all above 0.3 UMI per cell) XDS-M met the rule in {M['AtwoPthreeDecM']}, single-split XDS in {M['AtwoPthreeDecS']}, pseudobulk in {M['AtwoPthreeDecPB']}")
    if "AthreeHoldsXDSMHundred" in M: o3.append(f"(A3) XDS with permutation $p$-values held at every size; XDS-M held at 100 donors in {M['AthreeHoldsXDSMHundred']} subtypes, at 50 in {M['AthreeHoldsXDSMFifty']}, at 25 in {M['AthreeHoldsXDSMTwentyFive']}; LCP held at 25 and 50 and failed at 100")
    if "BtwoPBw" in M: o3.append("(B2) the frozen criterion that covariate adjustment does not fix pseudobulk was not met, because precision weighting left fewer than 100 mean false calls in two subtypes, although it made at least one false call in every run")
    if "BthreeXDSM" in M: o3.append(f"(B3) descriptive: XDS-M {M['BthreeXDSM']}, LCPd {M['BthreeLCPd']} runs with a false call")
    M["REVthreeOutcomes"] = "; ".join(o3) + "." if o3 else "[pending]"
except Exception as e:
    print("outcome3 failed", e)

# XDS-M: settings with >= 100 donors that hold the frozen run-level bound (REV2 A cells, REV3 A3 at 100, OneK1K E over 33 runs)
try:
    _cells = [S["A"][s][D]["XDSM"] for s in ("P1", "P2", "P3") if s in S["A"] for D in ("200", "400", "all")]
    if "A3" in dir() or os.path.exists(_a3):
        _ab = json.load(open(_a3))["A3"]; _cells += [_ab[s]["100"]["XDSM"] for s in ("P1", "P2", "P3") if _ab.get(s)]
    if S.get("E"): _cells += [S["E"][t_]["XDSM"] for t_ in ("B1", "B2", "B3") if t_ in S["E"]]
    _w = {i: w for i, w in enumerate("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen".split())}
    M["XDSMholdsLarge"] = f"{_w.get(sum(holds(c) for c in _cells), sum(holds(c) for c in _cells))} of {_w.get(len(_cells), len(_cells))}"
except Exception as e:
    print("XDSMholdsLarge failed", e)
M.setdefault("REVtwoOutcomes", S.get("outcomes2", "[pending]")); M.setdefault("REVthreeOutcomes", S.get("outcomes3", "[pending]"))
open(os.path.join(HERE, "rev2nums.tex"), "w").write("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
json.dump(M, open(os.path.join(HERE, "rev2_macros.json"), "w"), indent=1); print(len(M), "macros")
