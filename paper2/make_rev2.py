"""Run-level table and macros from REV2_SUMMARY.json (REV2_PREREG). python make_rev2.py <REV2_SUMMARY.json>"""
import sys, json, math, os
S = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "rev2_summary.json")))
HERE = os.path.dirname(os.path.abspath(__file__))
def holds(a): return a["any"] / a["n"] <= 0.05 + 2 * math.sqrt(0.05 * 0.95 / a["n"])
MS = [("PB", "PB"), ("LCP", "LCP"), ("LCPd", "LCPd"), ("XDS", "XDS"), ("XDSM", "XDS-M"), ("KX", "KX")]
cols = [(s, D) for s in ("P1", "P2", "P3") for D in ("200", "400", "all") if s in S["A"]]
ecols = [t for t in ("B1", "B2", "B3") if t in S.get("E", {})]
def cell(a):
    if a is None: return "--"
    x = f"{a['any']}"; return x if holds(a) else x + "$^\\dagger$"
L = ["\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{2.6pt}",
     "\\caption{Run-level error on the donor-centred real-cell nulls (frozen plan REV2\\_PREREG): number of reassignments with at least one BH call, which under a complete null equals the run-level FDR numerator. PsychAD: " + ", ".join(sorted({f'{s}: {S["A"][s]["all"]["PB"]["n"]}' for s in S["A"]})) + " reassignments per subtype at nested subsets of 200, 400 and all pseudo-donors. OneK1K: " + (", ".join(f'{t}: {S["E"][t]["PB"]["n"]}' for t in ecols) if ecols else "--") + " reassignments at full size. $^\\dagger$ fails the frozen rule (fraction $>0.05+2$\\,SE).}\\label{tab:runlevel}",
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
    M["OnekEfailXDSM"] = (" and ".join(fl) if fl else "none of the cell types"); M["OnekEholdXDSM"] = (" and ".join(ok) if ok else "none")

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
M.setdefault("REVtwoOutcomes", S.get("outcomes2", "[pending]")); M.setdefault("REVthreeOutcomes", S.get("outcomes3", "[pending]"))
open(os.path.join(HERE, "rev2nums.tex"), "w").write("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
json.dump(M, open(os.path.join(HERE, "rev2_macros.json"), "w"), indent=1); print(len(M), "macros")
