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
M.setdefault("REVtwoOutcomes", S.get("outcomes2", "[pending]")); M.setdefault("REVthreeOutcomes", S.get("outcomes3", "[pending]"))
open(os.path.join(HERE, "rev2nums.tex"), "w").write("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
json.dump(M, open(os.path.join(HERE, "rev2_macros.json"), "w"), indent=1); print(len(M), "macros")
