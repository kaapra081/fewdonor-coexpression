"""Generates nums.tex, table_sim.tex, table_null.tex, table_ablate.tex, nums_provenance.json and the compact JSON summaries
(r1_summary.json, onek_summary.json, dc_summary.json, dcg_summary.json, seanull_summary.json) from result files.

    python paper2/make_nums.py [--outdir DIR] [--strict] [--audit]

Input layouts. A path such as 'rev/out/R1_P1_k0.json' is looked up as <repo>/rev/out/... (the authors' working layout) and as
<repo>/results/rev/out/... (the committed layout). Inputs fall into three classes:
  1. Committed JSON results (results/...): required. A missing file stops the script and names both paths that were tried.
  2. Raw per-run arrays (*.npz, gitignored, not in the repository): used when present; the compact summary that stands in for
     them in a clean clone (paper2/*_summary.json) is then rewritten from them.
  3. Compact summaries (paper2/*_summary.json, committed): used when the raw arrays are absent (a NOTE is printed).
     If neither exists the block is SKIPPED with a WARNING naming the missing files: its macros keep the value found in the
     committed nums.tex (listed as STALE in nums_provenance.json) and its table file is not rewritten.
No block is skipped silently, and a committed summary is never overwritten unless it was recomputed from raw inputs.
--strict exits with status 2 if any block was skipped; --audit lists every input (found/missing, size) and exits;
--outdir writes all outputs elsewhere (used by tests/test_make_nums.py) while still reading the committed summaries from paper2/.
"""
import os, sys, json, glob, re, argparse, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..")
_ap = argparse.ArgumentParser(description="Regenerate nums.tex, tables and compact summaries.")
_ap.add_argument("--outdir", default=H, help="directory for nums.tex, table_*.tex, *_summary.json, nums_provenance.json (default: paper2/)")
_ap.add_argument("--strict", action="store_true", help="exit with status 2 if any block was skipped")
_ap.add_argument("--audit", action="store_true", help="list the inputs this script reads (found/missing, sizes) and exit")
ARGS = _ap.parse_args(); OUT = os.path.abspath(ARGS.outdir); os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(R, "dev")); from methods_v2 import bh  # noqa: E402
WARN, NOTES = [], []
def warn(block, msg): WARN.append((block, msg)); print(f"WARNING [{block}]: {msg}", file=sys.stderr, flush=True)
def note(block, msg): NOTES.append((block, msg)); print(f"NOTE [{block}]: {msg}", file=sys.stderr, flush=True)
def cands(rel): return [os.path.join(R, rel), os.path.join(R, "results", rel)]
def find(rel):
    """First existing path of rel in the working layout (<repo>/rel) or the committed layout (<repo>/results/rel); None if absent."""
    for c in cands(rel):
        if os.path.exists(c): return c
def files(rel_glob):
    """Files matching a glob in either layout (a file present in both counts once, working layout first), sorted by name."""
    seen = {}
    for base in (R, os.path.join(R, "results")):
        for f in sorted(glob.glob(os.path.join(base, rel_glob))): seen.setdefault(os.path.basename(f), f)
    return [seen[k] for k in sorted(seen)]
def J(p):
    """Required committed JSON (working layout <dir>/<file> or committed layout results/<dir>/<file>)."""
    c = find(p)
    if c is None: raise FileNotFoundError(f"required input '{p}' not found; looked in: " + "; ".join(cands(p)))
    return json.load(open(c))
def summary(name):
    """Committed compact summary paper2/<name> (always read from the source tree, never from --outdir); None if absent or empty."""
    f = os.path.join(H, name)
    if os.path.exists(f):
        d = json.load(open(f))
        if d: return d
def write_json(name, obj): json.dump(obj, open(os.path.join(OUT, name), "w"), indent=1)
# Every input read below: (block, class, pattern relative to the repository root or results/). Used by --audit.
INPUTS = [("core", "committed", p) for p in (
    "pop/pop_results.json", "dstd/d_results.json", "xsplit/xs_results.json", "fix/fix_results.json", "xu/xu_results.json", "xu/couple_results.json",
    "split/split_results.json", "diag/diag_real.json", "pop/pop_P[123]_info.json", "rev/pooling.json", "rev/rev_dc_r5.json", "rev/rev_sim_results.json",
    "rev/out/real_P[123].json", "onek/onek_B[123]_info.json", "rare/rare_prep_meta.json", "prereg/DSTD_PREREG.md", "prereg/DC_PREREG.md")] + [
    ("R1", "committed", "rev/out/R1_P[123]_k*.json"), ("R6", "committed", "rev/out/R6_P[123]_k*.json"),
    ("DCg", "raw, gitignored", "rev/out/DCg_P[123]_k*.npz"), ("DCg", "raw, gitignored", "rev/out/DCg_P[123]_real.npz"),
    ("table_ablate", "raw, gitignored", "rev/out/seanull_S[12]_k*.npz"),
    ("OneK1K", "raw, gitignored", "onek/out/ON_B[123]_k[0-9]_D*.npz"),
    ("DC", "raw, gitignored", "dc/out/DC_P[123]_k*.npz"), ("DC", "raw, gitignored", "dc/out/DC_P[123]_real.npz")] + [
    (b, "compact summary", "paper2/" + n) for b, n in (("R1", "r1_summary.json"), ("OneK1K", "onek_summary.json"), ("DC", "dc_summary.json"),
                                                     ("DCg", "dcg_summary.json"), ("table_ablate", "seanull_summary.json"))]
if ARGS.audit:
    print(f"{'block':13s} {'class':16s} {'files':>5s} {'MB':>8s}  pattern")
    for blk, cls, pat in INPUTS:
        fs = files(pat) if "*" in pat or "[" in pat else ([find(pat)] if find(pat) else [])
        print(f"{blk:13s} {cls:16s} {len(fs):5d} {sum(os.path.getsize(f) for f in fs) / 1e6:8.3f}  {pat}" + ("   <-- MISSING" if not fs else ""))
    sys.exit(0)
OLD = {}      # macro values in the committed nums.tex: what a skipped block keeps
if os.path.exists(os.path.join(H, "nums.tex")):
    for _l in open(os.path.join(H, "nums.tex")):
        _m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s*$", _l)
        if _m: OLD[_m.group(1)] = _m.group(2)
class Macros(dict):
    """Macro table that records which block set each macro. Assigning None means 'block skipped': keep the committed value."""
    block = "?"
    def __init__(self): super().__init__(); self.src = {}
    def __setitem__(self, k, v):
        if v is None:
            if k in OLD: v = OLD[k]; self.src[k] = f"STALE: kept from committed nums.tex, block '{self.block}' skipped"
            else: v = "\\textbf{[MISSING " + k + "]}"; self.src[k] = f"MISSING: block '{self.block}' skipped and no committed value"; warn(self.block, f"macro {k} has no value")
        else: self.src[k] = self.block
        super().__setitem__(k, v)
    def hand(self, k, v, why): dict.__setitem__(self, k, v); self.src[k] = "HAND-TYPED: " + why

pop = J("pop/pop_results.json"); dres = J("dstd/d_results.json"); xs = J("xsplit/xs_results.json"); fx = J("fix/fix_results.json")
xu = J("xu/xu_results.json"); sp = J("split/split_results.json"); dg = J("diag/diag_real.json")
info = {s: J(f"pop/pop_{s}_info.json") for s in ("P1", "P2", "P3")}
S = ("P1", "P2", "P3"); SUBN = {"P1": "Homeo\\_CECR2", "P2": "PVM\\_CD163", "P3": "ADAM\\_GPNMB"}
c = lambda x: f"{int(round(x)):,}".replace(",", "{,}")
f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"
def rng(vals, fmt=c):
    lo, hi = min(vals), max(vals); return fmt(lo) if fmt(lo) == fmt(hi) else f"{fmt(lo)}--{fmt(hi)}"
A = pop["A"]; M = Macros(); M.block = "committed results/ JSON (pop, dstd, xsplit, fix, xu, split, diag)"
M.hand("NnucPsy", "289{,}279", "n_obs of the PsychAD h5ad (CELLxGENE dataset 09642770-9a4c-4fe6-9f45-67e1249adc33); not recorded in any committed file")
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
B = pop["B"]; G = {s: info[s]["G"] * (info[s]["G"] - 1) // 2 for s in S}   # pairs per subtype = G(G-1)/2 from results/pop/pop_P*_info.json (was typed in)
M["RealXDSfrac"] = rng([100 * B[s]["LCPXd"] / G[s] for s in S], lambda x: f"{x:.0f}") + "\\%"
M["RealXDSrep"] = rng([100 * B[s]["rep_LCPXd"]["replicated"] / B[s]["rep_LCPXd"]["calls_h0"] for s in S], lambda x: f"{x:.0f}") + "\\%"
M["RealPBrep"] = rng([100 * B[s]["rep_PB"]["replicated"] / B[s]["rep_PB"]["calls_h0"] for s in S], lambda x: f"{x:.0f}") + "\\%"
M["RealPBonly"] = rng([B[s]["PB_only"] for s in S])
# coupling and transformation diagnostics
M["CouplingSEAAD"] = f3(dg["seaad"]["coupling_median"]); M["CouplingFrac"] = f"{100 * dg['seaad']['frac_abs_gt_0045']:.0f}\\%"
M["CouplingFracLee"] = f"{100 * dg['lee']['frac_abs_gt_0045']:.0f}\\%"
_ds = open(find("prereg/DSTD_PREREG.md")).read()   # disclosed exploratory diagnostic (P3, k = 0) stated in the frozen plan, lines 6-7
_r, _m = re.search(r"median over genes ([0-9.]+);", _ds), re.search(r"raw ratio means do not \((-?[0-9.]+)\)", _ds)
if _r and _m: M["Rcorrln"] = _r.group(1); M["Mcorrln"] = _m.group(1).lstrip("-") if float(_m.group(1)) == 0 else _m.group(1)
else:
    warn("core", "could not parse the Rcorrln/Mcorrln diagnostic from prereg/DSTD_PREREG.md; using typed values 0.155 / 0.000")
    M.hand("Rcorrln", "0.155", "prereg/DSTD_PREREG.md (disclosed exploratory diagnostic)"); M.hand("Mcorrln", "0.000", "prereg/DSTD_PREREG.md (disclosed exploratory diagnostic)")
M["PBposFrac"] = rng([100 * v for s in S for v in sp[s]["PB_false_frac_positive"]], lambda x: f"{x:.0f}") + "\\%"
cp = J("xu/couple_results.json")
fm = lambda m: ", ".join(f3(cp[d][m]["FDR"]) for d in ("LS200", "LS500", "LS1000"))
for m in ("PB", "LCP", "XDS", "XU"): M["Couple" + m] = fm(m)

# ---------------- review-response analyses (REV_PREREG, ONEK_PREREG)
import math
M.block = "R1 (rev/out/R1_*.json)"
r1 = {}
for s_ in S:
    Z = [json.load(open(f)) for f in files(f"rev/out/R1_{s_}_k*.json")]
    if not Z: continue
    r1[s_] = {m: dict(n=len(Z), runs_any=int(sum(z[m]["V"] > 0 for z in Z)), frac=float(np.mean([z[m]["V"] > 0 for z in Z])),
                      se=float(math.sqrt(0.05 * 0.95 / len(Z))), n5=float(np.mean([z[m]["n5"] for z in Z])), n4=float(np.mean([z[m]["n4"] for z in Z])),
                      t05=float(np.mean([z[m]["t05"] for z in Z])), V=float(np.mean([z[m]["V"] for z in Z]))) for m in ("PB", "XDS", "XDSM", "KX")}
    r1[s_]["P"] = Z[0]["P"]
    for m in ("XDS", "XDSM", "KX"): r1[s_][m]["holds"] = bool(r1[s_][m]["frac"] <= 0.05 + 2 * math.sqrt(0.05 * 0.95 / len(Z)))
if all(s_ in r1 for s_ in S): write_json("r1_summary.json", r1)
else:
    _r1s = summary("r1_summary.json")
    if _r1s and all(s_ in _r1s for s_ in S): r1 = _r1s; note("R1", "rev/out/R1_*.json not found; using the committed compact summary paper2/r1_summary.json")
    else: warn("R1", "rev/out/R1_P*_k*.json (committed under results/rev/out/) not found and paper2/r1_summary.json is absent or empty: R1 macros are STALE and table_null.tex / table_ablate.tex are not rewritten")
have_r1 = all(s_ in r1 for s_ in S)
if r1:
    M["RoneXDSruns"] = rng([r1[s_]["XDS"]["runs_any"] for s_ in r1]); M["RoneN"] = rng([r1[s_]["XDS"]["n"] for s_ in r1])
    M["RoneXDSMruns"] = rng([r1[s_]["XDSM"]["runs_any"] for s_ in r1]); M["RoneKXruns"] = rng([r1[s_]["KX"]["runs_any"] for s_ in r1])
    M["RoneXDSnfive"] = rng([r1[s_]["XDS"]["n5"] for s_ in r1], lambda x: f"{x:.2f}"); M["RoneExpnfive"] = rng([r1[s_]["P"] * 1e-5 for s_ in r1], lambda x: f"{x:.2f}")
    M["RoneXDSnfour"] = rng([r1[s_]["XDS"]["n4"] for s_ in r1], lambda x: f"{x:.1f}"); M["RoneExpnfour"] = rng([r1[s_]["P"] * 1e-4 for s_ in r1], lambda x: f"{x:.1f}")
else:
    for k_ in ("RoneXDSruns", "RoneN", "RoneXDSMruns", "RoneKXruns", "RoneXDSnfive", "RoneExpnfive", "RoneXDSnfour", "RoneExpnfour"): M[k_] = None
M.block = "R6 (rev/out/R6_*.json)"
r6 = {}
for s_ in S:
    Z = [json.load(open(f)) for f in files(f"rev/out/R6_{s_}_k*.json")]
    if Z: r6[s_] = Z
if r6:
    M["RsixPBcalls"] = rng([np.mean([z["PB"]["V"] for z in r6[s_]]) for s_ in r6]); M["RsixPBTI"] = rng([np.mean([z["PB"]["t05"] for z in r6[s_]]) for s_ in r6], f2)
    M["RsixXDScalls"] = rng([z["XDS"]["V"] for s_ in r6 for z in r6[s_]]); M["RsixXDSTI"] = rng([z["XDS"]["t05"] for s_ in r6 for z in r6[s_]], f3)
    M["RsixXDSTop"] = rng([z["XDS"]["top1"] for s_ in r6 for z in r6[s_]], f3); M["RsixK"] = rng([len(r6[s_]) for s_ in r6])
    M["RsixRw"] = rng([np.mean([z["median_abs_rw"] for z in r6[s_]]) for s_ in r6], f3)
else:
    warn("R6", "rev/out/R6_P*_k*.json not found: Rsix* macros are STALE")
    for k_ in ("RsixPBcalls", "RsixPBTI", "RsixXDScalls", "RsixXDSTI", "RsixXDSTop", "RsixK", "RsixRw"): M[k_] = None
M.block = "committed results/rev (pooling, rev_dc_r5, rev_sim_results, real_P*)"
pool = J("rev/pooling.json")
M["PoolRwReal"] = rng([pool[s_]["median_abs_real"] for s_ in S], f3); M["PoolRwPseudo"] = rng([pool[s_]["median_abs_pseudo"] for s_ in S], f3)
M["PoolRatio"] = rng([100 * (pool[s_]["ratio_median"] - 1) for s_ in S], lambda x: f"{x:.0f}") + "\\%"
r35 = J("rev/rev_dc_r5.json")
M["FalseRmed"] = rng([r35[s_]["false_r_median"] for s_ in S], f2); M["FalseRgtOne"] = rng([100 * r35[s_]["false_r_gt01"] for s_ in S], lambda x: f"{x:.0f}") + "\\%"
M["FalseRgtTwo"] = rng([100 * r35[s_]["false_r_gt02"] for s_ in S], lambda x: f"{x:.0f}") + "\\%"
M.block = "DCg (rev/out/DCg_*.npz or paper2/dcg_summary.json)"
dcg = {}
if all(files(f"rev/out/DCg_{s_}_k*.npz") and find(f"rev/out/DCg_{s_}_real.npz") for s_ in S):
    for s_ in S:
        nul = [np.load(f) for f in files(f"rev/out/DCg_{s_}_k*.npz")]
        real = np.load(find(f"rev/out/DCg_{s_}_real.npz"))
        dcg[s_] = dict(null_calls=[int(bh(z["XDSg"].astype(float)).sum()) for z in nul], null_typeI=[float((z["XDSg"] < 0.05).mean()) for z in nul],
                       real_calls=int(bh(real["XDSg"].astype(float)).sum()))
    write_json("dcg_summary.json", dcg)
else:
    dcg = summary("dcg_summary.json") or {}
    if dcg: note("DCg", "rev/out/DCg_*.npz (gitignored) not found; using the committed compact summary paper2/dcg_summary.json")
    else: warn("DCg", "rev/out/DCg_P*_k*.npz and DCg_P*_real.npz (gitignored *.npz) not found and paper2/dcg_summary.json is absent: macros DCgNull, DCgTI, DCgReal are STALE (kept from committed nums.tex). Run make_nums.py once where the .npz exist and commit paper2/dcg_summary.json")
if dcg:
    M["DCgNull"] = rng([v for s_ in S for v in dcg[s_]["null_calls"]]); M["DCgTI"] = rng([v for s_ in S for v in dcg[s_]["null_typeI"]], f3)
    M["DCgReal"] = ", ".join(c(dcg[s_]["real_calls"]) for s_ in S)
else: M["DCgNull"] = M["DCgTI"] = M["DCgReal"] = None
M.block = "committed results/rev (rev_dc_r5, real_P*, rev_sim_results)"
M["DCovP"] = ", ".join(f"{r35[s_]['overlap']:,} versus {r35[s_]['expected']:,.0f}".replace(",", "{,}") for s_ in S)
M["DCgHubs"] = ", ".join("\\emph{" + g + "}" for g, _ in r35["P1"]["XDSg_top"][:6])
sx = {s_: J(f"rev/out/real_{s_}.json") for s_ in S}
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

M.block = "OneK1K (onek/out/ON_*.npz or paper2/onek_summary.json)"
NUMW = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
on = {}; on_raw = True
for t_ in ("B1", "B2", "B3"):
    fi = find(f"onek/onek_{t_}_info.json")
    if fi is None: warn("OneK1K", f"results/onek/onek_{t_}_info.json not found"); on_raw = False; continue
    inf = json.load(open(fi)); res = {"info": inf}
    for D in ("all", "200"):
        Z = [np.load(f) for f in files(f"onek/out/ON_{t_}_k[0-9]_D{D}.npz")]
        if not Z: on_raw = False; continue
        res[D] = {}
        for m in ("PB", "LCP", "LCPXd", "XDSM"):
            if m not in Z[0]: continue
            calls = [int(bh(z[m].astype(float)).sum()) for z in Z]; t1 = float(np.mean([(z[m] < 0.05).mean() for z in Z]))
            top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Z]
            tt = np.concatenate([z[m][tp] < 0.05 for z, tp in zip(Z, top)]); se = float(np.sqrt(0.05 * 0.95 / len(tt)))
            res[D][m] = dict(K=len(Z), mean_calls=float(np.mean(calls)), calls=calls, typeI=t1, top1=float(tt.mean()), holds=bool(np.mean(calls) <= 1 and t1 <= 0.055 and tt.mean() <= 0.05 + 3 * se))
    on[t_] = res
if on_raw: write_json("onek_summary.json", on)
else:
    _ons = summary("onek_summary.json")
    if _ons and all(t_ in _ons and "all" in _ons[t_] and "200" in _ons[t_] for t_ in ("B1", "B2", "B3")): on = _ons; note("OneK1K", "onek/out/ON_*.npz (gitignored) not found; using the committed compact summary paper2/onek_summary.json")
    else: warn("OneK1K", "onek/out/ON_B*_k*_D*.npz (gitignored *.npz) not found and paper2/onek_summary.json is absent or incomplete: OnekText and OnekAbs are STALE (kept from committed nums.tex)")
ok_ = [t_ for t_ in on if "all" in on[t_] and "XDSM" in on[t_]["all"]]
if ok_:
    KW = rng([on[t_]["all"]["PB"]["K"] for t_ in ok_], lambda x: NUMW.get(int(x), str(int(x))))   # reassignments per cell type (was typed in as "three")
    nm = {"B1": "transitional B cells", "B2": "naive B cells", "B3": "naive CD4 T cells"}
    parts = [f"{nm[t_]} ({c(on[t_]['info']['donors'])} donors, median {c(on[t_]['info']['median_n'])} cells)" for t_ in ok_]
    M["OnekText"] = ("In OneK1K blood \\citep{Yazar2022}, queried from CELLxGENE Census release 2025-11-08, we repeated the null in " + ", ".join(parts) +
        " (frozen plan ONEK\\_PREREG, " + KW + " reassignments each). At full cohort size pseudobulk made " +
        rng([on[t_]["all"]["PB"]["mean_calls"] for t_ in ok_]) + " false calls per reassignment (type-I " + rng([on[t_]["all"]["PB"]["typeI"] for t_ in ok_], f2) +
        "), LCP " + rng([on[t_]["all"]["LCP"]["mean_calls"] for t_ in ok_]) + ", XDS " + rng([on[t_]["all"]["LCPXd"]["mean_calls"] for t_ in ok_], lambda x: f"{x:.1f}") +
        " and XDS-M " + rng([on[t_]["all"]["XDSM"]["mean_calls"] for t_ in ok_], lambda x: f"{x:.1f}") + " (type-I " + rng([on[t_]["all"]["XDSM"]["typeI"] for t_ in ok_], f3) +
        "). Pseudobulk's false calls at full size were " + ", ".join(f"{on[t_]['all']['PB']['mean_calls'] / on[t_]['200']['PB']['mean_calls']:.1f}" for t_ in ok_)
        + " times those at 200 donors, so the frozen growth criterion (at least 2) was met in " + str(sum(on[t_]['all']['PB']['mean_calls'] >= 2 * on[t_]['200']['PB']['mean_calls'] and on[t_]['all']['PB']['mean_calls'] >= 100 for t_ in ok_))
        + " of " + str(len(ok_)) + " cell types; XDS met its criterion in " + str(sum(on[t_]['all']['LCPXd']['holds'] and on[t_]['200']['LCPXd']['holds'] for t_ in ok_)) + " and XDS-M in "
        + str(sum(on[t_]['all']['XDSM']['holds'] and on[t_]['200']['XDSM']['holds'] for t_ in ok_))
        + "; XDS-M made no false call but exceeded the top-1\\% bound at 200 transitional B cell donors (" + f3(on['B1']['200']['XDSM']['top1']) + "), a reminder that its $t$ reference needs many donors. Even in naive CD4 T cells, with a median of " + c(on['B3']['info']['median_n']) + " cells per donor, pseudobulk failed at " + c(on['B3']['info']['donors']) + " donors.")
    M["OnekAbs"] = ("In three OneK1K blood cell types, pseudobulk made " + rng([on[t_]["all"]["PB"]["mean_calls"] for t_ in ok_]) + " false calls per reassignment, even with a median of " + c(on['B3']['info']['median_n']) + " cells per donor, and XDS-M none.")
else: M["OnekText"] = None; M["OnekAbs"] = None

# differential co-expression
M.block = "DC (dc/out/DC_*.npz or paper2/dc_summary.json)"
dc = {}
if all(files(f"dc/out/DC_{s}_k*.npz") and find(f"dc/out/DC_{s}_real.npz") for s in S):
    for s in S:
        nul = [np.load(f) for f in files(f"dc/out/DC_{s}_k*.npz")]
        real = np.load(find(f"dc/out/DC_{s}_real.npz"))
        dc[s] = {m: dict(null_calls=[int(bh(z[m].astype(float)).sum()) for z in nul], null_typeI=[float((z[m] < 0.05).mean()) for z in nul],
                         real_calls=int(bh(real[m].astype(float)).sum())) for m in ("PB", "LCP", "XDS")}
        dc[s]["n1"] = int(real["n1"][0]); dc[s]["n0"] = int(real["n0"][0])
    write_json("dc_summary.json", dc)
else:
    dc = summary("dc_summary.json") or {}
    if dc: note("DC", "dc/out/DC_*.npz (gitignored) not found; using the committed compact summary paper2/dc_summary.json")
    else: warn("DC", "dc/out/DC_P*_k*.npz and DC_P*_real.npz (gitignored *.npz) not found and paper2/dc_summary.json is absent or empty: DCtext is STALE (kept from committed nums.tex)")
pbm = {s: np.mean(dc[s]["PB"]["null_calls"]) for s in S} if dc else {}
dcyes = sum(pbm[s] >= 100 for s in S) >= 2 if dc else None
xdsok = all(np.mean(dc[s]["XDS"]["null_calls"]) <= 1 and max(dc[s]["XDS"]["null_typeI"]) <= 0.055 for s in S) if dc else None
_dp = open(find("prereg/DC_PREREG.md")).read()   # AD-vs-control median nucleus counts and Mann-Whitney p, stated in the frozen plan (lines 3-4); were typed into DCtext
_dm = {n_: re.search(r"(\d+) vs (\d+) " + n_, _dp) for n_ in ("Homeo_CECR2", "PVM_CD163", "ADAM_GPNMB")}; _dpv = re.search(r"Mann-Whitney p ([0-9.]+e-[0-9]+) to", _dp)
if all(_dm.values()) and _dpv: DCM = {n_: (m_.group(1), m_.group(2)) for n_, m_ in _dm.items()}; DCPs = sci(float(_dpv.group(1)))
else:
    warn("DC", "could not parse the nucleus-count medians / Mann-Whitney p from prereg/DC_PREREG.md; DCtext uses typed values 21/36, 23/18, 17/10 and 8e-5")
    DCM = {"Homeo_CECR2": (21, 36), "PVM_CD163": (23, 18), "ADAM_GPNMB": (17, 10)}; DCPs = sci(8e-5)
_a = lambda n_: f"{DCM[n_][0]} versus {DCM[n_][1]}"
M["DCtext"] = (
    f"Cell numbers differ between AD cases and controls in all three PsychAD subtypes (medians {_a('Homeo_CECR2')} nuclei in Homeo\\_CECR2, "
    f"{_a('PVM_CD163')} in PVM\\_CD163 and {_a('ADAM_GPNMB')} in ADAM\\_GPNMB; Mann--Whitney $p\\le{DCPs}$). " "By Propositions~\\ref{prop:leak} and~\\ref{prop:transform}, "
    "groups with different cell numbers have different biases, so comparing co-expression between groups should produce false "
    "differences. On the pooled null with each pseudo-donor carrying the AD label of the real donor whose cell number it has (three reassignments per subtype; frozen plan DC\\_PREREG), "
    f"pseudobulk differential co-expression (difference of Fisher-transformed correlations) made {rng([v for s in S for v in dc[s]['PB']['null_calls']])} false calls per reassignment "
    f"(per-pair type-I {rng([v for s in S for v in dc[s]['PB']['null_typeI']], f2)}), on average {100 * pbm['P1'] / dc['P1']['PB']['real_calls']:.0f}\\%, {100 * pbm['P2'] / dc['P2']['PB']['real_calls']:.0f}\\% and {100 * pbm['P3'] / dc['P3']['PB']['real_calls']:.0f}\\% "
    f"of the {c(dc['P1']['PB']['real_calls'])}, {c(dc['P2']['PB']['real_calls'])} and {c(dc['P3']['PB']['real_calls'])} AD-associated changes it reports on real donors, "
    f"whereas a Welch test of per-donor cross-split products between groups made {rng([v for s in S for v in dc[s]['XDS']['null_calls']])} (type-I {rng([v for s in S for v in dc[s]['XDS']['null_typeI']], f3)}). "
    ) if dc else None
M["SplitText"] = ("XDS uses one random split per donor. On real PsychAD donors we repeated XDS with three split seeds (post hoc). The numbers of BH calls were "
    + "; ".join(f"{SUBN[s]} {', '.join(c(v) for v in sp[s]['n_calls'])}" for s in S)
    + ", and the Jaccard index between call sets from different seeds was "
    + "; ".join(f"{rng(sp[s]['jaccard'], f2)} ({SUBN[s]})" for s in S)
    + ". Call counts are stable, but a fifth to a third of individual calls depend on the seed; XDS-M addresses this by averaging the per-donor products over five splits, which raised the Jaccard index between seeds to \\RealXDSMjac{} (Section~\\ref{sec:calib}).")

M.block = "?"
with open(os.path.join(OUT, "nums.tex"), "w") as fh:
    for k, v in M.items(): fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
TABLES = {"nums.tex": "regenerated"}
# ---------------- tables
def cell(d, m, key="FDR"):
    v = d[m]; return f"{v[key]:.3f}\\,({v[key + '_mcse']:.3f})"
rows = []
rm = J("rare/rare_prep_meta.json")   # SEA-AD design labels (donors, median nuclei) were typed in
dsim = dres["sim"]
def get(d, m):
    if m == "XDSM": return rs[d]["XDSM"]
    if m == "LCPd": return (E3[d]["LCPd"] if d in E3 else dsim[d]["LCPd"])
    return (E3[d] if d in E3 else E4[d])[m]
MS_ = ("PB", "LCP", "LCPd", "XDS", "XDSM")
for lab, d in (("16 donors, $\\sim$40 cells", "F_D16"), ("24 donors, $\\sim$40 cells", "F_D24"), ("40 donors, $\\sim$40 cells", "F_D40"),
               ("16 donors, $\\sim$100 cells", "F_D16_n100"), (f"{rm['S1']['D']} donors, median {int(rm['S1']['n_median'])} cells", "F_S1"), (f"{rm['S2']['D']} donors, median {int(rm['S2']['n_median'])} cells", "F_S2"),
               ("200 donors, median 20", "L200"), ("500 donors, median 20", "L500"), ("1,000 donors, median 20", "L1000")):
    rows.append(" & ".join([lab] + [f"{get(d, m)['FDR']:.3f}" for m in MS_] + [f"{get(d, m)['power']:.2f}" for m in MS_]) + " \\\\")
mcse = max(get(d, m)["FDR_mcse"] for d in ("F_D16", "F_D24", "F_D40", "F_D16_n100", "F_S1", "F_S2", "L200", "L500", "L1000") for m in MS_)
open(os.path.join(OUT, "table_sim.tex"), "w").write(
    "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3pt}\\caption{Simulations with exact truth: FDR at BH $q=0.05$ and power. $R=100$ replicates for designs with up to 83 donors, $R=50$ for larger designs; the largest FDR MCSE is "
    + f"{mcse:.3f}" + " (all MCSEs in the repository). XDS-M uses a normal reference and is intended for large cohorts. KX (not shown): power \\KXpowList{} at 200, 500 and 1{,}000 donors in the large designs (FDR \\KXfdrList).}\\label{tab:sim}\n"
    "\\begin{tabular}{lcccccccccc}\\toprule & \\multicolumn{5}{c}{FDR} & \\multicolumn{5}{c}{Power}\\\\\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}\n"
    "Design & PB & LCP & LCPd & XDS & XDS-M & PB & LCP & LCPd & XDS & XDS-M\\\\\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\\end{tabular}\\end{table}\n")
TABLES["table_sim.tex"] = "regenerated"
if have_r1:
    rows = []
    for s in S:
        a = full[s]; rp = pop["rep_null"][s]
        rows.append(f"{SUBN[s]} ({c(info[s]['donors'])}; {c(info[s]['median_nd'])}) & " + " & ".join(
            f"{c(a[m]['mean_calls'])} & {a[m]['typeI']:.3f} & {a[m]['typeI_top1']:.3f}" for m in ("PB", "LCP", "LCPXd")) +
            f" & {c(np.mean(rp['PB']['replicated_false']))} & {r1[s]['XDS']['runs_any']}/{r1[s]['XDS']['n']} & {r1[s]['XDSM']['runs_any']}/{r1[s]['XDSM']['n']} \\\\")
    open(os.path.join(OUT, "table_null.tex"), "w").write(
        "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3.5pt}\\caption{PsychAD pooled null at full cohort size (pseudo-donors; median nuclei), mean over three reassignments. "
        "Calls: BH false discoveries among $\\sim$89{,}000 pairs. TI: per-pair type-I error at 0.05, over all pairs and over the 1\\% most within-correlated pairs. Rep.: PB false calls replicated (same sign, $p<0.05$) in an independent half of the pseudo-donors. Runs: reassignments with at least one false call (equal to the FDR under this complete null) among further reassignments, for XDS and XDS-M.}\\label{tab:null}\n"
        "\\resizebox{\\linewidth}{!}{\\begin{tabular}{lcccccccccccc}\\toprule & \\multicolumn{3}{c}{PB} & \\multicolumn{3}{c}{LCP} & \\multicolumn{3}{c}{XDS} & Rep. & \\multicolumn{2}{c}{Runs $V\\ge1$}\\\\"
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){12-13}\n Subtype & Calls & TI & top 1\\% & Calls & TI & top 1\\% & Calls & TI & top 1\\% & PB & XDS & XDS-M\\\\\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\\end{tabular}}\\end{table}\n")
    TABLES["table_null.tex"] = "regenerated"
else: TABLES["table_null.tex"] = "SKIPPED: needs R1 summaries"; warn("table_null", "no R1 summaries: table_null.tex NOT rewritten (the committed file stays)")
M.block = "table_ablate (rev/out/seanull_*.npz or paper2/seanull_summary.json)"
sea_top1 = {}       # top-1% type-I of XDS-M and KX in the SEA-AD real-cell null, [S1, S2]
if all(files(f"rev/out/seanull_{st_}_k*.npz") for st_ in ("S1", "S2")):
    for key in ("XDSM", "KX"):
        sea_top1[key] = []
        for st_ in ("S1", "S2"):
            Zs = [np.load(f) for f in files(f"rev/out/seanull_{st_}_k*.npz")]
            top = [np.abs(z["r_w"]) >= np.nanquantile(np.abs(z["r_w"]), 0.99) for z in Zs]
            sea_top1[key].append(float(np.concatenate([z["p_" + key][tp] < 0.05 for z, tp in zip(Zs, top)]).mean()))
    write_json("seanull_summary.json", sea_top1)
else:
    sea_top1 = summary("seanull_summary.json") or {}
    if sea_top1: note("table_ablate", "rev/out/seanull_*.npz (gitignored) not found; using the committed compact summary paper2/seanull_summary.json")
    else: warn("table_ablate", "rev/out/seanull_S1_k*.npz, seanull_S2_k*.npz (gitignored *.npz) not found and paper2/seanull_summary.json is absent: the SEA-AD column for XDS-M and KX cannot be computed, so table_ablate.tex is NOT rewritten (the committed file stays). Run make_nums.py once where the .npz exist and commit paper2/seanull_summary.json")

if have_r1 and sea_top1:
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
        elif key in ("XDSM", "KX"): sea = rng(sea_top1[key], f3)
        rows.append(f"{nm} & {stdz} & {leak} & {left} & {calls} & {ti} & {sea} \\\\")
    open(os.path.join(OUT, "table_ablate.tex"), "w").write(
        "\\begin{table}[t]\\centering\\small\\setlength{\\tabcolsep}{3pt}\\caption{Each correction and the bias it leaves. Calls and TI: PsychAD pooled null at full size, mean calls per reassignment and per-pair type-I (range over subtypes; XU, XDS-M and KX on further reassignments). "
        "SEA-AD: type-I for the 1\\% most within-correlated pairs in the SEA-AD real-cell null (79 and 83 pseudo-donors). LCP2 (not shown) had SEA-AD top-1\\% type-I "
        + rng([fx["null"][s]["LCP2"]["typeI_top1"] for s in ("S1", "S2")], f3) + ".}\\label{tab:ablate}\n"
        "\\resizebox{\\linewidth}{!}{\\begin{tabular}{lllllll}\\toprule Method & Standardization & Leak handling & Bias left & Calls & TI & SEA-AD\\\\\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\\end{tabular}}\\end{table}\n")
    TABLES["table_ablate.tex"] = "regenerated"
else:
    TABLES["table_ablate.tex"] = "SKIPPED: needs R1 and SEA-AD null summaries"
    if not have_r1: warn("table_ablate", "no R1 summaries: table_ablate.tex NOT rewritten (the committed file stays)")
# ---------------- provenance, summary banner
_hand = [k for k, v in M.src.items() if v.startswith("HAND-TYPED")]; _stale = [k for k, v in M.src.items() if v.startswith(("STALE", "MISSING"))]
write_json("nums_provenance.json", dict(macros=M.src, tables=TABLES, hand_typed=_hand, stale=_stale, notes=[f"{b}: {m_}" for b, m_ in NOTES], warnings=[f"{b}: {m_}" for b, m_ in WARN]))
print(f"\n==== make_nums summary: {len(M)} macros, {len(M) - len(_stale)} regenerated (of which {len(_hand)} hand-typed), {len(_stale)} STALE; tables: " + ", ".join(f"{k} {v.split(':')[0]}" for k, v in TABLES.items()))
for k_ in _hand: print(f"  hand-typed  \\{k_}: {M.src[k_][12:]}")
for k_ in _stale: print(f"  STALE       \\{k_}: {M.src[k_].split(': ', 1)[1]}")
for b_, m_ in NOTES: print(f"  note        [{b_}] {m_}")
for b_, m_ in WARN: print(f"  WARNING     [{b_}] {m_}")
if ARGS.strict and WARN: sys.exit(2)
print(len(M), "macros"); print({k: v for k, v in M.items() if k not in ("DCtext", "SplitText")}); print(M["DCtext"][-700:])
