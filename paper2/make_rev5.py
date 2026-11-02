"""Macros for Plans 6, 6b, 6c, 6d (rev5/out/plan6*_summary.json) -> rev5nums.tex"""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(HERE, "..", "rev5", "out")
P6 = json.load(open(os.path.join(O, "plan6_summary.json"))); P6B = json.load(open(os.path.join(O, "plan6b_summary.json"))); P6C = json.load(open(os.path.join(O, "plan6c_summary.json")))
P6D = json.load(open(os.path.join(O, "plan6d_summary.json"))) if os.path.exists(os.path.join(O, "plan6d_summary.json")) else None
M = {}
def mac(k, v): M[k] = v
def f3(x): return f"{x:.3f}"
def f2(x): return f"{x:.2f}"
def rng(xs, f=f3):
    lo, hi = min(xs), max(xs); a, b = f(lo), f(hi); return a if a == b else f"{a}--{b}"
def ints(x): return f"{int(round(x)):,}".replace(",", "{,}")
S = P6["S"]; subs = ("P1", "P2", "P3")
mac("SpS", rng([S[f"{s}_full"]["s"] for s in subs], f2))
mac("SpPBfdr", rng([S[f"{s}_full"]["o"]["PB"]["FDR"] for s in subs], f2))
mac("SpXfdr", rng([S[f"{s}_full"]["o"]["XDSM"]["FDR"] for s in subs]))
mac("SpXpow", ", ".join(f2(S[f"{s}_full"]["o"]["XDSM"]["power"]) for s in subs))
mac("SpPBpow", ", ".join(f2(S[f"{s}_full"]["o"]["PB"]["power"]) for s in subs))
mac("SpKXpow", ", ".join(f2(S[f"{s}_full"]["o"]["KX"]["power"]) for s in subs))
mac("SpKXfdr", rng([S[f"{s}_full"]["o"]["KX"]["FDR"] for s in subs]))
mac("SpPBwfdr", ", ".join(f2(S[f"{s}_full"]["o"]["PBw"]["FDR"]) for s in subs))
mac("SpXrBias", f3(max(abs(x) for s in subs for x in S[f"{s}_full"]["e"]["XDSr"]["bias"])))
mac("SpPBrBias", rng([-x for s in subs for x in S[f"{s}_full"]["e"]["rPB"]["bias"]], f2))
mac("SpXfdrTwo", rng([S[f"{s}_200"]["o"]["XDSM"]["FDR"] for s in subs]))
mac("SpPBfdrTwo", rng([S[f"{s}_200"]["o"]["PB"]["FDR"] for s in subs], f2))
mac("SpXpowTwo", rng([S[f"{s}_200"]["o"]["XDSM"]["power"] for s in subs], f2))
mac("SpPBcalls", rng([S[f"{s}_full"]["o"]["PB"]["calls"] for s in subs], ints))
mac("SpXcalls", rng([S[f"{s}_full"]["o"]["XDSM"]["calls"] for s in subs], ints))
D = P6["D"]
mac("DmXany", ", ".join(f"{D[s]['XDSM']['any']}" for s in subs)); mac("DmXmean", ", ".join(f"{D[s]['XDSM']['meanV']:.1f}" for s in subs))
mac("DmXtI", ", ".join(f3(D[s]["XDSM"]["t05"]) for s in subs)); mac("DmPB", rng([D[s]["PB"]["meanV"] for s in subs], ints))
mac("DmKXany", ", ".join(f"{D[s]['KX']['any']}" for s in subs))
mac("DmXdany", ", ".join(f"{D[s]['XDSMd']['any']}" for s in subs)); mac("DmXdtI", rng([D[s]["XDSMd"]["t05"] for s in subs]))
mac("DmSDreal", rng([json.load(open(os.path.join(O, f"D_{s}_k0.json")))["sd_logrhat"]["real"] for s in subs], f2))
mac("DmSDthis", rng([json.load(open(os.path.join(O, f"D_{s}_k0.json")))["sd_logrhat"]["this"] for s in subs], f2))
mac("DmSDplain", rng([json.load(open(os.path.join(O, f"D_{s}_k0.json")))["sd_logrhat"]["plain_null"] for s in subs], f2))
C = P6B["C"]
mac("CfXcany", ", ".join(f"{C[s]['XDSMc']['any']}" for s in subs)); mac("CfXctI", ", ".join(f2(C[s]["XDSMc"]["t05p"]) for s in subs))
mac("CfXcaany", ", ".join(f"{C[s]['XDSMca']['any']}" for s in subs)); mac("CfXcatI", rng([C[s]["XDSMca"]["t05p"] for s in subs]))
mac("CfXcaADany", ", ".join(f"{C[s]['XDSMca_ad']['any']}" for s in subs))
mac("CfXplanted", rng([C[s]["XDSM"]["meanVp"] for s in subs], ints)); mac("CfPBc", rng([P6["C"][s]["PBc"]["meanV"] for s in subs], ints))
R6b = P6B["R"]
mac("CaKeep", rng([100 * R6b[s]["XDSMca"]["keep_of_XDSM"] for s in subs], lambda x: f"{x:.0f}") + "\\%")
mac("CaXISTY", ", ".join(str(R6b[s]["sex_XDSMca_all"]["XIST_Y"]["called"]) for s in subs))
mac("CaYY", ", ".join(str(R6b[s]["sex_XDSMca_all"]["Y_Y"]["called"]) for s in subs))
R = P6["R"]
mac("SexWithinXY", ", ".join(f"{R[s]['R4']['within_male']['XIST_Y']['called']}/{R[s]['R4']['within_female']['XIST_Y']['called']}" for s in subs))
mac("SexWithinMaleYY", ", ".join(str(R[s]["R4"]["within_male"]["Y_Y"]["called"]) for s in subs))
mac("SexWithinFemaleYY", ", ".join(str(R[s]["R4"]["within_female"]["Y_Y"]["called"]) for s in subs))
mac("SexLOY", ", ".join(f"{R[s]['R4']['LOY_spearman_age_Yscore_males']:.2f}" for s in subs))
mac("SexMismatch", ", ".join(str(R[s]["R4"]["mismatch"]) for s in subs))
mac("LibKeep", rng([100 * R[s]["R3"]["XDSMlib_keep_of_t1"] for s in subs], lambda x: f"{x:.0f}") + "\\%")
Lr = P6["L"]; Lb = P6B["L"]
mac("LnXany", ", ".join(str(Lr[s]["XDSM"]["any"]) for s in subs)); mac("LnXmean", ", ".join(f"{Lr[s]['XDSM']['meanV']:.0f}" for s in subs))
mac("LnSrcAny", ", ".join(str(Lb[s]["XDSMsrc"]["any"]) for s in subs)); mac("LnSrcMean", ", ".join(f"{Lb[s]['XDSMsrc']['meanV']:.1f}" for s in subs))
mac("LnOpp", ", ".join(f"{100 * Lb[s]['sign']['frac_opposite']:.0f}" for s in subs)); mac("LnOppBase", rng([100 * Lb[s]["sign"]["baseline"] for s in subs], lambda x: f"{x:.0f}"))
for t in ("B1", "B2", "B3"): mac(f"OnekReal{'ABC'['B1B2B3'.index(t)//2]}", f"{ints(R[t]['calls']['PB'])}/{ints(R[t]['calls']['XDSM'])}")
mac("OnekRealPB", ", ".join(ints(R[t]["calls"]["PB"]) for t in ("B1", "B2", "B3"))); mac("OnekRealX", ", ".join(ints(R[t]["calls"]["XDSM"]) for t in ("B1", "B2", "B3")))
O6 = P6C["O"]; mac("OnXdany", ", ".join(f"{O6[t]['XDSMd']['any']}" for t in ("B1", "B2", "B3"))); mac("OnXdaany", ", ".join(f"{O6[t]['XDSMda']['any']}" for t in ("B1", "B2", "B3")))
mac("XrSpear", rng([R[s]["R2"]["spearman_rPB_XDSr"] for s in subs], f2))
if P6D:
    for k, v in P6D.get("macros", {}).items(): mac(k, v)
with open(os.path.join(HERE, "rev5nums.tex"), "w") as fh:
    for k, v in M.items(): fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(len(M), "macros")
