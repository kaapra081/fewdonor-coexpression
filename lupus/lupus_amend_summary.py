"""Evaluate Amendments 1 and 2 against their rules. python lupus_amend_summary.py <dir A1 json> <dir A2 json> -> lupus_amend_summary.json"""
import sys, json, glob, os, numpy as np
D1 = sys.argv[1]; D2 = sys.argv[2]; SUBS = ("LM", "FC")
holds = lambda c, t: bool(sum(x > 0 for x in c) <= 2 and np.mean(t) <= 0.055)


def blk(d, mode, sub, meths, nested=None):
    R = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(d, f"{mode}_{sub}_k*.json")), key=lambda f: int(f.rsplit("_k", 1)[1][:-5]))]; o = {}
    for m in meths:
        rr = [(r[nested][m] if nested else r[m]) for r in R]; c = [x["calls"] for x in rr]; t = [x["t05"] for x in rr]
        o[m] = dict(n_runs=len(R), mean_calls=float(np.mean(c)), max_calls=int(max(c)), runs_ge1=int(sum(x > 0 for x in c)), mean_t05=float(np.mean(t)), holds=holds(c, t))
    if R and "D_lib" in R[0]: o["_D_lib_mean"] = float(np.mean([r["D_lib"] for r in R]))
    return o


S = dict(A1={}, A2={})
for sub in SUBS:
    S["A1"][sub] = dict(E1=blk(D1, "E1", sub, ("PB", "XDSM", "XDSMlib", "XDSMcoh", "XDSMpc", "XDSMlibcoh", "XDSMlibpc")), E2=blk(D1, "E2", sub, ("XDSM", "XDSMcoh", "XDSMpc")), ER=json.load(open(os.path.join(D1, f"ER_{sub}.json"))))
    S["A2"][sub] = dict(F1=blk(D2, "F1", sub, ("PB", "PBpool", "XDSM", "XDSMlib", "XDSMpool", "XDSMlibpool")), F1_DC=blk(D2, "F1", sub, ("PBDC", "XDSMgu", "XDSMgupool"), nested="DC"), F2=blk(D2, "F2", sub, ("XDSM", "XDSMpool", "PBpool")),
                        F2B=blk(D2, "F2B", sub, ("XDSMgu", "XDSMgupool"), nested="DC"), FR=json.load(open(os.path.join(D2, f"FR_{sub}.json"))))
r = {}
for v in ("XDSMcoh", "XDSMpc", "XDSMlibcoh", "XDSMlibpc"):
    r["E-a_" + v] = all(S["A1"][s]["E1"][v]["holds"] for s in SUBS)
    if v in ("XDSMcoh", "XDSMpc"): r["E-b_" + v] = all(S["A1"][s]["E2"][v]["holds"] for s in SUBS)
for v in ("XDSMpool", "XDSMlibpool"):
    r["F-a_" + v] = all(S["A2"][s]["F1"][v]["holds"] for s in SUBS)
    if v == "XDSMpool": r["F-b_" + v] = all(S["A2"][s]["F2"][v]["holds"] for s in SUBS)
r["F-c_XDSMgupool_restores"] = all(S["A2"][s]["F1_DC"]["XDSMgupool"]["holds"] for s in SUBS); r["F-c_XDSMgupool_valid"] = all(S["A2"][s]["F2B"]["XDSMgupool"]["holds"] for s in SUBS)
r["F-c_XDSMgu_raw_holds_F1"] = {s: S["A2"][s]["F1_DC"]["XDSMgu"]["holds"] for s in SUBS}
S["rules"] = r; json.dump(S, open("lupus_amend_summary.json", "w"), indent=1); print(json.dumps(r, indent=1))
