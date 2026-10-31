"""REV2_PREREG D2: true differential co-expression simulation. python rev2_D2.py"""
import os, sys, json, numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
from methods_v2 import bh  # noqa: E402
from dstd import design_R  # noqa: E402
from revm import stats_from_cells  # noqa: E402
pan = S.make_panel(seed=4112026)
c0 = dict(DESIGNS["D40"]); c0.update(D=300, nd_median=15, nd_min=4)
c1 = dict(DESIGNS["D40"]); c1.update(D=300, nd_median=40, nd_min=4, g=0.45)
st0, st1 = S.setting(**c0), S.setting(**c1); t0, t1 = S.truth(pan, st0), S.truth(pan, st1); iu = t0["iu"]; A_, B_ = iu
alt = np.abs(t1["sb"] - t0["sb"]) > 1e-12
od = os.path.join(HERE, "out", "D2"); os.makedirs(od, exist_ok=True)
def pear(Mx):
    Mc = Mx - Mx.mean(0); cov = Mc.T @ Mc / (len(Mx) - 1); sd = np.sqrt(np.diag(cov)); return np.clip((cov / np.outer(sd, sd))[iu], -0.999999, 0.999999)
def welch(q1, q0):
    z = (q1.mean(0) - q0.mean(0)) / np.sqrt(q1.var(0, ddof=1) / len(q1) + q0.var(0, ddof=1) / len(q0)); return 2 * norm.sf(np.abs(z))
res = {m: dict(fdp=[], pw=[]) for m in ("PB", "XDSg", "XDSMg", "XDSMgu")}
for rep in range(50):
    grp_out = []
    for g, st in ((0, st0), (1, st1)):
        sim = S.simulate(pan, st, 9_800_000 + 2 * rep + g); rows = [r for r in sim["donor_idx"] if len(r) >= 4]
        cells = [(sim["Y"][r].astype(np.float64), sim["s_hat"][r].astype(np.float64)) for r in rows]
        f_, hs = stats_from_cells(cells, 9_800_000 + 2 * rep + g)
        qs = []
        for h in hs:
            RA, _ = design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"]); RB, _ = design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"]); qs.append(RA[:, A_] * RB[:, B_])
        qu = [(h["mA"] - h["mA"].mean(0))[:, A_] * (h["mB"] - h["mB"].mean(0))[:, B_] for h in hs]   # REV2 Amendment 1
        grp_out.append(dict(L=np.log1p(1e6 * f_["m"]), q1=qs[0], qM=sum(qs) / len(qs), qU=sum(qu) / len(qu)))
    g0, g1 = grp_out; n0, n1 = len(g0["L"]), len(g1["L"])
    p = {"PB": 2 * norm.sf(np.abs((np.arctanh(pear(g1["L"])) - np.arctanh(pear(g0["L"]))) / np.sqrt(1 / (n1 - 3) + 1 / (n0 - 3)))),
         "XDSg": welch(g1["q1"], g0["q1"]), "XDSMg": welch(g1["qM"], g0["qM"]), "XDSMgu": welch(g1["qU"], g0["qU"])}
    for m, pv in p.items():
        rj = bh(pv); n = rj.sum(); res[m]["fdp"].append(float((rj & ~alt).sum() / max(n, 1))); res[m]["pw"].append(float((rj & alt).sum() / alt.sum()))
    print(rep, {m: (round(res[m]["fdp"][-1], 3), round(res[m]["pw"][-1], 3)) for m in res}, flush=True)
summ = {m: dict(FDR=float(np.mean(v["fdp"])), FDR_mcse=float(np.std(v["fdp"], ddof=1) / np.sqrt(len(v["fdp"]))), power=float(np.mean(v["pw"])),
                power_mcse=float(np.std(v["pw"], ddof=1) / np.sqrt(len(v["pw"]))), R=len(v["fdp"])) for m, v in res.items()}
summ["n_alt"] = int(alt.sum()); summ["n_pairs"] = int(len(alt))
json.dump(summ, open(os.path.join(HERE, "out", "D2_summary.json"), "w"), indent=1); print(summ)
