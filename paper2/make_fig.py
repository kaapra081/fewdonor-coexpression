"""Figure 2 (growth): (a) mean BH false calls vs pseudo-donors on the pooled null (log axis; exact zeros on a separate '0' strip);
(b) fraction of runs with >= 1 false call on the donor-centred null, split by reference distribution:
    left  = permutation p-values (REV4 R1, 2e6 permutations; D <= 100; LCP, LCPd, XDS one split),
    right = normal reference (REV3 A3 at 25-100 for PB, KX, XDS-M; REV2 A at 200-full for all methods).
All plotted numbers are read from the result JSONs (repo/results/pop, rev2, rev3, rev4); nothing is re-estimated.
Exported at the design width 6.5 in (= \\linewidth), so printed text size == nominal size (>= 7 pt).
Run from anywhere:  python paper2/make_fig.py   (or exec() it from a kernel started in the directory above repo/).
"""
import json, os, sys, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator, FixedLocator, FixedFormatter
from matplotlib.patches import ConnectionPatch
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.path.abspath("repo/paper2")
sys.path.insert(0, HERE)
from figstyle import apply_figure_style, panel_letter, overlap_report, audit_pdf
RR = os.path.dirname(HERE)
R = os.environ.get("RV2", "")          # optional second root holding the result trees
def J(p):
    for root in [RR] + ([R] if R else []):
        for sub in ("", "results"):
            c = os.path.join(root, sub, p)
            if os.path.exists(c): return json.load(open(c))
    raise FileNotFoundError(p)
pop = J("pop/pop_results.json"); info = {s: J(f"pop/pop_{s}_info.json") for s in ("P1", "P2", "P3")}
S2 = J("rev2/interim/REV2_SUMMARY_latest.json"); A3 = J("rev3/results/A3B3_SUMMARY.json")["A3"]; R1 = J("rev4/results/REV4_SUMMARY.json")["R1"]
apply_figure_style(sizes=(8, 7, 7))
W = 6.5
COL = {"PB": "#7f7f7f", "LCP": "#c98b2b", "LCPd": "#8c5fb3", "XDS": "#4f8fd1", "XDSM": "#0b2e5a", "KX": "#2a9d8f"}
LAB = {"PB": "Pseudobulk", "LCP": "Leak-corrected (LCP)", "LCPd": "LCPd", "XDS": "XDS, one split", "XDSM": "XDS-M", "KX": "KX"}
KEY = {"PB": "PB", "LCP": "LCP", "XDS": "LCPXd"}; SUBN = {"P1": "Homeo_CECR2", "P2": "PVM_CD163", "P3": "ADAM_GPNMB"}
# ---- geometry (inches) ----
Y_AX = 0.86; H_MAIN = 1.30; H_GAP = 0.05; H_STRIP = 0.17; H_B = H_MAIN + H_GAP + H_STRIP; H = Y_AX + H_B + 0.36
XA0 = 0.52; WA = 0.90; GA = 0.055                       # panel a: three sub-panels
XB1 = 3.98; WB1 = 0.80; XB2 = XB1 + WB1 + 0.07; WB2 = W - 0.04 - XB2
fig = plt.figure(figsize=(W, H))
def ax_in(x0, y0, w, h, **kw):   # axes from inch coordinates measured from the top-left corner of the figure
    return fig.add_axes([x0 / W, 1 - (y0 + h) / H, w / W, h / H], **kw)
def txt(x, y, s, **kw):          # figure text at inch coordinates from top-left
    return fig.text(x / W, 1 - y / H, s, **kw)
# ---------------- panel a ----------------
mains, strips = [], []
for i, s in enumerate(("P1", "P2", "P3")):
    x0 = XA0 + i * (WA + GA)
    am = ax_in(x0, Y_AX, WA, H_MAIN); az = ax_in(x0, Y_AX + H_MAIN + H_GAP, WA, H_STRIP)
    A = pop["A"][s]; Ds = sorted([int(k) for k in A if k != "all"]) + [info[s]["donors"]]; keys = [str(d) for d in Ds[:-1]] + ["all"]
    for m in ("PB", "LCP", "XDS"):   # LCP and XDS used 1e5 permutations at D <= 100 (too few for a BH call): plotted from 200 up only
        sel = [i_ for i_, d_ in enumerate(Ds) if m == "PB" or d_ >= 200]
        xs = np.array([Ds[i_] for i_ in sel], float); ys = np.array([A[keys[i_]][KEY[m]]["mean_calls"] for i_ in sel], float)
        pos = ys > 0
        am.plot(xs[pos], ys[pos], "-o", color=COL[m], ms=3.2, lw=1.1, zorder=2, clip_on=False)
        if (~pos).any():
            az.plot(xs[~pos], ys[~pos], "-o", color=COL[m], ms=3.2, lw=1.1, zorder=2, clip_on=False)
            for j in range(len(xs) - 1):   # dotted connector where the series passes between the log axis and the zero strip
                if pos[j] != pos[j + 1]:
                    (a_from, y_from), (a_to, y_to) = [((am, ys[q]) if pos[q] else (az, 0.0)) for q in (j, j + 1)]
                    fig.add_artist(ConnectionPatch((xs[j], y_from), (xs[j + 1], y_to), "data", "data", axesA=a_from, axesB=a_to, color=COL[m], lw=0.8, ls=(0, (1.5, 1.5)), zorder=1))
    for a_ in (am, az):
        a_.set_xscale("log"); a_.set_xlim(18, 1600); a_.xaxis.set_minor_locator(NullLocator())
    am.set_yscale("log"); am.set_ylim(0.2, 6e4); am.yaxis.set_minor_locator(NullLocator())
    am.set_yticks([1, 10, 100, 1000, 10000]); am.set_yticklabels(["1", "10", "100", "1k", "10k"] if i == 0 else [])
    am.set_xticks([]); am.spines["bottom"].set_visible(False)
    az.set_ylim(-1, 1); az.set_yticks([0]); az.set_yticklabels(["0"] if i == 0 else []); az.yaxis.set_minor_locator(NullLocator())
    az.set_xticks([25, 100, 400, 1000]); az.set_xticklabels(["25", "100", "400", "1k"])
    if i > 0: am.tick_params(axis="y", length=2.5); az.tick_params(axis="y", length=2.5)
    d = 0.035   # axis-break marks on the left spine between the log axis and the zero strip
    for a_, yy in ((am, 0.0), (az, 1.0)):
        a_.plot([-0.05, 0.05], [yy - 0.03, yy + 0.03], transform=a_.transAxes, color="black", lw=0.6, clip_on=False)
    txt(x0 + 0.02, Y_AX - 0.17, f"{SUBN[s]}", fontsize=7, ha="left", va="bottom"); txt(x0 + 0.02, Y_AX - 0.06, f"{info[s]['donors']:,} donors", fontsize=7, ha="left", va="bottom", color="0.35")
    mains.append(am); strips.append(az)
txt(XA0 + (3 * WA + 2 * GA) / 2, H - 0.02, "Pseudo-donors (pooled null, all calls false)", fontsize=8, ha="center", va="bottom")
fig.text(0.075 / W, 1 - (Y_AX + H_B / 2) / H, "Mean false calls (BH 5%)", rotation=90, fontsize=8, ha="center", va="center")
# ---------------- panel b ----------------
subs = ("P1", "P2", "P3")
def pooled(m, kind):
    """(percent, lo, hi, n) for method m: kind 'perm' = REV4 R1 permutation tests at D = 25, 50, 100; kind 'norm' = normal reference
    (REV3 A3 for PB, KX, XDS-M at 25-100; REV2 A at 200, 400, full for every method)."""
    k, n = [], []
    if kind == "perm":
        mm = "XDSperm" if m == "XDS" else m; sv = [s for s in subs if R1[s]["n"]]
        for D in ("25", "50", "100"):
            k.append(sum(R1[s][D][mm]["any"] for s in sv)); n.append(sum(R1[s][D][mm]["n"] for s in sv))
    else:
        if m in ("PB", "KX", "XDSM"):
            for D in ("25", "50", "100"):
                k.append(sum(A3[s][D][m]["any"] for s in subs)); n.append(sum(A3[s][D][m]["n"] for s in subs))
        for D in ("200", "400", "all"):
            k.append(sum(S2["A"][s][D][m]["any"] for s in subs)); n.append(sum(S2["A"][s][D][m]["n"] for s in subs))
    k, n = np.array(k), np.array(n); p = k / n
    lo = np.where(k == 0, 0, beta.ppf(0.025, k, n - k + 1)); hi = np.where(k == n, 1, beta.ppf(0.975, k + 1, n - k))
    return 100 * p, 100 * lo, 100 * hi, n, k
a1 = ax_in(XB1, Y_AX, WB1, H_B); a2 = ax_in(XB2, Y_AX, WB2, H_B)
XD = {"perm": [25, 50, 100], "norm": [25, 50, 100, 200, 400, 1000]}; plotted = []
def draw(ax, kind, ms_):
    for m in ms_:
        p, lo, hi, n, k = pooled(m, kind); xs = np.array(XD[kind][-len(p):], float); foc = m == "XDSM"
        ax.errorbar(xs, p, yerr=[p - lo, hi - p], fmt="-", color=COL[m], lw=1.9 if foc else 1.0, capsize=0, elinewidth=0.8, zorder=4 if foc else 2, clip_on=True)
        ax.plot(xs, p, "o", ms=3.8 if foc else 3.2, color=COL[m], zorder=5 if foc else 3, clip_on=False)
        for x_, p_, lo_, hi_, n_, k_ in zip(xs, p, lo, hi, n, k): plotted.append((kind, m, int(x_), int(k_), int(n_), p_, lo_, hi_))
draw(a1, "perm", ("LCP", "LCPd", "XDS")); draw(a2, "norm", ("PB", "LCP", "LCPd", "KX", "XDS", "XDSM"))
for ax in (a1, a2):
    ax.axhline(5, color="black", lw=0.6, ls=":", zorder=1); ax.set_yscale("log"); ax.set_xscale("log")
    ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator()); ax.set_ylim(0.4, 130)
a1.set_xlim(20, 130); a1.set_xticks([25, 50, 100]); a1.set_xticklabels(["25", "50", "100"])
a2.set_xlim(20, 1400); a2.set_xticks(XD["norm"]); a2.set_xticklabels(["25", "50", "100", "200", "400", "full"])
a1.set_yticks([1, 2, 5, 10, 20, 50, 100]); a1.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"])
a2.set_yticks([1, 2, 5, 10, 20, 50, 100]); a2.set_yticklabels([]); a2.tick_params(axis="y", length=2.5)
n1 = int(pooled("LCP", "perm")[3][0]); n2 = int(pooled("XDSM", "norm")[3][0])
txt(XB1 + 0.02, Y_AX - 0.17, "Permutation", fontsize=7, ha="left", va="bottom"); txt(XB1 + 0.02, Y_AX - 0.06, f"{n1} runs", fontsize=7, ha="left", va="bottom", color="0.35")
txt(XB2 + 0.02, Y_AX - 0.17, "Normal reference", fontsize=7, ha="left", va="bottom"); txt(XB2 + 0.02, Y_AX - 0.06, f"{n2} runs", fontsize=7, ha="left", va="bottom", color="0.35")
a2.text(21, 5.5, "5%", fontsize=7, ha="left", va="bottom")
txt(XB1 + (W - 0.04 - XB1) / 2 + 0.15, H - 0.02, "Pseudo-donors (donor-centred null)", fontsize=8, ha="center", va="bottom")
fig.text(3.60 / W, 1 - (Y_AX + H_B / 2) / H, "Runs with any false call (%)", rotation=90, fontsize=8, ha="center", va="center")
# ---------------- titles, per-panel legends ----------------
txt(0.02, 0.17, "a", fontsize=10, fontweight="bold", ha="left", va="baseline"); txt(0.20, 0.17, "Pseudobulk's false calls grow with the cohort", fontsize=8, ha="left", va="baseline")
txt(3.58, 0.17, "b", fontsize=10, fontweight="bold", ha="left", va="baseline"); txt(3.76, 0.17, "XDS-M: 3\u20138% of runs from 100 donors up", fontsize=8, ha="left", va="baseline")
def handles(ms_): return [mpl.lines.Line2D([], [], color=COL[m], marker="o", ms=3.2, lw=1.9 if m == "XDSM" else 1.0) for m in ms_]
la = fig.legend(handles(("PB", "LCP", "XDS")), ["Pseudobulk", "Leak-corrected (LCP)", "XDS, one split"], loc="upper left", bbox_to_anchor=(0.24 / W, 1 - 0.24 / H), bbox_transform=fig.transFigure, ncol=3, frameon=False,
                fontsize=7, handlelength=1.4, columnspacing=1.0, handletextpad=0.4, borderpad=0.1, borderaxespad=0)
bm = ("PB", "LCP", "LCPd", "KX", "XDS", "XDSM")
lb = fig.legend(handles(bm), ["Pseudobulk", "LCP", "LCPd", "KX", "XDS, one split", "XDS-M"], loc="upper left", bbox_to_anchor=(3.70 / W, 1 - 0.24 / H), bbox_transform=fig.transFigure, ncol=3, frameon=False,
                fontsize=7, handlelength=1.4, columnspacing=1.0, handletextpad=0.4, borderpad=0.1, borderaxespad=0)
out = os.path.join(RR, "paper2", "figures"); os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "fig_growth.pdf")); fig.savefig(os.path.join(out, "fig_growth.png"), dpi=300)
import csv
with open(os.path.join(out, "fig_growth_plotted.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["reference", "method", "pseudo_donors", "runs_with_any_false_call", "runs", "percent", "ci95_low_pct", "ci95_high_pct"]); w.writerows([(a, b, "full" if c == 1000 else c, d, e, round(p, 3), round(l, 3), round(h, 3)) for a, b, c, d, e, p, l, h in plotted])
fig.canvas.draw(); ov, outside = overlap_report(fig); print("overlaps", ov, "outside", outside)
audit_pdf(os.path.join(out, "fig_growth.pdf"))
print({(k_, m): [round(v, 1) for v in pooled(m, k_)[0]] for k_, ms_ in (("perm", ("LCP", "LCPd", "XDS")), ("norm", ("PB", "LCP", "LCPd", "KX", "XDS", "XDSM"))) for m in ms_})
