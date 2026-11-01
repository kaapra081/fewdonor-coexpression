"""Figure 2 from result JSONs: (a) mean false calls vs pseudo-donors (pooled null, POP_PREREG); (b) fraction of runs with >= 1 false call vs
pseudo-donors (donor-centred null, REV3 A3 at 25-100 and REV2 A at 200-full, pooled over three subtypes). Needs apply_figure_style, panel_letter."""
import json, os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from scipy.stats import beta
RR = os.path.dirname(os.path.dirname(os.path.abspath(__file__ if "__file__" in globals() else "repo/paper2/make_fig.py")))
R = os.environ.get("RV2", "/Users/kaavinprasanna/Documents/recomb27/recomb_v2")
def J(p):
    for c in (os.path.join(RR, p), os.path.join(RR, "results", p), os.path.join(R, p), os.path.join(R, "results", p)):
        if os.path.exists(c): return json.load(open(c))
    raise FileNotFoundError(p)
pop = J("pop/pop_results.json"); info = {s: J(f"pop/pop_{s}_info.json") for s in ("P1", "P2", "P3")}
S2 = J("rev2/interim/REV2_SUMMARY_latest.json"); A3 = J("rev3/results/A3B3_SUMMARY.json")["A3"]
apply_figure_style(sizes=(8, 7, 6))
COL = {"PB": "#7f7f7f", "LCP": "#c98b2b", "LCPd": "#8c5fb3", "XDS": "#4f8fd1", "XDSM": "#0b2e5a"}
LAB = {"PB": "Pseudobulk", "LCP": "Leak-corrected (LCP)", "LCPd": "LCPd", "XDS": "XDS, one split", "XDSM": "XDS-M"}
KEY = {"PB": "PB", "LCP": "LCP", "XDS": "LCPXd"}; SUBN = {"P1": "Homeo_CECR2", "P2": "PVM_CD163", "P3": "ADAM_GPNMB"}
fig = plt.figure(figsize=(7.0, 2.7))
outer = fig.add_gridspec(1, 2, width_ratios=[3, 1.45], wspace=0.30, left=0.085, right=0.985, bottom=0.19, top=0.76)
inner = outer[0].subgridspec(1, 3, wspace=0.07); axs = [fig.add_subplot(inner[0, i]) for i in range(3)]
for i, s in enumerate(("P1", "P2", "P3")):
    ax = axs[i]; A = pop["A"][s]; Ds = sorted([int(k) for k in A if k != "all"]) + [info[s]["donors"]]
    keys = [str(d) for d in Ds[:-1]] + ["all"]
    for m in ("PB", "LCP", "XDS"):
        ax.plot(Ds, [A[k][KEY[m]]["mean_calls"] for k in keys], "-o", color=COL[m], ms=3.2, lw=1.1, zorder=2, clip_on=False)
    ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1); ax.set_ylim(-0.25, 6e4)
    ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xticks([25, 100, 400, 1000]); ax.set_xticklabels(["25", "100", "400", "1k"])
    ax.set_yticks([0, 1, 10, 100, 1000, 10000]); ax.set_yticklabels(["0", "1", "10", "100", "1k", "10k"] if i == 0 else [])
    ax.set_title(f"{SUBN[s]} ({info[s]['donors']:,})", loc="left", fontsize=7); ax.set_xlim(18, 1600)
axs[0].set_ylabel("Mean false calls (BH 5%)", labelpad=2); axs[1].set_xlabel("Pseudo-donors (pooled null; every call is false)")
fig.text(0.085, 0.865, "Pseudobulk's false calls grow with the cohort", fontsize=8, ha="left")
# (b) run-level
ax = fig.add_subplot(outer[1]); subs = ("P1", "P2", "P3"); xs = [25, 50, 100, 200, 400, 1000]
def pooled(m):
    k, n = [], []
    for D in ("25", "50", "100"):
        mm = "XDSperm" if m == "XDS" else m
        k.append(sum(A3[s][D][mm]["any"] for s in subs)); n.append(sum(A3[s][D][mm]["n"] for s in subs))
    for D in ("200", "400", "all"):
        k.append(sum(S2["A"][s][D][m]["any"] for s in subs)); n.append(sum(S2["A"][s][D][m]["n"] for s in subs))
    k, n = np.array(k), np.array(n); p = k / n
    lo = np.where(k == 0, 0, beta.ppf(0.025, k, n - k + 1)); hi = np.where(k == n, 1, beta.ppf(0.975, k + 1, n - k))
    return 100 * p, 100 * lo, 100 * hi, n
for m in ("PB", "LCP", "LCPd", "XDS", "XDSM"):
    p, lo, hi, n = pooled(m); foc = m == "XDSM"
    ax.errorbar(xs, p, yerr=[p - lo, hi - p], fmt="-", color=COL[m], lw=1.9 if foc else 1.0, capsize=0, elinewidth=0.8, zorder=4 if foc else 2, clip_on=False)
    for j, (x_, y_) in enumerate(zip(xs, p)):
        openm = (m == "XDS" and j < 3)   # XDS with donor-permutation p-values at D <= 100 (normal-moment reference above)
        ax.plot([x_], [y_], "o", ms=3.6 if foc else 3.0, color=COL[m], mfc="white" if openm else COL[m], mec=COL[m], zorder=5 if foc else 3, clip_on=False)
ax.axhline(5, color="black", lw=0.6, ls=":"); ax.text(1300, 5.6, "5%", fontsize=6, ha="right", va="bottom")
ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1); ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator())
ax.set_xticks(xs); ax.set_xticklabels(["25", "50", "100", "200", "400", "full"]); ax.set_xlim(20, 1300); ax.set_ylim(-0.1, 130)
ax.set_yticks([0, 1, 5, 10, 50, 100]); ax.set_yticklabels(["0", "1", "5", "10", "50", "100"])
ax.set_xlabel("Pseudo-donors (donor-centred null)"); ax.set_ylabel("Runs with any false call (%)", labelpad=2)
fig.text(0.70, 0.865, "XDS-M rarely errs from 100 donors up", fontsize=8, ha="left")
hs = [mpl.lines.Line2D([], [], color=COL[m], marker="o", ms=3.2, lw=1.9 if m == "XDSM" else 1.0) for m in ("PB", "LCP", "LCPd", "XDS", "XDSM")]
fig.legend(hs, [LAB[m] for m in ("PB", "LCP", "LCPd", "XDS", "XDSM")], loc="upper left", bbox_to_anchor=(0.085, 1.0), ncol=5, frameon=False, fontsize=7,
           handlelength=1.6, columnspacing=1.4, handletextpad=0.5)
panel_letter(axs[0], "a"); panel_letter(ax, "b")
out = os.path.join(RR, "paper2", "figures"); fig.savefig(os.path.join(out, "fig_growth.pdf")); fig.savefig(os.path.join(out, "fig_growth.png"), dpi=300)
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print("overlaps", [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)])
print({m: [round(v, 1) for v in pooled(m)[0]] for m in ("PB", "LCP", "LCPd", "XDS", "XDSM")})
