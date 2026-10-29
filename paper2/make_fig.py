"""Figure 2 (growth) from result JSONs. Requires apply_figure_style and panel_letter (figure-style helpers) in the namespace."""
import json, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, os
from matplotlib.ticker import NullLocator
R = os.environ.get("RV2", "/Users/kaavinprasanna/Documents/recomb27/recomb_v2")
def J(p):
    for c in (os.path.join(R, p), os.path.join(R, "results", p)):
        if os.path.exists(c): return json.load(open(c))
pop = J("pop/pop_results.json"); xu = J("xu/xu_results.json"); rs = J("rev/rev_sim_results.json")
info = {s: J(f"pop/pop_{s}_info.json") for s in ("P1", "P2", "P3")}
apply_figure_style(sizes=(8, 7, 6))
COL = {"PB": "#7f7f7f", "LCP": "#c98b2b", "XDS": "#1f5fa8", "XDSM": "#0b2e5a"}
KEY = {"PB": "PB", "LCP": "LCP", "XDS": "LCPXd"}; SUBN = {"P1": "Homeo_CECR2", "P2": "PVM_CD163", "P3": "ADAM_GPNMB"}
fig = plt.figure(figsize=(7.0, 2.6))
outer = fig.add_gridspec(1, 2, width_ratios=[3, 1.35], wspace=0.28, left=0.085, right=0.985, bottom=0.19, top=0.83)
inner = outer[0].subgridspec(1, 3, wspace=0.07); axs = [fig.add_subplot(inner[0, i]) for i in range(3)]
for i, s in enumerate(("P1", "P2", "P3")):
    ax = axs[i]; A = pop["A"][s]; Ds = sorted([int(k) for k in A if k != "all"]) + [info[s]["donors"]]
    keys = [str(d) for d in Ds[:-1]] + ["all"]
    for m in ("PB", "LCP", "XDS"):
        ax.plot(Ds, [A[k][KEY[m]]["mean_calls"] for k in keys], "-o", color=COL[m], ms=3.2, lw=1.7 if m == "XDS" else 1.1, zorder=3 if m == "XDS" else 2, clip_on=False)
    ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1); ax.set_ylim(-0.25, 6e4)
    ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xticks([25, 100, 400, 1000]); ax.set_xticklabels(["25", "100", "400", "1k"])
    ax.set_yticks([0, 1, 10, 100, 1000, 10000]); ax.set_yticklabels(["0", "1", "10", "100", "1k", "10k"] if i == 0 else [])
    ax.set_title(f"{SUBN[s]}, {info[s]['donors']:,}", loc="left", fontsize=7); ax.set_xlim(18, 1600)
axs[0].set_ylabel("False discoveries (BH 5%)", labelpad=2); axs[1].set_xlabel("Pseudo-donors in the real-cell null (every call false)")
axs[0].text(0.04, 0.93, "Pseudobulk", color=COL["PB"], transform=axs[0].transAxes, fontsize=7, va="top")
axs[0].text(0.60, 0.44, "Leak-\ncorrected", color=COL["LCP"], transform=axs[0].transAxes, fontsize=7, va="top")
axs[0].text(0.04, 0.12, "XDS", color=COL["XDS"], transform=axs[0].transAxes, fontsize=7, fontweight="bold")
ax = fig.add_subplot(outer[1]); E3, E4 = xu["E3"], xu["E4"]
small = [("F_D16", 16), ("F_D24", 24), ("F_D40", 40)]; large = [("L200", 200), ("L500", 500), ("L1000", 1000)]
for m in ("PB", "LCP", "XDS", "XDSM"):
    for grp, mk, fill in ((small, "o", False), (large, "s", True)):
        get = (lambda k: rs[k]["XDSM"]) if m == "XDSM" else (lambda k, g=grp, m=m: (E4 if g is small else E3)[k][m])
        x = [d for _, d in grp]; y = [get(k)["FDR"] for k, _ in grp]; e = [2 * get(k)["FDR_mcse"] for k, _ in grp]
        ax.errorbar(x, y, yerr=e, fmt=mk + ("--" if m == "XDSM" else "-"), color=COL[m], ms=3.2, lw=1.7 if m.startswith("XDS") else 1.1, capsize=0,
                    mfc=COL[m] if fill else "white", mec=COL[m], zorder=3 if m.startswith("XDS") else 2)
ax.axhline(0.05, color="black", lw=0.6, ls=":"); ax.text(60, 0.068, "0.05", fontsize=6)
ax.set_xscale("log"); ax.xaxis.set_minor_locator(NullLocator()); ax.set_xticks([16, 40, 200, 1000]); ax.set_xticklabels(["16", "40", "200", "1k"]); ax.set_xlim(12, 1400)
ax.set_ylim(0, 0.72); ax.set_yticks([0, 0.2, 0.4, 0.6]); ax.set_xlabel("Donors (simulation)"); ax.set_ylabel("False discovery rate", labelpad=2)
ax.set_title("Simulation with exact truth", loc="left", fontsize=7)
ax.text(27, 0.25, "XDS-M", color=COL["XDSM"], fontsize=7, fontweight="bold")
panel_letter(axs[0], "a"); panel_letter(ax, "b")
out = os.path.join(R, "paper2", "figures"); fig.savefig(os.path.join(out, "fig_growth.pdf")); fig.savefig(os.path.join(out, "fig_growth.png"), dpi=300)
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print("overlaps", [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)])
