"""Figures for the paper. Inputs: dev_results/summary_all.csv (analyze_dev over runs_dev), real_summary/*.csv.
python make_figs.py"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V2 = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.join(V2, "paper", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "DejaVu Sans", "pdf.fonttype": 42})
COL = {"LCP": "#C0392B", "PB_logcpm": "#7F8C8D", "E5": "#A9CCE3", "CMB_naive": "#BDC3C7", "CMB_mom": "#2E86C1",
       "ORACLE": "#000000", "BOUND": "#E67E22"}
LAB = {"LCP": "LCP (this work)", "PB_logcpm": "Pseudobulk Pearson", "E5": "Standardized log means",
       "CMB_naive": "Count-scale means, uncorrected", "CMB_mom": "Moment-corrected + jackknife",
       "ORACLE": "Oracle: true donor means", "BOUND": "Noise-limited oracle (bound)"}
ORDER = ["ORACLE", "BOUND", "LCP", "PB_logcpm", "CMB_naive", "E5", "CMB_mom"]


def panel_letter(ax, s):
    ax.text(-0.18, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def fig_sim():
    s = pd.read_csv(os.path.join(V2, "results", "dev", "summary_all.csv"))
    designs = ["D8", "D12", "D16", "D24", "D40"]; xs = [8, 12, 16, 24, 40]
    s = s[s.design.isin(designs)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for m in ORDER:
        g = s[s.method == m].set_index("design").reindex(designs)
        lw = 2.0 if m == "LCP" else 1.0
        ls = "--" if m in ("ORACLE", "BOUND") else "-"
        kw = dict(color=COL[m], lw=lw, ls=ls, marker="o", ms=3 if m != "LCP" else 4,
                  zorder={"LCP": 3, "BOUND": 4}.get(m, 2))
        axes[0].errorbar(xs, g.FDR, yerr=2 * g.FDR_mcse, capsize=1.5, label=LAB[m], **kw)
        axes[1].errorbar(xs, g.power, yerr=2 * g.power_mcse, capsize=1.5, **kw)
        if m not in ("ORACLE", "BOUND"):
            axes[2].errorbar(xs, g.rej05_SS, yerr=2 * g.rej05_SS_mcse, capsize=1.5, **kw)
    axes[0].axhline(0.05, color="0.3", lw=0.6, ls=":")
    axes[0].set_ylim(-0.02, 0.92); axes[0].set_yticks([0, 0.2, 0.4, 0.6, 0.8]); axes[0].set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8"])
    axes[0].text(41.5, 0.05, "0.05", fontsize=6, va="center", ha="left", color="0.3")
    axes[0].set_ylabel("False discovery rate (BH, q = 0.05)")
    axes[0].set_title("LCP is calibrated from 24 donors", loc="left")
    axes[1].set_ylabel("Power (BH, q = 0.05)"); axes[1].set_title("LCP tracks the noise-limited bound", loc="left")
    axes[2].axhline(0.05, color="0.3", lw=0.6, ls=":")
    axes[2].set_ylabel("Type-I error at α = 0.05")
    axes[2].set_title("Leak: strong within-only pairs", loc="left")
    for ax in axes:
        ax.set_xscale("log"); ax.set_xticks(xs); ax.set_xticklabels([str(x) for x in xs]); ax.minorticks_off()
        ax.set_xlabel("Donors (≈40 cells per donor)"); ax.margins(x=0.06)
    for ax, l in zip(axes, "abc"):
        panel_letter(ax, l)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(os.path.join(OUT, "fig_sim.pdf")); fig.savefig(os.path.join(OUT, "fig_sim.png"), dpi=300)
    return fig


def fig_real():
    R = os.path.join(V2, "results", "real")
    sub = pd.read_csv(os.path.join(R, "subsample_vs_full.csv"))
    cap = pd.read_csv(os.path.join(R, "cap40_vs_full.csv"))
    goi = pd.read_csv(os.path.join(R, "genes_of_interest.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    mcol = {"LCP": COL["LCP"], "PB": COL["PB_logcpm"], "MC": COL["CMB_mom"], "naive": COL["CMB_naive"]}
    mlab = {"LCP": "LCP", "PB": "Pseudobulk", "MC": "Moment-corrected", "naive": "Uncorrected count means"}
    mk = {("seaad", 24): "o", ("seaad", 40): "s", ("lee", 24): "^", ("lee", 40): "D"}
    ax = axes[0]
    for m in ("MC", "PB", "LCP"):
        for (coh, D), mm in mk.items():
            g = sub[(sub.method == m) & (sub.cohort == coh) & (sub.D == D)]
            ax.scatter(g.recall_of_fullLCP, g.precision_vs_fullLCP, s=8, marker=mm, color=mcol[m], alpha=0.6, lw=0)
    for m in ("MC", "PB", "LCP"):
        ax.scatter([], [], color=mcol[m], s=10, label=mlab[m])
    ax.set_xlabel("Recall of full-cohort LCP calls"); ax.set_ylabel("Precision vs full-cohort LCP calls")
    ax.set_title("Donor subsamples (24 or 40 donors)", loc="left"); ax.legend(frameon=False, loc="lower left")
    ax.set_ylim(0.85, 1.005)
    ax = axes[1]
    meths = ["LCP", "MC", "PB", "naive"]
    for i, m in enumerate(meths):
        for j, coh in enumerate(["seaad", "lee"]):
            g = cap[(cap.method == m) & (cap.cohort == coh)]
            x = i + (j - 0.5) * 0.3
            ax.scatter(np.full(len(g), x), g.callrate_within_only, s=10, color=mcol[m], marker="o" if coh == "seaad" else "^", zorder=3)
            ax.hlines(g.callrate_within_only.mean(), x - 0.12, x + 0.12, color=mcol[m], lw=1.5)
    short = {"LCP": "LCP", "MC": "Moment-\ncorrected", "PB": "Pseudo-\nbulk", "naive": "Count\nmeans"}
    ax.set_xticks(range(len(meths))); ax.set_xticklabels([short[m] for m in meths], rotation=0)
    ax.set_ylabel("Call rate, within-dominant pairs")
    ax.set_title("Leak at 40 cells per donor", loc="left")
    ax.set_ylim(-0.005, None)
    ax = axes[2]
    for setname, c, lab in (("sex", "#8E44AD", "XIST / Y-linked / KDM6A"), ("heat_shock", "#16A085", "Heat-shock genes")):
        for coh, mm in (("seaad", "o"), ("lee", "^")):
            g = goi[(goi.set == setname) & (goi.cohort == coh)]
            ax.scatter(g.r_within, g.rhoB, s=10, color=c, marker=mm, alpha=0.8, lw=0, label=lab if coh == "seaad" else None)
    ax.axhline(0, color="0.6", lw=0.5); ax.axvline(0, color="0.6", lw=0.5)
    ax.set_xlabel("Within-donor correlation"); ax.set_ylabel("Between-donor correlation (corrected)")
    ax.set_title("Levels separate in both cohorts", loc="left")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], ls="", marker="o", color="#8E44AD", ms=4, label="XIST / Y-linked / KDM6A"),
                       Line2D([], [], ls="", marker="o", color="#16A085", ms=4, label="Heat-shock genes")],
              frameon=False, loc="center right")
    for ax, l in zip(axes, "abc"):
        panel_letter(ax, l); ax.margins(0.05)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_real.pdf")); fig.savefig(os.path.join(OUT, "fig_real.png"), dpi=300)
    return fig


if __name__ == "__main__":
    fig_sim(); fig_real()
