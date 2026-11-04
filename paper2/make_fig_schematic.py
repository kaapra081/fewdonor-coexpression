"""Figure 1 (schematic): pseudobulk (a) versus XDS (b). A drawing, no data. Re-created because the generator of the earlier fig_schematic.pdf
was not in the repository (README: 'has no generator'); layout, wording and symbols follow that figure. Changes: exported at the design width
6.5 in (= \\linewidth) with every text span >= 7 pt after printing (text and math at 9 pt; subscripts 7.7 pt and the radical sign 7.3 pt: mathtext sub/superscript shrink 0.85),
and both panels now show the *same* cells (same shades), as the caption says.  Run from anywhere:  python paper2/make_fig_schematic.py
"""
import os, sys, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.path.abspath("repo/paper2")
sys.path.insert(0, HERE)
from figstyle import apply_figure_style, overlap_report, audit_pdf, enlarge_mathtext_subscripts
RR = os.path.dirname(HERE)
apply_figure_style(sizes=(9, 9, 7)); enlarge_mathtext_subscripts(0.85)
W, H = 6.5, 2.88
RED, BLUE, PURP, GREY = "#a83c2a", "#1f5fa8", "#7a3fa0", "#555555"
fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
rng = np.random.default_rng(3); cmap = mpl.colormaps["Oranges"]
ns = [14, 6, 10]; shade = [rng.uniform(0.15, 0.95, n) for n in ns]          # one set of cells, drawn in both panels
BOXW, BOXH, TOPS = 1.42, 0.50, [0.38, 0.97, 1.56]
def panel(x0, split):
    for r, (n, top) in enumerate(zip(ns, TOPS)):
        ax.add_patch(FancyBboxPatch((x0, top), BOXW, BOXH, boxstyle="round,pad=0.01,rounding_size=0.07", fc="#f4f4f4", ec="#bbbbbb", lw=0.6))
        ax.text(x0 + 0.07, top + 0.11, f"donor {r + 1}, {n} cells", fontsize=7, va="center", color=GREY)
        xs = np.linspace(x0 + 0.14, x0 + BOXW - 0.14, n)
        for i, (x, s_) in enumerate(zip(xs, shade[r])):
            if split:
                ec = BLUE if i % 2 == 0 else PURP; ax.add_patch(Circle((x, top + 0.355), 0.042, fc=cmap(s_), ec=ec, lw=1.1))
            else:
                ax.add_patch(Circle((x, top + 0.355), 0.042, fc=cmap(s_), ec=GREY, lw=0.6))
        xl = x0 + BOXW + 0.10
        if split:
            ax.text(xl, top + 0.12, r"$m^{A}_{dj}$", color=BLUE, fontsize=9, va="center"); ax.text(xl + 0.50, top + 0.12, r"$m^{B}_{dk}$", color=PURP, fontsize=9, va="center")
            ax.text(xl, top + 0.39, r"$\div\sqrt{v^{*}(n_d,\,\mathrm{depth})}$", fontsize=9, va="center")
        else:
            ax.text(xl, top + 0.12, r"$m_{dj},\ m_{dk}$", fontsize=9, va="center")
            ax.text(xl, top + 0.39, r"log or $\div\sqrt{\hat{\psi}_d}$", fontsize=9, va="center")
XA, XB = 0.22, 3.47
panel(XA, False); panel(XB, True)
ax.text(0.0, 0.15, "a", fontsize=10, fontweight="bold", va="center"); ax.text(XA, 0.15, "Pseudobulk: both genes from the same cells", fontsize=9, va="center")
ax.text(XB - 0.22, 0.15, "b", fontsize=10, fontweight="bold", va="center"); ax.text(XB, 0.15, "XDS: gene $j$ from half A, gene $k$ from half B", fontsize=9, va="center")
ax.text(XA, 2.18, "cell-state programs leak in\n(fraction $\\kappa_d\\approx1/n_d$); the transform's\nbias depends on $n_d$, shared by all genes", fontsize=9, color=RED, va="top", linespacing=1.25)
ax.text(XB, 2.18, "no shared cells: no leak; design-only scale:\nno $n_d$-dependent bias; test", fontsize=9, color=BLUE, va="top", linespacing=1.25)
ax.text(XB, 2.54, "$\\sum_d R^{A}_{dj}R^{B}_{dk}$ by donor permutation", fontsize=9, color=BLUE, va="top")
out = os.path.join(RR, "paper2", "figures"); os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "fig_schematic.pdf")); fig.savefig(os.path.join(out, "fig_schematic.png"), dpi=300)
fig.canvas.draw(); ov, outside = overlap_report(fig); print("overlaps", ov, "outside", outside); audit_pdf(os.path.join(out, "fig_schematic.pdf"))
