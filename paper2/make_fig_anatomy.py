"""Figure 3 (anatomy of a pseudobulk network) from rev4/out/R6 tables (the coordinates uploaded to Mantis), the Mantis gene-map
clusters and rev4/out/A4/A5_summary.json.  Exported at 6.5 in (= \\linewidth), so printed text size == nominal size (>= 7 pt).
  row 1: a, b, c = donor-map thumbnails (c: real donors, true aspect ratio, colourbar underneath);
         d = first PC (left) and profile length (right) of the null profiles against cell number (Spearman rho from R6_summary.json);
  row 2: e = null false-edge degree by Mantis leaf cluster (clusters indexed 1..49 by rank of median); f = modularity z per edge set.
Outputs (paper2/figures): fig_anatomy.{pdf,png} = a-f;  fig_anatomy_top.{pdf,png} = a-d;  fig_anatomy_bottom.{pdf,png} = e-f
(the last two let e-f move to the appendix without touching the panel letters);  fig_anatomy_cluster_index.csv = index used on the x axis of e.
Every plotted number is read from the R6 / A5 result files; nothing is re-estimated.  Run from anywhere:  python paper2/make_fig_anatomy.py
"""
import os, sys, json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.path.abspath("repo/paper2")
sys.path.insert(0, HERE)
from figstyle import apply_figure_style, overlap_report, audit_pdf
RR = os.path.dirname(HERE)
R6 = os.path.join(RR, "rev4/out/R6"); dn = pd.read_csv(f"{R6}/donors_null.csv"); dr = pd.read_csv(f"{R6}/donors_real.csv")
gmap = pd.read_csv(os.path.join(RR, "rev4/mantis/gene_map_mantis.csv")); gn = pd.read_csv(f"{R6}/genes_numeric.csv")
g = gmap[["gene", "leaf_cluster_id", "leaf_cluster_label"]].merge(gn, on="gene")
mc = json.load(open(f"{R6}/R6_mantis_clusters.json")); HS = mc["top_null_cluster"]["label"]
R6S = json.load(open(f"{R6}/R6_summary.json")); A5 = json.load(open(os.path.join(RR, "rev4/out/A4/A5_summary.json")))
apply_figure_style(sizes=(8, 7, 7))
PBC, XC, HSC = "#7f7f7f", "#0b2e5a", "#d95f02"
W = 6.5; Y1 = 0.50; ROW1_END = 1.74; Y2 = 2.06; H2 = 1.42; H_FULL = Y2 + H2 + 0.40
norm = mpl.colors.Normalize(vmin=0.6, vmax=2.6); cmap = mpl.cm.viridis
lt = g.groupby("leaf_cluster_label").null_false_degree.median().sort_values(ascending=False, kind="stable"); order = list(lt.index); pos = {c: i for i, c in enumerate(order)}; nc = len(order)

def build(part):
    OFFS = (Y2 - 0.40) if part == "bottom" else 0.0
    H = {"full": H_FULL, "top": ROW1_END, "bottom": H_FULL - OFFS}[part]
    fig = plt.figure(figsize=(W, H))
    def ax_in(x0, y0, w, h, **kw): return fig.add_axes([x0 / W, 1 - (y0 - OFFS + h) / H, w / W, h / H], **kw)
    def txt(x, y, s, **kw): return fig.text(x / W, 1 - (y - OFFS) / H, s, **kw)
    def letter(x, y, s): txt(x, y, s, fontsize=10, fontweight="bold", ha="left", va="baseline")
    if part != "bottom":
        # ---------------- row 1, a-c: thumbnails
        W1, H1 = 1.02, 1.15; X_A, X_B, X_C = 0.20, 1.44, 2.68
        def emb(ax, df, x, y, title, tx):
            o = np.argsort(-df.log10_cells.values)   # few-cell donors drawn last
            ax.scatter(df[x].values[o], df[y].values[o], c=df.log10_cells.values[o], cmap=cmap, norm=norm, s=1.3, lw=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="datalim")
            xr = df[x].max() - df[x].min(); yr = df[y].max() - df[y].min(); xm = (df[x].max() + df[x].min()) / 2; ym = (df[y].max() + df[y].min()) / 2
            sc = max(xr / W1, yr / H1) * 1.04; ax.set_xlim(xm - sc * W1 / 2, xm + sc * W1 / 2); ax.set_ylim(ym - sc * H1 / 2, ym + sc * H1 / 2)
            for sp_ in ax.spines.values(): sp_.set_visible(False)
            txt(tx, 0.04, title, fontsize=7, ha="left", va="top", linespacing=1.15)
        emb(ax_in(X_A, Y1, W1, H1), dn, "x_L", "y_L", "Null pseudo-donors,\npseudobulk: ordered\nby cell count", X_A)
        emb(ax_in(X_B, Y1, W1, H1), dn, "x_R", "y_R", "Same pseudo-donors,\ndesign-standardized:\nno order", X_B)
        HC = 0.50; ax_c = ax_in(X_C, Y1, W1, HC); o = np.argsort(-dr.log10_cells.values)       # c: true aspect ratio, so two small clusters
        ax_c.scatter(dr.x_L.values[o], dr.y_L.values[o], c=dr.log10_cells.values[o], cmap=cmap, norm=norm, s=1.3, lw=0, rasterized=True)
        ax_c.set_xticks([]); ax_c.set_yticks([]); ax_c.set_aspect("equal", adjustable="datalim")
        xr = dr.x_L.max() - dr.x_L.min(); xm = (dr.x_L.max() + dr.x_L.min()) / 2; sc = xr / W1 * 1.04      # data units per inch
        ax_c.set_xlim(xm - sc * W1 / 2, xm + sc * W1 / 2); ax_c.set_ylim(dr.y_L.min() - 0.03 * sc, dr.y_L.min() - 0.03 * sc + sc * HC)
        for sp_ in ax_c.spines.values(): sp_.set_visible(False)
        for sx in ("female", "male"):
            sub = dr[dr.sex == sx]; ax_c.text(sub.x_L.median(), dr.y_L.max() + 0.02 * sc, sx, ha="center", va="bottom", fontsize=7, color="0.3")
        txt(X_C, 0.04, "Real donors,\npseudobulk: sex,\nthen cell count", fontsize=7, ha="left", va="top", linespacing=1.15)
        cax = ax_in(X_C + 0.02, Y1 + HC + 0.12, W1 - 0.04, 0.07); cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
        cb.set_ticks([1, 2]); cb.set_ticklabels(["10", "100"]); cb.ax.tick_params(labelsize=7, length=2, pad=1.5); cb.outline.set_linewidth(0.4)
        txt(X_C + W1 / 2, Y1 + HC + 0.12 + 0.07 + 0.22, "cells per donor", fontsize=7, ha="center", va="top")
        for x_, l_ in ((X_A, "a"), (X_B, "b"), (X_C, "c")): letter(x_ - 0.19, Y1 + 0.02, l_)
        # ---------------- row 1, d: first PC and profile length against cell number (null)
        XD1, WD, XD2, HD = 4.28, 0.80, 5.66, 0.85
        lx = np.log10(dn.cells.values); q = np.quantile(lx, np.linspace(0, 1, 11)); b = np.clip(np.searchsorted(q, lx, side="right") - 1, 0, 9)
        axd1 = ax_in(XD1, Y1, WD, HD); axd2 = ax_in(XD2, Y1, WD, HD)
        for ax_, cols, scal in ((axd1, ("pc1_L", "pc1_R"), lambda v: v / v.std()), (axd2, ("profile_norm_L", "profile_norm_R"), lambda v: v / np.median(v))):
            for col, zz in zip(cols, (PBC, XC)):
                y = scal(dn[col].values); ax_.scatter(lx, y, s=0.8, lw=0, color=zz, alpha=0.16, rasterized=True)
                ax_.plot([np.median(lx[b == i]) for i in range(10)], [np.median(y[b == i]) for i in range(10)], "-o", color=zz, ms=2.2, lw=1.1, zorder=3)
            ax_.set_xticks([1, 2]); ax_.set_xticklabels(["10", "100"]); ax_.set_xlim(0.45, 2.8)
        axd1.set_ylim(-2.5, 4.6); axd1.set_yticks([-2, 0, 2, 4]); axd1.set_yticklabels(["\u22122", "0", "2", "4"]); axd1.set_ylabel("First PC (SD)", labelpad=1)
        axd2.set_yscale("log"); axd2.set_ylim(0.3, 6); axd2.set_yticks([0.5, 1, 2, 4]); axd2.set_yticklabels(["0.5", "1", "2", "4"]); axd2.yaxis.set_minor_locator(mpl.ticker.NullLocator())
        axd2.set_ylabel("Profile length", labelpad=1)
        rL = R6S["null_L"]["spearman_norm_ncells"]; rR = R6S["null_R"]["spearman_norm_ncells"]; mn = lambda v: f"\u03c1 = {v:.2f}".replace("-", "\u2212")
        axd2.text(0.96, 0.98, mn(rL), transform=axd2.transAxes, fontsize=7, color=PBC, ha="right", va="top", fontweight="bold")
        axd2.text(0.96, 0.56, mn(rR), transform=axd2.transAxes, fontsize=7, color=XC, ha="right", va="bottom", fontweight="bold")
        txt((XD1 + XD2 + WD) / 2 + 0.03, Y1 + HD + 0.19, "Cells per pseudo-donor", fontsize=8, ha="center", va="top")
        txt(XD1 - 0.40, 0.04, "Null: PC1 and profile length track\ncell count for pseudobulk only", fontsize=7, ha="left", va="top", linespacing=1.15)
        hh = [mpl.lines.Line2D([], [], color=c_, marker="o", ms=2.6, lw=1.1) for c_ in (PBC, XC)]
        fig.legend(hh, ["pseudobulk", "design-standardized"], loc="upper left", bbox_to_anchor=((XD1 - 0.40) / W, 1 - (0.33 - OFFS) / H), bbox_transform=fig.transFigure, ncol=2, frameon=False,
                   fontsize=7, handlelength=1.3, handletextpad=0.35, columnspacing=0.9, borderpad=0.0, borderaxespad=0)
        letter(XD1 - 0.40, Y1 + 0.02, "d")
    if part != "top":
        # ---------------- row 2, e: null false-edge degree by Mantis leaf cluster, indexed by rank of the cluster median (1 = highest)
        axe = ax_in(0.64, Y2, 2.52, H2); rng_ = np.random.default_rng(0)
        for c_ in order:
            v = g[g.leaf_cluster_label == c_].null_false_degree.values; x0 = pos[c_]; hl = c_ == HS
            axe.scatter(x0 + rng_.uniform(-0.25, 0.25, len(v)), v, s=4 if hl else 2.2, lw=0, color=HSC if hl else "0.55", alpha=0.9 if hl else 0.6, zorder=3 if hl else 2)
            axe.plot([x0 - 0.38, x0 + 0.38], [np.median(v)] * 2, color=HSC if hl else "0.15", lw=1.3 if hl else 0.8, zorder=4)
        gmed = g.null_false_degree.median()
        axe.axhline(gmed, color="0.5", lw=0.6, ls=":"); axe.text(nc - 0.6, gmed + 4, "all genes", fontsize=7, ha="right", va="bottom", color="0.35", zorder=5, bbox=dict(fc="white", ec="none", pad=0.4))
        axe.set_xlim(-1.0, nc - 0.2); axe.set_ylim(-5, g.null_false_degree.max() * 1.12)
        axe.set_xticks([0, 9, 19, 29, 39, nc - 1]); axe.set_xticklabels(["1", "10", "20", "30", "40", str(nc)]); axe.set_yticks([0, 50, 100, 150])
        axe.set_xlabel(f"Gene cluster, ranked by median ({nc} Mantis clusters)", labelpad=2); axe.set_ylabel("Null false-edge degree", labelpad=1)
        axe.annotate("1  " + HS.replace("Chaperones And ", "Chaperones\n    And "), xy=(pos[HS] + 0.5, lt[HS]), xytext=(4.0, g.null_false_degree.max() * 1.02), fontsize=7, color=HSC, ha="left", va="center",
                     arrowprops=dict(arrowstyle="-", color=HSC, lw=0.6, shrinkA=1, shrinkB=1), linespacing=1.1)
        txt(0.64, Y2 - 0.12, "Homeo_CECR2 null: false-edge degree\nby function cluster", fontsize=7, ha="left", va="bottom", linespacing=1.15)
        # ---------------- row 2, f: modularity of the Mantis function partition on each edge set (REV4 Amendments 4-5; expression-stratified relabelling)
        ES = [("null20", "Null pseudobulk\n(20 runs)"), ("top5000_null_PB", "Null pseudobulk\ntop 5,000 (run 1)"), ("top5000_real_PB", "Real pseudobulk\ntop 5,000"),
              ("top5000_real_XDSM", "Real XDS-M\ntop 5,000"), ("PB_only", "Pseudobulk\nonly"), ("XDSM_only", "XDS-M\nonly")]
        MK = {"P1": ("o", "Homeo_CECR2"), "P2": ("s", "PVM_CD163"), "P3": ("^", "ADAM_GPNMB")}
        ECOL = {"null20": PBC, "top5000_null_PB": PBC, "top5000_real_PB": PBC, "top5000_real_XDSM": XC, "PB_only": PBC, "XDSM_only": XC}
        axf = ax_in(4.78, Y2, 1.68, H2)
        for j_, (es, lab_) in enumerate(ES):
            for k_, (sx, (mk, nm)) in enumerate(MK.items()):
                yy = j_ + (k_ - 1) * 0.24; filled = es in ("null20", "top5000_null_PB")
                if es == "null20":
                    z = A5[sx]["null20"]["z_mean"]; sd = A5[sx]["null20"]["z_sd"]; axf.plot([z - sd, z + sd], [yy, yy], color=ECOL[es], lw=0.8, zorder=2)
                else: z = A5[sx]["leaf"][es]["Q_z"]
                axf.scatter(z, yy, marker=mk, s=13, color=ECOL[es] if filled else "white", edgecolor=ECOL[es], lw=0.8, zorder=3)
        axf.axvline(0, color="0.5", lw=0.6); axf.axvline(2.33, color="0.5", lw=0.6, ls=":")
        axf.set_yticks(range(len(ES))); axf.set_yticklabels([l for _, l in ES], fontsize=7, linespacing=1.05); axf.set_xlim(-3, 6.8); axf.set_ylim(len(ES) - 0.45, -0.6)
        axf.set_xticks([-2, 0, 2, 4, 6]); axf.set_xticklabels(["\u22122", "0", "2", "4", "6"]); axf.set_xlabel("Function modularity (z)", labelpad=2)
        hh = [mpl.lines.Line2D([], [], marker=m_, ls="", mfc="white", mec="0.3", ms=3.6, mew=0.8) for m_, _ in MK.values()]
        axf.legend(hh, [nm for _, nm in MK.values()], frameon=False, fontsize=7, loc="lower right", handletextpad=0.2, borderaxespad=0.0, borderpad=0.1, labelspacing=0.25, ncol=1)
        txt(3.52, Y2 - 0.12, "False edges at least as coherent\nas real-donor edges", fontsize=7, ha="left", va="bottom", linespacing=1.15)
        letter(0.05, Y2 - 0.12, "e"); letter(3.36, Y2 - 0.12, "f")
    return fig

out = os.path.join(RR, "paper2", "figures"); os.makedirs(out, exist_ok=True)
for part, stem in (("full", "fig_anatomy"), ("top", "fig_anatomy_top"), ("bottom", "fig_anatomy_bottom")):
    fig = build(part); fig.savefig(os.path.join(out, stem + ".pdf")); fig.savefig(os.path.join(out, stem + ".png"), dpi=300)
    fig.canvas.draw(); ov, outside = overlap_report(fig); print(stem, "overlaps", ov, "outside", outside)
    audit_pdf(os.path.join(out, stem + ".pdf")); plt.close(fig)
ci = pd.DataFrame({"index": np.arange(1, nc + 1), "leaf_cluster_label": order}).merge(g.groupby("leaf_cluster_label").agg(n_genes=("gene", "size"), median_null_false_degree=("null_false_degree", "median")).reset_index(), on="leaf_cluster_label")
ci.to_csv(os.path.join(out, "fig_anatomy_cluster_index.csv"), index=False)
