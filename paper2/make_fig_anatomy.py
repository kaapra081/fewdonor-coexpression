"""Figure 3 (anatomy of a pseudobulk network) from rev4/out/R6 tables (the coordinates uploaded to Mantis) and the Mantis gene-map clusters.
Needs apply_figure_style, panel_letter (figure-style skill). Run from repo/."""
import os, json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
R6 = "rev4/out/R6"; dn = pd.read_csv(f"{R6}/donors_null.csv"); dr = pd.read_csv(f"{R6}/donors_real.csv")
gmap = pd.read_csv("rev4/mantis/gene_map_mantis.csv"); gn = pd.read_csv(f"{R6}/genes_numeric.csv"); g = gmap[["gene", "leaf_cluster_id", "leaf_cluster_label"]].merge(gn, on="gene")
mc = json.load(open(f"{R6}/R6_mantis_clusters.json")); HS = mc["top_null_cluster"]["label"]
apply_figure_style(sizes=(8, 7, 6))
PBC, XC, HSC = "#7f7f7f", "#0b2e5a", "#d95f02"
fig = plt.figure(figsize=(7.0, 4.15))
top = fig.add_gridspec(1, 3, wspace=0.06, left=0.03, right=0.86, bottom=0.55, top=0.9)
axs = [fig.add_subplot(top[0, i]) for i in range(3)]
bot = fig.add_gridspec(1, 2, width_ratios=[1, 2.25], wspace=0.22, left=0.075, right=0.99, bottom=0.075, top=0.43)
axd = fig.add_subplot(bot[0, 0]); axe = fig.add_subplot(bot[0, 1])
norm = mpl.colors.Normalize(vmin=0.6, vmax=2.6); cmap = mpl.cm.viridis
def emb(ax, df, x, y, title):
    o = np.argsort(-df.log10_cells.values)   # few-cell donors drawn last
    ax.scatter(df[x].values[o], df[y].values[o], c=df.log10_cells.values[o], cmap=cmap, norm=norm, s=2.0, lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.04); ax.set_title(title, loc="left", fontsize=7)
    for sp_ in ax.spines.values(): sp_.set_visible(False)
emb(axs[0], dn, "x_L", "y_L", "Null pseudo-donors, pseudobulk:\nordered by cell count")
emb(axs[1], dn, "x_R", "y_R", "Same pseudo-donors,\ndesign-standardized: no order")
emb(axs[2], dr, "x_L", "y_L", "Real donors, pseudobulk:\nsex, then cell count")
for sx in ("female", "male"):
    sub = dr[dr.sex == sx]; axs[2].text(sub.x_L.median(), sub.y_L.max() + 0.3, sx, ha="center", va="bottom", fontsize=6, color="0.3")
axs[2].set_ylim(None, dr.y_L.max() + 1.5)
cax = fig.add_axes([0.875, 0.60, 0.011, 0.26]); cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_ticks([1, 2]); cb.set_ticklabels(["10", "100"]); cb.set_label("cells per donor", fontsize=6, labelpad=2); cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.4)
# (d) first PC vs cell number, null
lx = np.log10(dn.cells.values); q = np.quantile(lx, np.linspace(0, 1, 11)); b = np.clip(np.searchsorted(q, lx, side="right") - 1, 0, 9)
for col, lab, zz in (("pc1_L", "pseudobulk", PBC), ("pc1_R", "design-standardized", XC)):
    y = dn[col].values / dn[col].std(); axd.scatter(lx, y, s=1.0, lw=0, color=zz, alpha=0.16, rasterized=True)
    axd.plot([np.median(lx[b == i]) for i in range(10)], [np.median(y[b == i]) for i in range(10)], "-o", color=zz, ms=2.4, lw=1.2, label=lab, zorder=3)
axd.set_xlabel("Cells per pseudo-donor", labelpad=3); axd.set_xticks([1, 2]); axd.set_xticklabels(["10", "100"]); axd.set_ylim(-2.5, 4.6)
axd.set_ylabel("First PC score (SD)", labelpad=1); axd.set_title("Null: first PC tracks cell count\nfor pseudobulk only", loc="left", fontsize=7)
axd.legend(frameon=False, fontsize=6, loc="upper right", handlelength=1.3, handletextpad=0.4, borderaxespad=0.0)
# (e) null false-edge degree by Mantis leaf cluster
lt = g.groupby("leaf_cluster_label").null_false_degree.median().sort_values(ascending=False); order = list(lt.index); pos = {c: i for i, c in enumerate(order)}
rng_ = np.random.default_rng(0)
for c_ in order:
    v = g[g.leaf_cluster_label == c_].null_false_degree.values; x0 = pos[c_]; hl = c_ == HS
    axe.scatter(x0 + rng_.uniform(-0.25, 0.25, len(v)), v, s=4 if hl else 2.2, lw=0, color=HSC if hl else "0.55", alpha=0.9 if hl else 0.6, zorder=3 if hl else 2)
    axe.plot([x0 - 0.38, x0 + 0.38], [np.median(v)] * 2, color=HSC if hl else "0.15", lw=1.3 if hl else 0.8, zorder=4)
axe.axhline(g.null_false_degree.median(), color="0.5", lw=0.6, ls=":"); axe.text(len(order) - 0.5, g.null_false_degree.median() - 4, "median, all genes", fontsize=6, ha="right", va="top", color="0.4", zorder=5, bbox=dict(fc="white", ec="none", pad=0.5))
axe.set_xlim(-0.8, len(order) - 0.2); axe.set_xticks([]); axe.set_xlabel(f"{len(order)} Mantis gene clusters (from function text only)", labelpad=3)
axe.set_ylabel("Null false-edge degree", labelpad=1); axe.set_ylim(-5, g.null_false_degree.max() * 1.12)
axe.annotate(HS.replace(" And ", " and "), xy=(pos[HS] + 0.3, lt[HS]), xytext=(7.0, g.null_false_degree.max() * 1.03), fontsize=6, color=HSC, ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=HSC, lw=0.6))
axe.set_title("Null false edges by gene function: highest in the heat-shock cluster", loc="left", fontsize=7)
for a_, l_ in zip(axs + [axd, axe], "abcde"): panel_letter(a_, l_)
os.makedirs("paper2/figures", exist_ok=True); fig.savefig("paper2/figures/fig_anatomy.pdf"); fig.savefig("paper2/figures/fig_anatomy.png", dpi=300)
