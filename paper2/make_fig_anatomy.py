"""Figure 3 (anatomy of a pseudobulk network) from rev4/out/R6 tables (the same coordinates uploaded to Mantis).
Needs apply_figure_style, panel_letter (figure-style skill). Run from repo/."""
import os, json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import spearmanr
R6 = "rev4/out/R6"; dn = pd.read_csv(f"{R6}/donors_null.csv"); dr = pd.read_csv(f"{R6}/donors_real.csv")
gm_path = "rev4/mantis/gene_map_mantis.csv"; GENE = os.path.exists(gm_path)
apply_figure_style(sizes=(8, 7, 6))
fig = plt.figure(figsize=(7.0, 2.35)); gs = fig.add_gridspec(1, 3, wspace=0.08, left=0.01, right=0.625, bottom=0.04, top=0.80)
gs2 = fig.add_gridspec(1, 1, left=0.775, right=0.99, bottom=0.17, top=0.80)
axs = [fig.add_subplot(gs[0, i]) for i in range(3)] + [fig.add_subplot(gs2[0, 0])]
norm = mpl.colors.Normalize(vmin=0.6, vmax=2.6); cmap = mpl.cm.viridis
def emb(ax, df, x, y, title):
    o = np.argsort(-df.log10_cells.values)   # few-cell donors drawn last
    ax.scatter(df[x].values[o], df[y].values[o], c=df.log10_cells.values[o], cmap=cmap, norm=norm, s=2.2, lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.04); ax.set_title(title, loc="left", fontsize=7)
    for sp_ in ax.spines.values(): sp_.set_visible(False)
emb(axs[0], dn, "x_L", "y_L", "Null: pseudobulk\nsorts by cell count")
emb(axs[1], dn, "x_R", "y_R", "Null: design-standardized\nno order")
emb(axs[2], dr, "x_L", "y_L", "Real donors: pseudobulk\nsex, then cell count")
# label the sex clusters in (c)
for sx, lab in (("female", "female"), ("male", "male")):
    sub = dr[dr.sex == sx]; axs[2].text(sub.x_L.median(), sub.y_L.max() + 0.35, lab, ha="center", va="bottom", fontsize=6, color="0.25")
axs[2].set_ylim(None, dr.y_L.max() + 1.6)
ax = axs[3]
if GENE:
    g = pd.read_csv(gm_path)
    v = g["null_false_degree"].values.astype(float); o = np.argsort(v)
    sc = ax.scatter(g.mantis_x.values[o], g.mantis_y.values[o], c=v[o], cmap="magma_r", s=4, lw=0)
    ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.04); ax.set_title("Genes by function (Mantis)\ncolour: false-edge degree", loc="left", fontsize=7)
else:
    lx = np.log10(dn.cells.values); q = np.quantile(lx, np.linspace(0, 1, 11)); b = np.clip(np.searchsorted(q, lx, side="right") - 1, 0, 9)
    for col, lab, zz in (("pc1_L", "pseudobulk", "#7f7f7f"), ("pc1_R", "design-standardized", "#0b2e5a")):
        y = dn[col].values / dn[col].std(); ax.scatter(lx, y, s=1.2, lw=0, color=zz, alpha=0.18, rasterized=True)
        xm = [np.median(lx[b == i]) for i in range(10)]; ym = [np.median(y[b == i]) for i in range(10)]
        ax.plot(xm, ym, "-o", color=zz, ms=2.6, lw=1.3 if lab != "pseudobulk" else 1.1, label=lab, zorder=3)
    ax.set_xlabel("Cells per pseudo-donor", labelpad=1); ax.set_xticks([1, 2]); ax.set_xticklabels(["10", "100"])
    ax.set_ylabel("First PC score (SD)", labelpad=1); ax.set_ylim(-2.5, 4.5); ax.margins(0.04)
    ax.set_title("Null: first PC is cell count\nfor pseudobulk only", loc="left", fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="upper right", handlelength=1.4, handletextpad=0.4, borderaxespad=0.1)
cax = fig.add_axes([0.635, 0.12, 0.009, 0.52]); cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_ticks([1, 2]); cb.set_ticklabels(["10", "100"]); cb.set_label("cells per donor", fontsize=6, labelpad=1); cb.ax.tick_params(labelsize=6)
for a_, l_ in zip(axs, "abcd"): panel_letter(a_, l_)
os.makedirs("paper2/figures", exist_ok=True); fig.savefig("paper2/figures/fig_anatomy.pdf"); fig.savefig("paper2/figures/fig_anatomy.png", dpi=300)
