"""REV4 Amendment 1, R6c: edge degrees across Mantis's automatic clusters of the gene map (function-text layout).
Inputs: rev4/mantis/gene_map_mantis.csv, gene_clusters.csv (read back from Mantis), rev4/out/R6/genes_numeric.csv. Run from repo/.
Pre-specified: Kruskal-Wallis of PB-only and XDS-M degree across clusters; clusters with the highest median PB-only degree.
Added (labelled): the same test for null false-edge degree and PB degree. Partition: Mantis leaf clusters (finest level; all >= 5 genes)."""
import json, pandas as pd
from scipy.stats import kruskal
gmap = pd.read_csv("rev4/mantis/gene_map_mantis.csv"); gcl = pd.read_csv("rev4/mantis/gene_clusters.csv"); gn = pd.read_csv("rev4/out/R6/genes_numeric.csv")
g = gmap[["gene", "leaf_cluster_id", "leaf_cluster_label"]].merge(gn, on="gene"); assert len(g) == len(gn) == 424
grp = [x for _, x in g.groupby("leaf_cluster_id")]
kw = {m: dict(zip(("H", "p"), map(float, kruskal(*[x[m].values for x in grp])))) for m in ("PB_only_degree", "XDSM_degree", "null_false_degree", "PB_degree")}
lt = g.groupby(["leaf_cluster_id", "leaf_cluster_label"]).agg(n=("gene", "size"), null_med=("null_false_degree", "median"), pbonly_med=("PB_only_degree", "median"), xdsm_med=("XDSM_degree", "median")).reset_index()
tp = lt.sort_values("pbonly_med", ascending=False).iloc[0]; tn = lt.sort_values("null_med", ascending=False).iloc[0]
out = dict(n_leaves=int(len(lt)), n_clusters_total=int(len(gcl)), n_levels=int(gcl.depth.max()), kw=kw,
           top_PBonly_cluster=dict(label=tp.leaf_cluster_label, n=int(tp.n), median=float(tp.pbonly_med), xdsm_median=float(tp.xdsm_med)),
           top_null_cluster=dict(label=tn.leaf_cluster_label, n=int(tn.n), median=float(tn.null_med)),
           gene_medians={m: float(g[m].median()) for m in ("null_false_degree", "PB_only_degree", "XDSM_degree")},
           hs_cluster_genes=sorted(g[g.leaf_cluster_label == tn.leaf_cluster_label].gene.tolist()))
lt.to_csv("rev4/out/R6/R6_mantis_cluster_table.csv", index=False); json.dump(out, open("rev4/out/R6/R6_mantis_clusters.json", "w"), indent=1); print(json.dumps(out)[:400])
