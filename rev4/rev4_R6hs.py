"""REV4 R6, POST HOC (labelled in paper): heat-shock and Y-linked gene pairs among null pseudobulk false calls (Homeo_CECR2)."""

import numpy as np, json, re
Z = np.load("rev4/out/R6_P1.npz", allow_pickle=True); sym = Z["sym"].astype(str); G = len(sym); a, b = np.triu_indices(G, 1)
c0 = Z["null_cPB"].astype(bool); r0 = Z["null_rPB"]; rw = Z["r_w"]
hs = np.array([bool(re.match(r"^(HSP|DNAJ|HSPH|HSPA|HSPB|HSPD|HSPE|BAG3|CHORDC1|STIP1|FKBP4|CACYBP|SERPINH1)", g)) for g in sym])
ylink = np.isin(sym, ["UTY","USP9Y","ZFY","DDX3Y","KDM5D","NLGN4Y","TTTY14","EIF1AY","RPS4Y1","XIST"])
pairHS = hs[a] & hs[b]; pairY = ylink[a] & ylink[b]
o = np.argsort(-np.abs(r0) * c0)[:25]
out = dict(n_hs_genes=int(hs.sum()), hs_genes=sym[hs].tolist(), n_calls=int(c0.sum()),
  hs_pairs=int(pairHS.sum()), hs_pairs_called=int((pairHS & c0).sum()), y_pairs=int(pairY.sum()), y_pairs_called=int((pairY & c0).sum()),
  frac_all_called=float(c0.mean()), median_r_hs_called=float(np.median(r0[pairHS & c0])) if (pairHS & c0).any() else None,
  top=[(sym[a[i]], sym[b[i]], round(float(r0[i]), 3), round(float(rw[i]), 4)) for i in o],
  median_rw_hs_pairs=float(np.median(rw[pairHS])), median_rw_all=float(np.median(rw)))
json.dump(out, open("rev4/out/R6/R6_hs.json", "w"), indent=0); print(json.dumps(out)[:200])
