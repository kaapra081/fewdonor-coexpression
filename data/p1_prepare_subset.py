#!/usr/bin/env python
"""
P1 preparation: extract a donor x cell subset of a large CELLxGENE h5ad without loading it.

Two chunked passes over raw/X (CSR) with h5py:
  pass 1: library size (all genes) and per-gene detection on the selected cells
  pass 2: counts for the selected cells restricted to the fixed gene universe
Writes a compact h5ad (cells x universe genes, raw counts in X and .raw) with obs columns
donor_id, lib_size (computed over ALL genes), the condition column, and a nested
within-donor random rank `sub_rank` so that smaller depths are subsets of larger ones.

Pre-registered P1 choices (pilot_P0_P1_preregistration.md): microglia only (drop subclass == PVM),
donors with >= 1600 microglia, 80 donors drawn at random stratified by disease, seed 20261022,
1600 cells per donor, universe = top 2000 genes by detection among those cells (>= 10%).
"""
import argparse, json, os, sys
import numpy as np, pandas as pd, h5py
from scipy import sparse

SEX_GENES = ["XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY"]  # as in p0_donor_null.py


def cat(h, col):
    g = h["obs"][col]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
        return cats[g["codes"][:]]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]])


def var_symbols(h, col="feature_name"):
    v = h["var"]
    if col in v:
        g = v[col]
        if isinstance(g, h5py.Group) and "categories" in g:
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
            return cats[g["codes"][:]]
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]])
    idx = v.attrs.get("_index", "_index")
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v[idx][:]])


def iter_rows(Xg, rows, chunk=5000):
    """Yield (row_ids, csr_block) for the sorted selected rows, reading contiguous row ranges."""
    indptr = Xg["indptr"]
    rows = np.sort(rows)
    for s in range(0, len(rows), chunk):
        r = rows[s:s + chunk]
        lo, hi = int(r[0]), int(r[-1]) + 1
        ip = indptr[lo:hi + 1].astype(np.int64)
        data = Xg["data"][ip[0]:ip[-1]]
        ind = Xg["indices"][ip[0]:ip[-1]]
        blk = sparse.csr_matrix((data, ind, ip - ip[0]), shape=(hi - lo, Xg.attrs["shape"][1]))
        yield r, blk[r - lo]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--donor-col", default="donor_id"); ap.add_argument("--state-col", default="subclass")
    ap.add_argument("--drop-state", default="PVM"); ap.add_argument("--cond-col", default="disease")
    ap.add_argument("--min-cells", type=int, default=1600); ap.add_argument("--n-donors", type=int, default=80)
    ap.add_argument("--cells-per-donor", type=int, default=1600); ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--min-detect", type=float, default=0.10); ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--chunk", type=int, default=5000, help="selected rows per read (memory)")
    a = ap.parse_args(argv)
    if os.path.exists(a.out):
        sys.exit(f"refusing to overwrite {a.out}")
    rng = np.random.default_rng(a.seed)
    with h5py.File(a.h5ad, "r") as h:
        Xg = h["raw/X"] if ("raw" in h and "X" in h["raw"]) else h["X"]
        donors = cat(h, a.donor_col); state = cat(h, a.state_col); cond = cat(h, a.cond_col)
        genes = var_symbols(h)
        vg = h["var"]; ids = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in vg[vg.attrs.get("_index", "_index")][:]])
        extra = {c: cat(h, c) for c in ("assay", "sex") if c in h["obs"]}
        ok = state != a.drop_state
        per = pd.Series(donors[ok]).value_counts()
        elig = per[per >= a.min_cells].index.values
        dcond = pd.Series(cond[ok], index=donors[ok]).groupby(level=0).first().loc[elig]
        # stratified draw of donors by condition, proportional allocation
        pick = []
        for c, grp in dcond.groupby(dcond):
            k = int(round(a.n_donors * len(grp) / len(dcond)))
            pick += list(rng.choice(grp.index.values, min(k, len(grp)), replace=False))
        pick = sorted(pick)[: a.n_donors]
        sel, rank = [], {}
        for d in pick:
            idx = np.nonzero(ok & (donors == d))[0]
            idx = rng.permutation(idx)[: a.cells_per_donor]
            sel.append(idx)
            for r, i in enumerate(idx):
                rank[i] = r
        rows = np.sort(np.concatenate(sel))
        G = Xg.attrs["shape"][1]
        assert len(genes) == G and len(ids) == G, "var does not match the count matrix"
        lib = np.zeros(len(rows)); det = np.zeros(G); pos = {r: j for j, r in enumerate(rows)}
        for r, blk in iter_rows(Xg, rows, a.chunk):             # pass 1
            j = np.array([pos[x] for x in r])
            lib[j] = np.asarray(blk.sum(axis=1)).ravel()
            det += np.asarray((blk > 0).sum(axis=0)).ravel()
        det /= len(rows)
        order = np.argsort(-det)
        cols = [g for g in order if det[g] >= a.min_detect][: a.n_genes]
        for sg in SEX_GENES:                                   # pre-registered: add if detected in >= 0.5%
            for g in np.nonzero(genes == sg)[0]:
                if det[g] > 0.005 and g not in cols:
                    cols.append(int(g))
        cols = np.array(sorted(cols))
        parts = []
        for r, blk in iter_rows(Xg, rows, a.chunk):             # pass 2
            parts.append(blk[:, cols])
        Xs = sparse.vstack(parts).tocsr()
    import anndata as ad
    obs = pd.DataFrame({"donor_id": pd.Categorical(donors[rows]), "lib_size": lib,
                        a.cond_col: pd.Categorical(cond[rows]), "sub_rank": [rank[i] for i in rows],
                        **{c: pd.Categorical(v[rows]) for c, v in extra.items()}},
                       index=[f"cell{i}" for i in rows])
    A = ad.AnnData(X=Xs.astype(np.float32), obs=obs,
                   var=pd.DataFrame({"feature_name": genes[cols]}, index=ids[cols]))
    A.raw = A
    A.write_h5ad(a.out, compression="gzip")
    meta = dict(source=os.path.abspath(a.h5ad), n_cells=int(len(rows)), n_donors=len(pick), donors=list(map(str, pick)),
                donor_conditions=pd.Series([str(dcond.loc[d]) for d in pick]).value_counts().to_dict(),
                n_genes=int(len(cols)), seed=a.seed,
                sex_genes_added=[g for g in SEX_GENES if g in set(genes[cols])],
                duplicated_symbols_in_universe=sorted(set(pd.Series(genes[cols])[pd.Series(genes[cols]).duplicated()])),
                donor_assay=(pd.Series(extra["assay"][rows]).groupby(donors[rows]).first().value_counts().to_dict() if "assay" in extra else None))
    open(a.out + ".universe_ids.txt", "w").write("\n".join(ids[cols]) + "\n")
    open(a.out + ".universe_symbols.txt", "w").write("\n".join(genes[cols]) + "\n")
    json.dump(meta, open(a.out + ".meta.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k != "donors"}))


if __name__ == "__main__":
    main()
