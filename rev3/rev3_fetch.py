"""REV3_PREREG A2 data: OneK1K B1/B2, genes ranked 1-2,000 by mean CP10k on a 5,000-cell subsample (seed 46,000,000). Runs on any host with internet."""
import os, json, urllib.parse, numpy as np, pandas as pd, scipy.sparse as sp, cellxgene_census, tiledbsoma as soma
HERE = os.path.dirname(os.path.abspath(__file__)); OD = os.path.join(HERE, "data"); os.makedirs(OD, exist_ok=True)
cfg = {"vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"}
pxs = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
if pxs:
    px = urllib.parse.urlparse(pxs); cfg.update({"vfs.s3.proxy_host": px.hostname, "vfs.s3.proxy_port": str(px.port), "vfs.s3.proxy_scheme": px.scheme or "http"})
ctx = soma.SOMATileDBContext(tiledb_config=cfg)
DS = "3faad104-2ab8-4434-816d-474d8d2641db"
TYPES = {"B1": "transitional stage B cell", "B2": "naive B cell"}
with cellxgene_census.open_soma(census_version="2025-11-08", context=ctx) as c:
    exp = c["census_data"]["homo_sapiens"]
    obs = exp.obs.read(value_filter=f"dataset_id == '{DS}'", column_names=["soma_joinid", "donor_id", "cell_type", "raw_sum"]).concat().to_pandas()
    for key, ct in TYPES.items():
        if os.path.exists(os.path.join(OD, f"a2_{key}_meta.npz")): continue
        o = obs[obs.cell_type == ct]; n = o.groupby("donor_id", observed=True).size(); o = o[o.donor_id.isin(n[n >= 4].index)]
        rng = np.random.default_rng(46_000_000); sub = np.sort(rng.choice(o.soma_joinid.values, size=min(5000, len(o)), replace=False))
        with exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(sub,))) as q:
            X = q.to_anndata(X_name="raw").X.tocsr(); ids = q.obs(column_names=["soma_joinid", "raw_sum"]).concat().to_pandas()
            vj = q.var(column_names=["soma_joinid"]).concat().to_pandas().soma_joinid.values
        cp = np.asarray(X.multiply(1e4 / ids.raw_sum.values[:, None]).mean(0)).ravel()
        top = np.sort(vj[np.argsort(-cp)[:2000]])
        with exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(np.sort(o.soma_joinid.values),)), var_query=soma.AxisQuery(coords=(top,))) as q:
            ad_ = q.to_anndata(X_name="raw", column_names={"obs": ["soma_joinid", "donor_id", "raw_sum"], "var": ["soma_joinid", "feature_name"]})
        Xf = sp.csr_matrix(ad_.X).astype(np.float32); umi = np.asarray(Xf.mean(0)).ravel()
        sp.save_npz(os.path.join(OD, f"a2_{key}_X.npz"), Xf)
        np.savez_compressed(os.path.join(OD, f"a2_{key}_meta.npz"), lib=ad_.obs.raw_sum.values.astype(float), donor=ad_.obs.donor_id.astype(str).values,
                            sym=ad_.var.feature_name.astype(str).values, umi=umi)
        info = dict(cell_type=ct, cells=int(ad_.n_obs), donors=int(ad_.obs.donor_id.nunique()), G=int(ad_.n_vars), census="2025-11-08",
                    strata=[int(((umi > lo) & (umi <= hi)).sum()) for lo, hi in ((1, 1e9), (0.3, 1), (0.1, 0.3), (0.03, 0.1), (-1, 0.03))])
        json.dump(info, open(os.path.join(OD, f"a2_{key}_info.json"), "w"), indent=1); print(info, flush=True)
