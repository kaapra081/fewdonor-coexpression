"""ONEK_PREREG: fetch three cell types (genes ranked on a 5,000-cell random subsample per type; disclosed deviation)."""
import os, json, urllib.parse, numpy as np, pandas as pd, scipy.sparse as sp, cellxgene_census, tiledbsoma as soma
px = urllib.parse.urlparse(os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy"))
ctx = soma.SOMATileDBContext(tiledb_config={"vfs.s3.proxy_host": px.hostname, "vfs.s3.proxy_port": str(px.port), "vfs.s3.proxy_scheme": px.scheme or "http",
                                            "vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"})
obs = pd.read_parquet("onek_obs.parquet")
TYPES = {"B1": "transitional stage B cell", "B2": "naive B cell", "B3": "naive thymus-derived CD4-positive, alpha-beta T cell"}
with cellxgene_census.open_soma(census_version="2025-11-08", context=ctx) as c:
    exp = c["census_data"]["homo_sapiens"]; var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
    for key, ct in TYPES.items():
        if os.path.exists(f"onek_{key}_info.json"): continue
        o = obs[obs.cell_type == ct]; n = o.groupby("donor_id", observed=True).size(); o = o[o.donor_id.isin(n[n >= 4].index)]
        rng = np.random.default_rng(46_000_000); sub = np.sort(rng.choice(o.soma_joinid.values, size=min(5000, len(o)), replace=False))
        with exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(sub,))) as q:
            X = q.to_anndata(X_name="raw").X.tocsr(); ids = q.obs(column_names=["soma_joinid", "raw_sum"]).concat().to_pandas()
            vj = q.var(column_names=["soma_joinid"]).concat().to_pandas().soma_joinid.values
        cp = np.asarray(X.multiply(1e4 / ids.raw_sum.values[:, None]).mean(0)).ravel()
        top = vj[np.argsort(-cp)[:400]]; top = np.sort(top)
        with exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(np.sort(o.soma_joinid.values),)), var_query=soma.AxisQuery(coords=(top,))) as q:
            ad_ = q.to_anndata(X_name="raw", column_names={"obs": ["soma_joinid", "donor_id", "raw_sum"], "var": ["soma_joinid", "feature_name"]})
        sp.save_npz(f"onek_{key}_X.npz", sp.csr_matrix(ad_.X).astype(np.float32))
        np.savez_compressed(f"onek_{key}_meta.npz", lib=ad_.obs.raw_sum.values.astype(float), donor=ad_.obs.donor_id.astype(str).values, sym=ad_.var.feature_name.astype(str).values)
        info = dict(cell_type=ct, cells=int(ad_.n_obs), donors=int(ad_.obs.donor_id.nunique()), median_n=float(ad_.obs.groupby("donor_id", observed=True).size().median()), G=int(ad_.n_vars), census="2025-11-08")
        json.dump(info, open(f"onek_{key}_info.json", "w"), indent=1); print(info, flush=True)
