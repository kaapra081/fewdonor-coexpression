"""ONEK_PREREG B3 fetch in chunks of cells with retries (same gene-ranking rule and seed as onek_fetch.py)."""
import os, json, time, urllib.parse, numpy as np, pandas as pd, scipy.sparse as sp, cellxgene_census, tiledbsoma as soma
px = urllib.parse.urlparse(os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy"))
ctx = soma.SOMATileDBContext(tiledb_config={"vfs.s3.proxy_host": px.hostname, "vfs.s3.proxy_port": str(px.port), "vfs.s3.proxy_scheme": px.scheme or "http",
                                            "vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true", "vfs.s3.connect_max_tries": "10", "sm.io_concurrency_level": "4"})
obs = pd.read_parquet("onek_obs.parquet"); key, ct = "B3", "naive thymus-derived CD4-positive, alpha-beta T cell"
o = obs[obs.cell_type == ct]; n = o.groupby("donor_id", observed=True).size(); o = o[o.donor_id.isin(n[n >= 4].index)]
def q_retry(fn, tries=6):
    for t in range(tries):
        try: return fn()
        except Exception as e:
            print("retry", t, str(e)[:120], flush=True); time.sleep(10 * (t + 1))
    raise RuntimeError("failed")
with cellxgene_census.open_soma(census_version="2025-11-08", context=ctx) as c:
    exp = c["census_data"]["homo_sapiens"]
    rng = np.random.default_rng(46_000_000); sub = np.sort(rng.choice(o.soma_joinid.values, size=min(5000, len(o)), replace=False))
    def rank():
        with exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(sub,))) as q:
            X = q.to_anndata(X_name="raw").X.tocsr(); ids = q.obs(column_names=["soma_joinid", "raw_sum"]).concat().to_pandas()
            vj = q.var(column_names=["soma_joinid"]).concat().to_pandas().soma_joinid.values
        cp = np.asarray(X.multiply(1e4 / ids.raw_sum.values[:, None]).mean(0)).ravel(); return np.sort(vj[np.argsort(-cp)[:400]])
    top = q_retry(rank)
    ids_all = np.sort(o.soma_joinid.values); blocks, obsb, var = [], [], None
    for i in range(0, len(ids_all), 25000):
        chunk = ids_all[i:i + 25000]
        def get():
            with exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(chunk,)), var_query=soma.AxisQuery(coords=(top,))) as q:
                return q.to_anndata(X_name="raw", column_names={"obs": ["soma_joinid", "donor_id", "raw_sum"], "var": ["soma_joinid", "feature_name"]})
        a_ = q_retry(get); blocks.append(sp.csr_matrix(a_.X).astype(np.float32)); obsb.append(a_.obs); var = a_.var
        print("chunk", i, a_.n_obs, flush=True)
X = sp.vstack(blocks).tocsr(); ob = pd.concat(obsb)
sp.save_npz(f"onek_{key}_X.npz", X)
np.savez_compressed(f"onek_{key}_meta.npz", lib=ob.raw_sum.values.astype(float), donor=ob.donor_id.astype(str).values, sym=var.feature_name.astype(str).values)
info = dict(cell_type=ct, cells=int(X.shape[0]), donors=int(ob.donor_id.nunique()), median_n=float(ob.groupby("donor_id", observed=True).size().median()), G=int(X.shape[1]), census="2025-11-08")
json.dump(info, open(f"onek_{key}_info.json", "w"), indent=1); print(info, flush=True)
