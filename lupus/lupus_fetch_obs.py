"""PLAN7 provenance: per-cell metadata of the Perez et al. 2022 lupus dataset (CELLxGENE dataset 218acb0f-9f2f-4f76-b90b-15a4b7c7f629,
version 4118e166-34f5-4c1f-9eed-c64b90a3dace, Census 2025-11-08). (1) Census obs columns; (2) h5ad obs columns library_uuid, sample_uuid,
Processing_Cohort, disease_state, read from the public h5ad with HTTP range requests (obs group only, no expression).
Writes lupus_obs_census.parquet and lupus_obs_h5ad_subset.parquet and verifies that the two have identical rows in identical order.
Run in the 'census' env: python lupus_fetch_obs.py"""
import os, io, time, collections, urllib.parse, numpy as np, pandas as pd, requests, h5py, cellxgene_census, tiledbsoma as soma
DID = "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"; VER = "4118e166-34f5-4c1f-9eed-c64b90a3dace"; CENSUS = "2025-11-08"
URL = f"https://datasets.cellxgene.cziscience.com/{VER}.h5ad"


class RangeFile(io.RawIOBase):
    def __init__(self, url, block=4 * 2 ** 20, maxblocks=64):
        self.url, self.block, self.pos, self.cache, self.maxblocks, self.s = url, block, 0, collections.OrderedDict(), maxblocks, requests.Session()
        self.size = int(self.s.head(url, allow_redirects=True, timeout=60).headers["Content-Length"])
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off); return self.pos
    def _blk(self, i):
        if i in self.cache: self.cache.move_to_end(i); return self.cache[i]
        a = i * self.block; b = min(a + self.block, self.size) - 1
        for att in range(5):
            try:
                r = self.s.get(self.url, headers={"Range": f"bytes={a}-{b}"}, timeout=120); r.raise_for_status(); break
            except Exception:
                if att == 4: raise
                time.sleep(2)
        self.cache[i] = r.content
        if len(self.cache) > self.maxblocks: self.cache.popitem(last=False)
        return r.content
    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        if n <= 0: return 0
        out = bytearray(); p = self.pos; end = p + n
        while p < end:
            i = p // self.block; d = self._blk(i); off = p - i * self.block; take = min(len(d) - off, end - p); out += d[off:off + take]; p += take
        buf[:n] = out; self.pos += n; return n
    def read(self, n=-1):
        if n < 0: n = self.size - self.pos
        b = bytearray(min(n, self.size - self.pos)); k = self.readinto(b); return bytes(b[:k])


def readcol(g, name):
    x = g[name]
    if isinstance(x, h5py.Group):
        codes = x["codes"][:]; cats = np.array([c.decode() if isinstance(c, bytes) else c for c in x["categories"][:]], dtype=object)
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
    return np.array([c.decode() if isinstance(c, bytes) else c for c in x[:]], dtype=object)


if __name__ == "__main__":
    px = urllib.parse.urlparse(os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy"))
    ctx = soma.SOMATileDBContext(tiledb_config={"vfs.s3.proxy_host": px.hostname, "vfs.s3.proxy_port": str(px.port), "vfs.s3.proxy_scheme": px.scheme or "http",
                                                "vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"})
    with cellxgene_census.open_soma(census_version=CENSUS, context=ctx) as c:
        exp = c["census_data"]["homo_sapiens"]
        cols = ["soma_joinid", "observation_joinid", "cell_type", "donor_id", "disease", "sex", "self_reported_ethnicity", "development_stage", "raw_sum", "nnz", "is_primary_data", "assay"]
        ob = exp.obs.read(value_filter=f"dataset_id == '{DID}'", column_names=cols).concat().to_pandas().sort_values("soma_joinid").reset_index(drop=True)
    assert (ob.soma_joinid.diff().dropna() == 1).all()
    for c_ in ob.columns:
        if str(ob[c_].dtype) == "category": ob[c_] = ob[c_].astype(str)
    ob.to_parquet("lupus_obs_census.parquet")
    h = h5py.File(RangeFile(URL), "r"); o = h["obs"]
    H = pd.DataFrame({k: readcol(o, k).astype(str) for k in ["library_uuid", "sample_uuid", "suspension_uuid", "Processing_Cohort", "donor_id", "cell_type", "disease", "sex", "disease_state", "observation_joinid"]})
    for k in ("donor_id", "cell_type", "disease", "sex", "observation_joinid"):
        assert (ob[k].values == H[k].values).all(), f"row mismatch in {k}"
    H.to_parquet("lupus_obs_h5ad_subset.parquet"); print("rows", len(ob), "aligned: donor, cell_type, disease, sex, observation_joinid all identical")
