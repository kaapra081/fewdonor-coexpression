"""PLAN7 section 0: extract one cell type of the Perez et al. 2022 lupus dataset from CELLxGENE Census 2025-11-08.
python lupus_prep.py LM|FC   (census env; needs lupus_obs_census.parquet and lupus_obs_h5ad_subset.parquet from lupus_fetch_obs.py)
Pass 1 streams the Census X 'raw' layer in cell blocks and accumulates cps_g = sum over included cells of count_g / raw_sum * 1e4 over ALL genes
(exact, no subsample). The panel is the PsychAD rule (pop_prep.py): top 400 by cps, plus forced genes absent from it with cps > 0, sorted by gene index.
Pass 2 extracts the panel columns for the included cells. Resumable: block files are skipped when present. Outputs lupus_<key>_X.npz (uint16 counts,
CSR), lupus_<key>_meta.npz, lupus_<key>_info.json."""
import os, sys, json, time, urllib.parse, numpy as np, pandas as pd, scipy.sparse as sp, cellxgene_census, tiledbsoma as soma
KEY = sys.argv[1]; assert KEY in ("LM", "FC")
CENSUS = "2025-11-08"; DID = "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"; VER = "4118e166-34f5-4c1f-9eed-c64b90a3dace"; COLL = "436154da-bcf1-4130-9c8b-120ff9a888f2"
BLOCK = 20000; G_TOP = 400; TMP = f"prep_tmp_{KEY}"; os.makedirs(TMP, exist_ok=True)
HS_SEX = ["HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1", "HSPD1", "HSPE1", "DNAJB1", "DNAJA1", "HSPA8", "HSPB1", "XIST", "UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14", "KDM5C", "KDM6A"]
ISG = ["ISG15", "IFI6", "IFI44L", "IFI44", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFI27", "RSAD2", "IFITM1", "IFITM3", "STAT1", "IRF7", "XAF1", "LY6E", "USP18",
       "EPSTI1", "SIGLEC1", "HERC5", "CMPK2", "SAMD9L", "PLSCR1", "IFI35", "BST2", "TRIM22"]
FORCE = HS_SEX + ISG

ob = pd.read_parquet("lupus_obs_census.parquet"); H = pd.read_parquet("lupus_obs_h5ad_subset.parquet")
for c_ in ob.columns:
    if str(ob[c_].dtype) == "category": ob[c_] = ob[c_].astype(str)
assert len(ob) == len(H) and (ob.donor_id.values == H.donor_id.values).all() and (ob.soma_joinid.diff().dropna() == 1).all()
BASE = int(ob.soma_joinid.min()); ob["library"] = H.library_uuid.values; ob["sample"] = H.sample_uuid.values; ob["cohort"] = H.Processing_Cohort.values
# few-cell type rule (PLAN7 section 0)
n = ob.groupby(["cell_type", "donor_id"]).size().rename("n").reset_index(); n4 = n[n.n >= 4]
tab = n4.groupby("cell_type").agg(donors=("donor_id", "nunique"), median_n=("n", "median")); tab = tab[tab.donors >= 200].sort_values("median_n")
FC_TYPE = tab.index[0]; CT = "classical monocyte" if KEY == "LM" else FC_TYPE
print("FC rule ->", FC_TYPE, tab.head(3).to_dict("index"), flush=True)
s = ob[ob.cell_type == CT]; cnt = s.groupby("donor_id").size(); s = s[s.donor_id.isin(cnt[cnt >= 4].index)].sort_values("soma_joinid").reset_index(drop=True)
ids = s.soma_joinid.values; lib = s.raw_sum.values.astype(float); NC = len(s); print(KEY, CT, "cells", NC, "donors", s.donor_id.nunique(), flush=True)
libfull = np.zeros(len(ob)); libfull[ids - BASE] = lib

px = urllib.parse.urlparse(os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy"))
ctx = soma.SOMATileDBContext(tiledb_config={"vfs.s3.proxy_host": px.hostname, "vfs.s3.proxy_port": str(px.port), "vfs.s3.proxy_scheme": px.scheme or "http", "vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"})
nb = (NC + BLOCK - 1) // BLOCK
with cellxgene_census.open_soma(census_version=CENSUS, context=ctx) as c:
    exp = c["census_data"]["homo_sapiens"]; var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas().sort_values("soma_joinid").reset_index(drop=True)
    NV = len(var); assert (var.soma_joinid.values == np.arange(NV)).all(); X = exp.ms["RNA"].X["raw"]
    # pass 1
    for b in range(nb):
        fn = f"{TMP}/p1_{b}.npy"
        if os.path.exists(fn): continue
        t0 = time.time(); bid = ids[b * BLOCK:(b + 1) * BLOCK]; acc = np.zeros(NV)
        for tb in X.read((bid, slice(None))).tables():
            d0 = tb.column("soma_dim_0").to_numpy(); d1 = tb.column("soma_dim_1").to_numpy(); v = tb.column("soma_data").to_numpy().astype(np.float64)
            acc += np.bincount(d1, weights=v / libfull[d0 - BASE] * 1e4, minlength=NV)
        np.save(fn, acc); print("p1", b + 1, "/", nb, round(time.time() - t0), "s", flush=True)
    cps = sum(np.load(f"{TMP}/p1_{b}.npy") for b in range(nb)) / NC
    order = np.argsort(-cps, kind="stable"); panel = list(order[:G_TOP]); sym = var.feature_name.astype(str).values; forced_added = []
    for f in FORCE:
        w = np.where(sym == f)[0]
        if len(w) and w[0] not in panel and cps[w[0]] > 0: panel.append(int(w[0])); forced_added.append(f)
    panel = np.array(sorted(int(p) for p in panel)); G = len(panel); print("panel", G, "forced added", forced_added, flush=True)
    # pass 2
    for b in range(nb):
        fn = f"{TMP}/p2_{b}.npz"
        if os.path.exists(fn): continue
        t0 = time.time(); bid = ids[b * BLOCK:(b + 1) * BLOCK]; rr, cc, vv = [], [], []
        for tb in X.read((bid, panel)).tables():
            rr.append(np.searchsorted(bid, tb.column("soma_dim_0").to_numpy())); cc.append(np.searchsorted(panel, tb.column("soma_dim_1").to_numpy())); vv.append(tb.column("soma_data").to_numpy())
        r_, c_, v_ = np.concatenate(rr), np.concatenate(cc), np.concatenate(vv)
        assert v_.max() < 65535 and np.all(v_ == np.round(v_)), "non-integer or too large counts"
        sp.save_npz(fn, sp.csr_matrix((v_.astype(np.uint16), (r_, c_)), shape=(len(bid), G)), compressed=False); print("p2", b + 1, "/", nb, round(time.time() - t0), "s", flush=True)
Xp = sp.vstack([sp.load_npz(f"{TMP}/p2_{b}.npz") for b in range(nb)]).tocsr(); Xp.sum_duplicates(); assert Xp.shape == (NC, G)
sp.save_npz(f"lupus_{KEY}_X.npz", Xp, compressed=True)
np.savez_compressed(f"lupus_{KEY}_meta.npz", lib=lib, donor=s.donor_id.astype(str).values, library=s.library.astype(str).values, sample=s["sample"].astype(str).values, cohort=s.cohort.astype(str).values,
                    disease=s.disease.astype(str).values, sex=s.sex.astype(str).values, soma_joinid=ids, panel=panel, sym=sym[panel], gene_id=var.feature_id.values[panel].astype(str), cps=cps[panel])
dd = s.drop_duplicates("donor_id"); nd = s.groupby("donor_id").size()
isg_in = [g for g in ISG if g in set(sym[panel])]; forced_in = [g for g in FORCE if g in set(sym[panel])]
info = dict(key=KEY, cell_type=CT, census=CENSUS, collection_id=COLL, dataset_id=DID, dataset_version_id=VER, cells=int(NC), donors=int(len(nd)), median_nd=float(nd.median()), iqr_nd=[float(nd.quantile(.25)), float(nd.quantile(.75))],
            donors_by_disease={k: int(v) for k, v in dd.disease.value_counts().items()}, median_nd_by_disease={k: float(v) for k, v in s.groupby(["donor_id", "disease"]).size().reset_index().groupby("disease")[0].median().items()},
            donors_by_sex={k: int(v) for k, v in dd.sex.value_counts().items()}, libraries=int(s.library.nunique()),
            multi_sample_donors=int((s.groupby("donor_id")["sample"].nunique() > 1).sum()), G=int(G), forced_added=forced_added, forced_in_panel=forced_in, isg_in_panel=isg_in, n_isg_in_panel=len(isg_in), nnz=int(Xp.nnz), nnz_per_cell=float(Xp.nnz / NC),
            fc_rule_table=tab.head(5).reset_index().to_dict("records"))
modal_c = s.groupby("donor_id").cohort.agg(lambda x: x.value_counts().index[0]); cd = pd.crosstab(modal_c.reindex(dd.donor_id).values, dd.disease.values)
info["cohort_by_disease_modal_cohort"] = {str(k): {str(a): int(b) for a, b in v.items()} for k, v in cd.to_dict("index").items()}
json.dump(info, open(f"lupus_{KEY}_info.json", "w"), indent=1, default=str); print(json.dumps(info, default=str)[:1500], flush=True)
