"""Plan 9 (breadth/PLAN9.md). python breadth_job.py select | extract <i> | null <i> <k0> <k1>   (writes to $OUT, default ./out)
select: scan Census obs, choose settings (9.1) -> out/selection.json, out/candidates.csv
extract: raw counts of the 400-gene panel for setting i (9.2) -> data/S<i>_X.npz, data/S<i>_meta.npz, out/S<i>_info.json
null: donor-centred null runs k0..k1-1 for setting i (9.3) -> out/b_S<i>_k<k>.json"""
import os, sys, json, time, urllib.parse, numpy as np, scipy.sparse as sp
OUT = os.environ.get("OUT", "out"); DATA = os.environ.get("BDATA", "data"); os.makedirs(OUT, exist_ok=True); os.makedirs(DATA, exist_ok=True)
CENSUS = "2025-11-08"; KEYS = ["psychad", "onek1k", "yazar", "lupus", "sea-ad", "seattle alzheimer", "cell type specific genetic control of autoimmune"]
USED_DIDS = ["09642770-9a4c-4fe6-9f45-67e1249adc33", "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"]   # PsychAD myeloid, Perez lupus (their whole collections are excluded)
USED_COLLS = ["84ce6837-548d-4a1f-919f-0bc0d9a3952f"] if os.environ.get("PLAN") == "10" else []   # Plan 10: PsychAD collection by id
MIN_DONORS = int(os.environ.get("MIN_DONORS", 200)); PER_COLL = int(os.environ.get("PER_COLL", 1)); MAX_SET = int(os.environ.get("MAX_SET", 8)); CAP = int(os.environ.get("CAP", 200)); SEED0 = int(os.environ.get("SEED0", 90_000_000))
stage = sys.argv[1]


def ctx():
    import tiledbsoma as soma
    p = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
    if not p: return None
    px = urllib.parse.urlparse(p)
    return soma.SOMATileDBContext(tiledb_config={"vfs.s3.proxy_host": px.hostname, "vfs.s3.proxy_port": str(px.port), "vfs.s3.proxy_scheme": px.scheme or "http", "vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"})


if stage == "select":
    import pandas as pd, cellxgene_census
    t0 = time.time()
    with cellxgene_census.open_soma(census_version=CENSUS, context=ctx()) as c:
        ds = c["census_info"]["datasets"].read().concat().to_pandas()
        obs = c["census_data"]["homo_sapiens"].obs.read(column_names=["dataset_id", "donor_id", "cell_type", "is_primary_data"]).concat().to_pandas()
    print("obs rows", len(obs), round(time.time() - t0), "s", flush=True)
    obs = obs[obs.is_primary_data]
    n = obs.groupby(["dataset_id", "cell_type", "donor_id"], observed=True).size().rename("n").reset_index()
    el = n[n.n >= 4].groupby(["dataset_id", "cell_type"], observed=True).agg(donors=("donor_id", "nunique"), cells=("n", "sum"), median_n=("n", "median")).reset_index()
    el = el[(el.donors >= MIN_DONORS) & (el.cell_type.astype(str).str.lower() != "unknown")]
    el = el.merge(ds[["dataset_id", "dataset_title", "collection_id", "collection_name", "collection_doi"]], on="dataset_id", how="left")
    txt = (el.collection_name.fillna("") + " | " + el.dataset_title.fillna("")).str.lower()
    used_coll = set(ds[ds.dataset_id.isin(USED_DIDS)].collection_id) | set(USED_COLLS)
    el["excluded"] = txt.apply(lambda s: any(k in s for k in KEYS)) | el.collection_id.isin(used_coll)
    el = el.sort_values("donors", ascending=False, kind="stable").reset_index(drop=True); el.to_csv(f"{OUT}/candidates.csv", index=False)
    keep = el[~el.excluded].groupby("collection_id", sort=False).head(PER_COLL).head(MAX_SET).reset_index(drop=True)
    sel = [dict(i=i + 1, **{k: (str(v) if not isinstance(v, (int, float, np.integer, np.floating)) else float(v)) for k, v in r.items()}) for i, r in keep.iterrows()]
    json.dump(dict(census=CENSUS, n_eligible_settings=int(len(el)), excluded_collections=sorted(set(el[el.excluded].collection_name.astype(str))), settings=sel), open(f"{OUT}/selection.json", "w"), indent=1)
    print(json.dumps(sel, indent=1)[:4000], flush=True)

elif stage == "extract":
    import pandas as pd, cellxgene_census
    i = int(sys.argv[2]); S = [s for s in json.load(open(f"{OUT}/selection.json"))["settings"] if s["i"] == i][0]; TMP = f"{DATA}/tmp_S{i}"; os.makedirs(TMP, exist_ok=True); BLOCK = 20000
    with cellxgene_census.open_soma(census_version=CENSUS, context=ctx()) as c:
        exp = c["census_data"]["homo_sapiens"]
        ob = exp.obs.read(value_filter=f"dataset_id == '{S['dataset_id']}' and is_primary_data == True", column_names=["soma_joinid", "donor_id", "cell_type", "raw_sum", "assay", "disease"]).concat().to_pandas()
        ob = ob[ob.cell_type.astype(str) == S["cell_type"]]; cnt = ob.groupby("donor_id", observed=True).size(); ob = ob[ob.donor_id.isin(cnt[cnt >= 4].index)]
        rng = np.random.default_rng(91_000_000 + i); keep = []
        for d, g in ob.groupby("donor_id", observed=True):
            idx = g.index.values; keep.append(idx if (CAP <= 0 or len(idx) <= CAP) else rng.choice(idx, CAP, replace=False))
        s = ob.loc[np.concatenate(keep)].sort_values("soma_joinid").reset_index(drop=True); ids = s.soma_joinid.values.astype(np.int64); lib = s.raw_sum.values.astype(float); NC = len(s)
        print("S", i, S["cell_type"], "cells", NC, "donors", s.donor_id.nunique(), flush=True)
        var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas().sort_values("soma_joinid").reset_index(drop=True); NV = len(var); X = exp.ms["RNA"].X["raw"]
        libmap = dict(zip(ids, lib)); nb = (NC + BLOCK - 1) // BLOCK; acc = np.zeros(NV)
        for b in range(nb):
            bid = ids[b * BLOCK:(b + 1) * BLOCK]; bl = lib[b * BLOCK:(b + 1) * BLOCK]
            for tb in X.read((bid, slice(None))).tables():
                d0 = tb.column("soma_dim_0").to_numpy(); d1 = tb.column("soma_dim_1").to_numpy(); v = tb.column("soma_data").to_numpy().astype(np.float64)
                acc += np.bincount(d1, weights=v / bl[np.searchsorted(bid, d0)] * 1e4, minlength=NV)
            print("p1", i, b + 1, "/", nb, flush=True)
        cps = acc / NC; panel = np.sort(np.argsort(-cps, kind="stable")[:400]); G = len(panel); sym = var.feature_name.astype(str).values
        mats = []
        for b in range(nb):
            bid = ids[b * BLOCK:(b + 1) * BLOCK]; rr, cc, vv = [], [], []
            for tb in X.read((bid, panel)).tables():
                rr.append(np.searchsorted(bid, tb.column("soma_dim_0").to_numpy())); cc.append(np.searchsorted(panel, tb.column("soma_dim_1").to_numpy())); vv.append(tb.column("soma_data").to_numpy())
            r_, c_, v_ = np.concatenate(rr), np.concatenate(cc), np.concatenate(vv)
            assert np.all(v_ == np.round(v_)), "non-integer counts"
            mats.append(sp.csr_matrix((v_.astype(np.float32), (r_, c_)), shape=(len(bid), G))); print("p2", i, b + 1, "/", nb, flush=True)
    Xp = sp.vstack(mats).tocsr(); Xp.sum_duplicates(); assert Xp.shape == (NC, G)
    sp.save_npz(f"{DATA}/S{i}_X.npz", Xp, compressed=True)
    np.savez_compressed(f"{DATA}/S{i}_meta.npz", lib=lib, donor=s.donor_id.astype(str).values, sym=sym[panel], gene_id=var.feature_id.values[panel].astype(str), soma_joinid=ids, assay=s.assay.astype(str).values, disease=s.disease.astype(str).values)
    nd = s.groupby("donor_id", observed=True).size()
    json.dump(dict(i=i, dataset_id=S["dataset_id"], dataset_title=S["dataset_title"], collection_name=S["collection_name"], collection_doi=S["collection_doi"], cell_type=S["cell_type"],
                   cells=int(NC), donors=int(len(nd)), median_nd=float(nd.median()), min_nd=int(nd.min()), max_nd=int(nd.max()), capped_donors=(int((cnt[cnt >= 4] > CAP).sum()) if CAP > 0 else 0), cap=CAP, G=int(G),
                   assays={k: int(v) for k, v in s.assay.astype(str).value_counts().items()}, median_lib=float(np.median(lib)), nnz_per_cell=float(Xp.nnz / NC)), open(f"{OUT}/S{i}_info.json", "w"), indent=1)

elif stage == "null":
    from scipy.stats import norm
    i, k0, k1 = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); REPO = os.environ.get("REPO", os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    for p in ("dev", "dstd", "pop", "rev", "xsplit"): sys.path.insert(0, os.path.join(REPO, p))
    from methods_v2 import bh, paule_mandel  # noqa: E402
    from dstd import _ms, design_R  # noqa: E402
    from revm import stats_from_cells, methods as rmethods  # noqa: E402
    X = sp.load_npz(f"{DATA}/S{i}_X.npz").tocsr().astype(np.float32); M = np.load(f"{DATA}/S{i}_meta.npz", allow_pickle=True); lib = M["lib"].astype(float); don = M["donor"].astype(str)
    G = X.shape[1]; iu = np.triu_indices(G, 1)
    ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]; st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
    m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)   # Plan 7 mode A centring (lupus_run.centre, exact=False)
    th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for j, r in enumerate(rows): idx[r] = j
    Xs = sp.csr_matrix(X.multiply(f[idx]))
    def tt(q): return 2 * norm.sf(np.abs(q.mean(0) / (q.std(0, ddof=1) / np.sqrt(q.shape[0]))))
    def xdsm(Rs):
        A_, B_ = iu; p = np.empty(len(A_))
        for c0 in range(0, len(A_), 10000):
            aa, bb = A_[c0:c0 + 10000], B_[c0:c0 + 10000]; q = sum(RA[:, aa] * RB[:, bb] for RA, RB in Rs) / len(Rs); p[c0:c0 + 10000] = tt(q)
        return p
    def resid(Z, W): return Z - W @ np.linalg.lstsq(W, Z, rcond=None)[0]
    for k in range(k0, k1):
        fn = f"{OUT}/b_S{i}_k{k}.json"
        if os.path.exists(fn): continue
        t0 = time.time(); seed = SEED0 + 10 * k + i; rng = np.random.default_rng(seed); lab = rng.permutation(don); pu = np.unique(lab)
        rr = [np.where(lab == u)[0] for u in pu]; cl = [(Xs[r].toarray().astype(np.float64), lib[r]) for r in rr]; nd = np.array([len(r) for r in rr], float)
        fs, hs = stats_from_cells(cl, seed); o = rmethods(fs, hs, iu)
        Rs = [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in hs]; nA = np.floor(nd / 2); nB = nd - nA; Rh = []
        for (RA, RB), h in zip(Rs, hs):
            WA = np.column_stack([np.ones(len(nd)), np.log(h["kappaA"] * h["SA"]), np.log(nA)]); WB = np.column_stack([np.ones(len(nd)), np.log(h["kappaB"] * h["SB"]), np.log(nB)]); Rh.append((resid(RA, WA), resid(RB, WB)))
        pv = dict(PB=o["PB"], KX=o["KX"], XDSM=o["XDSM"], XDSMh=xdsm(Rh), XDSMcheck=xdsm(Rs))
        res = {mm: dict(calls=int(bh(np.asarray(v, float)).sum()), t05=float((np.asarray(v, float) < 0.05).mean())) for mm, v in pv.items()}
        res.update(D=int(len(pu)), sec=round(time.time() - t0)); json.dump(res, open(fn, "w")); print("b S", i, k, {mm: v["calls"] for mm, v in res.items() if isinstance(v, dict)}, res["sec"], "s", flush=True)
