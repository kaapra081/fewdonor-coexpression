"""REV4 Amendment 4 (post hoc): is each edge set organized by the Mantis function partition?
Reads rev4/out/A4/A4_<s>.npz (compact call vectors extracted from rev4_R6.py output) and rev4/mantis/gene_map_mantis.csv;
writes rev4/out/A4/A4_summary.json.  Run from repo/: python rev4/rev4_A4.py"""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OD = os.path.join(HERE, "out", "A4")
gm = pd.read_csv(os.path.join(HERE, "mantis", "gene_map_mantis.csv"))
gm["L4"] = gm.L4_cluster_label.fillna(gm.L3_cluster_label).fillna(gm.L2_cluster_label)   # D1: no level-4 node -> level-3 cluster (11 genes); 8 lncRNAs that are a leaf at level 2 keep it (deviation from D1 wording)
PART = {"leaf": dict(zip(gm.gene, gm.leaf_cluster_label)), "L4": dict(zip(gm.gene, gm.L4))}
NPERM = 10_000


def unpack(Z, k): return np.unpackbits(Z[k])[: int(Z["n"])].astype(bool)


def stats(A, B, lab, K):
    m = len(A); w = int((lab[A] == lab[B]).sum()); d = np.bincount(lab[A], minlength=K) + np.bincount(lab[B], minlength=K)
    return w, w / m - float(((d / (2 * m)) ** 2).sum())


def perm_stats(A, B, lab, K, perms):
    W = np.empty(len(perms)); Q = np.empty(len(perms)); m = len(A)
    for i, p in enumerate(perms):
        L = lab[p]; w = (L[A] == L[B]).sum(); d = np.bincount(L[A], minlength=K) + np.bincount(L[B], minlength=K)
        W[i] = w; Q[i] = w / m - ((d / (2 * m)) ** 2).sum()
    return W, Q


out = {}
for s in ("P1", "P2", "P3"):
    fn = os.path.join(OD, f"A4_{s}.npz")
    if not os.path.exists(fn): continue
    Z = np.load(fn, allow_pickle=True); sym = Z["sym"].astype(str); G = len(sym); iu = np.triu_indices(G, 1)
    c0, cP, cX = unpack(Z, "null_cPB"), unpack(Z, "real_cPB"), unpack(Z, "real_cXDSM")
    sets = {"null_PB": c0, "real_PB": cP, "real_XDSM": cX, "PB_only": cP & ~cX, "XDSM_only": cX & ~cP}
    keep = np.array([g in PART["leaf"] for g in sym]); pair_ok = keep[iu[0]] & keep[iu[1]]
    out[s] = dict(n_genes=int(keep.sum()), n_pairs=int(pair_ok.sum()))
    for pn, pmap in PART.items():
        names = sorted(set(pmap[g] for g in sym[keep])); code = {n: i for i, n in enumerate(names)}
        lab = np.array([code[pmap[g]] if k else -1 for g, k in zip(sym, keep)]); K = len(names)
        gi = np.where(keep)[0]; rng = np.random.default_rng(20_261_103)
        perms = []
        for _ in range(NPERM):                      # permute labels among kept genes only; dropped genes keep -1 and carry no edges
            p = np.arange(G); p[gi] = rng.permutation(gi); perms.append(p)
        res = {}
        allA, allB = iu[0][pair_ok], iu[1][pair_ok]; base_within = float((lab[allA] == lab[allB]).mean())
        for en, cs in sets.items():
            sel = cs & pair_ok; A, B = iu[0][sel], iu[1][sel]
            if len(A) == 0: res[en] = dict(m=0); continue
            w, q = stats(A, B, lab, K); W, Q = perm_stats(A, B, lab, K, perms)
            res[en] = dict(m=int(len(A)), within_frac=w / len(A), within_frac_perm_median=float(np.median(W) / len(A)),
                           p_within=float((1 + (W >= w).sum()) / (NPERM + 1)), Q=q, Q_perm_median=float(np.median(Q)),
                           Q_z=float((q - Q.mean()) / Q.std()), p_Q=float((1 + (Q >= q).sum()) / (NPERM + 1)))
        out[s][pn] = dict(n_clusters=K, all_pairs_within_frac=base_within, sets=res)
        print(s, pn, K, {en: (r.get("m"), round(r.get("within_frac", 0), 4), round(r.get("Q_z", 0), 1), r.get("p_Q")) for en, r in res.items()}, flush=True)
# D4: expression-adjusted eta^2 (P1)
gn = pd.read_csv(os.path.join(HERE, "out", "R6", "genes_numeric.csv")).merge(gm[["gene", "leaf_cluster_label"]], on="gene")
Xd = np.column_stack([np.ones(len(gn)), np.log(gn.mean_umi_per_cell), gn.detection_fraction]); labs = gn.leaf_cluster_label.values
rng = np.random.default_rng(20_261_104)


def eta2(y, lab):
    mu = y.mean(); ss = ((y - mu) ** 2).sum(); return sum(len(y[lab == l]) * (y[lab == l].mean() - mu) ** 2 for l in np.unique(lab)) / ss


PL = [rng.permutation(labs) for _ in range(NPERM)]; d4 = {}
for col in ("null_false_degree", "PB_only_degree", "XDSM_degree"):
    y = gn[col].values.astype(float); r = y - Xd @ np.linalg.lstsq(Xd, y, rcond=None)[0]; e0 = eta2(r, labs); nul = np.array([eta2(r, p) for p in PL])
    d4[col] = dict(eta2_resid=float(e0), perm_median=float(np.median(nul)), p=float((1 + (nul >= e0).sum()) / (NPERM + 1)),
                   r2_expression=float(1 - (r ** 2).sum() / ((y - y.mean()) ** 2).sum()))
out["D4_P1"] = d4; print("D4", d4, flush=True)
json.dump(out, open(os.path.join(OD, "A4_summary.json"), "w"), indent=1)
