"""REV4 Amendment 5 E1 and E3 (post hoc), run locally from repo/: python rev4/rev4_A5local.py
Inputs: rev4/out/A4/A4_<s>.npz, A5null20_<s>.npz, A5top_<s>.npz, rev4/mantis/gene_map_mantis.csv, rev4/mantis/uniprot_go_bp_raw.tsv
Output: rev4/out/A4/A5_summary.json"""
import os, re, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OD = os.path.join(HERE, "out", "A4")
gm = pd.read_csv(os.path.join(HERE, "mantis", "gene_map_mantis.csv"))
gm["L4"] = gm.L4_cluster_label.fillna(gm.L3_cluster_label).fillna(gm.L2_cluster_label)
PART = {"leaf": dict(zip(gm.gene, gm.leaf_cluster_label)), "L4": dict(zip(gm.gene, gm.L4))}
go = pd.read_csv(os.path.join(HERE, "mantis", "uniprot_go_bp_raw.tsv"), sep="\t"); go.columns = ["acc", "gene", "gobp"]
go = go[go.gene.isin(gm.gene)]; go["t"] = go.gobp.fillna("").apply(lambda s: set(re.findall(r"GO:\d{7}", s)))
GOT = go.groupby("gene").t.apply(lambda x: set().union(*x)).to_dict()


def unpack(v, n): return np.unpackbits(v)[:n].astype(bool)


def strat_perms(G, keep_idx, bins, n, seed):
    rng = np.random.default_rng(seed); out = np.empty((n, G), dtype=np.int64); base = np.arange(G)
    groups = [keep_idx[bins == b] for b in np.unique(bins)]
    for i in range(n):
        p = base.copy()
        for gi in groups: p[gi] = rng.permutation(gi)
        out[i] = p
    return out


def q_stats(A, B, lab, K):
    m = len(A); w = (lab[A] == lab[B]).sum(); d = np.bincount(lab[A], minlength=K) + np.bincount(lab[B], minlength=K)
    return w / m, w / m - float(((d / (2 * m)) ** 2).sum())


def mod_test(A, B, lab, K, perms):
    f0, q0 = q_stats(A, B, lab, K); F = np.empty(len(perms)); Q = np.empty(len(perms))
    for i, p in enumerate(perms): F[i], Q[i] = q_stats(A, B, lab[p], K)
    return dict(m=int(len(A)), within_frac=float(f0), within_ratio=float(f0 / np.median(F)), Q=float(q0), Q_z=float((q0 - Q.mean()) / Q.std()),
                p_Q=float((1 + (Q >= q0).sum()) / (len(perms) + 1)), p_within=float((1 + (F >= f0).sum()) / (len(perms) + 1)))


def go_test(A, B, Sh, perms):
    s0 = Sh[A, B].mean(); P = np.array([Sh[p[A], p[B]].mean() for p in perms])
    return dict(m=int(len(A)), shared_frac=float(s0), ratio=float(s0 / np.median(P)), p=float((1 + (P >= s0).sum()) / (len(perms) + 1)))


out = {}
for si, s in enumerate(("P1", "P2", "P3")):
    Z = np.load(os.path.join(OD, f"A4_{s}.npz"), allow_pickle=True); sym = Z["sym"].astype(str); G = len(sym); n = int(Z["n"]); iu = np.triu_indices(G, 1)
    keep = np.array([g in PART["leaf"] for g in sym]); ki = np.where(keep)[0]; pair_ok = keep[iu[0]] & keep[iu[1]]
    lm = np.log(np.maximum(Z["mean_umi"][ki], 1e-12)); bins = np.minimum((np.argsort(np.argsort(lm)) * 10) // len(ki), 9)
    c0, cP, cX = (unpack(Z[k], n) for k in ("null_cPB", "real_cPB", "real_cXDSM"))
    sets = {"null_PB": c0, "real_PB": cP, "real_XDSM": cX, "PB_only": cP & ~cX, "XDSM_only": cX & ~cP}
    T = np.load(os.path.join(OD, f"A5top_{s}.npz"), allow_pickle=True)
    for k in ("null_PB", "real_PB", "real_XDSM"):
        v = np.zeros(n, bool); v[T[k]] = True; sets[f"top5000_{k}"] = v
    R = out[s] = dict(n_genes=int(keep.sum()))
    P10k = strat_perms(G, ki, bins, 10_000, 20_261_105)
    for pn, pmap in PART.items():                                     # E1(a), E1(c)
        names = sorted(set(pmap[g] for g in sym[keep])); code = {x: i for i, x in enumerate(names)}
        lab = np.array([code[pmap[g]] if kk else -1 for g, kk in zip(sym, keep)]); K = len(names); res = {}
        for en, cs in sets.items():
            sel = cs & pair_ok; res[en] = mod_test(iu[0][sel], iu[1][sel], lab, K, P10k)
        R[pn] = res; print(s, pn, {e: (r["m"], round(r["within_ratio"], 3), round(r["Q_z"], 1), r["p_Q"]) for e, r in res.items()}, flush=True)
        if pn == "leaf": labL, KL = lab, K
    # E1(d): GO BP terms annotated to 2..50 kept genes
    cnt = {}
    for g in sym[keep]:
        for t in GOT.get(g, ()): cnt[t] = cnt.get(t, 0) + 1
    use = {t for t, v in cnt.items() if 2 <= v <= 50}; TS = [GOT.get(g, set()) & use if kk else set() for g, kk in zip(sym, keep)]
    tl = sorted(use); ti = {t: i for i, t in enumerate(tl)}; M = np.zeros((G, len(tl)), np.float32)
    for gi_, ts in enumerate(TS):
        for t in ts: M[gi_, ti[t]] = 1
    Sh = (M @ M.T) > 0; np.fill_diagonal(Sh, False)
    P10kg = strat_perms(G, ki, bins, 10_000, 20_261_107); gres = {}
    for en, cs in sets.items():
        sel = cs & pair_ok; gres[en] = go_test(iu[0][sel], iu[1][sel], Sh, P10kg)
    R["GO"] = dict(n_terms=len(tl), genes_with_terms=int(sum(len(t) > 0 for t in TS)), all_pairs_shared=float(Sh[iu[0][pair_ok], iu[1][pair_ok]].mean()), sets=gres)
    print(s, "GO", len(tl), {e: (round(r["ratio"], 3), r["p"]) for e, r in gres.items()}, flush=True)
    # E1(e): are Mantis leaves GO-coherent? (all pairs, not edges)
    allA, allB = iu[0][pair_ok], iu[1][pair_ok]; wl = labL[allA] == labL[allB]; f_in = float(Sh[allA[wl], allB[wl]].mean()); f_all = float(Sh[allA, allB].mean())
    Pe = strat_perms(G, ki, bins, 10_000, 20_261_108); nul = np.array([Sh[allA[labL[p][allA] == labL[p][allB]], allB[labL[p][allA] == labL[p][allB]]].mean() for p in Pe])
    R["leaf_GO_check"] = dict(within_leaf_shared=f_in, all_pairs_shared=f_all, perm_median=float(np.median(nul)), p=float((1 + (nul >= f_in).sum()) / (len(Pe) + 1)))
    print(s, "leafGO", R["leaf_GO_check"], flush=True)
    # E1(b): 20 null reassignments, stratified Q z (1,000 relabellings each)
    N20 = np.load(os.path.join(OD, f"A5null20_{s}.npz"), allow_pickle=True); zs = []; ps = []; ms = []
    for k in range(20):
        ck = unpack(N20["calls"][k], n) & pair_ok; r_ = mod_test(iu[0][ck], iu[1][ck], labL, KL, strat_perms(G, ki, bins, 1_000, 20_261_106 + k)); zs.append(r_["Q_z"]); ps.append(r_["p_Q"]); ms.append(r_["m"])
    R["null20"] = dict(z=zs, p=ps, m=ms, z_mean=float(np.mean(zs)), z_sd=float(np.std(zs, ddof=1)), n_p_lt_001=int(sum(p < 0.01 for p in ps)))
    print(s, "null20", round(np.mean(zs), 2), round(np.std(zs, ddof=1), 2), R["null20"]["n_p_lt_001"], flush=True)
    # E3: odds ratio, pseudobulk-only real calls vs null false calls (first reassignment), all pairs of the subtype
    pbo = cP & ~cX; t = [[int((pbo & c0).sum()), int((pbo & ~c0).sum())], [int((~pbo & c0).sum()), int((~pbo & ~c0).sum())]]
    R["PBonly_null_OR"] = dict(table=t, OR=float(t[0][0] * t[1][1] / max(t[0][1] * t[1][0], 1)))
L = out; L["rule"] = dict(a=sum(out[s]["leaf"]["null_PB"]["p_Q"] < 0.01 for s in ("P1", "P2", "P3")), b=sum(out[s]["null20"]["z_mean"] > 2.33 for s in ("P1", "P2", "P3")),
                          d=sum(out[s]["GO"]["sets"]["null_PB"]["p"] < 0.01 for s in ("P1", "P2", "P3")))
print("rule", L["rule"]); json.dump(out, open(os.path.join(OD, "A5_summary.json"), "w"), indent=1)
