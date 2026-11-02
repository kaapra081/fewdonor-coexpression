exec(open('crit_common.py').read())

# ---------------------------------------------------------------- FIX
plan('FIX', 'FIX_PREREG (real-cell pseudo-donor null in SEA-AD; LCP2)', ['prereg/FIX_PREREG.md', 'prereg/FIX_RESULTS.md'],
     'prereg/FIX_PREREG.sha256 matches', 'Item 1: K = 5 reassignments x 2 supertypes; Item 2: R = 100 x 6 designs.')
F = J('results/fix/fix_results.json'); N = F['null']; S = F['sim']
se = N['S1']['PB']['typeI_top1_se']; bt = 0.05 + 3 * se
leak = all(N[s]['PB']['mean_calls'] >= 50 and N[s]['PB']['mean_calls'] >= 10 * max(N[s]['LCP']['mean_calls'], 1) for s in ('S1', 'S2'))
assert leak and all(N[s]['decision_leak_real'] == 'YES' for s in ('S1', 'S2'))
add('FIX', 'FIX-1a', 'Leak produces spurious pseudobulk discoveries in real cells: PB mean calls \u226550 and \u226510\u00d7 LCP\u2019s (LCP counted as \u22651).',
    'PASS', f"PB mean false calls {N['S1']['PB']['mean_calls']:,.1f} (S1), {N['S2']['PB']['mean_calls']:,.1f} (S2); LCP {fx(N['S1']['LCP']['mean_calls'],1)}, {fx(N['S2']['LCP']['mean_calls'],1)}. Decision YES in both.",
    'results/fix/fix_results.json: null.{S1,S2}.PB.mean_calls, .LCP.mean_calls, .decision_leak_real')
for meth, rid in (('LCP', 'FIX-1b'), ('LCP2', 'FIX-1c')):
    a = [N[s][meth]['typeI'] for s in ('S1', 'S2')]; t = [N[s][meth]['typeI_top1'] for s in ('S1', 'S2')]
    holds = [x <= 0.055 and y <= bt for x, y in zip(a, t)]
    assert holds == [N[s][meth + '_holds'] for s in ('S1', 'S2')] and not any(holds)
    add('FIX', rid, f"{meth} holds the real-cell null: all-pair type-I \u22640.055 and top-1% within-correlated type-I \u22640.05+3 SE.", 'FAIL',
        f"All-pair type-I {fx(a[0])} / {fx(a[1])} (S1/S2; {'met' if all(x<=0.055 for x in a) else 'not met in S1' if a[0]>0.055 and a[1]<=0.055 else 'not met'}); top-1% {fx(t[0])} / {fx(t[1])} vs bound {fx(bt,4)} (SE {fx(se,4)}): not met. 0 BH calls in all {N['S1']['K']} reassignments.",
        f"results/fix/fix_results.json: null.{{S1,S2}}.{meth}.typeI, .typeI_top1, .typeI_top1_se, .{meth}_holds")
# Item 2
ok = lambda x: x['FDR'] <= 0.05 + 2 * x['FDR_mcse']
a_ = [d for d in ('F_D24', 'F_D40', 'F_S2', 'F_D16_n100') if ok(S[d]['LCP2'])]
b_ = [d for d in ('F_D16', 'F_S1') if ok(S[d]['LCP2'])]
r24 = S['F_D24']['LCP2']['power'] / S['F_D24']['LCP']['power']; r40 = S['F_D40']['LCP2']['power'] / S['F_D40']['LCP']['power']
dec = F['LCP2_decision']
assert (len(a_) == 4) == dec['a'] and (len(b_) >= 1) == dec['b'] and (r24 >= .9 and r40 >= .9) == dec['c'] and dec['adopt'] is False
add('FIX', 'FIX-2', 'Adopt LCP2 if (a) FDR \u22640.05+2 MCSE in F_D24, F_D40, F_S2, F_D16_n100; (b) also in F_D16 or F_S1; (c) power \u22650.9\u00d7 LCP (F_D24, F_D40).',
    'FAIL', f"(a) met: LCP2 FDR {fx(S['F_D24']['LCP2']['FDR'])}, {fx(S['F_D40']['LCP2']['FDR'])}, {fx(S['F_S2']['LCP2']['FDR'])}, {fx(S['F_D16_n100']['LCP2']['FDR'])}. (b) met (F_D16 {fx(S['F_D16']['LCP2']['FDR'])}, F_S1 {fx(S['F_S1']['LCP2']['FDR'])}). (c) not met: power ratio {fx(r24,2)} at F_D24, {fx(r40,2)} at F_D40. Not adopted.",
    'results/fix/fix_results.json: sim.<design>.{LCP,LCP2}.{FDR,FDR_mcse,power}; LCP2_decision.{a,b,c,adopt}')
xcheck('(appendix FIX text) power 0.87 of LCP at 24 donors', fx(r24, 2))
note('FIX', 'No deviation from the plan (FIX_RESULTS.md, Deviations). The pooling caveat (null_caveat.json) is a descriptive addition.')

# ---------------------------------------------------------------- XS
plan('XS', 'XS_PREREG (LCPX, cross-split LCP) with Amendment 1', ['prereg/XS_PREREG.md', 'prereg/POP_RESULTS.md (XS section)'],
     'prereg/XS_PREREG.sha256 matches the amended file; XS_PREREG_original.sha256 hashes the pre-amendment file', 'R = 100 x 6 designs; SEA-AD null K = 5 x 2.')
X = J('results/xsplit/xs_results.json'); XS = X['sim']; XN = X['null']
des = ['F_D16', 'F_D16_n100', 'F_D24', 'F_D40', 'F_S1', 'F_S2']
fail_a = [d for d in des if not ok(XS[d]['LCPX'])]
rat = {d: XS[d]['LCPX']['power'] / XS[d]['LCP']['power'] for d in ('F_D24', 'F_D40')}
dx = X['decision']
assert dx['a'] is False and dx['b'] is False and dx['c'] is False and all(r < .8 for r in rat.values()) and fail_a == ['F_S1', 'F_S2']
add('XS', 'XS-1', 'Adopt LCPX if (a) FDR \u22640.05+2 MCSE in all six designs; (b) holds SEA-AD null in S1, S2; (c) power \u22650.8\u00d7 LCP in F_D24, F_D40.',
    'FAIL', f"(a) not met: LCPX FDR {fx(XS['F_S1']['LCPX']['FDR'])} in F_S1 (bound {fx(.05+2*XS['F_S1']['LCPX']['FDR_mcse'])}) and {fx(XS['F_S2']['LCPX']['FDR'])} in F_S2 ({fx(.05+2*XS['F_S2']['LCPX']['FDR_mcse'])}). (b) not met: type-I {fx(XN['S1']['LCPX']['typeI'])}/{fx(XN['S2']['LCPX']['typeI'])}, top-1% {fx(XN['S1']['LCPX']['typeI_top1'])}/{fx(XN['S2']['LCPX']['typeI_top1'])}. (c) not met: power ratio {fx(rat['F_D24'],2)} (F_D24), {fx(rat['F_D40'],2)} (F_D40).",
    'results/xsplit/xs_results.json: sim.<design>.{LCPX,LCP}; null.{S1,S2}.LCPX; decision.{a,b,c,adopt}')
w = {d: XS[d]['LCPXS'] for d in ('F_D24', 'F_D40', 'F_S1', 'F_S2')}
wn_fail = [d for d in w if w[d]['WN'] > 0.05 + 2 * w[d]['WN_mcse']]; fdr_fail = [d for d in w if not ok(w[d])]
assert dx['weak_null_claim'] is False and wn_fail == ['F_S1'] and fdr_fail == ['F_S1', 'F_S2']
add('XS', 'XS-2', 'Claim LCPX-S tests \u03a3_B=0 only if weak-null type-I and FDR \u22640.05+2 MCSE in F_D24, F_D40, F_S1, F_S2 (weak null = QQ_cross).',
    'FAIL', f"LCPX-S F_S1: FDR {fx(w['F_S1']['FDR'])} (bound {fx(.05+2*w['F_S1']['FDR_mcse'])}), weak-null {fx(w['F_S1']['WN'])} ({fx(.05+2*w['F_S1']['WN_mcse'])}); F_S2 FDR {fx(w['F_S2']['FDR'],4)} (bound {fx(.05+2*w['F_S2']['FDR_mcse'],4)}). F_D24, F_D40 within bounds. Claim not made.",
    'results/xsplit/xs_results.json: sim.<design>.LCPXS.{FDR,FDR_mcse,WN,WN_mcse}; decision.weak_null_claim')
note('XS', 'Amendment 1 (before any full run; only a B = 2000 smoke test had been run) restricted the weak null to QQ_cross; the table uses the amended definition. F_S2 LCPX-S FDR exceeds its bound by 0.00005 (0.0557 vs 0.0556).')

# ---------------------------------------------------------------- POP
plan('POP', 'POP_PREREG (PsychAD population-scale null and real donors)', ['prereg/POP_PREREG.md', 'prereg/POP_RESULTS.md'],
     'prereg/POP_PREREG.sha256 matches (path pop/POP_PREREG.md)', 'Part A: K = 3 per subtype at D = 25, 50, 100 (K = 1 only), 200, 400, 800, all.')
P = J('results/pop/pop_results.json'); A = P['A']
ratio_m = {}; ratio_k = {}
for s in A:
    ratio_m[s] = A[s]['all']['PB']['mean_calls'] / A[s]['100']['PB']['mean_calls']
    ratio_k[s] = A[s]['all']['PB']['calls'][0] / A[s]['100']['PB']['calls'][0]
assert all(r >= 3 and A[s]['all']['PB']['mean_calls'] >= 100 for s, r in ratio_m.items())
add('POP', 'POP-1', 'PB false discoveries grow with cohort size: PB mean calls at D=all \u22653\u00d7 those at D=100 and \u2265100 (per subtype).',
    'PASS', f"YES in P1, P2, P3. Calls at D=all {ic(A['P1']['all']['PB']['mean_calls'])}/{ic(A['P2']['all']['PB']['mean_calls'])}/{ic(A['P3']['all']['PB']['mean_calls'])}; ratio to D=100 {fx(ratio_k['P1'],1)}/{fx(ratio_k['P2'],1)}/{fx(ratio_k['P3'],1)} (same reassignment; with K=3 means {fx(ratio_m['P1'],1)}/{fx(ratio_m['P2'],1)}/{fx(ratio_m['P3'],1)}). D=100 point has K=1.",
    'results/pop/pop_results.json: A.<P>.all.PB.{calls,mean_calls}, A.<P>.100.PB.calls; prereg/POP_RESULTS.md')
Dlist = ['25', '50', '100', '200', '400', '800', 'all']
def lcpx_holds(s, D):
    x = A[s][D]['LCPX']; return x['typeI'] <= 0.055 and x['typeI_top1'] <= 0.05 + 3 * x['se_top1'] and x['mean_calls'] <= 1
holds_by = {s: [D for D in Dlist if D in A[s] and lcpx_holds(s, D)] for s in A}
assert all(len(v) == 0 for v in holds_by.values())
add('POP', 'POP-2', 'LCPX holds the real-cell null at population scale: at every D, all-pair type-I \u22640.055, top-1% \u22640.05+3 SE, mean calls \u22641.',
    'FAIL', f"NO in all three subtypes; criterion fails at every D. LCPX mean false calls at D=all {ic(A['P1']['all']['LCPX']['mean_calls'])}/{ic(A['P2']['all']['LCPX']['mean_calls'])}/{ic(A['P3']['all']['LCPX']['mean_calls'])}; all-pair type-I {fx(A['P1']['all']['LCPX']['typeI'])}/{fx(A['P2']['all']['LCPX']['typeI'])}/{fx(A['P3']['all']['LCPX']['typeI'])}; already {fx(A['P1']['25']['LCPX']['typeI'])}/{fx(A['P2']['25']['LCPX']['typeI'])}/{fx(A['P3']['25']['LCPX']['typeI'])} at D=25. The stated prediction (LCPX near 0 at every D) was not borne out.",
    'results/pop/pop_results.json: A.<P>.<D>.LCPX.{mean_calls,typeI,typeI_top1,se_top1,holds}')
rn = P['rep_null']
assert all(min(rn[s]['PB']['replicated_false']) >= 50 for s in rn)
add('POP', 'POP-3', 'Replication does not protect pseudobulk: PB replicated false calls \u226550 in the Part A null at D=half (per subtype).',
    'PASS', f"YES in all three. Mean replicated false calls PB {ic(rn['P1']['PB']['mean'])}/{ic(rn['P2']['PB']['mean'])}/{ic(rn['P3']['PB']['mean'])} (minimum single run {ic(min(rn['P1']['PB']['replicated_false']))}); LCP {ic(rn['P1']['LCP']['mean'])}/{ic(rn['P2']['LCP']['mean'])}/{ic(rn['P3']['LCP']['mean'])}; XDS 0 in all 9.",
    'results/pop/pop_results.json: rep_null.<P>.{PB,LCP,LCPXd}.{replicated_false,mean}')
B = P['B']
rep = lambda s, m: B[s][m]['replicated'] / B[s][m]['calls_h0']
add('POP', 'POP-4', 'Part B, real donors at full D: calls, overlaps, replication in donor halves. \u201cResults are reported whatever they are\u201d (no threshold).',
    'DESCRIPTIVE', f"PB calls {ic(B['P1']['PB'])}/{ic(B['P2']['PB'])}/{ic(B['P3']['PB'])}; XDS {ic(B['P1']['LCPXd'])}/{ic(B['P2']['LCPXd'])}/{ic(B['P3']['LCPXd'])}; PB-only {ic(B['P1']['PB_only'])}/{ic(B['P2']['PB_only'])}/{ic(B['P3']['PB_only'])}. Half-to-half replication XDS {fx(100*rep('P1','rep_LCPXd'),0)}/{fx(100*rep('P2','rep_LCPXd'),0)}/{fx(100*rep('P3','rep_LCPXd'),0)}%, PB {fx(100*rep('P1','rep_PB'),0)}/{fx(100*rep('P2','rep_PB'),0)}/{fx(100*rep('P3','rep_PB'),0)}%.",
    'results/pop/pop_results.json: B.<P>.{PB,LCPXd,PB_only,rep_PB,rep_LCPXd}')
xcheck('PBnullCalls', f"{ic(min(A[s]['all']['PB']['mean_calls'] for s in A))}\u2013{ic(max(A[s]['all']['PB']['mean_calls'] for s in A))}")
xcheck('PBrepFalse', f"{ic(min(min(rn[s]['PB']['replicated_false']) for s in rn))}\u2013{ic(max(max(rn[s]['PB']['replicated_false']) for s in rn))}")
xcheck('RealPBrep', f"{fx(100*min(rep(s,'rep_PB') for s in B),0)}\u2013{fx(100*max(rep(s,'rep_PB') for s in B),0)}%")
xcheck('RealXDSrep', f"{fx(100*min(rep(s,'rep_LCPXd') for s in B),0)}\u2013{fx(100*max(rep(s,'rep_LCPXd') for s in B),0)}%")
note('POP', 'POP_RESULTS.md lists the real-donor XDS call counts as 32,790 / 49,852 / 55,021, which is not P1/P2/P3 order; the table uses pop_results.json (P1 49,852; P2 55,021; P3 32,790). Deviation (stated in appendix): the D<=100 points were run after the larger ones and with one reassignment per subtype (K = 1, not 3). P3 uses reassignments k = 1..3 (k = 0 was the exploratory reassignment excluded in DSTD_PREREG). POP-1 ratios use the same reassignment (k = 0 for P1, P2; k = 1 for P3), as in POP_RESULTS.md.')

# ---------------------------------------------------------------- DSTD
plan('DSTD', 'DSTD_PREREG (LCPd, LCPXd = XDS)', ['prereg/DSTD_PREREG.md', 'prereg/POP_RESULTS.md (DSTD section)'],
     'prereg/DSTD_PREREG.sha256 matches', 'E1: R = 100 x 6 designs; E2: K = 5 x 2; E3: K = 3 x 3 (P1, P2 k = 0-2; P3 k = 1-3).')
D = J('results/dstd/d_results.json'); DS = D['sim']; DN = D['null']
def dstd(meth):
    a_ok = [d for d in ('F_D24', 'F_D40', 'F_S2', 'F_D16_n100') if ok(DS[d][meth])]
    b_ok = [s for s in ('S1', 'S2') if DN[s][meth]['holds']]
    r24 = DS['F_D24'][meth]['power'] / DS['F_D24']['LCP']['power']; r40 = DS['F_D40'][meth]['power'] / DS['F_D40']['LCP']['power']
    e3 = {}
    for s in A:
        x = A[s]['all'][meth]; e3[s] = (x['mean_calls'], x['typeI'], x['typeI_top1'], x['mean_calls'] <= 1 and x['typeI'] <= 0.055 and x['typeI_top1'] <= 0.05 + 3 * x['se_top1'])
    return dict(a=len(a_ok) == 4, b=len(b_ok) == 2, c=(r24 >= .9 and r40 >= .9), d=all(v[3] for v in e3.values()), r24=r24, r40=r40, e3=e3)
dl = dstd('LCPd'); dx_ = dstd('LCPXd')
assert dl['a'] == D['decision_LCPd']['a'] and dl['b'] == D['decision_LCPd']['b'] and dl['c'] == D['decision_LCPd']['c']
assert dx_['a'] == D['decision_LCPXd']['a'] and dx_['b'] == D['decision_LCPXd']['b'] and dx_['c'] == D['decision_LCPXd']['c']
assert (dl['a'], dl['b'], dl['c'], dl['d']) == (True, False, False, False) and (dx_['a'], dx_['b'], dx_['c'], dx_['d']) == (True, True, False, True)
e = dl['e3']
add('DSTD', 'DSTD-1', 'Adopt LCPd over LCP if (a) FDR bound in four designs; (b) holds SEA-AD null; (c) power \u22650.9\u00d7 LCP; (d) holds PsychAD null at D=all.',
    'FAIL', f"(a) met (FDR {fx(min(DS[d]['LCPd']['FDR'] for d in ('F_D24','F_D40','F_S2','F_D16_n100')))}\u2013{fx(max(DS[d]['LCPd']['FDR'] for d in ('F_D24','F_D40','F_S2','F_D16_n100')))}). (b) not met: top-1% type-I {fx(DN['S1']['LCPd']['typeI_top1'])}/{fx(DN['S2']['LCPd']['typeI_top1'])}. (c) not met: power ratio {fx(dl['r24'],2)} (D24), {fx(dl['r40'],2)} (D40). (d) not met: PsychAD mean false calls {fx(e['P1'][0],1)}/{fx(e['P2'][0],1)}/{ic(e['P3'][0])}; type-I {fx(e['P1'][1])}/{fx(e['P2'][1])}/{fx(e['P3'][1])}.",
    'results/dstd/d_results.json: sim, null, decision_LCPd.{a,b,c}; (d) computed from results/pop/pop_results.json A.<P>.all.LCPd')
e = dx_['e3']
add('DSTD', 'DSTD-2', 'Same adoption rule (a)\u2013(d) for LCPXd (= XDS) as secondary method.',
    'FAIL', f"(a) met (FDR {fx(min(DS[d]['LCPXd']['FDR'] for d in des))}\u2013{fx(max(DS[d]['LCPXd']['FDR'] for d in des))}, all six designs). (b) met: SEA-AD type-I {fx(DN['S1']['LCPXd']['typeI'])}/{fx(DN['S2']['LCPXd']['typeI'])}, 0 calls. (c) not met: power ratio {fx(dx_['r24'],2)} (D24), {fx(dx_['r40'],2)} (D40). (d) met: mean calls {fx(e['P1'][0],1)}/{fx(e['P2'][0],2)}/{fx(e['P3'][0],1)}, type-I {fx(e['P1'][1])}/{fx(e['P2'][1])}/{fx(e['P3'][1])}. XDS does not replace LCP by the frozen rule.",
    'results/dstd/d_results.json: sim.<design>.LCPXd, null.<S>.LCPXd, decision_LCPXd.{a,b,c}; (d) from results/pop/pop_results.json A.<P>.all.LCPXd')
lc = {s: A[s]['all']['LCP']['typeI'] for s in A}; ld = {s: A[s]['all']['LCPd']['typeI'] for s in A}
assert sum(v > .055 for v in lc.values()) >= 2 and not (sum(v > .055 for v in ld.values()) == 0)
add('DSTD', 'DSTD-3', 'Claim \u201cdesign-induced bias is a distinct artifact\u201d if LCP all-pair type-I >0.055 in \u22652 subtypes at D=all while LCPd\u2019s does not.',
    'FAIL', f"LCP type-I {fx(lc['P1'])}/{fx(lc['P2'])}/{fx(lc['P3'])} (>0.055 in 3 of 3) but LCPd {fx(ld['P1'])}/{fx(ld['P2'])}/{fx(ld['P3'])} also >0.055 in 3 of 3. Claim not supported as frozen.",
    'results/pop/pop_results.json: A.<P>.all.{LCP,LCPd}.typeI')
xcheck('XDSnullTI', f"{fx(min(A[s]['all']['LCPXd']['typeI'] for s in A))}\u2013{fx(max(A[s]['all']['LCPXd']['typeI'] for s in A))}")
xcheck('LCPnullTI', f"{fx(min(lc.values()),2)}\u2013{fx(max(lc.values()),2)}")
xcheck('LCPnullCalls', f"{ic(min(A[s]['all']['LCP']['mean_calls'] for s in A))}\u2013{ic(max(A[s]['all']['LCP']['mean_calls'] for s in A))}")
note('DSTD', 'Written after one exploratory reassignment (PsychAD P3, k = 0), excluded from every LCPd/LCPXd evaluation; interim XS simulation results (~45/100 replicates) were also seen before the plan. Criterion (d) is not stored as a flag in d_results.json; it was recomputed here from pop_results.json (same K = 3 reassignments).')

# ---------------------------------------------------------------- XU
plan('XU', 'XU_PREREG (exact within-donor U-statistic, studentized across donors)', ['prereg/XU_PREREG.md'],
     'prereg/XU_PREREG.sha256 matches', 'E1: K = 3 x 3 subtypes, D in {100, 200, 400, all}; E2: K = 5 x 2; E3: R = 50 x 3 designs; E4: R = 100 x 6 designs.')
U = J('results/xu/xu_results.json'); E1 = U['E1']; E3 = U['E3']
cells = [(s, d) for s in E1 for d in ('200', '400', 'all')]
fails = [(s, d) for s, d in cells if not E1[s][d]['XU']['holds']]
assert len(fails) == 9 and U['decision']['holds_pop'] is False
add('XU', 'XU-1', 'XU holds the real-cell null at population scale: at every D\u2265200, every subtype, mean calls \u22641, all-pair type-I \u22640.055, top-1% \u22640.05+3 SE.',
    'FAIL', f"Fails in all 9 subtype\u00d7D cells. XU mean false calls at D=200: {fx(E1['P1']['200']['XU']['mean_calls'],1)}/{fx(E1['P2']['200']['XU']['mean_calls'],1)}/{fx(E1['P3']['200']['XU']['mean_calls'],1)}; at D=all {ic(E1['P1']['all']['XU']['mean_calls'])}/{ic(E1['P2']['all']['XU']['mean_calls'])}/{ic(E1['P3']['all']['XU']['mean_calls'])}. All-pair type-I at D=200 {fx(E1['P1']['200']['XU']['typeI'])}/{fx(E1['P2']['200']['XU']['typeI'])}/{fx(E1['P3']['200']['XU']['typeI'])}.",
    'results/xu/xu_results.json: E1.<P>.<D>.XU.{mean_calls,typeI,typeI_top1,holds}; decision.holds_pop')
fdr_ok = {L: E3[L]['XU']['FDR'] <= 0.05 + 2 * E3[L]['XU']['FDR_mcse'] for L in E3}
pr = {L: E3[L]['XU']['power'] / E3[L]['XDS']['power'] for L in E3}
assert all(fdr_ok.values()) and pr['L500'] < 1.5 and pr['L1000'] < 1.5 and U['decision']['recommended'] is False
add('XU', 'XU-2', 'XU is the recommended large-cohort test: XU-1 AND E3 FDR \u22640.05+2 MCSE in L200, L500, L1000 AND power \u22651.5\u00d7 XDS in L500, L1000.',
    'FAIL', f"XU-1 not met. E3 FDR within bound ({fx(E3['L200']['XU']['FDR'])}/{fx(E3['L500']['XU']['FDR'])}/{fx(E3['L1000']['XU']['FDR'])}). Power ratio XU/XDS {fx(pr['L500'],2)} (L500), {fx(pr['L1000'],2)} (L1000), below 1.5 (XDS power {fx(E3['L500']['XDS']['power'],2)}, {fx(E3['L1000']['XDS']['power'],2)}). XDS remains the large-cohort test, as the plan prescribes.",
    'results/xu/xu_results.json: E3.<L>.{XU,XDS}.{FDR,FDR_mcse,power}; decision.recommended')
wn_ok = {L: E3[L]['XU']['WN'] <= 0.05 + 2 * E3[L]['XU']['WN_mcse'] for L in E3}
assert all(wn_ok.values()) and U['decision']['weak_null'] is True
add('XU', 'XU-3', 'Weak-null validity: XU weak-null type-I \u22640.05+2 MCSE in L200, L500, L1000.',
    'PASS', f"XU weak-null type-I {fx(E3['L200']['XU']['WN'])}/{fx(E3['L500']['XU']['WN'])}/{fx(E3['L1000']['XU']['WN'])} vs bounds {fx(.05+2*E3['L200']['XU']['WN_mcse'])}/{fx(.05+2*E3['L500']['XU']['WN_mcse'])}/{fx(.05+2*E3['L1000']['XU']['WN_mcse'])}.",
    'results/xu/xu_results.json: E3.<L>.XU.{WN,WN_mcse}; decision.weak_null')
E2 = U['E2']
add('XU', 'XU-4', 'E2 (SEA-AD null, XU with t_{D-1}) and E4 (six small designs): \u201creported, no decision.\u201d',
    'DESCRIPTIVE', f"E2 XU mean false calls {fx(E2['S1']['mean_calls'],1)} (S1), {fx(E2['S2']['mean_calls'],1)} (S2); all-pair type-I {fx(E2['S1']['typeI'])}/{fx(E2['S2']['typeI'])}; top-1% {fx(E2['S1']['typeI_top1'])}/{fx(E2['S2']['typeI_top1'])}. E4 results are in the bundle (E4.<design>.XU), not tabulated here.",
    'results/xu/xu_results.json: E2.{S1,S2}; E4')
xcheck('XUnullCalls', f"{ic(min(E1[s]['all']['XU']['mean_calls'] for s in E1))}\u2013{ic(max(E1[s]['all']['XU']['mean_calls'] for s in E1))}")
note('XU', 'Designed after all XS/DSTD/POP results were seen; evaluated on fresh seeds only. XU-3 passes but does not make XU recommended.')

# ---------------------------------------------------------------- COUPLE
plan('COUPLE', 'COUPLE_PREREG (depth-state coupling stress test)', ['prereg/COUPLE_PREREG.md'], 'prereg/COUPLE_PREREG.sha256 matches', 'R = 50 x 3 designs (LS200, LS500, LS1000).')
C = J('results/xu/couple_results.json'); L3 = ['LS200', 'LS500', 'LS1000']
g = lambda m: '/'.join(fx(C[L][m]['FDR']) for L in L3)
pred = all(C[L]['PB']['FDR'] > .2 and C[L]['LCP']['FDR'] > .1 and C[L]['XDS']['FDR'] <= C[L]['LCP']['FDR'] for L in L3) and C['LS500']['XDS']['FDR'] < C['LS500']['XU']['FDR'] and C['LS1000']['XDS']['FDR'] < C['LS1000']['XU']['FDR']
assert pred
add('COUPLE', 'COUPLE-1', 'No adoption decision. Prediction: PB and LCP FDR far above 0.05; XDS FDR inflated less than LCP and XU.',
    'DESCRIPTIVE', f"Prediction consistent with results. FDR at LS200/LS500/LS1000: PB {g('PB')}; LCP {g('LCP')}; XU {g('XU')}; XDS {g('XDS')} (below LCP in all three, below XU at LS500, LS1000, equal at LS200: {fx(C['LS200']['XDS']['FDR'],4)} vs {fx(C['LS200']['XU']['FDR'],4)}).",
    'results/xu/couple_results.json: <LS>.<method>.FDR')
xcheck('CouplePB', ', '.join(fx(C[L]['PB']['FDR']) for L in L3)); xcheck('CoupleLCP', ', '.join(fx(C[L]['LCP']['FDR']) for L in L3))
xcheck('CoupleXDS', ', '.join(fx(C[L]['XDS']['FDR']) for L in L3)); xcheck('CoupleXU', ', '.join(fx(C[L]['XU']['FDR']) for L in L3))

# ---------------------------------------------------------------- DC
plan('DC', 'DC_PREREG (AD differential co-expression)', ['prereg/DC_PREREG.md'], 'prereg/DC_PREREG.sha256 matches', 'K = 3 null reassignments x 3 subtypes (k = 10, 11, 12).')
Dc = J('results/dc/dc_summary.json')
pbm = {s: sum(Dc[s]['PB']['null_calls']) / 3 for s in Dc}
assert sum(v >= 100 for v in pbm.values()) >= 2
add('DC', 'DC-1', 'Pseudobulk DC is confounded by cell number: PB-DC mean false calls \u2265100 in at least two subtypes.',
    'PASS', f"YES in 3 of 3: mean false calls {ic(pbm['P1'])}/{ic(pbm['P2'])}/{ic(pbm['P3'])} (range {ic(min(min(Dc[s]['PB']['null_calls']) for s in Dc))}\u2013{ic(max(max(Dc[s]['PB']['null_calls']) for s in Dc))}); per-pair type-I {fx(min(min(Dc[s]['PB']['null_typeI']) for s in Dc),2)}\u2013{fx(max(max(Dc[s]['PB']['null_typeI']) for s in Dc),2)}.",
    'results/dc/dc_summary.json: <P>.PB.null_calls, .null_typeI')
xm = {s: sum(Dc[s]['XDS']['null_calls']) / 3 for s in Dc}; xt = [t for s in Dc for t in Dc[s]['XDS']['null_typeI']]
assert all(v <= 1 for v in xm.values()) and max(xt) <= 0.055
add('DC', 'DC-2', 'XDS-DC holds: in every subtype mean false calls \u22641 and all-pair type-I \u22640.055.',
    'PASS', f"YES. XDS-DC false calls 0 in all 9 runs; type-I {fx(min(xt),3)}\u2013{fx(max(xt),3)} (max {fx(max(xt),4)}). LCP-DC (not decided): 0\u20136 calls, type-I {fx(min(t for s in Dc for t in Dc[s]['LCP']['null_typeI']),3)}\u2013{fx(max(t for s in Dc for t in Dc[s]['LCP']['null_typeI']),3)}.",
    'results/dc/dc_summary.json: <P>.XDS.null_calls, .null_typeI; <P>.LCP')
add('DC', 'DC-3', 'Real-donor results reported descriptively (number of DC calls per method, overlap).',
    'DESCRIPTIVE', f"Real-donor DC calls P1/P2/P3: PB {ic(Dc['P1']['PB']['real_calls'])}/{ic(Dc['P2']['PB']['real_calls'])}/{ic(Dc['P3']['PB']['real_calls'])}; LCP {ic(Dc['P1']['LCP']['real_calls'])}/{ic(Dc['P2']['LCP']['real_calls'])}/{ic(Dc['P3']['LCP']['real_calls'])}; XDS {ic(Dc['P1']['XDS']['real_calls'])}/{ic(Dc['P2']['XDS']['real_calls'])}/{ic(Dc['P3']['XDS']['real_calls'])}. Overlap between methods not tabulated in the bundle.",
    'results/dc/dc_summary.json: <P>.{PB,LCP,XDS}.real_calls')
xcheck('XDSnull (DC type-I range, nums.tex DCtext)', f"{fx(min(xt),3)}\u2013{fx(max(xt),3)}")
note('DC', 'The plan states only that real-donor results are reported descriptively; group sizes n1/n0 per subtype are in dc_summary.json.')

# ---------------------------------------------------------------- REV
plan('REV', 'REV_PREREG (internal review round 1; R1\u2013R5, amendment R6)', ['prereg/REV_PREREG.md'], 'prereg/REV_PREREG.sha256 matches', 'R1: 100 runs x 3 subtypes; R2: SEA-AD k = 10..14; R3: K = 3; R6: K = 3 x 3.')
R1 = {}
for s in ('P1', 'P2', 'P3'):
    R1[s] = [json.load(open(f)) for f in sorted(glob.glob(ROOT + f'results/rev/out/R1_{s}_k*.json'))]
    assert len(R1[s]) == 100
cnt = lambda m: {s: sum(r[m]['V'] >= 1 for r in R1[s]) for s in R1}
cx, cm, ck, cp = cnt('XDS'), cnt('XDSM'), cnt('KX'), cnt('PB')
mx = max_runs(100)
assert mx == 9
trip = lambda c: ', '.join(str(c[s]) for s in ('P1', 'P2', 'P3'))
add('REV', 'REV-R1a', 'XDS controls FDR at population scale (per subtype): fraction of 100 runs with \u22651 BH call \u22640.05+2 SE (0.0936; \u22649 runs).',
    'FAIL', f"Fails in 3 of 3: runs with a call {trip(cx)} of 100 (P1, P2, P3). Context: PB {trip(cp)} of 100 (mean calls {ic(sum(r['PB']['V'] for r in R1['P1'])/100)}/{ic(sum(r['PB']['V'] for r in R1['P2'])/100)}/{ic(sum(r['PB']['V'] for r in R1['P3'])/100)}).",
    'results/rev/out/R1_<P>_k0..k99.json: XDS.V, PB.V (count of runs with V>=1)')
assert all(v <= mx for v in cm.values())
add('REV', 'REV-R1b', 'Same decision for XDS-M (\u22649 of 100 runs with a call, per subtype).',
    'PASS', f"Holds in 3 of 3: runs with a call {trip(cm)} of 100 (P1, P2, P3).", 'results/rev/out/R1_<P>_k*.json: XDSM.V')
assert [cx_ for cx_ in ck.values()] == [4, 8, 36]
add('REV', 'REV-R1c', 'Same decision for KX (\u22649 of 100 runs with a call, per subtype).',
    'MIXED', f"Holds in P1, P2 ({ck['P1']}, {ck['P2']} of 100); fails in P3 ({ck['P3']} of 100; max {max(r['KX']['V'] for r in R1['P3'])} calls in one run).", 'results/rev/out/R1_<P>_k*.json: KX.V')
xcheck('RoneXDSruns', rng(list(cx.values()))); xcheck('RoneXDSMruns', rng(list(cm.values()))); xcheck('RoneKXruns', rng(list(ck.values())))
mn = lambda m, k: [sum(r[m][k] for r in R1[s]) / 100 for s in R1]
add('REV', 'REV-R1d', 'Per run, counts of p<1e-3, 1e-4, 1e-5 recorded (expected 89, 8.9, 0.89 for ~89,000 pairs); no threshold.',
    'DESCRIPTIVE', f"XDS-M mean counts per run: p<1e-3 {rf(mn('XDSM','n3'),1)}, p<1e-4 {rf(mn('XDSM','n4'),1)}, p<1e-5 {rf(mn('XDSM','n5'),2)} (range over subtypes; pairs = {R1['P1'][0]['P']:,}). XDS: {rf(mn('XDS','n3'),1)}, {rf(mn('XDS','n4'),1)}, {rf(mn('XDS','n5'),2)}.",
    'results/rev/out/R1_<P>_k*.json: <method>.{n3,n4,n5}, P')
xcheck('RoneXDSMnfive (nums.tex)', rng(mn('XDSM', 'n5'), lambda x: fx(x, 2)))
xcheck('RoneXDSnfive (nums.tex)', rng(mn('XDS', 'n5'), lambda x: fx(x, 2)))
sn = J('paper2/seanull_summary.json'); Sx = J('results/rev/rev_sim_results.json')
LG = ('L200', 'L500', 'L1000', 'LS200', 'LS500', 'LS1000')
xs_small = [Sx[d]['XDSM']['FDR'] for d in ('F_D16', 'F_D24', 'F_D16_n100')]; kx_small = [Sx[d]['KX']['FDR'] for d in ('F_D16', 'F_D24', 'F_D40', 'F_D16_n100')]
add('REV', 'REV-R2', 'R2 defines XDS-M (S=5) and KX; evaluated in L/LS and six small designs and the SEA-AD null (k=10..14). No decision rule.',
    'DESCRIPTIVE', f"XDS-M FDR {rf([Sx[d]['XDSM']['FDR'] for d in LG],3)} in the six large designs; {rf(xs_small,3)} in F_D16, F_D24, F_D16_n100. KX FDR {rf(kx_small,3)} in F_D16, F_D24, F_D40, F_D16_n100. SEA-AD top-1% type-I (S1/S2): XDS-M {fx(sn['XDSM'][0])}/{fx(sn['XDSM'][1])}, KX {fx(sn['KX'][0])}/{fx(sn['KX'][1])}.",
    'results/rev/rev_sim_results.json: <design>.{XDSM,KX}.FDR; paper2/seanull_summary.json: XDSM, KX')
xcheck('XDSMsmallFDR', rf(xs_small, 3)); xcheck('KXsmallFDR', rf(kx_small, 2))
xcheck('SeaXDSMTop', f"{fx(min(sn['XDSM']),3)}\u2013{fx(max(sn['XDSM']),3)}"); xcheck('SeaKXTop', f"{fx(min(sn['KX']),3)}\u2013{fx(max(sn['KX']),3)}")
xcheck('XDSMLFDR', f"{fx(min(Sx[d]['XDSM']['FDR'] for d in ('L200','L500','L1000','LS200','LS500','LS1000')),3)}\u2013{fx(max(Sx[d]['XDSM']['FDR'] for d in ('L200','L500','L1000','LS200','LS500','LS1000')),3)}")
Dg = J('paper2/dcg_summary.json'); RD = J('results/rev/rev_dc_r5.json')
dgt = [t for s in Dg for t in Dg[s]['null_typeI']]; dgm = {s: sum(Dg[s]['null_calls']) / 3 for s in Dg}
assert all(v <= 1 for v in dgm.values()) and max(dgt) <= .055
add('REV', 'REV-R3', 'XDS-DCg holds \u201cas in DC_PREREG\u201d: per subtype mean false calls \u22641 and all-pair type-I \u22640.055.',
    'PASS', f"YES. False calls in the 3 runs per subtype: P1 {', '.join(map(str, Dg['P1']['null_calls']))}; P2 {', '.join(map(str, Dg['P2']['null_calls']))}; P3 {', '.join(map(str, Dg['P3']['null_calls']))}; type-I {fx(min(dgt),3)}\u2013{fx(max(dgt),3)}. Real donors: {Dg['P1']['real_calls']}/{Dg['P2']['real_calls']}/{Dg['P3']['real_calls']} calls; overlap with pooled-centre XDS-DC calls not determinable from bundle.",
    'paper2/dcg_summary.json: <P>.null_calls, .null_typeI, .real_calls; results/rev/rev_dc_r5.json: <P>.XDSg_real, .XDSg_in_pooled')
xcheck('DCgTI', f"{fx(min(dgt),3)}\u2013{fx(max(dgt),3)}"); xcheck('DCgReal', ', '.join(str(Dg[s]['real_calls']) for s in Dg))
rr = {s: J(f'results/rev/out/real_{s}.json') for s in ('P1', 'P2', 'P3')}
xy = [x for s in rr for x in rr[s]['sexpairs'] if 'XIST' in x['pair']]; yy = [x for s in rr for x in rr[s]['sexpairs'] if 'XIST' not in x['pair']]
assert len(xy) == 21 and all(x['sign_XDS'] == -1 for x in xy) and len(yy) == 63 and all(x['sign_XDS'] == 1 for x in yy)
add('REV', 'REV-R4', 'Positive control (descriptive): XDS on real donors, XIST\u2013Y-linked pairs; \u201cExpected: negative between-donor association (sex), detected by XDS.\u201d',
    'DESCRIPTIVE', f"Expectation met: {len(xy)} of {len(xy)} XIST\u2013Y pairs negative with p \u2264 {max(x['p_XDS'] for x in xy):.0e} (3 subtypes \u00d7 7 pairs). Y\u2013Y pairs (reported in the paper, not in this plan): {len(yy)} of {len(yy)} positive, max p {max(x['p_XDS'] for x in yy):.0e}.",
    'results/rev/out/real_<P>.json: sexpairs[].{pair,p_XDS,sign_XDS}')
xcheck('PosXYmaxp', f"{max(x['p_XDS'] for x in xy):.0e}"); xcheck('PosYYmaxp', f"{max(x['p_XDS'] for x in yy):.0e}")
add('REV', 'REV-R5', 'Size of pseudobulk false edges (descriptive): absolute r of PB false calls in the POP null at full D; median and counts above 0.1, 0.2.',
    'DESCRIPTIVE', f"Median absolute r {fx(min(RD[s]['false_r_median'] for s in RD),2)}\u2013{fx(max(RD[s]['false_r_median'] for s in RD),2)}; share >0.1: {fx(100*min(RD[s]['false_r_gt01'] for s in RD),0)}\u2013{fx(100*max(RD[s]['false_r_gt01'] for s in RD),0)}%; share >0.2: {fx(100*min(RD[s]['false_r_gt02'] for s in RD),0)}\u2013{fx(100*max(RD[s]['false_r_gt02'] for s in RD),0)}% (P1/P2/P3, n = {ic(RD['P1']['false_r_n'])}/{ic(RD['P2']['false_r_n'])}/{ic(RD['P3']['false_r_n'])} false calls).",
    'results/rev/rev_dc_r5.json: <P>.false_r_median, .false_r_gt01, .false_r_gt02, .false_r_n')
xcheck('FalseRmed', f"{fx(min(RD[s]['false_r_median'] for s in RD),2)}\u2013{fx(max(RD[s]['false_r_median'] for s in RD),2)}")
xcheck('FalseRgtOne', f"{fx(100*min(RD[s]['false_r_gt01'] for s in RD),0)}\u2013{fx(100*max(RD[s]['false_r_gt01'] for s in RD),0)}%")
xcheck('FalseRgtTwo', f"{fx(100*min(RD[s]['false_r_gt02'] for s in RD),0)}\u2013{fx(100*max(RD[s]['false_r_gt02'] for s in RD),0)}%")
R6 = {s: [json.load(open(f)) for f in sorted(glob.glob(ROOT + f'results/rev/out/R6_{s}_k*.json'))] for s in ('P1', 'P2', 'P3')}
assert all(len(v) == 3 for v in R6.values())
pbm6 = [sum(r['PB']['V'] for r in R6[s]) / 3 for s in R6]; xv = [r['XDS']['V'] for s in R6 for r in R6[s]]
pool = J('results/rev/pooling.json')
add('REV', 'REV-R6', 'Amendment R6: donor-centred null to bound pooling inflation (K=3, full D); BH calls, type-I, median within-donor r vs real. Descriptive; no decision.',
    'DESCRIPTIVE', f"PB mean calls {ic(min(pbm6))}\u2013{ic(max(pbm6))} per subtype; XDS {min(xv)}\u2013{max(xv)} per run (type-I {fx(min(r['XDS']['t05'] for s in R6 for r in R6[s]),3)}\u2013{fx(max(r['XDS']['t05'] for s in R6 for r in R6[s]),3)}). Pooling raised the median absolute within-pseudo-donor r by {fx(100*(min(pool[s]['ratio_median'] for s in pool)-1),0)}\u2013{fx(100*(max(pool[s]['ratio_median'] for s in pool)-1),0)}% over real donors.",
    'results/rev/out/R6_<P>_k0..k2.json; results/rev/pooling.json: <P>.ratio_median')
xcheck('RsixPBcalls', f"{ic(min(pbm6))}\u2013{ic(max(pbm6))}"); xcheck('PoolRatio', f"{fx(100*(min(pool[s]['ratio_median'] for s in pool)-1),0)}\u2013{fx(100*(max(pool[s]['ratio_median'] for s in pool)-1),0)}%")
note('REV', 'R1 decision rule: SE is binomial at p = 0.05 (confirmed by the later plans, \u201cSE at p = 0.05\u201d; the rule gives \u22649 of 100). Under the alternative reading (SE at the observed fraction) XDS in P1 and P3 (11 of 100) would be borderline; the paper\u2019s reading is the one used. REV_PREREG states no explicit decision for PB, for R2 or for R4\u2013R6.')

# ---------------------------------------------------------------- ONEK
plan('ONEK', 'ONEK_PREREG (OneK1K blood) with Amendment 1', ['prereg/ONEK_PREREG.md'], 'prereg/ONEK_PREREG.sha256 matches (hashes the file including the amendment)', 'K = 3 reassignments x 3 cell types, D = 200 and all.')
O = J('paper2/onek_summary.json')
gr = {b: O[b]['all']['PB']['mean_calls'] / O[b]['200']['PB']['mean_calls'] for b in O}
gok = {b: (gr[b] >= 2 and O[b]['all']['PB']['mean_calls'] >= 100) for b in O}
assert gok == {'B1': True, 'B2': False, 'B3': True}
nm = {'B1': 'transitional B', 'B2': 'naive B', 'B3': 'naive CD4 T'}
add('ONEK', 'ONEK-1', 'PB false discoveries grow with cohort size (per cell type): PB mean calls at D=all \u22652\u00d7 those at D=200 and \u2265100.',
    'MIXED', f"YES in 2 of 3. Ratio all/200: {fx(gr['B1'],2)} (transitional B; {ic(O['B1']['all']['PB']['mean_calls'])} calls at D=all), {fx(gr['B2'],2)} (naive B; not met), {fx(gr['B3'],2)} (naive CD4 T). Median cells per donor {int(O['B1']['info']['median_n'])}/{int(O['B2']['info']['median_n'])}/{int(O['B3']['info']['median_n'])}.",
    'paper2/onek_summary.json: <B>.all.PB.mean_calls, <B>.200.PB.mean_calls; results/onek/onek_<B>_info.json')
# SE for top-1%: 798 pairs x K = 3
G = 400; npairs = G * (G - 1) // 2; ntop = int(math.ceil(0.01 * npairs)); se_o = math.sqrt(0.05 * 0.95 / (ntop * 3)); bo = 0.05 + 3 * se_o
def okh(x): return x['mean_calls'] <= 1 and x['typeI'] <= 0.055 and x['top1'] <= bo
for b in O:
    for D_ in ('all', '200'):
        for m in ('LCPXd', 'XDSM'):
            assert okh(O[b][D_][m]) == O[b][D_][m]['holds'], (b, D_, m)
xd = {b: all(O[b][D_]['LCPXd']['holds'] for D_ in ('all', '200')) for b in O}; xm_ = {b: all(O[b][D_]['XDSM']['holds'] for D_ in ('all', '200')) for b in O}
assert all(xd.values()) and xm_ == {'B1': False, 'B2': True, 'B3': True}
add('ONEK', 'ONEK-2', 'XDS holds (per cell type) at both D: mean calls \u22641, all-pair type-I \u22640.055, top-1% type-I \u22640.05+3 SE.',
    'PASS', f"Holds in 3 of 3. XDS mean calls \u2264{fx(max(O[b][D_]['LCPXd']['mean_calls'] for b in O for D_ in ('all','200')),1)}; all-pair type-I {fx(min(O[b][D_]['LCPXd']['typeI'] for b in O for D_ in ('all','200')),3)}\u2013{fx(max(O[b][D_]['LCPXd']['typeI'] for b in O for D_ in ('all','200')),3)}; top-1% \u2264{fx(max(O[b][D_]['LCPXd']['top1'] for b in O for D_ in ('all','200')),3)} (bound {fx(bo,3)}).",
    'paper2/onek_summary.json: <B>.{all,200}.LCPXd.{mean_calls,typeI,top1,holds}')
add('ONEK', 'ONEK-3', 'Amendment 1: XDS-M holds (per cell type), same rule at both D.',
    'MIXED', f"Holds in 2 of 3 (naive B, naive CD4 T). Fails in transitional B at D=200: top-1% type-I {fx(O['B1']['200']['XDSM']['top1'],3)} > {fx(bo,3)} (all-pair {fx(O['B1']['200']['XDSM']['typeI'],3)}, 0 calls). Run-level frequencies over K=3 (too few for a decision): {sum(c >= 1 for b in O for D_ in ('all','200') for c in O[b][D_]['XDSM']['calls'])} of {sum(len(O[b][D_]['XDSM']['calls']) for b in O for D_ in ('all','200'))} XDS-M runs had a call across all cells.",
    'paper2/onek_summary.json: <B>.{all,200}.XDSM.{typeI,top1,holds,calls}')
add('ONEK', 'ONEK-4', 'Real donors: XDS and PB calls, PB-only count (descriptive).',
    'NOT DETERMINABLE', 'No OneK1K real-donor results (XDS/PB calls or PB-only counts) are in the bundle.', '(none in bundle)')
xcheck('OnekText ratios (2.2, 1.7, 8.5)', ', '.join(fx(gr[b], 1) for b in ('B1', 'B2', 'B3')))
xcheck('OnekText PB all calls range', f"{ic(min(O[b]['all']['PB']['mean_calls'] for b in O))}\u2013{ic(max(O[b]['all']['PB']['mean_calls'] for b in O))}")
note('ONEK', f"Top-1% SE for ONEK-2/3 is not stored in the summary; recomputed as binomial over {ntop} top pairs \u00d7 K = 3 reassignments (SE {fx(se_o,4)}, bound {fx(bo,4)}); the recomputed verdicts agree with the stored holds flags in all 12 cells.")
