from itertools import product
SUBS = ('P1', 'P2', 'P3')

# ---------------------------------------------------------------- REV2
plan('REV2', 'REV2_PREREG (internal review round 2) with Amendment 1', ['rev2/REV2_PREREG.md', 'rev2/REV2_AMENDMENT1.md'],
     f"no .sha256 file in bundle; sha256 of bundle files: PREREG {sha('rev2/REV2_PREREG.md')[:12]}\u2026, AMENDMENT1 {sha('rev2/REV2_AMENDMENT1.md')[:12]}\u2026", 'A: 100 runs x 3 subtypes x D in {200, 400, all}; E: 30 runs x 3 cell types.')
R2d = J('rev2/interim/REV2_SUMMARY_latest.json'); RA = R2d['A']
DD = ('200', '400', 'all')
def cellsA(m): return [(s, D, RA[s][D][m]['any'], RA[s][D][m]['n']) for s in SUBS for D in DD]
def holdsA(m): return sum(ok_rl(k, n) for _, _, k, n in cellsA(m))
assert all(n == 100 for m in ('PB', 'LCP', 'LCPd', 'XDS', 'XDSM', 'KX') for _, _, _, n in cellsA(m))
hx = holdsA('XDSM'); failx = [(s, D, k) for s, D, k, n in cellsA('XDSM') if not ok_rl(k, n)]
assert hx == 8 and failx == [('P3', '400', 10)]
cr = lambda m: rf([k for _, _, k, _ in cellsA(m)], 0)
add('REV2', 'REV2-A1', 'Part A, per method/subtype/D (200, 400, all): \u201cholds\u201d if fraction of 100 runs with V\u22651 \u22640.05+2 SE (\u22649 runs). XDS-M.',
    'MIXED', f"Holds in {hx} of 9 cells; fails P3 at D=400 ({failx[0][2]} of 100 runs). Runs with a call across cells: {cr('XDSM')} of 100.",
    'rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.XDSM.{any,n}')
assert holdsA('KX') == 9
add('REV2', 'REV2-A2', 'Same rule, KX.', 'PASS', f"Holds in 9 of 9 cells; runs with a call {cr('KX')} of 100.", 'rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.KX.{any,n}')
assert holdsA('XDS') == 0
add('REV2', 'REV2-A3', 'Same rule, XDS (single split, normal-moment reference).', 'FAIL', f"Holds in 0 of 9 cells; runs with a call {cr('XDS')} of 100.", 'rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.XDS.{any,n}')
assert holdsA('PB') == 0 and holdsA('LCP') == 0
add('REV2', 'REV2-A4', 'Same rule, PB and LCP.', 'FAIL', f"PB and LCP hold in 0 of 9 cells (100 of 100 runs with a call in every cell). Mean BH calls PB {ic(min(RA[s][D]['PB']['meanV'] for s in SUBS for D in DD))}\u2013{ic(max(RA[s][D]['PB']['meanV'] for s in SUBS for D in DD))}; LCP {ic(min(RA[s][D]['LCP']['meanV'] for s in SUBS for D in DD))}\u2013{ic(max(RA[s][D]['LCP']['meanV'] for s in SUBS for D in DD))}.",
    'rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.{PB,LCP}.{any,meanV}')
hl = holdsA('LCPd'); hl_cells = [(s, D) for s, D, k, n in cellsA('LCPd') if ok_rl(k, n)]
both = [(s, D) for s, D, k, n in cellsA('LCPd') if ok_rl(k, n) and ok_rl(*[x for x in cellsA('XDSM') if x[0] == s and x[1] == D][0][2:])]
assert hl == 3 and all(s == 'P1' for s, _ in hl_cells)
add('REV2', 'REV2-A5', 'Primary question: does LCPd hold where XDS-M holds? If LCPd holds at every D and subtype it is reported as co-recommended.',
    'FAIL', f"LCPd holds in {hl} of 9 cells (all in P1: {', '.join(str(RA['P1'][D]['LCPd']['any']) for D in DD)} of 100 runs); fails P2 ({', '.join(str(RA['P2'][D]['LCPd']['any']) for D in DD)}) and P3 ({', '.join(str(RA['P3'][D]['LCPd']['any']) for D in DD)}). Of XDS-M\u2019s {hx} holding cells LCPd holds in {len(both)}. Not co-recommended.",
    'rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.LCPd.{any,n}')
xcheck('RtwoHoldsXDSM', f"{hx} of 9"); xcheck('RtwoHoldsKX', f"{holdsA('KX')} of 9"); xcheck('RtwoHoldsLCPd', f"{hl} of 9"); xcheck('RtwoHoldsXDS', f"{holdsA('XDS')} of 9")
xcheck('RtwoAnyXDSM', f"{sum(k for _,_,k,_ in cellsA('XDSM'))}/900"); xcheck('RtwoAnyXDS', f"{sum(k for _,_,k,_ in cellsA('XDS'))}/900"); xcheck('RtwoAnyLCPd', f"{sum(k for _,_,k,_ in cellsA('LCPd'))}/900"); xcheck('RtwoAnyKX', f"{sum(k for _,_,k,_ in cellsA('KX'))}/900")
xcheck('RtwoRangeXDS', cr('XDS')); xcheck('RtwoRangeXDSM', cr('XDSM'))
# B
def spear(x, y):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
        for pos, i in enumerate(o): r[i] = pos + 1
        return r
    a, b = rk(x), rk(y); n = len(x); ma = sum(a) / n
    num = sum((p - ma) * (q - ma) for p, q in zip(a, b)); den = math.sqrt(sum((p - ma) ** 2 for p in a) * sum((q - ma) ** 2 for q in b)); return num / den
a5 = {s: R2d['C'][s]['A5']['mean_R2'] for s in SUBS}
sp = {s: spear(list(range(10)), a5[s]) for s in SUBS}
add('REV2', 'REV2-B', 'Part B, descriptive: mean over genes of squared design-standardized half-donor values (R^A) by decile of n_d (A5 diagnostic).',
    'DESCRIPTIVE', f"Mean (R^A)^2 by decile {rf([v for s in SUBS for v in a5[s]],2)} across subtypes (P1 {rf(a5['P1'],3)}, P2 {rf(a5['P2'],3)}, P3 {rf(a5['P3'],3)}).",
    'rev2/interim/REV2_SUMMARY_latest.json: C.<P>.A5.{decile_median_n,mean_R2}')
xcheck('AfiveRtwo', rf([v for s in SUBS for v in a5[s]], 2))
# C
C2 = {s: R2d['C'][s] for s in SUBS}
add('REV2', 'REV2-C', 'Part C, descriptive: real donors with XDS-M (calls, split-half replication, three split seeds with Jaccard, XIST\u2013Y and Y\u2013Y pairs).',
    'DESCRIPTIVE', f"XDS-M calls {ic(C2['P1']['seeds']['n_calls'][0])}/{ic(C2['P2']['seeds']['n_calls'][0])}/{ic(C2['P3']['seeds']['n_calls'][0])} ({fx(100*C2['P1']['seeds']['frac_pairs'],0)}/{fx(100*C2['P2']['seeds']['frac_pairs'],0)}/{fx(100*C2['P3']['seeds']['frac_pairs'],0)}% of pairs); Jaccard between split seeds {rf([x for s in SUBS for x in C2[s]['seeds']['jaccard']],2)}; half-to-half replication {rf([100*C2[s]['halves']['replicated']/C2[s]['halves']['calls_h0'] for s in SUBS],0)}%. Sex pairs: XDS-M XIST\u2013Y {sum(C2[s]['sex_XDSM']['xist_y']['sig_right_sign'] for s in SUBS)} of {sum(C2[s]['sex_XDSM']['xist_y']['n'] for s in SUBS)}, Y\u2013Y {sum(C2[s]['sex_XDSM']['y_y']['sig_right_sign'] for s in SUBS)} of {sum(C2[s]['sex_XDSM']['y_y']['n'] for s in SUBS)} significant with the expected sign.",
    'rev2/interim/REV2_SUMMARY_latest.json: C.<P>.{seeds,halves,sex_XDSM}')
xcheck('RealXDSMjac', rf([x for s in SUBS for x in C2[s]['seeds']['jaccard']], 2)); xcheck('RealXDSMrep', rf([100 * C2[s]['halves']['replicated'] / C2[s]['halves']['calls_h0'] for s in SUBS], 0) + '%')
xcheck('RealXDSMcalls', ', '.join(ic(C2[s]['seeds']['n_calls'][0]) for s in SUBS))
# D1
D1 = R2d['D1']; ks = ('k10', 'k11', 'k12')
pbn = [D1[s][k]['PB']['V'] for s in SUBS for k in ks]
nul = {m: [D1[s][k][m]['V'] for s in SUBS for k in ks] for m in ('XDSg', 'XDSMg', 'XDSMgu')}
ti = [D1[s][k][m]['t05'] for s in SUBS for k in ks for m in ('XDSg', 'XDSMg', 'XDSMgu')]
add('REV2', 'REV2-D1', 'Part D1, donor-centred DC null (k=10\u201312): PB-DC, XDS-DCg, XDS-Mg, plus XDS-Mgu (Amendment 1: \u201csame decision rules as XDS-Mg\u201d). No threshold stated.',
    'DESCRIPTIVE', f"PB-DC false calls {ic(min(pbn))}\u2013{ic(max(pbn))} per run (9 runs). XDS-DCg {max(nul['XDSg'])} calls in all 9 runs; XDS-Mg {sum(nul['XDSMg'])} call in total (max {max(nul['XDSMg'])}); XDS-Mgu {max(nul['XDSMgu'])} in all 9; type-I {fx(min(ti),3)}\u2013{fx(max(ti),3)}. Real donors (P1/P2/P3): PB {ic(D1['P1']['real']['PB']['V'])}/{ic(D1['P2']['real']['PB']['V'])}/{ic(D1['P3']['real']['PB']['V'])}; XDS-DCg {D1['P1']['real']['XDSg']['V']}/{D1['P2']['real']['XDSg']['V']}/{D1['P3']['real']['XDSg']['V']}; XDS-Mg {D1['P1']['real']['XDSMg']['V']}/{D1['P2']['real']['XDSMg']['V']}/{D1['P3']['real']['XDSMg']['V']}; XDS-Mgu {D1['P1']['real']['XDSMgu']['V']}/{D1['P2']['real']['XDSMgu']['V']}/{D1['P3']['real']['XDSMgu']['V']}.",
    'rev2/interim/REV2_SUMMARY_latest.json: D1.<P>.{k10,k11,k12,real}.<method>.{V,t05}')
xcheck('DCdcPB', f"{ic(min(pbn))}\u2013{ic(max(pbn))}"); xcheck('DCdcPBTI', rf([D1[s][k]['PB']['t05'] for s in SUBS for k in ks], 2)); xcheck('DCdcXDSTI', f"{fx(min(ti),3)}\u2013{fx(max(ti),3)}")
xcheck('DCrealPB', ', '.join(ic(D1[s]['real']['PB']['V']) for s in SUBS)); xcheck('DCrealM', ', '.join(str(D1[s]['real']['XDSMg']['V']) for s in SUBS)); xcheck('DCrealU', ', '.join(str(D1[s]['real']['XDSMgu']['V']) for s in SUBS)); xcheck('DCrealG', ', '.join(str(D1[s]['real']['XDSg']['V']) for s in SUBS))
# D2, D3
D2 = R2d['D2']
add('REV2', 'REV2-D2', 'Part D2, true-DC simulation (R=50; 101 true-DC pairs of 3,160): FDR and power of PB-DC, XDS-DCg, XDS-Mg (+XDS-Mgu). Reported; no threshold.',
    'DESCRIPTIVE', f"FDR/power: PB-DC {fx(D2['PB']['FDR'])}/{fx(D2['PB']['power'])}; XDS-DCg {fx(D2['XDSg']['FDR'])}/{fx(D2['XDSg']['power'])}; XDS-Mg {fx(D2['XDSMg']['FDR'])}/{fx(D2['XDSMg']['power'])}; XDS-Mgu {fx(D2['XDSMgu']['FDR'])}/{fx(D2['XDSMgu']['power'])} (n_alt {D2['n_alt']}, n_pairs {ic(D2['n_pairs'])}, R = {D2['PB']['R']}).",
    'rev2/interim/REV2_SUMMARY_latest.json: D2.<method>.{FDR,power}, n_alt, n_pairs')
xcheck('DCtwoFDRPB', fx(D2['PB']['FDR'], 3)); xcheck('DCtwoPowPB', fx(D2['PB']['power'], 3)); xcheck('DCtwoFDRM', fx(D2['XDSMg']['FDR'], 3)); xcheck('DCtwoPowM', fx(D2['XDSMg']['power'], 3)); xcheck('DCtwoFDRU', fx(D2['XDSMgu']['FDR'], 3)); xcheck('DCtwoPowU', fx(D2['XDSMgu']['power'], 3)); xcheck('DCtwoFDRG', fx(D2['XDSg']['FDR'], 3)); xcheck('DCtwoPowG', fx(D2['XDSg']['power'], 3))
D3 = R2d['D3']['P1']
add('REV2', 'REV2-D3', 'Part D3, descriptive: split-half replication of real-donor XDS-DCg, XDS-Mg and XDS-Mgu calls (P1).',
    'DESCRIPTIVE', f"Replicated / calls in half 1 (P1): PB-DC {D3['PB']['replicated']}/{ic(D3['PB']['calls_h0'])}; XDS-DCg {D3['XDSg']['replicated']}/{D3['XDSg']['calls_h0']}; XDS-Mg {D3['XDSMg']['replicated']}/{D3['XDSMg']['calls_h0']}; XDS-Mgu {D3['XDSMgu']['replicated']}/{D3['XDSMgu']['calls_h0']}.",
    'rev2/interim/REV2_SUMMARY_latest.json: D3.P1.<method>.{calls_h0,replicated}')
# E
E_ = R2d['E']; nE = {b: 3 + E_[b]['XDSM']['n'] for b in E_}
# XDS-M runs with a call over all 33 runs: 30 (REV2 E) + 3 (ONEK); take from rev2nums macro (any/33) only for the 3 ONEK runs
onek_calls = {b: sum(c >= 1 for D_ in ('all',) for c in J('paper2/onek_summary.json')[b][D_]['XDSM']['calls']) for b in E_}
anyx = {b: E_[b]['XDSM']['any'] + onek_calls[b] for b in E_}
b33 = max_runs(33)
add('REV2', 'REV2-E', 'Part E: 30 further OneK1K reassignments per cell type; run-level fraction for XDS-M over all 33 runs; PBlin vs PB (descriptive). No threshold stated.',
    'DESCRIPTIVE', f"XDS-M runs with a call: {anyx['B1']}/33, {anyx['B2']}/33, {anyx['B3']}/33 (transitional B, naive B, naive CD4 T); all within the Part-A rule for n=33 (\u2264{b33} of 33) \u2014 the paper\u2019s \u201cwithin the bound\u201d. PB {E_['B1']['PB']['any']}/30 in all; PBlin mean calls {ic(E_['B1']['PBlin']['meanV'])}/{ic(E_['B2']['PBlin']['meanV'])}/{ic(E_['B3']['PBlin']['meanV'])} vs PB {ic(E_['B1']['PB']['meanV'])}/{ic(E_['B2']['PB']['meanV'])}/{ic(E_['B3']['PB']['meanV'])}.",
    'rev2/interim/REV2_SUMMARY_latest.json: E.<B>.{PB,PBlin,LCP,XDSM}.{any,n,meanV}; paper2/onek_summary.json: <B>.all.XDSM.calls (3 earlier runs)')
xcheck('OnekEanyXDSM', ', '.join(f"{anyx[b]}/33" for b in ('B1', 'B2', 'B3'))); xcheck('OnekEPBlin', ', '.join(ic(E_[b]['PBlin']['meanV']) for b in ('B1', 'B2', 'B3'))); xcheck('OnekEPB', ', '.join(ic(E_[b]['PB']['meanV']) for b in ('B1', 'B2', 'B3')))
# F
F2 = R2d['F']
add('REV2', 'REV2-F', 'Part F: A5 overdispersion stress simulation (L200, L500; R=50): FDR and power. \u201cNo decision.\u201d',
    'DESCRIPTIVE', f"FDR at L200/L500: XDS-M {fx(F2['L200']['XDSM']['FDR'])}/{fx(F2['L500']['XDSM']['FDR'])}; KX {fx(F2['L200']['KX']['FDR'])}/{fx(F2['L500']['KX']['FDR'])}; LCPd {fx(F2['L200']['LCPd']['FDR'])}/{fx(F2['L500']['LCPd']['FDR'])}; XDS {fx(F2['L200']['XDS']['FDR'])}/{fx(F2['L500']['XDS']['FDR'])}; LCP {fx(F2['L200']['LCP']['FDR'])}/{fx(F2['L500']['LCP']['FDR'])}; PB {fx(F2['L200']['PB']['FDR'])}/{fx(F2['L500']['PB']['FDR'])}. Power XDS-M {fx(F2['L200']['XDSM']['power'],2)}/{fx(F2['L500']['XDSM']['power'],2)}, KX {fx(F2['L200']['KX']['power'],2)}/{fx(F2['L500']['KX']['power'],2)}.",
    'rev2/interim/REV2_SUMMARY_latest.json: F.<L>.<method>.{FDR,power}')
xcheck('FtwoFDRXDSM', fx(F2['L200']['XDSM']['FDR'], 3)); xcheck('FfiveFDRXDSM', fx(F2['L500']['XDSM']['FDR'], 3)); xcheck('FtwoFDRLCPd', fx(F2['L200']['LCPd']['FDR'], 3)); xcheck('FtwoFDRXDS', fx(F2['L200']['XDS']['FDR'], 3)); xcheck('FtwoPowKX', fx(F2['L200']['KX']['power'], 2)); xcheck('FfivePowKX', fx(F2['L500']['KX']['power'], 2))
add('REV2', 'REV2-G', 'Part G: oracle (matched-FDP) power from existing simulation p-values, per replicate (descriptive).',
    'NOT DETERMINABLE', 'No matched-FDP / oracle-power result is in the bundle (the \u201coracle\u201d hits in results/dev and results/heldout refer to the DEV/HELDOUT ORACLE test, not this analysis).', '(none in bundle)')
note('REV2', f"Spearman correlation of decile rank with mean (R^A)^2, computed here (n = 10 deciles, not in the bundle): {fx(sp['P1'],2)}/{fx(sp['P2'],2)}/{fx(sp['P3'],2)} for P1/P2/P3 \u2014 a small drift in magnitude (range 0.895\u20130.996); the paper\u2019s statement of no trend refers to magnitude.")
note('REV2', 'The only REV2 summary in the bundle is rev2/interim/REV2_SUMMARY_latest.json (file name says \u201cinterim\u201d); it holds complete counts for every part (A: 100 runs in each of 9 cells; E: 30 per cell type; D1 k = 10\u201312; D2/F R = 50) and reproduces the rev2nums.tex macros, but the bundle does not show whether a later \u201cfinal\u201d file differs.')
note('REV2', 'REV2_PREREG states no decision rule for D1, D2, D3, E or F; Amendment 1 says XDS-Mgu follows \u201cthe same decision rules as those frozen for XDS-Mg\u201d, but no rule for XDS-Mg is written in REV2_PREREG (REV_PREREG R3 refers to the DC_PREREG rule: mean false calls \u22641, type-I \u22640.055, which the D1 null runs would meet for XDS-DCg, XDS-Mg and XDS-Mgu). Part E: the paper judges XDS-M \u201cwithin the bound\u201d over 33 runs using the Part-A rule (0.05+2 SE, \u22644 of 33); REV2_PREREG does not state this threshold for E.')

# ---------------------------------------------------------------- REV3
plan('REV3', 'REV3_PREREG (internal review round 3) with Amendment 1', ['rev3/REV3_PREREG.md', 'rev3/REV3_AMENDMENT1.md'],
     f"no .sha256 file in bundle; sha256 of bundle files: PREREG {sha('rev3/REV3_PREREG.md')[:12]}\u2026, AMENDMENT1 {sha('rev3/REV3_AMENDMENT1.md')[:12]}\u2026", 'A2: 20 runs per null (B1, B2, P3); A2x: 10; A3: 100 x 3 subtypes; B2: 20 x 3; B3: 40.')
A2s = J('rev3/results/A2_SUMMARY.json')
METH = {'PB': 'PB', 'XDS': 'XDS', 'XDSM1': 'XDSM1', 'XDSM5': 'XDSM5', 'XDSM10': 'XDSM10'}
def tally(m, sel):
    c = [0, 0]
    for ds in A2s:
        for D in ('200', 'all'):
            for st, x in A2s[ds][D][m].items():
                if st in ('Vtotal_mean', 'anyV_total'): continue
                if sel(ds, int(st)): c[1] += 1; c[0] += bool(x['calibrated'])
    return c
selhi = lambda ds, st: ds[:2] in ('B1', 'B2') and not ds.startswith('B1x') and st <= 2
selb2s4 = lambda ds, st: ds.startswith('B2') and st == 3
selp3 = lambda ds, st: ds.startswith('P3')
selx3 = lambda ds, st: ds.startswith('B1x') and st == 2
selx4 = lambda ds, st: ds.startswith('B1x') and st == 3
tl = lambda m, sel: '%d/%d' % tuple(tally(m, sel))
runs_ok = all(A2s[ds]['invalid'] == 0 for ds in A2s)
mc = max(A2s[ds]['maxcheck'] for ds in A2s)
assert runs_ok and mc < 0.05
add('REV3', 'REV3-A2v', 'A2 validity: a run is invalid if the largest GPU\u2013CPU difference in log10 p exceeds 0.05 for p>1e-8 (2,000 random pairs).',
    'PASS', f"0 invalid runs in all 8 data set\u00d7null blocks ({sum(A2s[ds]['runs'] for ds in A2s)} runs); largest check discrepancy {mc:.1e}.",
    'rev3/results/A2_SUMMARY.json: <block>.{runs,invalid,maxcheck}')
cell_txt = lambda m: f"B1/B2 genes >0.1 UMI/cell {tl(m, selhi)}; B2 0.03\u20130.1 UMI/cell (50 genes; 1e-5 not assessable) {tl(m, selb2s4)}; P3 {tl(m, selp3)}"
CRIT_A2 = '\u201cTail-calibrated\u201d per data set/null/D/stratum: obs/exp at p\u22641e-4 and \u22641e-5 both \u22641.5 AND run-level fraction \u22640.147 (n=20). '
assert tally('PB', lambda ds, st: True)[0] == 0
add('REV3', 'REV3-A2a', 'A2 tail calibration rule (obs/exp \u22641.5 at p\u22641e-4 and 1e-5; run-level \u22640.147), PB.', 'FAIL', f"PB tail-calibrated in 0 of 36 combinations ({cell_txt('PB')}).", 'rev3/results/A2_SUMMARY.json: <block>.<D>.PB.<stratum>.{calibrated,ratio,run_ok}')
add('REV3', 'REV3-A2b', 'Same rule, XDS (single split 0, Hoeffding-moment normal reference).', 'MIXED', f"{cell_txt('XDS')}.", 'rev3/results/A2_SUMMARY.json: <block>.<D>.XDS.<stratum>.calibrated')
assert tally('XDSM5', selhi) == [23, 24] and tally('XDSM5', selp3) == [8, 8]
fail5 = [(ds, D, st) for ds in A2s if ds[:2] in ('B1', 'B2') and not ds.startswith('B1x') for D in ('200', 'all') for st, x in A2s[ds][D]['XDSM5'].items() if st not in ('Vtotal_mean', 'anyV_total') and int(st) <= 2 and not x['calibrated']]
fx5 = A2s['B1_n1']['200']['XDSM5']['0']
add('REV3', 'REV3-A2c', 'Same rule, XDS-M with S=5 splits (the validated-range method).', 'MIXED',
    f"{cell_txt('XDSM5')}. One failure: B1 donor-centred null, D=200, >1 UMI/cell ({fx5['anyV']} of {fx5['runs']} runs with a call; bound 2).",
    'rev3/results/A2_SUMMARY.json: <block>.<D>.XDSM5.<stratum>.{calibrated,run_ok,anyV,runs,ratio,exp}')
add('REV3', 'REV3-A2d', 'Same rule, XDS-M with S=1 and S=10 splits.', 'MIXED', f"S=1: {cell_txt('XDSM1')}. S=10: {cell_txt('XDSM10')}.", 'rev3/results/A2_SUMMARY.json: <block>.<D>.{XDSM1,XDSM10}.<stratum>.calibrated')
cal5 = {st: [0, 0] for st in range(5)}
for ds in A2s:
    for D in ('200', 'all'):
        for st, x in A2s[ds][D]['XDSM5'].items():
            if st in ('Vtotal_mean', 'anyV_total'): continue
            cal5[int(st)][1] += 1; cal5[int(st)][0] += bool(x['calibrated'])
lab = {0: '>1', 1: '0.3\u20131', 2: '0.1\u20130.3', 3: '0.03\u20130.1'}
assert cal5[1] == [12, 12] and cal5[2] == [12, 12] and cal5[0] == [11, 12] and cal5[3] == [5, 8]
add('REV3', 'REV3-A2r', 'Validated gene range for XDS-M (S=5): strata tail-calibrated in every data set, null and D.',
    'DESCRIPTIVE', f"Derived here from the flags (UMI/cell): calibrated in every combination at {lab[1]} ({cal5[1][0]}/{cal5[1][1]}) and {lab[2]} ({cal5[2][0]}/{cal5[2][1]}); not at >1 ({cal5[0][0]}/{cal5[0][1]}) or {lab[3]} ({cal5[3][0]}/{cal5[3][1]}). Strict range: 0.1\u20131 UMI/cell; the appendix gives counts, not the range.",
    'rev3/results/A2_SUMMARY.json: <block>.<D>.XDSM5.<stratum>.calibrated (stratum 0..3 = >1, 0.3\u20131, 0.1\u20130.3, 0.03\u20130.1 UMI/cell)')
xcheck('AtwoHiDecM', '%d of %d' % tuple(tally('XDSM5', selhi))); xcheck('AtwoHiDecS', '%d of %d' % tuple(tally('XDS', selhi))); xcheck('AtwoHiDecMone', '%d of %d' % tuple(tally('XDSM1', selhi))); xcheck('AtwoHiDecMten', '%d of %d' % tuple(tally('XDSM10', selhi))); xcheck('AtwoHiDecPB', '%d of %d' % tuple(tally('PB', selhi)))
xcheck('AtwoLoDecM', '%d of %d' % tuple(tally('XDSM5', selx4))); xcheck('AtwoLoDecS', '%d of %d' % tuple(tally('XDS', selx4))); xcheck('AtwoLoDecMone', '%d of %d' % tuple(tally('XDSM1', selx4))); xcheck('AtwoLoDecMten', '%d of %d' % tuple(tally('XDSM10', selx4))); xcheck('AtwoLoDecPB', '%d of %d' % tuple(tally('PB', selx4)))
xcheck('AtwoPthreeDecM', '%d of %d' % tuple(tally('XDSM5', selp3))); xcheck('AtwoPthreeDecS', '%d of %d' % tuple(tally('XDS', selp3))); xcheck('AtwoPthreeDecPB', '%d of %d' % tuple(tally('PB', selp3)))
b10 = 0.05 + 2 * math.sqrt(0.05 * 0.95 / 10)
add('REV3', 'REV3-A2x', 'Amendment 1: same decision rule on OneK1K B1 genes ranked 2,001\u20133,000, n=10 runs (bound 0.188, \u22641 of 10 runs).',
    'MIXED', f"Genes 0.1\u20130.3 UMI/cell (640 genes): XDS-M S=5 {tl('XDSM5', selx3)}, S=1 {tl('XDSM1', selx3)}, S=10 {tl('XDSM10', selx3)}, XDS {tl('XDS', selx3)}, PB {tl('PB', selx3)}. Genes 0.03\u20130.1 UMI/cell (360 genes): XDS-M S=5 {tl('XDSM5', selx4)}, S=1 {tl('XDSM1', selx4)}, S=10 {tl('XDSM10', selx4)}, XDS {tl('XDS', selx4)}, PB {tl('PB', selx4)}. (Combinations = null \u00d7 D.)",
    'rev3/results/A2_SUMMARY.json: B1x_n0, B1x_n1 .<D>.<method>.<stratum>.calibrated')
# A3
A3 = J('rev3/results/A3B3_SUMMARY.json')['A3']
assert all(A3[s]['n'] == 100 for s in SUBS)
def a3cells(m): return [(s, D, A3[s][D][m]['any'], A3[s][D][m]['n']) for D in ('25', '50', '100') for s in SUBS]
def a3h(m): return sum(ok_rl(k, n) for _, _, k, n in a3cells(m))
def a3txt(m): return '; '.join(f"D={D}: " + ', '.join(str(A3[s][D][m]['any']) for s in SUBS) for D in ('25', '50', '100'))
A3CRIT = 'A3, per subtype and D (25, 50, 100), donor-centred null: \u201cholds\u201d if \u22640.0936 of 100 runs (\u22649) have a BH call. '
assert a3h('PB') == 0 and a3h('LCP') == 6 and a3h('LCPd') == 7 and a3h('XDSperm') == 9 and a3h('XDSM') == 4 and a3h('KX') == 0
add('REV3', 'REV3-A3a', 'A3 run-level rule (\u22649 of 100 runs with a call; D=25, 50, 100; per subtype), PB.', 'FAIL', f"Holds in 0 of 9 cells; 100 of 100 runs with a call in every cell.", 'rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.PB.{any,n}')
add('REV3', 'REV3-A3b', 'Same rule, LCP.', 'MIXED', f"Holds in {a3h('LCP')} of 9 cells (D=25, 50 all subtypes); fails at D=100 in all three. Runs with a call (P1, P2, P3): {a3txt('LCP')}.", 'rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.LCP.{any,n}')
add('REV3', 'REV3-A3c', 'Same rule, LCPd.', 'MIXED', f"Holds in {a3h('LCPd')} of 9 cells; fails at D=25 in P2, P3. Runs with a call: {a3txt('LCPd')}.", 'rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.LCPd.{any,n}')
add('REV3', 'REV3-A3d', 'Same rule, XDS with permutation p-values (LCPXd, as in DSTD; B=1e5).', 'PASS', f"Holds in 9 of 9 cells; 0 of 900 runs with a call ({a3txt('XDSperm')}). Caveat (REV4 plan): with B=1e5 the smallest p (~1e-5) cannot reach the BH threshold for ~88,000 pairs. Single-split XDS with normal reference (also in file) holds in {a3h('XDS')} of 9.",
    'rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.XDSperm.{any,n} (and .XDS)')
add('REV3', 'REV3-A3e', 'Same rule, XDS-M (S=5).', 'MIXED', f"Holds in {a3h('XDSM')} of 9 cells: D=25 none of 3, D=50 one of 3 (P1), D=100 three of 3. Runs with a call: {a3txt('XDSM')}.", 'rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.XDSM.{any,n}')
add('REV3', 'REV3-A3f', 'Same rule, KX.', 'FAIL', f"Holds in 0 of 9 cells. Runs with a call: {a3txt('KX')}.", 'rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.KX.{any,n}')
xcheck('AthreeXDSMtwentyfive', '%d%%' % A3['P1']['25']['XDSM']['any']); xcheck('AthreeXDSMfifty', rf([A3[s]['50']['XDSM']['any'] for s in SUBS], 0) + '%'); xcheck('AthreeXDSMhundred', rf([A3[s]['100']['XDSM']['any'] for s in SUBS], 0) + '%')
xcheck('AthreeLCPhundred', rf([A3[s]['100']['LCP']['any'] for s in SUBS], 0) + '%'); xcheck('AthreeLCPdtwentyfive', rf([A3[s]['25']['LCPd']['any'] for s in SUBS], 0) + '%'); xcheck('AthreeLCPsmall', rf([A3[s][D]['LCP']['any'] for s in SUBS for D in ('25', '50')], 0) + '%')
B2s = J('rev3/results/B2_summary.json')
b2cond = all(B2s[s]['PBres']['meanV'] >= 100 and B2s[s]['PBw']['meanV'] >= 100 for s in SUBS)
assert not b2cond and all(B2s[s]['n'] == 20 if 'n' in B2s[s] else B2s[s]['PB']['n'] == 20 for s in SUBS)
add('REV3', 'REV3-B2', 'B2: \u201ccovariate adjustment does not fix pseudobulk\u201d if mean calls \u2265100 for both PBres and PBw in all subtypes (20 runs, D=all).',
    'FAIL', f"Criterion not met (claim not supported). PBres mean calls {ic(B2s['P1']['PBres']['meanV'])}/{ic(B2s['P2']['PBres']['meanV'])}/{ic(B2s['P3']['PBres']['meanV'])} (\u2265100 in all); PBw {fx(B2s['P1']['PBw']['meanV'],1)}/{ic(B2s['P2']['PBw']['meanV'])}/{fx(B2s['P3']['PBw']['meanV'],1)} (<100 in P1, P3). Both made \u22651 call in 20 of 20 runs in every subtype; PBw all-pair type-I {fx(B2s['P1']['PBw']['t05'],3)}/{fx(B2s['P2']['PBw']['t05'],3)}/{fx(B2s['P3']['PBw']['t05'],3)}.",
    'rev3/results/B2_summary.json: <P>.{PB,PBres,PBw}.{meanV,any,n,t05}')
xcheck('BtwoPBw', f"{ic(min(B2s[s]['PBw']['meanV'] for s in SUBS))}\u2013{ic(max(B2s[s]['PBw']['meanV'] for s in SUBS))}"); xcheck('BtwoPBres', f"{ic(min(B2s[s]['PBres']['meanV'] for s in SUBS))}\u2013{ic(max(B2s[s]['PBres']['meanV'] for s in SUBS))}")
B3 = J('rev3/results/A3B3_SUMMARY.json')['B3']['P3']
add('REV3', 'REV3-B3', 'B3: donor-centring sensitivity (clip [0.5, 2]; P3; 40 reassignments; D=200, 400, all). Descriptive comparison; no decision.',
    'DESCRIPTIVE', f"Runs with a call (of {B3['n']}) at D=200/400/all: XDS-M {B3['200']['XDSM']['any']}/{B3['400']['XDSM']['any']}/{B3['all']['XDSM']['any']}; LCPd {B3['200']['LCPd']['any']}/{B3['400']['LCPd']['any']}/{B3['all']['LCPd']['any']}; XDS {B3['200']['XDS']['any']}/{B3['400']['XDS']['any']}/{B3['all']['XDS']['any']}; KX {B3['200']['KX']['any']}/{B3['400']['KX']['any']}/{B3['all']['KX']['any']}; PB and LCP 40/40/40.",
    'rev3/results/A3B3_SUMMARY.json: B3.P3.<D>.<method>.{any,n}')
xcheck('BthreeXDSM', ', '.join(f"{B3[D]['XDSM']['any']}/40" for D in ('200', '400', 'all'))); xcheck('BthreeLCPd', ', '.join(f"{B3[D]['LCPd']['any']}/40" for D in ('200', '400', 'all')))
note('REV3', 'The appendix\u2019s \u201c1 of 4 for 0.03\u20130.1 UMI per cell\u201d (macro AtwoLoDec*) counts only the A2x block (B1, genes 2,001\u20133,000). The original A2 block B2 also has a 0.03\u20130.1 stratum (50 genes), where XDS-M (S=5) is flagged calibrated in 4 of 4 combinations, but its expected count at 1e-5 is 19.7 (<20, not assessable at that threshold); pooling both blocks gives 5 of 8. The A2 criterion for the macro \u201cabove 0.1 UMI per cell\u201d (23 of 24) covers B1 and B2 strata >0.1 only; A2x 0.1\u20130.3 UMI/cell (4 of 4) and P3 (8 of 8) are counted separately.')
note('REV3', 'A3 \u201cXDS\u201d in the frozen plan is the permutation version (LCPXd); the summary file also carries a single-split XDS with normal reference (reported above for context).')

# ---------------------------------------------------------------- REV4
plan('REV4', 'REV4_PREREG (external review) with Amendments 1 and 2', ['rev4/REV4_PREREG.md', 'rev4/REV4_AMENDMENT1.md', 'rev4/REV4_AMENDMENT2.md'],
     f"no .sha256 file in bundle; sha256 of bundle files: PREREG {sha('rev4/REV4_PREREG.md')[:12]}\u2026, AMENDMENT1 {sha('rev4/REV4_AMENDMENT1.md')[:12]}\u2026, AMENDMENT2 {sha('rev4/REV4_AMENDMENT2.md')[:12]}\u2026",
     'R1 PARTIAL (13, 12, 0 of 20 runs); R2 40 per subtype (R2c P1: 36); R3 20/10; R5 20; R6 descriptive.')
R4 = J('rev4/results/REV4_SUMMARY.json'); R1_4 = R4['R1']
nR1 = {s: R1_4[s]['n'] for s in SUBS}
assert nR1 == {'P1': 13, 'P2': 12, 'P3': 0}
MAXR = 2   # plan: at most 2 of 20 runs (0.147); also the bound's whole-run count at n = 12, 13
assert max_runs(12) == 2 and max_runs(13) == 2 and max_runs(20) == 2
def r1cells(m): return [(s, D, R1_4[s][D][m]['any'], R1_4[s][D][m]['n']) for D in ('25', '50', '100') for s in ('P1', 'P2')]
def r1h(m): return sum(k <= MAXR for _, _, k, _ in r1cells(m))
def r1tot(m, D): return sum(R1_4[s][D][m]['any'] for s in ('P1', 'P2')), sum(R1_4[s][D]['PB']['n'] for s in ('P1', 'P2'))
def r1txt(m): return '; '.join(f"D={D}: {r1tot(m, D)[0]}/{r1tot(m, D)[1]}" for D in ('25', '50', '100'))
R1CRIT = 'R1, per subtype and D (25, 50, 100): \u201cholds\u201d if at most 2 of 20 runs have \u22651 BH call (\u22640.147). '
firm = lambda m: all(k >= 3 for _, _, k, _ in r1cells(m) if k > MAXR)
pt = lambda m: f"Holds in {r1h(m)} of 6 evaluated cells (P1, P2 \u00d7 D=25, 50, 100); runs with a call (pooled P1+P2, 25 runs): {r1txt(m)}."
assert r1h('PB') == 0 and r1h('LCP') == 0 and r1h('LCPd') == 3 and r1h('XDSperm') == 6 and r1h('XDS') == 5 and r1h('XDSM') == 4 and r1h('KX') == 2
assert firm('PB') and firm('LCP') and firm('LCPd') and firm('XDS') and firm('XDSM') and firm('KX')
add('REV4', 'REV4-R1a', 'R1 (all R1 rows): \u201cholds\u201d if \u22642 of 20 runs have a BH call, per subtype and D=25, 50, 100. PB and LCP.',
    'FAIL*', f"PB: holds in 0 of 6 cells ({', '.join(str(r1tot('PB', D)[0]) + ' of ' + str(r1tot('PB', D)[1]) for D in ('25','50','100'))} runs with a call at D=25, 50, 100). LCP: 0 of 6; {r1txt('LCP')}. Failures cannot reverse with more runs.",
    'rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.{PB,LCP}.{any,n}')
add('REV4', 'REV4-R1b', 'R1, LCPd.', 'MIXED*', f"{pt('LCPd')} Fails at D=25 (both) and at D=50 in P2 (3 of 12); holds at D=100 (1 of 13, 2 of 12).", 'rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.LCPd.{any,n}')
add('REV4', 'REV4-R1c', 'R1, XDS with permutation p-values (B=2,000,000 in the plan).', 'PASS*', f"{pt('XDSperm')} Provisional: P3 not run; B actually used is not recorded in the summary.", 'rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.XDSperm.{any,n}')
add('REV4', 'REV4-R1d', 'R1, XDS (single split, normal-moment reference, recomputed on the same runs).', 'MIXED*', f"{pt('XDS')} Fails at D=100 in P2 (4 of 12).", 'rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.XDS.{any,n}')
add('REV4', 'REV4-R1e', 'R1, XDS-M (S=5).', 'MIXED*', f"{pt('XDSM')} Fails at D=25 (6 of 13, 7 of 12); holds at D=50, 100.", 'rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.XDSM.{any,n}')
add('REV4', 'REV4-R1f', 'R1, KX.', 'MIXED*', f"{pt('KX')} Fails at D=25 and D=50; holds at D=100 (2 of 13, 1 of 12).", 'rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.KX.{any,n}')
add('REV4', 'REV4-R1g', 'R1 in P3 (ADAM_GPNMB): 20 runs planned at D=25, 50, 100.', 'NOT EVALUATED', f"0 of 20 runs in the bundle (R1.P3.n = {nR1['P3']}). P1: {nR1['P1']} of 20 runs, P2: {nR1['P2']} of 20 (partial). Plan: \u201creport the runs completed (stated as such)\u201d if compute does not allow 20.", 'rev4/results/REV4_SUMMARY.json: R1.<P>.n')
for nm_, mac, comp in (('LCP', 'RfLCPTwentyFive', r1tot('LCP', '25')), ('LCPd', 'RfLCPdTwentyFive', r1tot('LCPd', '25')), ('XDSperm', 'RfXDSpermTwentyFive', r1tot('XDSperm', '25')), ('XDSM', 'RfXDSMTwentyFive', r1tot('XDSM', '25')), ('KX', 'RfKXTwentyFive', r1tot('KX', '25')),
                       ('LCP', 'RfLCPFifty', r1tot('LCP', '50')), ('LCPd', 'RfLCPdFifty', r1tot('LCPd', '50')), ('XDSperm', 'RfXDSpermFifty', r1tot('XDSperm', '50')), ('XDSM', 'RfXDSMFifty', r1tot('XDSM', '50')), ('KX', 'RfKXFifty', r1tot('KX', '50')),
                       ('LCP', 'RfLCPHundred', r1tot('LCP', '100')), ('LCPd', 'RfLCPdHundred', r1tot('LCPd', '100')), ('XDSperm', 'RfXDSpermHundred', r1tot('XDSperm', '100')), ('XDSM', 'RfXDSMHundred', r1tot('XDSM', '100')), ('KX', 'RfKXHundred', r1tot('KX', '100'))):
    xcheck(mac, f"{comp[0]}/{comp[1]}")
# R2
R2 = R4['R2']
def r2c(m): return [R2[s][m]['any'] for s in SUBS]
assert all(R2[s]['XDSM']['n'] == 40 for s in SUBS) and max_runs(40) == 4
hxm = sum(ok_rl(k, 40) for k in r2c('XDSM')); assert hxm == 1 and r2c('XDSM') == [15, 8, 2]
oth = lambda m: ', '.join(f"{k}/40" for k in r2c(m))
add('REV4', 'REV4-R2', 'R2 decision: XDS-M holds by the run-level rule (n=40: \u22644 of 40 runs with a call) in each subtype (brain-bank-stratified null).',
    'MIXED', f"Holds in 1 of 3 (P3: 2/40); fails P1 15/40 and P2 8/40. Other methods (P1, P2, P3; same rule, not a named decision): PB {oth('PB')}; LCP {oth('LCP')}; LCPd {oth('LCPd')}; XDS {oth('XDS')}; KX {oth('KX')}.",
    'rev4/results/REV4_SUMMARY.json: R2.<P>.<method>.{any,n}')
xcheck('RfBankXDSM', oth('XDSM')); xcheck('RfBankPB', oth('PB')); xcheck('RfBankLCP', oth('LCP')); xcheck('RfBankLCPd', oth('LCPd')); xcheck('RfBankXDS', oth('XDS')); xcheck('RfBankKX', oth('KX'))
# R2c
R2c = R4['R2c']
xb = [(R2c[s]['XDSMb']['any'], R2c[s]['XDSMb']['n']) for s in SUBS]; xm_c = [(R2c[s]['XDSM']['any'], R2c[s]['XDSM']['n']) for s in SUBS]
rule = [ok_rl(k, n) for k, n in xb if True]
assert xb == [(4, 36), (1, 40), (2, 40)] and all(rule)
add('REV4', 'REV4-R2c', 'Amendment 2 (post hoc; rule fixed before R2c): if XDS-Mb fails \u2264 the bound\u2019s runs where XDS-M fails, attribute R2 failures to residual bank structure.',
    'PASS*', f"Rule met. XDS-M fails in P1 ({xm_c[0][0]}/{xm_c[0][1]} runs) and P2 ({xm_c[1][0]}/{xm_c[1][1]}); XDS-Mb {xb[0][0]}/{xb[0][1]} (P1; exactly at the 4-run bound), {xb[1][0]}/{xb[1][1]} (P2), {xb[2][0]}/{xb[2][1]} (P3). Provisional: P1 has 36 of 40 runs. PB and bank-centred PB (PBb): 40/40 runs with a call in P2, P3 and 36/36 in P1.",
    'rev4/results/REV4_SUMMARY.json: R2c.<P>.{XDSM,XDSMb,PB,PBb}.{any,n}')
xcheck('RfBankXDSMb', ', '.join(f"{k}/{n}" for k, n in xb)); xcheck('RfBankXDSMc', ', '.join(f"{k}/{n}" for k, n in xm_c)); xcheck('RfBankXDSMbHolds', str(sum(rule)))
# R3
R3 = R4['R3']; R3o = R4['R3o']
assert all(R3[s][m]['any'] == R3[s][m]['n'] for s in SUBS for m in R3[s]) and all(R3o[b][m]['any'] == R3o[b][m]['n'] for b in R3o for m in R3o[b])
add('REV4', 'REV4-R3', 'R3: PC-adjusted (k=5, 10) and weighted pseudobulk baselines, PsychAD (20 runs) and OneK1K (10 runs). Descriptive; no decision.',
    'DESCRIPTIVE', f"Every variant made \u22651 call in every run (PsychAD 20/20; OneK1K 10/10). Mean calls PsychAD: PBres {ic(min(R3[s]['PBres']['meanV'] for s in SUBS))}\u2013{ic(max(R3[s]['PBres']['meanV'] for s in SUBS))}; PBw {ic(min(R3[s]['PBw']['meanV'] for s in SUBS))}\u2013{ic(max(R3[s]['PBw']['meanV'] for s in SUBS))}; PC5/PC10 {ic(min(R3[s][m]['meanV'] for s in SUBS for m in ('PC5','PC10')))}\u2013{ic(max(R3[s][m]['meanV'] for s in SUBS for m in ('PC5','PC10')))}. OneK1K: PB {ic(min(R3o[b]['PB']['meanV'] for b in R3o))}\u2013{ic(max(R3o[b]['PB']['meanV'] for b in R3o))}; PBw {ic(min(R3o[b]['PBw']['meanV'] for b in R3o))}\u2013{ic(max(R3o[b]['PBw']['meanV'] for b in R3o))}; PC {ic(min(R3o[b][m]['meanV'] for b in R3o for m in ('PC5','PC10')))}\u2013{ic(max(R3o[b][m]['meanV'] for b in R3o for m in ('PC5','PC10')))}.",
    'rev4/results/REV4_SUMMARY.json: R3.<P>.<method>.{any,n,meanV}; R3o.<B>.<method>')
xcheck('RfPCadj', f"{ic(min(R3[s][m]['meanV'] for s in SUBS for m in ('PC5','PC10')))}\u2013{ic(max(R3[s][m]['meanV'] for s in SUBS for m in ('PC5','PC10')))}"); xcheck('RfOnekPB', f"{ic(min(R3o[b]['PB']['meanV'] for b in R3o))}\u2013{ic(max(R3o[b]['PB']['meanV'] for b in R3o))}"); xcheck('RfOnekPBw', f"{ic(min(R3o[b]['PBw']['meanV'] for b in R3o))}\u2013{ic(max(R3o[b]['PBw']['meanV'] for b in R3o))}"); xcheck('RfOnekPC', f"{ic(min(R3o[b][m]['meanV'] for b in R3o for m in ('PC5','PC10')))}\u2013{ic(max(R3o[b][m]['meanV'] for b in R3o for m in ('PC5','PC10')))}")
# R4
Rr = R4['R4']
add('REV4', 'REV4-R4', 'R4: real-donor consequences on PsychAD: PB vs XDS-M calls and overlap, r_w and (r_PB \u2212 rho) by class, hub overlap. Descriptive; no decision.',
    'DESCRIPTIVE', f"Calls PB {ic(Rr['P1']['calls']['PB'])}/{ic(Rr['P2']['calls']['PB'])}/{ic(Rr['P3']['calls']['PB'])}; XDS-M {ic(Rr['P1']['calls']['XDSM'])}/{ic(Rr['P2']['calls']['XDSM'])}/{ic(Rr['P3']['calls']['XDSM'])}. PB-only {ic(Rr['P1']['classes']['PBonly']['n'])}/{ic(Rr['P2']['classes']['PBonly']['n'])}/{ic(Rr['P3']['classes']['PBonly']['n'])}; XDS-M-only {ic(Rr['P1']['classes']['XDSMonly']['n'])}/{ic(Rr['P2']['classes']['XDSMonly']['n'])}/{ic(Rr['P3']['classes']['XDSMonly']['n'])}; shared {ic(Rr['P1']['classes']['shared']['n'])}/{ic(Rr['P2']['classes']['shared']['n'])}/{ic(Rr['P3']['classes']['shared']['n'])}. Top-20 hub overlap {rf([Rr[s]['hub_overlap_top20'] for s in SUBS],0)}; top-5,000-pair overlap {rf([100*Rr[s]['top5000_pair_overlap']/5000 for s in SUBS],0)}%; Spearman(r_PB\u2212rho, r_w) {fx(min(Rr[s]['spearman_delta_rw'] for s in SUBS),2)} to {fx(max(Rr[s]['spearman_delta_rw'] for s in SUBS),2)}.",
    'rev4/results/REV4_SUMMARY.json: R4.<P>.{calls,classes,hub_overlap_top20,top5000_pair_overlap,spearman_delta_rw}')
xcheck('RfRealPB', ', '.join(ic(Rr[s]['calls']['PB']) for s in SUBS)); xcheck('RfRealXDSM', ', '.join(ic(Rr[s]['calls']['XDSM']) for s in SUBS)); xcheck('RfRealPBonly', ', '.join(ic(Rr[s]['classes']['PBonly']['n']) for s in SUBS)); xcheck('RfRealXonly', ', '.join(ic(Rr[s]['classes']['XDSMonly']['n']) for s in SUBS)); xcheck('RfRealShared', ', '.join(ic(Rr[s]['classes']['shared']['n']) for s in SUBS))
xcheck('RfHubOverlap', rf([Rr[s]['hub_overlap_top20'] for s in SUBS], 0)); xcheck('RfTopOverlap', rf([100 * Rr[s]['top5000_pair_overlap'] / 5000 for s in SUBS], 0) + '%')
# R5
R5 = R4['R5']
assert all(R5[s]['XDSMgu']['n'] == 20 and R5[s]['XDSMgu']['any'] <= 2 for s in SUBS) and max_runs(20) == 2
hg = sum(R5[s]['XDSMg']['any'] <= 2 for s in SUBS)
add('REV4', 'REV4-R5', 'R5 decision: XDS-Mgu holds by the run-level rule (n=20: \u22642 of 20 runs with a call) in each subtype (DC donor-centred null).',
    'PASS', f"Holds in 3 of 3: XDS-Mgu {R5['P1']['XDSMgu']['any']}/20, {R5['P2']['XDSMgu']['any']}/20, {R5['P3']['XDSMgu']['any']}/20 runs with a call. Same rule on others (not named decisions): XDS-Mg {R5['P1']['XDSMg']['any']}, {R5['P2']['XDSMg']['any']}, {R5['P3']['XDSMg']['any']} of 20 (holds in {hg} of 3; fails P2); XDS-g {R5['P1']['XDSg']['any']}, {R5['P2']['XDSg']['any']}, {R5['P3']['XDSg']['any']}; PB 20, 20, 20 (mean calls {ic(R5['P1']['PB']['meanV'])}/{ic(R5['P2']['PB']['meanV'])}/{ic(R5['P3']['PB']['meanV'])}).",
    'rev4/results/REV4_SUMMARY.json: R5.<P>.{XDSMgu,XDSMg,XDSg,PB}.{any,n,meanV}')
xcheck('RfDCguHolds', str(sum(R5[s]['XDSMgu']['any'] <= 2 for s in SUBS))); xcheck('RfDCgHolds', str(hg)); xcheck('RfDCgu', ', '.join(str(R5[s]['XDSMgu']['any']) for s in SUBS)); xcheck('RfDCg', ', '.join(str(R5[s]['XDSMg']['any']) for s in SUBS))
xcheck('RfDCPB', f"{ic(min(R5[s]['PB']['meanV'] for s in SUBS))}\u2013{ic(max(R5[s]['PB']['meanV'] for s in SUBS))}")
# R6
R6_ = R4['R6']
add('REV4', 'REV4-R6', 'Amendment 1 (R6a\u2013c, Mantis maps): donor maps, pair-level and gene-level Spearman statistics. \u201cDescriptive only; no pass/fail decision.\u201d',
    'DESCRIPTIVE', f"P1 profile norm vs n_d (Spearman): null PB {fx(R6_['null_L']['spearman_norm_ncells'],2)}, design-standardized {fx(R6_['null_R']['spearman_norm_ncells'],2)}; real PB {fx(R6_['real_L']['spearman_norm_ncells'],2)}, standardized {fx(R6_['real_R']['spearman_norm_ncells'],2)}. Pairs vs r_w: null PB {fx(R6_['pairs']['spearman_nullPB_rw'],2)}, real PB {fx(R6_['pairs']['spearman_realPB_rw'],2)}, real XDS-M {fx(R6_['pairs']['spearman_realXDSMz_rw'],2)}. Genes: null false degree vs PB-only degree {fx(R6_['genes']['spearman_nullfalse_PBonly'],2)}. Mantis-cluster Kruskal\u2013Wallis tests: not determinable from bundle.",
    'rev4/results/REV4_SUMMARY.json: R6.{null_L,null_R,real_L,real_R,pairs,genes}')
xcheck('RsixNormPB', fx(R6_['null_L']['spearman_norm_ncells'], 2)); xcheck('RsixNormR', fx(R6_['null_R']['spearman_norm_ncells'], 2)); xcheck('RsixRealNormPB', fx(R6_['real_L']['spearman_norm_ncells'], 2)); xcheck('RsixRealNormR', fx(R6_['real_R']['spearman_norm_ncells'], 2))
xcheck('RsixPairNull', fx(R6_['pairs']['spearman_nullPB_rw'], 2)); xcheck('RsixPairRealPB', fx(R6_['pairs']['spearman_realPB_rw'], 2)); xcheck('RsixPairRealX', fx(R6_['pairs']['spearman_realXDSMz_rw'], 2)); xcheck('RsixPairs', ic(R6_['pairs']['n'])); xcheck('RsixNullCalls', ic(R6_['pairs']['null_calls']))
note('REV4', 'R1 is PARTIAL: 13 (P1), 12 (P2) and 0 (P3) of the planned 20 runs per subtype are in the summary (25 of 60). The appendix list item for REV4 does not say so. Verdicts marked partial (table: *; Markdown: (partial)) rest on these runs. Because the bound is at most 2 runs for n = 12, 13 and 20 alike, a cell with \u22653 runs with a call has failed for good (FAIL partial is final for those cells); \u201cholds\u201d cells (PASS and MIXED partial) are provisional.')
note('REV4', 'Discrepancy for the authors: the appendix says \u201cLCP and LCPd fail at 25\u2013100 donors\u201d. In the R1 data LCPd holds at D=100 in both evaluated subtypes (3/25 runs) and in P1 at D=50; it fails at D=25 (20/25 runs). The statement is accurate for LCP only.')
note('REV4', 'R1 plan: permutation tests with B = 2,000,000; the summary does not record B, so it cannot be confirmed from the bundle. R2 plan item \u201cmedian library size and cells per donor by bank\u201d is not in the summary; R2c bank shifts are in R2c.<P>.bank_shift_median. R6c Mantis-cluster results are referenced by a macro (GeneMapAppendix, in appendix.tex) that is not defined in the supplied macro files.')
note('REV4', 'R2c P1 has 36 of 40 reassignments in the bundle; the XDS-M comparison in R2c P1 (14/36) therefore differs slightly from R2 P1 (15/40).')
