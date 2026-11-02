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

# ---------------------------------------------------------------- REV4 (data-driven; final runs)
plan('REV4', 'REV4_PREREG (external review) with Amendments 1-5', ['rev4/REV4_PREREG.md'] + [f'rev4/REV4_AMENDMENT{i}.md' for i in range(1, 6)],
     'sha256 (first 12): ' + ', '.join(f"{n} {sha(f)[:12]}" for n, f in [('PREREG', 'rev4/REV4_PREREG.md')] + [(f'A{i}', f'rev4/REV4_AMENDMENT{i}.md') for i in range(1, 6)]),
     'R1 20 runs per subtype; R2/R2c 40; R3 20/10; R5 20; C1 20; C2 40; C4 20; R6, C3, C5, C6, A4, A5 descriptive or single-run.')
R4 = J('rev4/results/REV4_SUMMARY.json'); R1_4 = R4['R1']
nR1 = {s: R1_4[s]['n'] for s in SUBS}
PART = any(n < 20 for n in nR1.values())
def r1cells(m): return [(s, D, R1_4[s][D][m]['any'], R1_4[s][D][m]['n']) for D in ('25', '50', '100') for s in SUBS if nR1[s]]
def r1h(m): return sum(ok_rl(k, n) for _, _, k, n in r1cells(m))
def r1tot(m, D): return sum(R1_4[s][D][m]['any'] for s in SUBS if nR1[s]), sum(R1_4[s][D][m]['n'] for s in SUBS if nR1[s])
def r1txt(m): return '; '.join(f"D={D}: {r1tot(m, D)[0]}/{r1tot(m, D)[1]}" for D in ('25', '50', '100'))
def r1fail(m): return ', '.join(f"{s} D={D} ({k}/{n})" for s, D, k, n in r1cells(m) if not ok_rl(k, n)) or 'none'
def r1out(m):
    h, t = r1h(m), len(r1cells(m)); o = 'PASS' if h == t else ('FAIL' if h == 0 else 'MIXED'); return o + ('*' if PART else '')
for rid, m, lab in (('REV4-R1a', 'PB', 'pseudobulk'), ('REV4-R1b', 'LCP', 'LCP'), ('REV4-R1c', 'LCPd', 'LCPd'), ('REV4-R1d', 'XDSperm', 'XDS, permutation p-values (B = 2,000,000)'),
                    ('REV4-R1e', 'XDS', 'XDS, normal-moment reference'), ('REV4-R1f', 'XDSM', 'XDS-M (S = 5)'), ('REV4-R1g', 'KX', 'KX')):
    add('REV4', rid, f"R1, {lab}: at most 2 of 20 runs with a BH call, per subtype and D = 25, 50, 100.", r1out(m),
        f"Holds in {r1h(m)} of {len(r1cells(m))} cells; runs with a call, pooled over subtypes: {r1txt(m)}. Failing cells: {r1fail(m)}.", f"rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.{m}.{{any,n}}")
for D, w in (('25', 'TwentyFive'), ('50', 'Fifty'), ('100', 'Hundred')):
    for m in ('LCP', 'LCPd', 'XDSperm', 'XDSM', 'KX'):
        k, n = r1tot(m, D); xcheck(f'Rf{m}{w}', f"{k}/{n}")
# R2
R2 = R4['R2']
def r2c(m): return [(R2[s][m]['any'], R2[s][m]['n']) for s in SUBS]
oth = lambda m: ', '.join(f"{k}/{n}" for k, n in r2c(m))
hxm = sum(ok_rl(k, n) for k, n in r2c('XDSM'))
add('REV4', 'REV4-R2', 'R2 decision: XDS-M within the run-level bound (n = 40: at most 4 runs with a call) in each subtype, brain-bank-stratified null.',
    'PASS' if hxm == 3 else ('FAIL' if hxm == 0 else 'MIXED'), f"Holds in {hxm} of 3 (XDS-M {oth('XDSM')}). Other methods (same rule, not named decisions): PB {oth('PB')}; LCP {oth('LCP')}; LCPd {oth('LCPd')}; XDS {oth('XDS')}; KX {oth('KX')}.",
    'rev4/results/REV4_SUMMARY.json: R2.<P>.<method>.{any,n}')
xcheck('RfBankXDSM', ', '.join(f"{k}/{n}" for k, n in r2c('XDSM'))); xcheck('RfBankKX', ', '.join(f"{k}/{n}" for k, n in r2c('KX')))
# R2c (Amendment 2)
R2c = R4['R2c']; xb = [(R2c[s]['XDSMb']['any'], R2c[s]['XDSMb']['n']) for s in SUBS]; hb = sum(ok_rl(k, n) for k, n in xb)
add('REV4', 'REV4-R2c', 'Amendment 2 (post hoc; rule written before R2c ran): bank-centred XDS-Mb within the bound in every subtype attributes R2 failures to residual bank means.',
    'PASS' if hb == 3 else 'MIXED', f"XDS-Mb {', '.join(f'{k}/{n}' for k, n in xb)} (holds in {hb} of 3); bank-centred PB {', '.join(str(R2c[s]['PBb']['any']) + '/' + str(R2c[s]['PBb']['n']) for s in SUBS)} runs with a call.",
    'rev4/results/REV4_SUMMARY.json: R2c.<P>.{XDSMb,PBb}.{any,n}')
xcheck('RfBankXDSMb', ', '.join(f"{k}/{n}" for k, n in xb)); xcheck('RfBankXDSMbHolds', str(hb))
# R3, R4, R5 (unchanged logic)
R3 = R4['R3']; R3o = R4['R3o']
add('REV4', 'REV4-R3', 'R3: PC-adjusted (k = 5, 10), residualized and weighted pseudobulk, PsychAD (20 runs) and OneK1K (10 runs). Descriptive.',
    'DESCRIPTIVE', f"Every variant made at least one call in every run. PsychAD mean calls: PBres {ic(min(R3[s]['PBres']['meanV'] for s in SUBS))}\u2013{ic(max(R3[s]['PBres']['meanV'] for s in SUBS))}; PBw {ic(min(R3[s]['PBw']['meanV'] for s in SUBS))}\u2013{ic(max(R3[s]['PBw']['meanV'] for s in SUBS))}; PC5/PC10 {ic(min(R3[s][m]['meanV'] for s in SUBS for m in ('PC5','PC10')))}\u2013{ic(max(R3[s][m]['meanV'] for s in SUBS for m in ('PC5','PC10')))}. OneK1K PB {ic(min(R3o[b]['PB']['meanV'] for b in R3o))}\u2013{ic(max(R3o[b]['PB']['meanV'] for b in R3o))}; PBw {ic(min(R3o[b]['PBw']['meanV'] for b in R3o))}\u2013{ic(max(R3o[b]['PBw']['meanV'] for b in R3o))}.",
    'rev4/results/REV4_SUMMARY.json: R3.<P>.<method>; R3o.<B>.<method>')
Rr = R4['R4']
add('REV4', 'REV4-R4', 'R4: real-donor pseudobulk and XDS-M calls, overlap and hubs, PsychAD. Descriptive.', 'DESCRIPTIVE',
    f"Calls PB {'/'.join(ic(Rr[s]['calls']['PB']) for s in SUBS)}; XDS-M {'/'.join(ic(Rr[s]['calls']['XDSM']) for s in SUBS)}; PB-only {'/'.join(ic(Rr[s]['classes']['PBonly']['n']) for s in SUBS)}; XDS-M-only {'/'.join(ic(Rr[s]['classes']['XDSMonly']['n']) for s in SUBS)}.",
    'rev4/results/REV4_SUMMARY.json: R4.<P>.{calls,classes}')
xcheck('RfRealPB', ', '.join(ic(Rr[s]['calls']['PB']) for s in SUBS)); xcheck('RfRealXDSM', ', '.join(ic(Rr[s]['calls']['XDSM']) for s in SUBS))
R5 = R4['R5']; hg = sum(ok_rl(R5[s]['XDSMg']['any'], R5[s]['XDSMg']['n']) for s in SUBS); hgu = sum(ok_rl(R5[s]['XDSMgu']['any'], R5[s]['XDSMgu']['n']) for s in SUBS)
add('REV4', 'REV4-R5', 'R5 decision: XDS-Mgu within the run-level bound (n = 20: at most 2 runs with a call) in each subtype, differential co-expression null.',
    'PASS' if hgu == 3 else 'MIXED', f"XDS-Mgu {', '.join(str(R5[s]['XDSMgu']['any']) for s in SUBS)} of 20 (holds in {hgu} of 3); XDS-Mg {', '.join(str(R5[s]['XDSMg']['any']) for s in SUBS)} (holds in {hg}); PB 20 of 20 in each.",
    'rev4/results/REV4_SUMMARY.json: R5.<P>.{XDSMgu,XDSMg,PB}')
xcheck('RfDCgu', ', '.join(str(R5[s]['XDSMgu']['any']) for s in SUBS)); xcheck('RfDCg', ', '.join(str(R5[s]['XDSMg']['any']) for s in SUBS))
# R6 (Amendment 1)
R6_ = R4['R6']; mc = J('rev4/out/R6/R6_mantis_clusters.json')
add('REV4', 'REV4-R6', 'Amendment 1 (Mantis maps): donor and gene maps; Kruskal-Wallis of real pseudobulk-only and XDS-M degree across leaf clusters. Descriptive.', 'DESCRIPTIVE',
    f"Null PB profile length vs cells, Spearman {fx(R6_['null_L']['spearman_norm_ncells'],2)} (design-standardized {fx(R6_['null_R']['spearman_norm_ncells'],2)}). Kruskal-Wallis p: PB-only degree {mc['kw']['PB_only_degree']['p']:.2g}, XDS-M degree {mc['kw']['XDSM_degree']['p']:.2g}, null false degree (added) {mc['kw']['null_false_degree']['p']:.2g}; top PB-only leaf: {mc['top_PBonly_cluster']['label']}.",
    'rev4/results/REV4_SUMMARY.json: R6; rev4/out/R6/R6_mantis_clusters.json: kw, top_PBonly_cluster')
# Amendment 3
C = R4['C']
def cell3(blk, m): return [(blk[s][m]['any'], blk[s][m]['n']) for s in SUBS if s in blk]
c1 = cell3(C['C1'], 'XDSMc'); h1 = sum(ok_rl(k, n) for k, n in c1)
add('REV4', 'REV4-C1', 'Amendment 3 C1: covariate-adjusted XDS-Mc on a covariate-carrying donor-centred null; at most 2 of 20 runs with a call, each subtype.',
    'PASS' if h1 == 3 else 'MIXED', f"XDS-Mc {', '.join(f'{k}/{n}' for k, n in c1)} (holds in {h1} of 3); covariate-adjusted PB {', '.join(str(C['C1'][s]['PBc']['any']) + '/' + str(C['C1'][s]['PBc']['n']) for s in SUBS)}. Real donors: XDS-Mc {'/'.join(ic(C['C1real'][s]['XDSMc']) for s in SUBS)} calls, XDS-M {'/'.join(ic(C['C1real'][s]['XDSM']) for s in SUBS)}.",
    'rev4/results/REV4_SUMMARY.json: C.C1.<P>.{XDSMc,PBc}; C.C1real.<P>')
c2 = cell3(C['C2'], 'XDSM'); h2 = sum(ok_rl(k, n) for k, n in c2)
add('REV4', 'REV4-C2', 'Amendment 3 C2: bank-mean-removed bank-stratified null; XDS-M at most 4 of 40 in every subtype attributes R2 failures to residual bank means.',
    'PASS' if h2 == 3 else 'MIXED', f"XDS-M {', '.join(f'{k}/{n}' for k, n in c2)} (holds in {h2} of 3); PB {', '.join(str(C['C2'][s]['PB']['any']) + '/' + str(C['C2'][s]['PB']['n']) for s in SUBS)}.",
    'rev4/results/REV4_SUMMARY.json: C.C2.<P>.{XDSM,PB}')
xcheck('BmX', ', '.join(f"{k}/{n}" for k, n in c2)); xcheck('CoNullXc', ', '.join(f"{k}/{n}" for k, n in c1))
C3 = C['C3']
add('REV4', 'REV4-C3', 'Amendment 3 C3: noise baselines for real-donor overlaps (top 5,000 pairs; 20 hubs). Descriptive.', 'DESCRIPTIVE',
    f"Top-5,000 shared: XDS-M seed vs seed {rf([C3[s]['XDSM_seed_vs_seed']['top5000'] for s in SUBS],0)}; PB vs XDS-M {rf([C3[s]['PB_vs_XDSM']['top5000'] for s in SUBS],0)}; XDS-M halves {rf([C3[s]['XDSM_half_vs_half']['top5000'] for s in SUBS],0)}; PB halves {rf([C3[s]['PB_half_vs_half']['top5000'] for s in SUBS],0)}.",
    'rev4/results/REV4_SUMMARY.json: C.C3.<P>')
c4 = C.get('C4', {}).get('P1')
if c4:
    add('REV4', 'REV4-C4', 'Amendment 3 C4: permutation XDS (2,000,000 permutations) at 200 donors, Homeo_CECR2; at most 2 of 20 runs with a call.',
        'PASS' if ok_rl(c4['XDS_perm']['any'], c4['XDS_perm']['n']) else 'FAIL', f"XDS (perm.) {c4['XDS_perm']['any']}/{c4['XDS_perm']['n']}; XDS (normal) {c4['XDS_normal']['any']}/{c4['XDS_normal']['n']}; LCPd (perm.) {c4['LCPd_perm']['any']}/{c4['LCPd_perm']['n']}.",
        'rev4/results/REV4_SUMMARY.json: C.C4.P1')
else:
    add('REV4', 'REV4-C4', 'Amendment 3 C4: permutation XDS at 200 donors, at most 2 of 20 runs with a call.', 'NOT EVALUATED', 'No C4 results in the bundle.', 'rev4/out/C4_P1_k*.json')
C5 = C['C5']
add('REV4', 'REV4-C5', 'Amendment 3 C5: depth-state coupling statistic (median over genes; fraction with |c| > 0.045). Descriptive.', 'DESCRIPTIVE',
    '; '.join(f"{s} {C5[s]['coupling_median']:.3f} ({100*C5[s]['frac_abs_gt_0045']:.0f}%)" for s in ('P1', 'P2', 'P3', 'B1', 'B2', 'B3') if s in C5) + '.', 'rev4/results/REV4_SUMMARY.json: C.C5.<sub>')
rep_ = C.get('C6_rep', {}); perm = C.get('C6_perm', {})
add('REV4', 'REV4-C6', 'Amendment 3 C6: eta-squared of gene degree across Mantis leaves vs random partitions (P1); heat-shock leaf rank in P2, P3. Descriptive.', 'DESCRIPTIVE',
    f"Null false degree eta2 {perm['null_false_degree']['eta2']:.2f} (permutation median {perm['null_false_degree']['null_median']:.2f}, p {perm['null_false_degree']['perm_p']:.3f}); heat-shock leaf rank {', '.join(str(v['hs_rank']) + ' of ' + str(v['n_leaves']) for v in rep_.values())} in P2, P3.",
    'rev4/results/REV4_SUMMARY.json: C.C6_perm, C.C6_rep')
add('REV4', 'REV4-C7', 'Amendment 3 C7 (computed shortly before the amendment): overlap of REV3 runs with an XDS-M call at 25 donors across subtypes. Descriptive.', 'DESCRIPTIVE',
    'Runs with a call share 24, 22 and 27 indices between subtype pairs and 13 across all three (49 of 100 in each); independence predicts 24 and 12.', 'computed from rev3 per-run files; stated in rev4/REV4_AMENDMENT3.md')
# Amendment 4
A4 = J('rev4/out/A4/A4_summary.json'); h4 = sum(A4[s]['leaf']['sets']['null_PB']['p_Q'] < 0.01 for s in SUBS)
z4 = ', '.join(format(A4[s]['leaf']['sets']['null_PB']['Q_z'], '.1f') for s in SUBS)
add('REV4', 'REV4-A4', 'Amendment 4 (post hoc): null pseudobulk edges function-coherent if leaf-partition modularity p < 0.01 in at least 2 of 3 subtypes.',
    'PASS' if h4 >= 2 else 'FAIL', f"Met in {h4} of 3: z {z4}.", 'rev4/out/A4/A4_summary.json: <P>.leaf.sets.null_PB.{Q_z,p_Q}')
# Amendment 5
A5 = J('rev4/out/A4/A5_summary.json'); ru = A5['rule']
add('REV4', 'REV4-A5a', 'Amendment 5 (post hoc) E1(a): expression-stratified relabelling, null PB leaf modularity p < 0.01 in at least 2 of 3 subtypes.',
    'PASS' if ru['a'] >= 2 else 'FAIL', f"Met in {ru['a']} of 3: z {', '.join(format(A5[s]['leaf']['null_PB']['Q_z'], '.1f') for s in SUBS)}.", 'rev4/out/A4/A5_summary.json: <P>.leaf.null_PB')
add('REV4', 'REV4-A5b', 'E1(b): mean stratified z over 20 null reassignments above 2.33 in at least 2 of 3 subtypes.',
    'PASS' if ru['b'] >= 2 else 'FAIL', f"Met in {ru['b']} of 3: mean z {', '.join(format(A5[s]['null20']['z_mean'], '.1f') for s in SUBS)}; p < 0.01 in {', '.join(str(A5[s]['null20']['n_p_lt_001']) for s in SUBS)} of 20.", 'rev4/out/A4/A5_summary.json: <P>.null20')
add('REV4', 'REV4-A5d', 'E1(d): GO named only if shared-term enrichment of null PB edges has p < 0.01 in at least 2 of 3 subtypes.',
    'PASS' if ru['d'] >= 2 else 'FAIL', f"Met in {ru['d']} of 3: ratio {', '.join(format(A5[s]['GO']['sets']['null_PB']['ratio'], '.2f') for s in SUBS)}, p {', '.join(format(A5[s]['GO']['sets']['null_PB']['p'], '.3f') for s in SUBS)}. GO is therefore not named in the abstract.", 'rev4/out/A4/A5_summary.json: <P>.GO.sets.null_PB')
sx = {s: J(f'rev4/out/A5sex_{s}.json') for s in SUBS}
add('REV4', 'REV4-A5e2', 'E2: sex positive control; XIST-Y and Y-Y pairs under XDS-M and covariate-adjusted XDS-Mc on real donors. Descriptive.', 'DESCRIPTIVE',
    f"XIST-Y called (of 7): XDS-M {', '.join(str(sx[s]['XIST_Y']['called_XDSM']) for s in SUBS)}, XDS-Mc {', '.join(str(sx[s]['XIST_Y']['called_XDSMc']) for s in SUBS)}; Y-Y (of 21): {', '.join(str(sx[s]['Y_Y']['called_XDSM']) for s in SUBS)} and {', '.join(str(sx[s]['Y_Y']['called_XDSMc']) for s in SUBS)}.", 'rev4/out/A5sex_<P>.json')
add('REV4', 'REV4-A5e3', 'E3: odds ratio that a real PB-only call is a false call in the first null reassignment. Descriptive.', 'DESCRIPTIVE',
    f"OR {', '.join(format(A5[s]['PBonly_null_OR']['OR'], '.2f') for s in SUBS)}.", 'rev4/out/A4/A5_summary.json: <P>.PBonly_null_OR')
note('REV4', 'R1 used B = 2,000,000 donor permutations (code default; the override variable was unset in the running processes, recorded in Amendment 5).')
