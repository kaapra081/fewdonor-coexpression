import math, re, collections

# ---------- LaTeX conversion
UNI = {'\u2264': r'$\le$', '\u2265': r'$\ge$', '\u00d7': r'$\times$', '\u2212': r'$-$', '\u2013': '--', '\u2014': '---', '\u2019': "'", '\u2018': '`',
       '\u201c': '``', '\u201d': "''", '\u2026': r'\ldots{}', '\u2020': r'$^\dagger$', '\u2248': r'$\approx$', '\u2192': r'$\to$', '\u03a3': r'$\Sigma$', '\u03c1': r'$\rho$', '\u00b1': r'$\pm$'}
MATHTOK = [('(R^A)^2', r'$(R^{A})^{2}$'), ('|r_w|', r'$|r_w|$'), ('r_w', r'$r_w$'), ('n_d', r'$n_d$'), ('r_PB', r'$r_{\mathrm{PB}}$'), ('\u03a3_B', r'$\Sigma_{B}$'), ('t_{D-1}', r'$t_{D-1}$'), ('R^A', r'$R^{A}$')]


def tex(s):
    ph = []

    def stash(t):
        ph.append(t); return f'@@{len(ph)-1}@@'
    for a, b in MATHTOK:
        if a in s:
            s = s.replace(a, stash(b))
    for a, b in UNI.items():
        if a in s:
            s = s.replace(a, stash(b))
    assert '\\' not in s, s
    s = s.replace('&', r'\&').replace('%', r'\%').replace('#', r'\#').replace('_', r'\_')
    s = s.replace('<', '@@LT@@').replace('>', '@@GT@@').replace('~', r'\textasciitilde{}')
    s = s.replace('@@LT@@', '$<$').replace('@@GT@@', '$>$')
    s = re.sub(r'@@(\d+)@@', lambda m: ph[int(m.group(1))], s)
    assert '@@' not in s
    assert all(ord(ch) < 128 for ch in s), [ch for ch in s if ord(ch) >= 128]
    return s


OUT_TEX = {'PASS': 'PASS', 'FAIL': 'FAIL', 'MIXED': 'MIXED', 'DESCRIPTIVE': 'DESCR.', 'NOT EVALUATED': 'NOT EVAL.', 'NOT DETERMINABLE': 'n/d',
           'PASS*': r'PASS$^{*}$', 'FAIL*': r'FAIL$^{*}$', 'MIXED*': r'MIXED$^{*}$'}
PLAN_SHORT = {
    'FIX': 'FIX\\_PREREG (real-cell null in SEA-AD; LCP2)',
    'XS': 'XS\\_PREREG (LCPX), with Amendment 1',
    'POP': 'POP\\_PREREG (PsychAD population scale)',
    'DSTD': 'DSTD\\_PREREG (LCPd, LCPXd = XDS)',
    'XU': 'XU\\_PREREG (XU)',
    'COUPLE': 'COUPLE\\_PREREG (depth-state coupling)',
    'DC': 'DC\\_PREREG (AD differential co-expression)',
    'REV': 'REV\\_PREREG (internal review 1; R1--R5, amendment R6)',
    'ONEK': 'ONEK\\_PREREG (OneK1K), with Amendment 1',
    'REV2': 'REV2\\_PREREG (internal review 2), with Amendment 1',
    'REV3': 'REV3\\_PREREG (internal review 3), with Amendment 1',
    'REV4': 'REV4\\_PREREG (external review), with Amendments 1--5',
}
W = dict(id=0.10, crit=0.285, out=0.085, ev=0.45)
CPL_C, CPL_E, LH = 29.0, 43.0, 8.3   # chars per line (crit, evidence) at \scriptsize, \linewidth = 347 pt (LNCS); line height
LEGEND = (r"Every frozen decision criterion and its outcome. "
          r"PASS/FAIL: the frozen rule was met/not met; for rules evaluated per subtype, cohort size or method cell, PASS means every cell met it and MIXED means some did (counts in Evidence). "
          r"DESCR.: reported without a frozen threshold. n/d: not determinable from the supplied files. NOT EVAL.: planned, no results available. "
          r""
          r"Run-level bound: fraction of runs with a BH call $\le 0.05+2\,\mathrm{SE}$, SE binomial at $p=0.05$ (at most 9 of 100, 4 of 40, 2 of 20 runs). "
          r"Numbers are taken from the result files; sources are listed in CRITERIA.md.")


def row_height(r):
    nl = max(math.ceil(len(r['crit']) / CPL_C), math.ceil(len(r['ev']) / CPL_E), 1)
    return nl * LH + 3.0


def _pack(rows, budget):
    groups = collections.OrderedDict()
    for r in rows: groups.setdefault(r['plan'], []).append(r)
    tables, cur, h = [], [], 0.0
    first_extra = 70.0
    for pk, rs in groups.items():
        hdr_h = 14.0
        for i, r in enumerate(rs):
            need = row_height(r) + (hdr_h if (i == 0 or not cur or cur[-1][0] != pk) else 0)
            lim = budget - (first_extra if not tables else 0)
            if cur and h + need > lim:
                tables.append(cur); cur, h = [], 0.0
                need = row_height(r) + hdr_h
            cur.append((pk, r)); h += need
    if cur: tables.append(cur)
    return tables


def build_tex(rows, budget_max=485.0):
    k_min = len(_pack(rows, budget_max)); tables = None
    for b_ in [x * 5.0 for x in range(60, int(budget_max / 5) + 1)]:
        t_ = _pack(rows, b_)
        if len(t_) == k_min: tables = t_; break
    out = ["% criteria_table.tex -- generated; every number is read from the files in criteria_inputs.tgz (see CRITERIA.md for file/key sources).",
           "% Requires: \\usepackage{booktabs,array}. Uses no paper macros. Each table is sized to fit one page at \\linewidth ~ 347 pt (single-column LNCS).",
           "% Not compiled in the build sandbox (no TeX installation); checked by a static lint only (balanced braces, 4 cells per row, ASCII only, whitelisted commands)."]
    for ti, tb in enumerate(tables):
        lab = 'tab:criteria' if ti == 0 else f'tab:criteria{ti+1}'
        cap = LEGEND if ti == 0 else r"Every frozen decision criterion and its outcome (continued; part %d of %d)." % (ti + 1, len(tables))
        out += ["", r"\begin{table}[p]", r"\centering\scriptsize", r"\setlength{\tabcolsep}{2.5pt}\renewcommand{\arraystretch}{1.1}",
                r"\caption{%s}\label{%s}" % (cap, lab),
                r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{%.3f\linewidth}>{\raggedright\arraybackslash}p{%.3f\linewidth}>{\raggedright\arraybackslash}p{%.3f\linewidth}>{\raggedright\arraybackslash}p{%.3f\linewidth}@{}}" % (W['id'], W['crit'], W['out'], W['ev']),
                r"\toprule", r"Plan & Criterion & Outcome & Evidence \\", r"\midrule"]
        lastp = None
        for pk, r in tb:
            if pk != lastp:
                if lastp is not None: out.append(r"\addlinespace[3pt]")
                cont = ' (continued)' if (pk in [x[0] for t in tables[:ti] for x in t]) and tb[0][0] == pk else ''
                out.append(r"\multicolumn{4}{@{}>{\raggedright\arraybackslash}p{0.90\linewidth}}{\textbf{%s}%s} \\" % (PLAN_SHORT[pk], cont))
                out.append(r"\cmidrule(lr){1-4}")
                lastp = pk
            out.append(r"%s & %s & %s & %s \\" % (tex(r['id']), tex(r['crit']), OUT_TEX[r['out']], tex(r['ev'])))
        out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return '\n'.join(out) + '\n', tables


def lint_tex(txt):
    probs = []
    body = '\n'.join(l for l in txt.split('\n') if not l.startswith('%'))
    if body.count('{') != body.count('}'): probs.append(('brace count', body.count('{'), body.count('}')))
    # ignore escaped braces (none expected)
    b = len(re.findall(r'\\begin\{(\w+)\}', body)); e = len(re.findall(r'\\end\{(\w+)\}', body))
    if b != e: probs.append(('begin/end', b, e))
    allowed = {'begin', 'end', 'centering', 'scriptsize', 'setlength', 'tabcolsep', 'renewcommand', 'arraystretch', 'caption', 'label', 'tabular', 'raggedright', 'arraybackslash', 'linewidth',
               'toprule', 'midrule', 'bottomrule', 'cmidrule', 'addlinespace', 'multicolumn', 'textbf', 'le', 'ge', 'times', 'ldots', 'dagger', 'approx', 'to', 'Sigma', 'rho', 'pm', 'mathrm', 'textasciitilde',
               'table', 'p', 'quad', ',', '&', '%', '_', '#', '\\', '$', ' ', 'it'}
    used = set(re.findall(r'\\([A-Za-z]+|[^A-Za-z])', body))
    if used - allowed: probs.append(('commands', sorted(used - allowed)))
    # row structure: every body row ends with \\ and has exactly 3 unescaped ampersands (except multicolumn rows with none)
    for ln in body.split('\n'):
        if ln.endswith(r'\\') and not ln.startswith(r'\multicolumn') and ln not in (r'Plan & Criterion & Outcome & Evidence \\',):
            amp = len(re.findall(r'(?<!\\)&', ln))
            if amp != 3: probs.append(('ampersands', amp, ln[:80]))
    if any(ord(ch) > 127 for ch in body): probs.append(('non-ascii',))
    dollars = body.count('$') - body.count('\\$')
    if dollars % 2: probs.append(('odd dollars', dollars))
    for ln in body.split('\n'):
        if (ln.count('$') % 2): probs.append(('odd $ in line', ln[:90]))
    return probs


OUT_MD = {'PASS': 'PASS', 'FAIL': 'FAIL', 'MIXED': 'MIXED', 'DESCRIPTIVE': 'DESCRIPTIVE', 'NOT EVALUATED': 'NOT EVALUATED', 'NOT DETERMINABLE': 'NOT DETERMINABLE FROM BUNDLE',
          'PASS*': 'PASS (partial)', 'FAIL*': 'FAIL (partial)', 'MIXED*': 'MIXED (partial)'}


def md_cell(s):
    return s.replace('|', '\\|').replace('<', '&lt;').replace('>', '&gt;').replace('\n', ' ')


def build_md(rows, plans, notes, xcheck_rows, extra_sections):
    L = []
    L.append('# CRITERIA: every frozen decision criterion and its outcome\n')
    L.append('Built from `criteria_inputs.tgz` (prereg/, rev2/, rev3/, rev4/, results/, paper2/). Every number below is read from a file in that bundle; the file and key are given in the Source column. '
             'Where a result file stores a decision flag, the verdict was recomputed from the underlying numbers and asserted equal to the flag. Plans are in the order FIX, XS, POP, DSTD, XU, COUPLE, DC, REV, ONEK, REV2, REV3, REV4.\n')
    L.append('## Legend\n')
    L.append('- **PASS / FAIL**: the frozen rule was met / not met. For rules evaluated per subtype, cohort size or method cell, PASS means every cell met it.')
    L.append('- **MIXED**: the rule was evaluated in several cells and some met it, some did not (counts in Evidence). This label is an addition to the PASS / FAIL / NOT EVALUATED / DESCRIPTIVE set, used only where a single PASS or FAIL would hide the split.')
    L.append('- **DESCRIPTIVE**: the plan reports the quantity without a frozen threshold (no criterion). **NOT EVALUATED**: planned, no results in the bundle. **NOT DETERMINABLE FROM BUNDLE**: the outcome cannot be read from the supplied files.')
    pass
    L.append('- **Run-level rule**: a method holds in a cell if the fraction of runs with at least one BH call (q = 0.05) is at most 0.05 + 2 SE, SE binomial at p = 0.05 (REV3 and REV4 plans state this; it gives at most 9 of 100, 4 of 40, 2 of 20, 1 of 10, 4 of 33 runs).')
    L.append('- **MCSE**: Monte Carlo standard error as stored in the result files; **top-1%**: type-I over the 1% of pairs most correlated within pseudo-donors; **SE** in top-1% bounds is the stored binomial SE (pairs x reassignments).\n')
    cnt = collections.Counter(r['out'] for r in rows)
    L.append('## Summary\n')
    L.append(f"{len(rows)} criteria rows. " + '; '.join(f"{OUT_MD[k]}: {v}" for k, v in sorted(cnt.items(), key=lambda kv: -kv[1])) + '.\n')
    per = collections.OrderedDict()
    for r in rows: per.setdefault(r['plan'], collections.Counter())[r['out']] += 1
    L.append('| Plan | Rows | ' + ' | '.join(OUT_MD[k] for k in ('PASS', 'FAIL', 'MIXED', 'DESCRIPTIVE', 'NOT EVALUATED', 'NOT DETERMINABLE', 'PASS*', 'FAIL*', 'MIXED*')) + ' |')
    L.append('|---|---|' + '---|' * 9)
    for pk, c in per.items():
        L.append(f"| {pk} | {sum(c.values())} | " + ' | '.join(str(c.get(k, 0) or '') for k in ('PASS', 'FAIL', 'MIXED', 'DESCRIPTIVE', 'NOT EVALUATED', 'NOT DETERMINABLE', 'PASS*', 'FAIL*', 'MIXED*')) + ' |')
    L.append('')
    for pk, pl in plans.items():
        L.append(f"## {pl['title']}\n")
        L.append(f"- Files: {', '.join('`'+f+'`' for f in pl['files'])}")
        L.append(f"- Hash: {pl['hashnote']}")
        L.append(f"- Scope in bundle: {pl['scope']}\n")
        L.append('| ID | Criterion | Outcome | Evidence | Source (file: key) |')
        L.append('|---|---|---|---|---|')
        for r in rows:
            if r['plan'] == pk:
                L.append(f"| {r['id']} | {md_cell(r['crit'])} | **{OUT_MD[r['out']]}** | {md_cell(r['ev'])} | `{r['src']}` |")
        L.append('')
        ns = [t for p_, t in notes if p_ == pk]
        if ns:
            L.append('Notes and deviations:\n')
            for t in ns: L.append(f"- {t}")
            L.append('')
    for h, body in extra_sections:
        L.append(f"## {h}\n"); L.append(body); L.append('')
    L.append('## Cross-check of the paper macros against the bundle\n')
    ok_ = [x for x in xcheck_rows if x[3]]; bad_ = [x for x in xcheck_rows if not x[3]]
    L.append(f"{len(xcheck_rows)} macro values (nums.tex, rev2nums.tex, rev4nums.tex) were compared with the value recomputed from the result files for this table: {len(ok_)} agree after normalising formatting, {len(bad_)} differ.\n")
    if bad_:
        L.append('| Macro | In macro file | Recomputed | Status |'); L.append('|---|---|---|---|')
        for m, a, b_, k in bad_: L.append(f"| `{m}` | `{a}` | `{b_}` | {k} |")
        L.append('')
    L.append('<details><summary>All compared macros</summary>\n')
    L.append('| Macro | In macro file | Recomputed |'); L.append('|---|---|---|')
    for m, a, b_, k in xcheck_rows: L.append(f"| `{m}` | `{a}` | `{b_}` |")
    L.append('\n</details>\n')
    return '\n'.join(L)
