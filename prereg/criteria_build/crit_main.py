"""Rebuild criteria_table.tex and CRITERIA.md from criteria_inputs.tgz.
Usage (from an empty directory containing these scripts and the .tgz):  python crit_main.py criteria_inputs.tgz
Extracts the bundle to work/, runs crit_rows_a.py and crit_rows_b.py (row definitions; every number is computed from the bundle files),
then writes out/criteria_table.tex and out/CRITERIA.md. Extra prose sections come from crit_extra.json."""
import sys, os, json, re, tarfile, collections
tgz = sys.argv[1] if len(sys.argv) > 1 else 'criteria_inputs.tgz'
if not os.path.isdir('work'):
    os.makedirs('work'); tarfile.open(tgz).extractall('work')
os.makedirs('out', exist_ok=True)
exec(open('crit_common.py').read()); exec(open('crit_rows_a.py').read()); exec(open('crit_rows_b.py').read()); exec(open('crit_emit.py').read())
assert not LONG, LONG
MAC = read_macros()
def nz2(s):
    s = norm_macro(s).replace('\u2212', '-').replace('$', '').replace('\\times', 'x').replace('^{', 'e').replace('}', '')
    m = re.fullmatch(r'(\d)x10e(-?\d+)', s)
    if m: s = f"{m.group(1)}e{int(m.group(2))}"
    m = re.fullmatch(r'(\d)e(-?)0*(\d+)', s)
    if m: s = f"{m.group(1)}e{m.group(2)}{m.group(3)}"
    return s
xrows = [(m, MAC[m], v, nz2(MAC[m]) == nz2(v)) for m, v in XCHECK if m in MAC]
assert all(x[3] for x in xrows), [x for x in xrows if not x[3]]
extra = [tuple(x) for x in json.load(open('crit_extra.json'))]
tex_txt, tables = build_tex(ROWS); probs = lint_tex(tex_txt); assert not probs, probs
open('out/criteria_table.tex', 'w').write(tex_txt)
open('out/CRITERIA.md', 'w').write(build_md(ROWS, PLANS, NOTES, xrows, extra))
print(len(ROWS), 'rows;', len(tables), 'tables;', len(xrows), 'macro checks all agree')
