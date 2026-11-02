"""Shared helpers for building the frozen-criteria table from the input bundle (work/)."""
import json, math, glob, os, re, hashlib
from collections import OrderedDict

ROOT = 'work/'


def J(p):
    with open(ROOT + p) as fh:
        return json.load(fh)


def bnd(n):
    """Run-level bound 0.05 + 2 SE, SE binomial at p = 0.05 (REV3/REV4 plans: 'SE at p = 0.05')."""
    return 0.05 + 2 * math.sqrt(0.05 * 0.95 / n)


def ok_rl(k, n):
    return k / n <= bnd(n) + 1e-12


def max_runs(n):
    return int(math.floor(bnd(n) * n + 1e-9))


def fx(x, n=3):
    return f"{x:.{n}f}".replace('-', '\u2212')


def ic(x):
    return f"{int(round(x)):,}"


def sha(p):
    with open(ROOT + p, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


ROWS = []
LONG = []
PLANS = OrderedDict()
NOTES = []        # (plan, text) deviations / flags
XCHECK = []       # (macro, computed string)


def plan(key, title, files, hashnote, scope):
    PLANS[key] = dict(key=key, title=title, files=files, hashnote=hashnote, scope=scope)


def add(pl, rid, crit, out, ev, src):
    if len(crit.split()) > 25: LONG.append((rid, len(crit.split()), crit))
    assert out in ('PASS', 'FAIL', 'MIXED', 'DESCRIPTIVE', 'NOT EVALUATED', 'NOT DETERMINABLE', 'PASS*', 'FAIL*', 'MIXED*'), out
    ROWS.append(dict(plan=pl, id=rid, crit=crit, out=out, ev=ev, src=src))


def note(pl, text):
    NOTES.append((pl, text))


def xcheck(macro, computed):
    XCHECK.append((macro, computed))


def rng(vals, f=lambda x: x):
    lo, hi = min(vals), max(vals)
    return f(lo) if f(lo) == f(hi) else f"{f(lo)}\u2013{f(hi)}"


def read_macros():
    mac = {}
    for fn in ('paper2/nums.tex', 'paper2/rev2nums.tex', 'paper2/rev4nums.tex'):
        txt = open(ROOT + fn).read()
        for m in re.finditer(r'\\newcommand\{\\(\w+)\}\{(.*)\}\s*$', txt, flags=re.M):
            mac[m.group(1)] = m.group(2)
    return mac


def norm_macro(s):
    s = s.replace('{,}', ',').replace('\\%', '%').replace('--', '\u2013').replace('$-$', '-').replace('\\ ', ' ')
    s = s.replace('\u2212', '-')
    return re.sub(r'\s+', ' ', s).strip()


def rf(vals, nd):
    lo, hi = min(vals), max(vals)
    a, b = f"{lo:.{nd}f}", f"{hi:.{nd}f}"
    return a if a == b else f"{a}\u2013{b}"
