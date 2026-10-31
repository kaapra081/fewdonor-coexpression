"""Convert the committed interim snapshot rev2/interim/A_interim_0615UTC.txt into the JSON layout that paper2/make_rev2.py reads.

    python rev2/interim_to_summary.py            # writes rev2/interim/REV2_SUMMARY_interim.json
    python paper2/make_rev2.py rev2/interim/REV2_SUMMARY_interim.json

The snapshot lists, per subtype and cohort size, "any/n" (runs with at least one BH call / runs), mean calls and mean per-pair type-I
for every method; only 'any', 'n', 'meanV' and 'mean_t05' are available from it (one decimal for meanV). The final summary comes
from rev2/rev2_summary.py on the full per-run outputs (rev2/out/*.json), which are not in the repository (see repro_notes in README).
"""
import os, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); src = os.path.join(HERE, "interim", "A_interim_0615UTC.txt")
S = {"A": {}, "E": {}}
for line in open(src):
    f = line.split()
    if len(f) < 3 or not re.match(r"P[123]$", f[0]): continue
    cell = S["A"].setdefault(f[0], {}).setdefault(f[1], {})
    for tok in f[2:]:
        m, rest = tok.split(":"); k_n, mv, t = rest.split(",")
        k, n = k_n.split("/"); cell[m] = dict(any=int(k), n=int(n), meanV=float(mv.split("=")[1]), mean_t05=float(t.split("=")[1]))
out = os.path.join(HERE, "interim", "REV2_SUMMARY_interim.json"); json.dump(S, open(out, "w"), indent=1); print("wrote", out)
