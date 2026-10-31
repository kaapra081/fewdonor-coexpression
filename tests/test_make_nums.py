"""Checks that paper2/make_nums.py runs from the committed layout and never degrades silently (synthetic inputs, ~10 s).

test_clean_clone_regenerates_committed_files compares the regenerated nums.tex and tables with the committed ones, so it
fails if nums.tex or a table is edited by hand instead of by make_nums.py (re-run `python paper2/make_nums.py`, commit the output).
"""
import os, sys, shutil, subprocess, filecmp
import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COPY = ("paper2", "dev", "prereg", "results", "rev", "onek", "dc")      # what make_nums.py reads, in either layout
GENERATED = ("nums.tex", "table_sim.tex", "table_null.tex", "r1_summary.json")


def copy_repo(dst):
    """Copy the files make_nums.py reads (skipping the 14 MB rare/*.csv.gz and __pycache__)."""
    for d in COPY:
        shutil.copytree(os.path.join(REPO, d), os.path.join(dst, d), ignore=shutil.ignore_patterns("*.csv.gz", "__pycache__", "*.pdf", "*.png"))
    return dst


def add_synthetic_raw(root, P=3000, seed=0):
    """Write the gitignored per-run arrays make_nums.py would read on the authors' machine (random p-values; shapes and keys as in the scripts)."""
    rng = np.random.default_rng(seed)

    def pv(frac):
        p = rng.uniform(size=P); k = int(frac * P); p[:k] = 10 ** -rng.uniform(2, 9, size=k); rng.shuffle(p); return p.astype(np.float32)

    def save(rel, **a):
        f = os.path.join(root, rel); os.makedirs(os.path.dirname(f), exist_ok=True); np.savez_compressed(f, **a)

    for t in ("B1", "B2", "B3"):
        for k in range(3):
            for D in ("all", "200"):
                save(f"onek/out/ON_{t}_k{k}_D{D}.npz", PB=pv(.3), LCP=pv(.1), LCPXd=pv(.001), XDSM=pv(0), r_w=rng.normal(size=P).astype(np.float32))
    for s in ("P1", "P2", "P3"):
        for k in range(3):
            save(f"dc/out/DC_{s}_k{k}.npz", PB=pv(.2), LCP=pv(.01), XDS=pv(0)); save(f"rev/out/DCg_{s}_k{k}.npz", XDSg=pv(0))
        save(f"dc/out/DC_{s}_real.npz", PB=pv(.2), LCP=pv(.1), XDS=pv(.05), n1=np.array([380]), n0=np.array([470])); save(f"rev/out/DCg_{s}_real.npz", XDSg=pv(.01))
    for st in ("S1", "S2"):
        for k in range(3): save(f"rev/out/seanull_{st}_k{k}.npz", p_XDSM=pv(0), p_KX=pv(0), r_w=rng.normal(size=P).astype(np.float32))


def run(root, outdir, *args):
    return subprocess.run([sys.executable, os.path.join(root, "paper2", "make_nums.py"), "--outdir", outdir, *args], capture_output=True, text=True, timeout=300)


def test_clean_clone_regenerates_committed_files(tmp_path):
    root = copy_repo(str(tmp_path / "repo")); out = str(tmp_path / "out")
    r = run(root, out)
    assert r.returncode == 0, r.stderr[-2000:]
    for f in GENERATED:
        assert filecmp.cmp(os.path.join(out, f), os.path.join(REPO, "paper2", f), shallow=False), f"{f} differs from the committed file"
    # blocks whose raw .npz inputs are not in the repository must announce themselves and must not touch their outputs
    assert "WARNING [DCg]" in r.stderr and "WARNING [table_ablate]" in r.stderr
    assert not os.path.exists(os.path.join(out, "table_ablate.tex")) and not os.path.exists(os.path.join(out, "dcg_summary.json"))
    assert "NOTE [OneK1K]" in r.stderr and "NOTE [DC]" in r.stderr           # served from the committed compact summaries
    assert r.stdout.count("STALE") >= 3 and "\\DCgReal" in r.stdout
    assert run(root, out, "--strict").returncode == 2                         # --strict turns the skips into a failure


def test_raw_arrays_then_compact_summaries_give_identical_output(tmp_path):
    a = copy_repo(str(tmp_path / "a")); add_synthetic_raw(a)
    out_a = str(tmp_path / "out_a"); ra = run(a, out_a)
    assert ra.returncode == 0, ra.stderr[-2000:]
    assert "WARNING" not in ra.stderr
    for f in ("dc_summary.json", "dcg_summary.json", "onek_summary.json", "seanull_summary.json", "r1_summary.json"):
        assert os.path.exists(os.path.join(out_a, f)), f"{f} not written although the raw inputs were present"
    b = copy_repo(str(tmp_path / "b"))             # a clean clone in which the authors committed the summaries and no .npz
    for f in os.listdir(out_a):
        if f.endswith("_summary.json"): shutil.copy(os.path.join(out_a, f), os.path.join(b, "paper2", f))
    out_b = str(tmp_path / "out_b"); rb = run(b, out_b, "--strict")
    assert rb.returncode == 0, rb.stderr[-2000:]
    for f in ("nums.tex", "table_sim.tex", "table_null.tex", "table_ablate.tex"):
        assert filecmp.cmp(os.path.join(out_a, f), os.path.join(out_b, f), shallow=False), f"{f}: summary fallback differs from raw computation"


def test_missing_required_input_stops_with_both_paths(tmp_path):
    root = copy_repo(str(tmp_path / "repo")); os.remove(os.path.join(root, "results", "pop", "pop_results.json"))
    r = run(root, str(tmp_path / "out"))
    assert r.returncode != 0 and "pop/pop_results.json" in r.stderr and "results" in r.stderr


def test_absent_summary_keeps_committed_macro_and_warns(tmp_path):
    root = copy_repo(str(tmp_path / "repo")); os.remove(os.path.join(root, "paper2", "dc_summary.json")); out = str(tmp_path / "out")
    r = run(root, out)
    assert r.returncode == 0 and "WARNING [DC]" in r.stderr
    assert not os.path.exists(os.path.join(out, "dc_summary.json"))           # nothing is written in place of the missing summary
    kept = [l for l in open(os.path.join(out, "nums.tex")) if l.startswith("\\newcommand{\\DCtext}")]
    assert kept == [l for l in open(os.path.join(REPO, "paper2", "nums.tex")) if l.startswith("\\newcommand{\\DCtext}")]
