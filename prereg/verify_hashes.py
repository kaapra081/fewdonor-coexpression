"""Check every frozen-plan hash in the repository against the file it names.   python prereg/verify_hashes.py [--strict]

Reads prereg/*.sha256 and rev3/*.sha256 (lines 'sha256  filename'; the first line is the plan, further lines in DEV/HELDOUT files
are code files frozen with it) and reports MATCH / DIFFERS / MISSING for each. Known, documented differences are labelled and do
not fail --strict: DEV_PREREG.sha256 and XS_PREREG_original.sha256 hash plan versions before an amendment recorded in the later
hash files (DEV_PREREG_change1.sha256, XS_PREREG.sha256), and the code files hashed with DEV/HELDOUT plans were edited after the
plan (the amendments are listed in prereg/GO_NOGO.md and in the plans themselves).
The hashes are self-generated text files: they show that a plan was not edited after hashing, not when it was written. The only
external time stamp is the git commit date, and plans and results often share a commit (see README.md, section Frozen plans).
"""
import os, sys, glob, hashlib

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
KNOWN = {"DEV_PREREG.sha256": "plan amended afterwards; current file is hashed in DEV_PREREG_change1.sha256",
         "XS_PREREG_original.sha256": "plan amended afterwards; current file is hashed in XS_PREREG.sha256"}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def resolve(name, sha_dir):
    for c in (os.path.join(ROOT, name), os.path.join(sha_dir, os.path.basename(name)), os.path.join(ROOT, "prereg", os.path.basename(name)), os.path.join(ROOT, "dev", os.path.basename(name))):
        if os.path.exists(c): return c


def main():
    bad = 0
    for sf in sorted(glob.glob(os.path.join(ROOT, "prereg", "*.sha256")) + glob.glob(os.path.join(ROOT, "rev3", "*.sha256"))):
        lines = [l.split() for l in open(sf) if len(l.split()) == 2]
        for i, (h, name) in enumerate(lines):
            f = resolve(name, os.path.dirname(sf)); kind = "plan" if i == 0 else "code"
            status = "MISSING" if f is None else ("MATCH" if sha(f) == h else "DIFFERS")
            why = KNOWN.get(os.path.basename(sf), "") if (status == "DIFFERS" and i == 0) else ("edited after the plan (documented amendment)" if status == "DIFFERS" else "")
            if status != "MATCH" and not why: bad += 1
            print(f"{os.path.relpath(sf, ROOT):38s} {kind:4s} {name:22s} {h[:12]}  {status:8s} {why}")
    print(f"\n{bad} unexplained mismatch(es)")
    return 1 if (bad and "--strict" in sys.argv) else 0


if __name__ == "__main__":
    sys.exit(main())
