"""Extract compact call vectors for Amendments 4-5 from rev4_R6.py output: python rev4_A4extract.py P1 P2 P3 -> out/A4/A4_<s>.npz
(the exact command used on the compute host on 2 Nov 2026, written out as a script)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); os.makedirs(os.path.join(HERE, "out", "A4"), exist_ok=True)
for s in sys.argv[1:]:
    Z = np.load(os.path.join(HERE, "out", f"R6_{s}.npz"), allow_pickle=True)
    np.savez_compressed(os.path.join(HERE, "out", "A4", f"A4_{s}.npz"), null_cPB=np.packbits(Z["null_cPB"].astype(bool)), real_cPB=np.packbits(Z["real_cPB"].astype(bool)),
                        real_cXDSM=np.packbits(Z["real_cXDSM"].astype(bool)), n=len(Z["null_cPB"]), sym=Z["sym"].astype(str), mean_umi=Z["mean_umi"], detect=Z["detect"])
    print(s, int(Z["null_cPB"].sum()), int(Z["real_cPB"].sum()), int(Z["real_cXDSM"].sum()))
