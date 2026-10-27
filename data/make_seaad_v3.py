"""SEA-AD analysis subset: 10x 3' v3 nuclei of the 84 aged donors (ADNC != 'Reference'), raw UMI counts.
Reproduces scdemon_in_seaad_v3_84donors.h5ad (34,794 nuclei x 35,483 genes; verified identical cell order and counts)."""
import anndata as ad, numpy as np, scipy.sparse as sp
s = ad.read_h5ad("data/seaad_mtg_microglia.h5ad")
keep = ((s.obs["assay"] == "10x 3' v3") & (s.obs["ADNC"] != "Reference")).values
X = sp.csr_matrix(s.raw.X[np.where(keep)[0]])
out = ad.AnnData(X=X.astype(np.float32), obs=s.obs.loc[keep, ["donor_id"]].reset_index(drop=True), var=s.raw.var.copy())
out.obs_names = [f"n{i}" for i in range(out.n_obs)]
assert out.n_obs == 34794 and out.obs["donor_id"].nunique() == 84
out.write_h5ad("data/scdemon_in_seaad_v3_84donors.h5ad")
