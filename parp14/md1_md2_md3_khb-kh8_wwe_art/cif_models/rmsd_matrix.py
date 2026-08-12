"""Pairwise CA RMSD (after optimal Kabsch superposition) across all model.cif in this folder."""
import numpy as np
from Bio.PDB import MMCIFParser

N = 25
parser = MMCIFParser(QUIET=True)


def ca_coords(path):
    structure = parser.get_structure("m", path)
    return np.array(
        [res["CA"].coord for res in structure.get_residues() if "CA" in res],
        dtype=float,
    )


coords = [ca_coords(f"{i}.cif") for i in range(1, N + 1)]
n_ca = [len(c) for c in coords]
assert len(set(n_ca)) == 1, f"CA count mismatch across models: {n_ca}"
print(f"{N} models, {n_ca[0]} CA atoms each")


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    P_rot = (R @ Pc.T).T
    diff = P_rot - Qc
    return np.sqrt((diff ** 2).sum(axis=1).mean())


mat = np.zeros((N, N))
for i in range(N):
    for j in range(i + 1, N):
        r = kabsch_rmsd(coords[i], coords[j])
        mat[i, j] = mat[j, i] = r

np.savetxt("rmsd_matrix.csv", mat, delimiter=",", fmt="%.3f")

iu = np.triu_indices(N, k=1)
vals = mat[iu]
print(f"pairwise CA RMSD (A): mean={vals.mean():.2f}  std={vals.std():.2f}  min={vals.min():.2f}  max={vals.max():.2f}")

i_max, j_max = np.unravel_index(np.argmax(mat), mat.shape)
print(f"most divergent pair: {i_max+1}.cif vs {j_max+1}.cif = {mat[i_max, j_max]:.2f} A")

mean_per_model = mat.sum(axis=0) / (N - 1)
best = np.argmin(mean_per_model) + 1
print(f"most 'central' model (lowest avg RMSD to all others): {best}.cif ({mean_per_model[best-1]:.2f} A avg)")
