import numpy as np
import pandas as pd
import MDAnalysis as mda
from calvados.analysis import calc_rg, calc_ete, calc_cmap, fit_scaling_exp

# Load trajectory
u = mda.Universe('top.pdb', 'proteins_only.dcd')
ag = u.select_atoms('all')

# Load residue parameters
residues = pd.read_csv('../input/residues_CALVADOS3.csv').set_index('three')

# Radius of gyration
rgs = calc_rg(u, ag, ag.resnames.tolist(), residues, start=100)
print(f"Rg = {np.mean(rgs):.2f} ± {np.std(rgs):.2f} nm")

# End-to-end distance
rees, ree_mean, ree_sem = calc_ete(u, ag, start=100)
print(f"Ree = {ree_mean:.2f} ± {ree_sem:.2f} nm")

# Scaling exponent (for IDRs)
ij, dij, r0, nu, nu_err = fit_scaling_exp(u, ag, start=100)
print(f"ν = {nu:.3f} ± {nu_err:.3f}")

# Contact map
cmap = calc_cmap(ag, ag, cutoff=1.0)
np.save('contact_map.npy', cmap)
