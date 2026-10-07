"""Add hydrogens to the all-atom PDB before at2hyres.

HyRes keeps the backbone amide H explicitly (it drives the hydrogen-bond term,
eps_hb), so the HyRes PSF expects one per non-proline residue. at2hyres only
MAPS atoms that exist -- it does not build them -- so an input without
hydrogens yields a PDB 1,515 particles short of its own PSF.

This is the step CHARMM-GUI performs in the documented workflow. Doing it with
OpenMM's Modeller against CHARMM36 gives the same result without the web form.
"""
import sys
from openmm.app import PDBFile, Modeller, ForceField
inp, out = sys.argv[1], sys.argv[2]
pdb = PDBFile(inp)
m = Modeller(pdb.topology, pdb.positions)
m.addHydrogens(ForceField('charmm36.xml'), pH=7.0)
with open(out, 'w') as g:
    PDBFile.writeFile(m.topology, m.positions, g, keepIds=True)
n = sum(1 for a in m.topology.atoms())
nh = sum(1 for a in m.topology.atoms() if a.element is not None and a.element.symbol == 'H')
print(f'{out}: {n} atoms ({nh} H) from {pdb.topology.getNumAtoms()}')
