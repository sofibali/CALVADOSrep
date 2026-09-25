#!/usr/bin/env python
"""
Put an experimental reference structure into the same representation as a
CALVADOS model, so RMSD and Rg comparisons are like-for-like.

WHY THIS IS NEEDED
------------------
CALVADOS runs here use `use_com: true` (its shipped default), so
`build.geometry_from_pdb` places each bead at its **residue centre of mass**,
not at the CA. Comparing such a model against CA coordinates from a crystal
structure therefore carries a floor that has nothing to do with simulation
quality: for PARP14 the whole-chain Kabsch RMSD between the CA and COM
representations of the *same* structure is **0.135 nm (~1.3 A)**, and per-domain
0.126-0.147 nm. Domain Rg differs by 0.4-2.3% the same way.

See `docs/GEOMETRY_AUDIT.md` for how that was established.

Converting the *reference* to residue COM removes the floor without touching
the simulations, which is the point: the trajectories stay exactly as they were
run and only the yardstick changes.

MATCHING CALVADOS EXACTLY
-------------------------
Two things matter:

* **The same COM definition.** CALVADOS uses MDAnalysis
  `residue.atoms.center_of_mass()`, mass-weighted over whatever atoms the file
  contains. This module calls the same thing.

* **The same atom set.** The AF2/AF3 structures these simulations were built
  from contain **no hydrogens**, so their residue COM is a heavy-atom COM.
  Most crystal files match that, but NMR entries do not -- `1x4r.pdb` is 50%
  hydrogen. Including those H would shift its COM relative to every other
  reference and quietly re-introduce an offset for the WWE comparison alone.
  So `heavy_only=True` is the default here.

USAGE
-----
    from ref_geometry import residue_com

    coords, resids = residue_com('input/xtal_refs/3q6z.pdb', (792, 978))

    python ref_geometry.py --report      # show the floor this removes
"""
import argparse
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')


def _pick_chain(u, chain):
    """
    Residues of ONE chain.

    Several of these references are multi-copy: 3GOY has four chains A-D, every
    one numbered 1532-1720. Filtering on residue number alone silently returns
    all four copies stacked together -- 441 "residues" for a 118-residue domain,
    and an Rg of 3.2 nm instead of 1.5. Always pin the chain.
    """
    segs = [s for s in u.segments if len(s.atoms.select_atoms('name CA'))]
    if chain is not None:
        for s in segs:
            if (s.segid or '').strip() == chain:
                return s.residues
        raise ValueError(f'chain {chain!r} not found; have '
                         f'{[s.segid for s in segs]}')
    return segs[0].residues if segs else u.residues


def residue_com(pdb_path, resid_range=None, heavy_only=True, chain=None):
    """
    Per-residue centre of mass, in nm, matching CALVADOS's `use_com: true`.

    Returns (coords[n,3] in nm, resids[n]). `resid_range` is an inclusive
    (first, last) filter on author residue numbers. `chain` picks one copy in a
    multi-chain file; the first chain with CA atoms is used by default.
    """
    import MDAnalysis as mda

    u = mda.Universe(pdb_path)
    sel = u.select_atoms('not (type H or name H*)') if heavy_only else u.atoms
    coords, resids = [], []
    for res in _pick_chain(u, chain):
        atoms = res.atoms.intersection(sel)
        if len(atoms) == 0:
            continue
        if resid_range and not (resid_range[0] <= res.resid <= resid_range[1]):
            continue
        coords.append(atoms.center_of_mass())
        resids.append(int(res.resid))
    if not coords:
        return None, []
    return np.array(coords) / 10.0, resids


def residue_ca(pdb_path, resid_range=None, chain=None):
    """Per-residue CA, in nm -- the representation this replaces."""
    import MDAnalysis as mda

    u = mda.Universe(pdb_path)
    coords, resids = [], []
    for res in _pick_chain(u, chain):
        ca = res.atoms.select_atoms('name CA')
        if len(ca) == 0:
            continue
        if resid_range and not (resid_range[0] <= res.resid <= resid_range[1]):
            continue
        coords.append(ca.positions[0])
        resids.append(int(res.resid))
    if not coords:
        return None, []
    return np.array(coords) / 10.0, resids


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition, nm."""
    P = P - P.mean(0)
    Q = Q - Q.mean(0)
    V, _, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1, 1, d]) @ W
    return float(np.sqrt(((P @ R - Q) ** 2).sum(axis=1).mean()))


def rg(X):
    c = X.mean(axis=0)
    return float(np.sqrt(((X - c) ** 2).sum(axis=1).mean()))


def _report():
    """Show, per crystal reference, the floor that switching to COM removes."""
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    xtal = os.path.join(root, 'input', 'xtal_refs')
    refs = {
        'MD1 (3Q6Z)':  ('3q6z.pdb', (792, 978)),
        'MD2 (3VFQ)':  ('3vfq.pdb', (1005, 1191)),
        'WWE (3GOY)':  ('3goy.pdb', (1534, 1602)),
        'ART (3GOY)':  ('3goy.pdb', (1603, 1720)),
        'WWE (1X4R)':  ('1x4r.pdb', None),
    }
    print('Floor removed by comparing against residue COM instead of CA')
    print('(same structure, two representations -- nothing to do with the simulation)\n')
    print(f"{'reference':<14}{'res':>5}{'CA-vs-COM RMSD':>17}{'Rg CA':>9}{'Rg COM':>9}{'dRg':>8}")
    print('-' * 64)
    for name, (f, rng) in refs.items():
        p = os.path.join(xtal, f)
        if not os.path.isfile(p):
            print(f'{name:<14}  missing {f}')
            continue
        com, r1 = residue_com(p, rng)
        ca, r2 = residue_ca(p, rng)
        if com is None or ca is None:
            print(f'{name:<14}  no residues in range')
            continue
        keep = [i for i, r in enumerate(r2) if r in set(r1)]
        idx = {r: i for i, r in enumerate(r1)}
        ca = ca[keep]
        com = com[[idx[r2[i]] for i in keep]]
        print(f'{name:<14}{len(ca):>5}{kabsch_rmsd(ca, com):>16.3f}n'
              f'{rg(ca):>9.3f}{rg(com):>9.3f}{100 * (rg(com) / rg(ca) - 1):>+7.1f}%')
    print('\nNote: 1x4r is an NMR entry and is ~50% hydrogen. heavy_only=True keeps its')
    print('COM on the same footing as the H-free AF2/AF3 structures the sims were built')
    print('from -- otherwise the WWE comparison alone would carry a different offset.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--report', action='store_true')
    a = ap.parse_args()
    if a.report:
        _report()
    else:
        ap.print_help()
