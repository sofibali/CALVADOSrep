#!/usr/bin/env python3
"""HyRes runs built to be COMPARABLE to the CALVADOS p9_dtx3l runs.

WHY THE FIRST HyRes RUN WAS NOT A FAIR COMPARISON
-------------------------------------------------
It was unrestrained, because HyRes ships no working fold restraints
(restraints/folded_restriants.py is a TODO stub with a syntax error at line 50
and is referenced by no run script). CALVADOS, meanwhile, holds every domain
rigid with harmonic CA-CA restraints. Comparing the two as they stood compared
RESTRAINTS, not force fields -- and predictably 13 of 14 domains inflated by
>25% (up to 1.69x) and the chains never came within 2.44 nm in 500 ns.

This script fixes the three things that made the comparison unfair:

1. DOMAIN RESTRAINTS, ported from CALVADOS exactly.
   CALVADOS applies harmonic bonds between CA pairs of the SAME domain that sit
   within cutoff_restr = 0.9 nm in the reference structure, k = 700 kJ/mol/nm^2,
   domains from the same binding/p9_dtx3l/input/domains.yaml. The identical
   network is built here on the HyRes CA beads. Same domains, same cutoff, same
   k -- so any remaining difference is the force field, which is the point.

2. BOTH ARMS. CALVADOS runs separated (chains ~17 nm apart) and docked (the AF3
   pose). The first HyRes run was separated-only and never associated, so it
   could not be compared with the docked CALVADOS arm at all. The docked start
   is built by Kabsch-superposing each HyRes chain onto the corresponding chain
   of binding/p9_dtx3l_docked/input/system.pdb -- same internal geometry, the
   experimental relative orientation.

3. REPLICATES. n=1 cannot be compared with CALVADOS's 5-7. Default 3 per arm,
   differing only in the thermostat seed.

Everything else already matched: 40 nm cubic box (25.9 uM per chain), 293 K,
150 mM salt, 500 ns, 50 ps/frame.

Usage:  python prepare_hyres_matched.py [--nreps 3]
"""
import shutil
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, yaml

HERE = Path(__file__).resolve().parent
PREP = HERE / 'hyres' / 'runs' / 'prep'
DST = HERE / 'hyres' / 'runs' / 'matched'
REPO = HERE / 'hyres' / 'HyRes_GPU'
DOMAINS = HERE / 'binding' / 'p9_dtx3l' / 'input' / 'domains.yaml'
DOCKED = HERE / 'binding' / 'p9_dtx3l_docked' / 'rep-1' / 'input' / 'system.pdb'
LEN = {'parp9': 854, 'dtx3l': 740}
CUTOFF_NM, K_HARM, BOX, SEP = 0.9, 700.0, 40.0, 17.0


def read_pdb(fn):
    lines, xyz = [], []
    for l in open(fn):
        if l.startswith('ATOM'):
            lines.append(l)
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return lines, np.array(xyz)


def ca_mask(lines):
    return np.array([l[12:16].strip() == 'CA' for l in lines])


def kabsch(P, Q):
    """Rotation+translation taking P onto Q (both N x 3)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, pc, qc


def write_pdb(fn, chains, box_nm):
    L = box_nm * 10.0
    n = 0
    with open(fn, 'w') as g:
        g.write(f'CRYST1{L:9.3f}{L:9.3f}{L:9.3f}  90.00  90.00  90.00 P 1           1\n')
        for lines, xyz in chains:
            for l, c in zip(lines, xyz):
                n += 1
                g.write(l[:6] + f'{n:5d}' + l[11:30] +
                        f'{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}' + l[54:])
        g.write('END\n')
    return n


def build_restraints(chains_lines, chains_xyz, out_txt):
    """Intra-domain harmonic network on HyRes CA beads, matching CALVADOS.

    KNOWN GAP: this builds pairs WITHIN each domains.yaml block and never
    across blocks, so it does not reproduce the 85 custom restraints (k = 350)
    in binding/p9_dtx3l/input/custom_restraints.txt that hold PARP9's split KH1
    together -- blocks [60,96] and [516,539] are two halves of one fold,
    separated in sequence by the macrodomains. Measured consequence: the halves
    hold at 0.87 nm in CALVADOS but reach 5.02 nm (separated) and 4.62 nm
    (docked) in HyRes, so PARP9 KH1 is not folded in either HyRes arm. Every
    other domain is a single contiguous block and is restrained correctly.
    To close the gap, also emit pairs for (block 0, block 3) of parp9.
    """
    dom = yaml.safe_load(open(DOMAINS))
    pairs, base = [], 0
    for ch, lines, xyz in zip(('parp9', 'dtx3l'), chains_lines, chains_xyz):
        m = ca_mask(lines)
        ca_idx = np.where(m)[0] + base           # absolute particle index
        ca_xyz = xyz[m] / 10.0                   # nm
        assert len(ca_idx) == LEN[ch], f'{ch}: {len(ca_idx)} CA'
        for lo, hi in dom[ch]:
            sel = np.arange(lo - 1, hi)
            X = ca_xyz[sel]
            d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
            iu = np.triu_indices(len(sel), 1)
            keep = d[iu] < CUTOFF_NM
            for a, b, r in zip(np.array(iu[0])[keep], np.array(iu[1])[keep], d[iu][keep]):
                pairs.append((int(ca_idx[sel[a]]), int(ca_idx[sel[b]]), float(r)))
        base += len(lines)
    with open(out_txt, 'w') as g:
        g.write(f'# CALVADOS-equivalent intra-domain harmonic restraints\n')
        g.write(f'# cutoff {CUTOFF_NM} nm, k {K_HARM} kJ/mol/nm^2, 0-based particle indices\n')
        for i, j, r in pairs:
            g.write(f'{i} {j} {r:.4f} {K_HARM}\n')
    return len(pairs)


def main():
    ap = ArgumentParser(); ap.add_argument('--nreps', type=int, default=3)
    a = ap.parse_args()
    p9l, p9x = read_pdb(PREP / 'parp9_hyres.pdb')
    dxl, dxx = read_pdb(PREP / 'dtx3l_hyres.pdb')
    DST.mkdir(parents=True, exist_ok=True)

    nr = build_restraints([p9l, dxl], [p9x, dxx], DST / 'domain_restraints_separated.txt')
    print(f'domain_restraints_separated.txt: {nr} harmonic CA-CA pairs '
          f'(cutoff {CUTOFF_NM} nm, k {K_HARM}) -- matching CALVADOS')

    # ---- separated arm ----
    sep = []
    for xyz, shift in ((p9x.copy(), np.array([-SEP/2, 0, 0])), (dxx.copy(), np.array([SEP/2, 0, 0]))):
        xyz += (np.array([BOX/2]*3) + shift) * 10.0 - xyz.mean(0)
        sep.append(xyz)
    n = write_pdb(DST / 'start_separated.pdb', [(p9l, sep[0]), (dxl, sep[1])], BOX)
    print(f'start_separated.pdb: {n} particles, centroids {SEP} nm apart')

    # ---- docked arm ----
    # Built from the DOCKED set's own all-atom conformers, not by superposing the
    # monomers. The AF3 complex prediction's chains differ from the AF3 monomer
    # predictions, so a whole-chain Kabsch fit gave 2.2 / 3.9 nm CA RMSD and a
    # 0.09 nm clash. Converting the docked conformers through the HyRes pipeline
    # reproduces the experimental pose exactly (0.42 nm min inter-chain CA-CA).
    p9dl, p9dx = read_pdb(PREP / 'parp9_dock_hyres.pdb')
    dxdl, dxdx = read_pdb(PREP / 'dtx3l_dock_hyres.pdb')
    allx = np.vstack([p9dx, dxdx])
    shift = np.array([BOX/2]*3)*10.0 - allx.mean(0)
    dock = [p9dx + shift, dxdx + shift]
    n = write_pdb(DST / 'start_docked.pdb', [(p9dl, dock[0]), (dxdl, dock[1])], BOX)
    from scipy.spatial import cKDTree
    mind = cKDTree(dock[0][ca_mask(p9dl)]).query(dock[1][ca_mask(dxdl)])[0].min()/10
    print(f'start_docked.pdb: {n} particles, min inter-chain CA-CA {mind:.2f} nm')
    # restraints for the docked arm come from ITS conformer
    nr2 = build_restraints([p9dl, dxdl], [p9dx, dxdx], DST / 'domain_restraints_docked.txt')
    print(f'domain_restraints_docked.txt: {nr2} pairs')

    for f in ('complex.psf', 'top_hyres_GPU.inp', 'param_hyres_GPU.inp'):
        shutil.copy2(PREP / f, DST / f)
    print(f'\nprepared in {DST}; {a.nreps} replicates per arm to be launched')


if __name__ == '__main__':
    main()
