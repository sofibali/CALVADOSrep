#!/usr/bin/env python
"""
Make representative CG states full-atom AND remove the clashes that back-mapping
introduces.

WHY CLASHES HAPPEN: back-mapping grafts an all-atom reference domain onto each
domain's CG (CA-bead) positions (per-domain Kabsch). CALVADOS one-bead-per-residue
contacts pack tighter than all-atom domains can, so when two domains touch in the CG
frame their reconstructed side chains/backbones overlap (steric clashes).

CLASH REMOVAL (recommended recipe, implemented here):
  Restrained energy minimization in OpenMM.
    - add hydrogens (Modeller) and build an amber14 system (vacuum, NoCutoff),
    - restrain CA atoms to their back-mapped positions (harmonic, k=1000 kJ/mol/nm^2)
      so the *domain arrangement is preserved* while side chains + backbone relax,
    - L-BFGS minimize (minimizeEnergy); the steep LJ repulsion pushes overlapping
      atoms apart, the CA restraints stop the fold/arrangement from drifting,
    - strip hydrogens, write the heavy-atom structure.
  This clears local/side-chain and mild backbone clashes (the common case).

  For DEEP inter-domain interpenetration (domains that overlap by whole residues),
  CA-restrained minimization can't separate interpenetrating backbones. Then:
    (a) rigid-body declash first (translate domains apart along their COM axis to
        van-der-Waals contact, preserving orientation) -- see bindcraft_md/
        prepare_targets.py:declash() -- then minimize; or
    (b) re-run with a weaker CA restraint (--k 200) and/or 2 passes (--passes 2),
        which the script does automatically if residual clashes remain.

Usage:
  python fullatom_minimize_states.py \
      --states-dir representative_frames/2026-06-25/md_full_pose_tica \
      --reference md_full/input/ref_allatom.pdb \
      --fdomains md_full/input/domains.yaml
  # FL: --reference input/parp14.pdb --fdomains input/domains.yaml
"""
import os, sys, glob, argparse, warnings
warnings.simplefilter('ignore')
import numpy as np
import yaml
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from backmap_states import backmap_state, parse_pdb_atoms

import openmm
from openmm import app, unit, CustomExternalForce, LangevinIntegrator, Platform
from pdbfixer import PDBFixer


def count_clashes(pdb, cut=2.0, seq_excl=1):
    """Heavy-atom pairs closer than cut (A), excluding same/near-sequential residues."""
    st = PDBParser(QUIET=True).get_structure('x', pdb)
    atoms = [a for a in st.get_atoms() if a.element != 'H']
    ns = NeighborSearch(atoms)
    n = 0
    for a, b in ns.search_all(cut, 'A'):
        ra, rb = a.get_parent(), b.get_parent()
        if ra.get_parent().id == rb.get_parent().id and abs(ra.id[1] - rb.id[1]) <= seq_excl:
            continue
        n += 1
    return n


def declash_allatom(pdb_path, domains, clash_cut=2.6, step=0.3, max_iter=300):
    """Rigid-body separate interpenetrating domains to vdW contact (preserves
    orientation). Partition residues into domain blocks at the midpoints between
    cores, fix the central block, translate the others outward along their
    COM->center axis until no heavy atom is within clash_cut of another block.
    Overwrites pdb_path. For the deep-interpenetration case minimization can't fix.
    """
    st = PDBParser(QUIET=True).get_structure('x', pdb_path)
    model = next(st.get_models())
    res = [r for ch in model for r in ch if r.id[0] == ' ']
    cores = [tuple(d) for d in domains]
    # block boundaries at midpoints between consecutive cores
    bnds = [cores[0][0]]
    for i in range(len(cores) - 1):
        bnds.append((cores[i][1] + cores[i + 1][0]) // 2)
    bnds.append(cores[-1][1] + 10**6)
    def block_of(rid):
        for i in range(len(cores)):
            if bnds[i] <= rid < bnds[i + 1]:
                return i
        return len(cores) - 1
    blocks = {i: [] for i in range(len(cores))}
    for r in res:
        blocks[block_of(r.id[1])].append(r)
    def atoms(i): return [a for r in blocks[i] for a in r if a.element != 'H']
    def com(i):
        x = np.array([a.get_coord() for a in atoms(i)]); return x.mean(0)
    fixed = len(cores) // 2
    fcom = com(fixed)
    for i in range(len(cores)):
        if i == fixed:
            continue
        axis = com(i) - fcom; n = np.linalg.norm(axis)
        axis = axis / n if n > 1e-6 else np.array([1.0, 0, 0])
        others = [j for j in range(len(cores)) if j != i]
        oa = [a for j in others for a in atoms(j)]
        for _ in range(max_iter):
            ns = NeighborSearch(atoms(i) + oa)
            clash = any(a.element != 'H' and b.element != 'H'
                        for a, b in ns.search_all(clash_cut, 'A'))
            if not clash:
                break
            for a in atoms(i):
                a.set_coord(a.get_coord() + axis * step)
    io = __import__('Bio.PDB', fromlist=['PDBIO']).PDBIO()
    io.set_structure(st); io.save(pdb_path)


def minimize(in_pdb, out_pdb, k=1000.0, max_iter=3000):
    """Restrained amber14 minimization; writes heavy-atom-only output."""
    # PDBFixer: add OXT / terminal atoms + any missing heavy atoms (NOT missing
    # residues -- we keep linkers as back-mapped), so amber templates match.
    fixer = PDBFixer(filename=in_pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    ff = app.ForceField('amber14-all.xml')
    mod = app.Modeller(fixer.topology, fixer.positions)
    mod.addHydrogens(ff)
    system = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    # CA positional restraints (preserve arrangement)
    restr = CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restr.addGlobalParameter('k', k * unit.kilojoule_per_mole / unit.nanometer**2)
    for p in ('x0', 'y0', 'z0'):
        restr.addPerParticleParameter(p)
    pos = mod.positions
    for atom in mod.topology.atoms():
        if atom.name == 'CA':
            xyz = pos[atom.index].value_in_unit(unit.nanometer)
            restr.addParticle(atom.index, xyz)
    system.addForce(restr)
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               0.001 * unit.picoseconds)
    sim = app.Simulation(mod.topology, system, integ,
                         Platform.getPlatformByName('CPU'))
    sim.context.setPositions(mod.positions)
    sim.minimizeEnergy(maxIterations=max_iter)
    state = sim.context.getState(getPositions=True)
    # write heavy atoms only
    with open(out_pdb, 'w') as fh:
        app.PDBFile.writeFile(mod.topology, state.getPositions(),
                              fh, keepIds=True)
    # strip H lines for a clean heavy-atom file
    lines = [l for l in open(out_pdb)
             if not (l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() == 'H')]
    open(out_pdb, 'w').writelines(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--states-dir', required=True, type=Path)
    ap.add_argument('--reference', '-r', required=True, type=Path)
    ap.add_argument('--fdomains', required=True, type=Path)
    ap.add_argument('--out-dir', type=Path, default=None)
    ap.add_argument('--k', type=float, default=1000.0, help='CA restraint k (kJ/mol/nm^2)')
    ap.add_argument('--passes', type=int, default=2, help='Re-minimize (weaker k) if clashes remain')
    args = ap.parse_args()

    out_dir = args.out_dir or (args.states_dir.parent / (args.states_dir.name + '_allatom'))
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_atoms = parse_pdb_atoms(args.reference)
    dom = yaml.safe_load(open(args.fdomains))
    domains = list(dom.values())[0] if isinstance(dom, dict) else dom  # [[s,e],...]

    states = sorted(glob.glob(str(args.states_dir / 'state_*.pdb')))
    print(f"{len(states)} states | ref {args.reference.name} | {len(domains)} domains\n")
    print(f"{'state':10} {'clash_backmap':>13} {'clash_min':>10} {'status':>8}")
    for sp in states:
        sp = Path(sp)
        bm = out_dir / (sp.stem + '_allatom.pdb')
        backmap_state(sp, ref_atoms, domains, bm, pulchra=False, verbose=False)
        c0 = count_clashes(bm)
        minimize(str(bm), str(bm), k=args.k)
        c1 = count_clashes(bm)
        kk = args.k
        p = 1
        while c1 > 0 and p < args.passes:        # weaker restraint if residual clashes
            kk *= 0.3
            minimize(str(bm), str(bm), k=kk)
            c1 = count_clashes(bm); p += 1
        declashed = False
        if c1 > 3:                               # deep interpenetration: declash then re-min
            backmap_state(sp, ref_atoms, domains, bm, pulchra=False, verbose=False)
            declash_allatom(bm, domains)
            minimize(str(bm), str(bm), k=args.k)
            c1 = count_clashes(bm); declashed = True
        status = ('clean' if c1 == 0 else ('ok' if c1 <= 3 else 'CHECK'))
        status += '+declash' if declashed else ''
        print(f"{sp.stem:10} {c0:>13} {c1:>10} {status:>10}")
    print(f"\nWrote full-atom, minimized states -> {out_dir}")
    print("If any row says CHECK (deep interpenetration), rigid-body declash first "
          "(bindcraft_md/prepare_targets.py:declash) then re-run.")


if __name__ == '__main__':
    main()
