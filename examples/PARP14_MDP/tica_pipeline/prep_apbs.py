#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare a real Poisson-Boltzmann electrostatic potential map for one
all-atom state PDB, via pdb2pqr (charge/radius assignment) + APBS
(SBgrid installs of both), instead of PyMOL's built-in
util.protein_vacuum_esp -- which has a reproducible internal bug that
makes it unusable in a scripted batch context (see make_surface_gallery.py
history / tica_pipeline/README notes).

Steps:
  1. Recenter the PDB near the origin. CALVADOS sim boxes place coordinates
     at large absolute values (often >1000 A); pdb2pqr's fixed-width PQR
     columns overflow at that range and APBS's PQR parser then reads
     adjacent coordinates as one unparseable run-together token. A plain
     translation fixes this with no effect on shape/electrostatics.
  2. pdb2pqr --ff=amber --with-ph=7.0  (AMBER charges/radii, protonate at pH 7)
  3. Auto-size an APBS mg-auto grid from the recentered structure's actual
     extent (pdb2pqr's own --apbs-input auto-sizer crashes on this
     structure with a TypeError; sized manually instead).
  4. Run apbs -> writes a .dx potential map (kT/e units, single-level
     mg-auto, no ionic strength focusing refinement -- adequate for
     qualitative whole-molecule/pocket surface coloring, not a rigorous
     binding free-energy calculation).

Requires SBgrid (`source /programs/sbgrid.shrc`) for pdb2pqr + apbs.

Usage:
    python prep_apbs.py --pdb state_1_allatom.pdb --out-dir OUTDIR
Produces:
    OUTDIR/centered.pdb   (recentered input, matches the .dx map's frame)
    OUTDIR/esp.dx         (APBS potential map, kT/e)
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

SBGRID_SOURCE = 'source /programs/sbgrid.shrc >/dev/null 2>&1'
ALLOWED_DIME = [65, 97, 129, 161, 193, 225]
MAX_GRID_SPACING = 1.4  # Angstrom; coarser -> faster, still fine at this scale
PADDING = 15.0  # Angstrom padding around actual extent for the fine grid


def recenter_pdb(in_pdb, out_pdb):
    """Translate all ATOM/HETATM coords so their centroid is at the origin.
    Returns the (dx, dy, dz) extent of the recentered structure."""
    lines = []
    coords = []
    with open(in_pdb) as f:
        for line in f:
            lines.append(line)
            if line.startswith(('ATOM', 'HETATM')):
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    coords = np.array(coords)
    center = coords.mean(axis=0)

    with open(out_pdb, 'w') as f:
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                x = float(line[30:38]) - center[0]
                y = float(line[38:46]) - center[1]
                z = float(line[46:54]) - center[2]
                f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
            else:
                f.write(line)

    recentered = coords - center
    extent = recentered.max(axis=0) - recentered.min(axis=0)
    return extent


def pick_dime(length):
    """Smallest allowed APBS dime giving grid spacing <= MAX_GRID_SPACING
    over `length`, capped at the largest allowed value."""
    for d in ALLOWED_DIME:
        if length / (d - 1) <= MAX_GRID_SPACING:
            return d
    return ALLOWED_DIME[-1]


def write_apbs_input(pqr_path, dx_stem, extent, apbs_in_path):
    fglen = extent + 2 * PADDING
    cglen = fglen * 1.3
    dime = [pick_dime(fl) for fl in fglen]
    apbs_in_path.write_text(f"""read
    mol pqr {pqr_path.name}
end
elec
    mg-auto
    dime {dime[0]} {dime[1]} {dime[2]}
    cglen {cglen[0]:.1f} {cglen[1]:.1f} {cglen[2]:.1f}
    fglen {fglen[0]:.1f} {fglen[1]:.1f} {fglen[2]:.1f}
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl mdh
    pdie 2.0
    sdie 78.54
    chgm spl2
    srfm smol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy no
    calcforce no
    write pot dx {dx_stem}
end
quit
""")


def run(cmd, cwd, log_path):
    with open(log_path, 'w') as log:
        result = subprocess.run(
            f"{SBGRID_SOURCE} && {cmd}", shell=True, cwd=cwd,
            stdout=log, stderr=subprocess.STDOUT, executable='/bin/bash')
    return result.returncode


def prepare(pdb_path, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    centered_pdb = out_dir / 'centered.pdb'
    pqr_path = out_dir / 'state.pqr'
    apbs_in = out_dir / 'apbs.in'
    dx_path = out_dir / 'esp.dx'

    extent = recenter_pdb(pdb_path, centered_pdb)

    rc = run(f"pdb2pqr --ff=amber --with-ph=7.0 {centered_pdb.name} {pqr_path.name}",
             out_dir, out_dir / 'pdb2pqr.log')
    if rc != 0 or not pqr_path.exists():
        return None, (out_dir / 'pdb2pqr.log').read_text()[-2000:]

    write_apbs_input(pqr_path, 'esp', extent, apbs_in)
    rc = run(f"apbs {apbs_in.name}", out_dir, out_dir / 'apbs.log')
    if rc != 0 or not dx_path.exists():
        return None, (out_dir / 'apbs.log').read_text()[-2000:]

    return {'centered_pdb': str(centered_pdb), 'dx': str(dx_path)}, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdb', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    result, err = prepare(args.pdb, args.out_dir)
    if result is None:
        print(f"APBS prep FAILED:\n{err}", file=sys.stderr)
        return 1
    print(f"centered_pdb={result['centered_pdb']}")
    print(f"dx={result['dx']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
