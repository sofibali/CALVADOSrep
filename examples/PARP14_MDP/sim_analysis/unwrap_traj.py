#!/usr/bin/env python
"""
Remove periodic-boundary jumps from a monomer trajectory, for visualisation.

WHAT IS AND IS NOT A PROBLEM
----------------------------
CALVADOS writes the trajectory through OpenMM's DCDReporter with
`enforcePeriodicBox` left at its default, which for a periodic system wraps
whole MOLECULES -- it translates each bonded molecule as a unit rather than
splitting it. Verified on these trajectories: the largest consecutive CA-CA
distance is 0.9 nm against boxes of 70-300 nm, i.e. no chain is ever broken.

So **the analyses are fine**. Rg, inter-domain distances, contact maps, RMSF
and the active-site metrics are all computed on intact chains and need no
unwrapping.

What does break is **watching** the trajectory. When the molecule's centre of
mass crosses a boundary the whole chain is translated by one box length, so a
viewer sees it teleport. In the 1 us extension runs this happens often enough
to be distracting:

    md_full       80 nm box, 100k frames, 191 wraps
    md1_md2       70 nm                   206
    md2_md3       70 nm                   187
    mka_full     100 nm                    45
    core_full_go 120 nm                    29

The original 3-10k-frame runs in 120-300 nm boxes never wrap.

Anything that uses ABSOLUTE position or centre-of-mass displacement -- a mean
squared displacement or a diffusion coefficient -- would also be wrong on the
raw trajectory. Nothing in this project currently computes those, but if you
add one, use `--mode unwrap` output rather than the raw DCD.

MODES
-----
  center   (default) subtract the COM every frame. The molecule sits still and
           you see only conformational change. Best for visualisation.
  unwrap   follow the COM continuously across boundaries instead of removing
           it. Keeps real diffusion, just without the jumps. Use this if you
           want MSD/diffusion.

USAGE
-----
    python unwrap_traj.py <run_dir>                 # writes {sysname}_viz.dcd
    python unwrap_traj.py <run_dir> --mode unwrap
    python unwrap_traj.py --all                     # every run that wraps
    python unwrap_traj.py <run_dir> --check         # report only, write nothing
"""
import argparse
import glob
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)


def wrap_events(xyz, box):
    """(n_wraps, max_jump_nm) for the molecule's centre of mass."""
    com = xyz.mean(axis=1)
    jump = np.linalg.norm(np.diff(com, axis=0), axis=1)
    return int((jump > box / 4).sum()), float(jump.max() if len(jump) else 0.0)


def fix(xyz, lengths, mode):
    """Return coordinates with PBC jumps removed."""
    com = xyz.mean(axis=1)
    if mode == 'center':
        return xyz - com[:, None, :]
    # unwrap: accumulate the box-length steps the COM took, and undo them
    shift = np.zeros_like(com)
    acc = np.zeros(3)
    for i in range(1, len(com)):
        d = com[i] - com[i - 1]
        acc -= np.round(d / lengths[i]) * lengths[i]
        shift[i] = acc
    return xyz + shift[:, None, :]


def run_one(d, mode, check, stride):
    import mdtraj as md
    import yaml

    top = os.path.join(d, 'top.pdb')
    dcds = [f for f in os.listdir(d)
            if f.endswith('.dcd') and not f.startswith(('backup', 'equilibration'))
            and not f.endswith(('_viz.dcd', '_centered.dcd'))]
    if not dcds or not os.path.isfile(top):
        return None
    name = os.path.basename(os.path.normpath(d))
    t = md.load(os.path.join(d, dcds[0]), top=top, stride=stride)
    box = float(t.unitcell_lengths[0, 0])
    n, mx = wrap_events(t.xyz, box)
    print(f"  {name:<28} box {box:>5.0f} nm  {t.n_frames:>6} frames  "
          f"{n:>4} wraps (max jump {mx:.1f} nm)")
    if check:
        return n
    if n == 0:
        print(f"      no wraps -- nothing to fix")
        return 0
    t.xyz = fix(t.xyz, t.unitcell_lengths, mode)
    out = os.path.join(d, dcds[0].replace('.dcd', '_viz.dcd'))
    t.save_dcd(out)
    n2, _ = wrap_events(t.xyz, box)
    print(f"      wrote {os.path.basename(out)} ({mode}); wraps after: {n2}")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', nargs='?')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--mode', choices=('center', 'unwrap'), default='center')
    ap.add_argument('--check', action='store_true', help='report only')
    ap.add_argument('--stride', type=int, default=1)
    args = ap.parse_args()

    if args.all:
        import sim_registry as reg
        targets = []
        for k in reg.SETS:
            try:
                d = reg.get_sim_dir(k, 1, 0)
            except Exception:
                continue
            if os.path.isdir(d):
                targets.append(d)
    else:
        if not args.run_dir:
            sys.exit('give a run directory or --all')
        targets = [os.path.abspath(args.run_dir)]

    print(f"{'checking' if args.check else args.mode} {len(targets)} run(s)\n")
    total = 0
    for d in targets:
        r = run_one(d, args.mode, args.check, args.stride)
        if r:
            total += r
    print(f"\n{total} wrap event(s) across the runs examined")
    if args.check and total:
        print("re-run without --check to write *_viz.dcd")


if __name__ == '__main__':
    main()
