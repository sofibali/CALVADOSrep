#!/usr/bin/env python
"""
Extend the full-length PARP14 runs (fl_optimized/) toward a TARGET total sampling,
keeping the domain restraints AND the KH7a-KHb custom restraints already configured.

CALVADOS checkpoint restart runs `steps` ADDITIONAL steps and APPENDS to the DCD.
This script is TARGET-DRIVEN and idempotent: it reads each replicate's CURRENT
trajectory length and sets `steps` to exactly the remaining amount needed to reach
the target, so re-running toward the same goal is a no-op and any lagging replicate
gets topped up. Then launch `fl_optimized/run_extend.sh`.

Target (pick one):
  --target-us  N        aggregate microseconds across all replicates (DEFAULT 10)
  --target-ns-per-rep N  nanoseconds per replicate (overrides --target-us)

  # extend the 25 reps from 100 ns each (2.5 us) to 10 us aggregate = 400 ns/rep:
  python setup_extend_fl.py --target-us 10
  # or an explicit per-replicate length:
  python setup_extend_fl.py --target-ns-per-rep 400
  bash fl_optimized/run_extend.sh parallel     # launch when ready (does NOT auto-launch)

CHAIN-SPLIT GUARD: the "chain A not continuous at residue X" artifact comes from an
INPUT pdb with >1 MDAnalysis segment (a TER / chainID / segid change) -> seq_from_pdb
sets an internal terminus -> Protein.bond_check drops that backbone bond. This script
refuses to extend any replicate whose input parp14.pdb is multi-segment or whose
top.pdb has an internal TER (fix chain continuity first).

NOTE: FL is ~3x the beads of md_full; 25 x +300 ns (10 us aggregate) is a heavy,
multi-day CPU job. Scale replicates/target to the machine.
"""
import os, glob, argparse, warnings
warnings.simplefilter('ignore')
import numpy as np
from MDAnalysis import Universe

DT_PS = 0.01                       # CALVADOS timestep (ps); 1 ns = 1e5 steps
NS_PER_STEP = DT_PS / 1000.0
STEPS_PER_NS = 1.0 / NS_PER_STEP   # = 100000
CWD = os.path.dirname(os.path.abspath(__file__))
SET = os.path.join(CWD, 'fl_optimized')


def n_segments(pdb):
    return len(Universe(pdb).atoms.segments)


def internal_ter(pdb):
    lines = open(pdb).read().splitlines()
    atom_res = [l[22:26].strip() for l in lines if l.startswith('ATOM')]
    last = atom_res[-1] if atom_res else None
    return [l[22:26].strip() for l in lines
            if l.startswith('TER') and l[22:26].strip() and l[22:26].strip() != last]


def current_ns(rep):
    """Current trajectory length of a replicate in ns (0 if no dcd)."""
    dcd = glob.glob(os.path.join(rep, '*.dcd'))
    top = os.path.join(rep, 'top.pdb')
    if not dcd or not os.path.isfile(top):
        return 0.0
    try:
        u = Universe(top, dcd[0])
        return len(u.trajectory) * 0.01     # 0.01 ns/frame (wfreq 1000 * dt 0.01ps)
    except Exception:
        return 0.0


def set_steps(cfg, steps):
    lines = open(cfg).read().splitlines()
    out = [f'steps: {steps}' if (l.strip().startswith('steps:') and 'steps_eq' not in l) else l
           for l in lines]
    open(cfg, 'w').write('\n'.join(out) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target-us', type=float, default=10.0,
                    help='Aggregate target across all replicates, microseconds (default 10).')
    ap.add_argument('--target-ns-per-rep', type=float, default=None,
                    help='Per-replicate target in ns (overrides --target-us).')
    args = ap.parse_args()

    reps = sorted(glob.glob(os.path.join(SET, 'seed-*_sample-*')))
    reps = [r for r in reps if os.path.isfile(os.path.join(r, 'config.yaml'))]
    n_reps = len(reps)
    if n_reps == 0:
        print("No replicates found in fl_optimized/"); return

    target_ns = (args.target_ns_per_rep if args.target_ns_per_rep is not None
                 else args.target_us * 1000.0 / n_reps)
    print(f"fl_optimized: {n_reps} replicates")
    print(f"Target: {target_ns:.0f} ns/rep  "
          f"(= {target_ns * n_reps / 1000:.1f} us aggregate)\n")

    n_ok = n_skip = n_split = n_done = 0
    total_cur = total_extra = 0.0
    print(f"{'replicate':22} {'now(ns)':>8} {'extra(ns)':>10} {'steps':>12}  status")
    for rep in reps:
        name = os.path.basename(rep)
        chk = os.path.join(rep, 'restart.chk')
        inp = os.path.join(rep, 'input', 'parp14.pdb')
        top = os.path.join(rep, 'top.pdb')
        if not os.path.isfile(chk):
            print(f"{name:22} {'-':>8} {'-':>10} {'-':>12}  SKIP (no checkpoint)"); n_skip += 1; continue
        problem = None
        if os.path.isfile(inp) and n_segments(inp) > 1:
            problem = f"{n_segments(inp)} segments in input pdb"
        elif os.path.isfile(top) and internal_ter(top):
            problem = f"internal TER {internal_ter(top)}"
        if problem:
            print(f"{name:22} {'-':>8} {'-':>10} {'-':>12}  SPLIT-RISK: {problem}"); n_split += 1; continue

        cur = current_ns(rep); total_cur += cur
        extra = max(0.0, target_ns - cur)
        steps = int(round(extra * STEPS_PER_NS))
        if steps <= 0:
            print(f"{name:22} {cur:8.0f} {0:10.0f} {0:>12}  already at target"); n_done += 1; continue
        set_steps(os.path.join(rep, 'config.yaml'), steps)
        total_extra += extra
        print(f"{name:22} {cur:8.0f} {extra:10.0f} {steps:>12}  set"); n_ok += 1

    print(f"\n{n_ok} set to extend, {n_done} already at target, "
          f"{n_skip} no-checkpoint, {n_split} split-risk.")
    print(f"Current aggregate: {total_cur/1000:.2f} us  ->  target {target_ns*n_reps/1000:.1f} us "
          f"(+{total_extra/1000:.2f} us to run).")

    if n_ok:
        run_sh = os.path.join(SET, 'run_extend.sh')
        with open(run_sh, 'w') as f:
            f.write('#!/bin/bash\n# Extend fl_optimized from checkpoint (domain + custom restraints retained).\n')
            f.write('# Each run appends to the DCD until this round\'s target. Arg: "parallel" or "serial".\n')
            f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\nPIDS=()\n')
            f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
            f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
            f.write('  [ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; continue; }\n')
            f.write('  echo "EXTEND $d"\n')
            f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & PIDS+=($!); else ( cd "$d" && $PY run.py ); fi\n')
            f.write('done\n[ "$MODE" = parallel ] && wait\necho "fl_optimized extension done."\n')
        os.chmod(run_sh, 0o755)
        print(f"\nLauncher: {run_sh}")
        print("  bash fl_optimized/run_extend.sh parallel   # heavy, multi-day CPU job")
        print("Re-run this script anytime to top up toward the target (idempotent).")


if __name__ == '__main__':
    main()
