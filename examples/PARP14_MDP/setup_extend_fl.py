#!/usr/bin/env python
"""
Extend the full-length PARP14 runs (fl_optimized/) to long trajectories, keeping the
domain restraints AND the KH7a-KHb custom restraints already configured there.

fl_optimized already has 25 replicates at 20 ns with domain restraints + custom
restraints (input/custom_restraints.txt) and a continuous single chain. CALVADOS
checkpoint restart runs `steps` ADDITIONAL steps and appends to the DCD, so we just
raise `steps` and re-run run.py (same mechanism as setup_extend_md_full.py).

Target: 100 ns/replicate (= +80 ns now). FL is ~3x the beads of md_full, so this is a
heavy job — adjust with --extra-ns and launch when ready (this script does NOT launch).

CHAIN-SPLIT GUARD (the bug you saw):
  The "chain A not continuous at residue X" artifact comes from the INPUT pdb having
  >1 MDAnalysis segment (a TER / chainID / segid change). calvados.sequence.seq_from_pdb
  then sets an internal N/C-terminus at the segment break, and Protein.bond_check drops
  the backbone bond there (components.py), so top.pdb shows a chain break. It is NOT
  caused by the custom restraints themselves. This script verifies every replicate's
  input parp14.pdb is a SINGLE segment and that top.pdb (if present) has no internal TER,
  and refuses to set up any replicate that would split. (To repair a multi-segment pdb:
  rewrite all ATOM records to chain A / segid A and keep a single trailing TER.)

    python setup_extend_fl.py              # +80 ns -> 100 ns total, with guard
    python setup_extend_fl.py --extra-ns 180   # -> 200 ns total
    bash fl_optimized/run_extend.sh parallel   # launch when ready
"""
import os, glob, argparse, warnings
warnings.simplefilter('ignore')
from MDAnalysis import Universe

DT_PS = 0.01
CWD = os.path.dirname(os.path.abspath(__file__))
SET = os.path.join(CWD, 'fl_optimized')

def n_segments(pdb):
    return len(Universe(pdb).atoms.segments)

def internal_ter(pdb):
    """Return resSeq of any TER that is not the final residue (an internal chain break)."""
    lines = open(pdb).read().splitlines()
    last_atom_res = None; bad = []
    atom_res = [l[22:26].strip() for l in lines if l.startswith('ATOM')]
    last = atom_res[-1] if atom_res else None
    for l in lines:
        if l.startswith('TER'):
            r = l[22:26].strip()
            if r and r != last:
                bad.append(r)
    return bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--extra-ns', type=float, default=80.0,
                    help='Additional ns per replicate this round (default 80 -> 100 ns total).')
    args = ap.parse_args()
    extra_steps = int(round(args.extra_ns / (DT_PS / 1000)))

    reps = sorted(glob.glob(os.path.join(SET, 'seed-*_sample-*')))
    n_ok = n_skip = n_split = 0
    for rep in reps:
        cfg = os.path.join(rep, 'config.yaml')
        chk = os.path.join(rep, 'restart.chk')
        inp = os.path.join(rep, 'input', 'parp14.pdb')
        top = os.path.join(rep, 'top.pdb')
        name = os.path.basename(rep)
        if not (os.path.isfile(cfg) and os.path.isfile(chk)):
            print(f"  SKIP {name}: missing config/checkpoint"); n_skip += 1; continue

        # --- chain-split guard ---
        problem = None
        if os.path.isfile(inp) and n_segments(inp) > 1:
            problem = f"input parp14.pdb has {n_segments(inp)} segments"
        elif os.path.isfile(top) and internal_ter(top):
            problem = f"top.pdb has internal TER at residue(s) {internal_ter(top)}"
        if problem:
            print(f"  SPLIT-RISK {name}: {problem} -> NOT extending (fix chain continuity first)")
            n_split += 1; continue

        lines = open(cfg).read().splitlines()
        out = [f'steps: {extra_steps}' if (l.strip().startswith('steps:') and 'steps_eq' not in l) else l
               for l in lines]
        open(cfg, 'w').write('\n'.join(out) + '\n')
        n_ok += 1
    print(f"\nSet steps = {extra_steps} (+{args.extra_ns:g} ns) in {n_ok} replicates; "
          f"{n_skip} missing, {n_split} split-risk (skipped).")

    if n_ok:
        run_sh = os.path.join(SET, 'run_extend.sh')
        with open(run_sh, 'w') as f:
            f.write('#!/bin/bash\n# Extend fl_optimized from checkpoint (domain + custom restraints retained).\n')
            f.write('# Each run appends to the DCD. Arg: "parallel" or "serial" (default).\n')
            f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\nPIDS=()\n')
            f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
            f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
            f.write('  [ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; continue; }\n')
            f.write('  echo "EXTEND $d"\n')
            f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & PIDS+=($!); else ( cd "$d" && $PY run.py ); fi\n')
            f.write('done\n[ "$MODE" = parallel ] && wait\necho "fl_optimized extension done."\n')
        os.chmod(run_sh, 0o755)
        print(f"Launcher: {run_sh}  ->  bash fl_optimized/run_extend.sh parallel")
        print("NOTE: FL is ~3x md_full beads; 25 x +%g ns is a heavy CPU job." % args.extra_ns)

if __name__ == '__main__':
    main()
