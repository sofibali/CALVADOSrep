#!/usr/bin/env python
"""
Set up the md_full extension. CALVADOS checkpoint restart runs `steps` ADDITIONAL
steps and APPENDS to the DCD (calvados/sim.py: loadCheckpoint + simulation.step(steps);
the DCD is only backed up when NOT restarting from checkpoint).

So to extend each replicate we just raise `steps` in its config.yaml to the desired
EXTRA length and re-run run.py (restart='checkpoint', frestart='restart.chk' already set).

Why extend: the 5 ns runs are too short for the 10-20+ ns slow processes -> MSM ITS
never plateau (see sweep_its.py results). Target 100 ns/replicate (= +95 ns now) so
trajectories exceed the slow timescale several-fold and the kinetic analysis converges.

    python setup_extend_md_full.py          # set steps for +95 ns (-> 100 ns total)
    python setup_extend_md_full.py --extra-ns 45   # smaller round (-> 50 ns total)
    bash md_full/run_extend.sh parallel     # then launch
Re-running run_extend.sh adds another EXTRA_NS each time (rounds are cumulative).
"""
import os, glob, argparse

DT_PS = 0.01                      # CALVADOS timestep (ps); 1 ns = 100000 steps
CWD = os.path.dirname(os.path.abspath(__file__))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--extra-ns', type=float, default=95.0,
                    help='Additional ns to run per replicate this round (default 95 -> 100 ns total).')
    args = ap.parse_args()
    extra_steps = int(round(args.extra_ns / (DT_PS / 1000)))   # ns -> steps

    configs = sorted(glob.glob(os.path.join(CWD, 'md_full', 'seed-*_sample-*', 'config.yaml')))
    n = 0
    for cfg in configs:
        rep_dir = os.path.dirname(cfg)
        if not os.path.isfile(os.path.join(rep_dir, 'restart.chk')):
            print(f"  SKIP {os.path.basename(rep_dir)}: no restart.chk"); continue
        lines = open(cfg).read().splitlines()
        out = []
        for ln in lines:
            if ln.strip().startswith('steps:') and 'steps_eq' not in ln:
                out.append(f'steps: {extra_steps}')
            else:
                out.append(ln)
        open(cfg, 'w').write('\n'.join(out) + '\n')
        n += 1
    print(f"Set steps = {extra_steps} (+{args.extra_ns:g} ns) in {n} md_full replicates.")

    # launcher (mirrors run_extend.sh)
    run_sh = os.path.join(CWD, 'md_full', 'run_extend.sh')
    with open(run_sh, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# Extend md_full from checkpoint (+EXTRA_NS set by setup_extend_md_full.py).\n')
        f.write('# Each run appends to the DCD. Arg: "parallel" (all at once) or "serial" (default).\n')
        f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\nPIDS=()\n')
        f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
        f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
        f.write('  [ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; continue; }\n')
        f.write('  echo "EXTEND $d"\n')
        f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & PIDS+=($!); else ( cd "$d" && $PY run.py ); fi\n')
        f.write('done\n[ "$MODE" = parallel ] && wait\necho "md_full extension done."\n')
    os.chmod(run_sh, 0o755)
    print(f"Launcher: {run_sh}  ->  bash md_full/run_extend.sh parallel")
    print("After it finishes, re-run: python sweep_its.py --set md_full --vamp --cv "
          "--its-lags 80 200 400 800 1600 --skip inter_exposed")

if __name__ == '__main__':
    main()
