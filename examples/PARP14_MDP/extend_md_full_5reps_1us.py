#!/usr/bin/env python
"""
Extend 5 of md_full's 25 existing replicates (seed-{1..5}_sample-0 -- one per
seed) from their current checkpoint to a cumulative 1 us total, to pair with
prepare_mka_full.py's fresh 5 x 1us MD1-ART replicates.

Every md_full replicate has already run to 100 ns (10,000,000 steps; verified
via each replicate's .log). CALVADOS checkpoint restart runs `steps` ADDITIONAL
steps from restart.chk and appends to the DCD (calvados/sim.py), so reaching
1000 ns total from here means +900 ns = 90,000,000 more steps.

Only touches the 5 chosen replicates' config.yaml (steps, threads) -- the
other 20 md_full replicates are left exactly as prepare_md_full.py /
setup_extend_md_full.py last left them.

    python extend_md_full_5reps_1us.py       # set steps for the 5 chosen reps
    bash md_full/run_5reps_1us.sh parallel   # then launch
"""
import os

CWD = os.path.dirname(os.path.abspath(__file__))
DT_PS = 0.01  # CALVADOS timestep (ps); 1 ns = 100,000 steps

CHOSEN_REPS = [f'seed-{s}_sample-0' for s in range(1, 6)]
CURRENT_NS = 100.0
TARGET_NS = 1000.0
EXTRA_STEPS = int(round((TARGET_NS - CURRENT_NS) / (DT_PS / 1000)))
THREADS = 16


def main():
    n = 0
    for rep in CHOSEN_REPS:
        rep_dir = os.path.join(CWD, 'md_full', rep)
        cfg_path = os.path.join(rep_dir, 'config.yaml')
        chk_path = os.path.join(rep_dir, 'restart.chk')
        if not os.path.isfile(chk_path):
            print(f"  SKIP {rep}: no restart.chk")
            continue
        lines = open(cfg_path).read().splitlines()
        out = []
        for ln in lines:
            if ln.strip().startswith('steps:') and 'steps_eq' not in ln:
                out.append(f'steps: {EXTRA_STEPS}')
            elif ln.strip().startswith('threads:'):
                out.append(f'threads: {THREADS}')
            else:
                out.append(ln)
        open(cfg_path, 'w').write('\n'.join(out) + '\n')
        print(f"  {rep}: steps -> {EXTRA_STEPS} (+{TARGET_NS - CURRENT_NS:.0f} ns, "
              f"-> {TARGET_NS:.0f} ns total), threads -> {THREADS}")
        n += 1

    print(f"\nSet {n}/{len(CHOSEN_REPS)} replicates for the 1us extension.")

    run_sh = os.path.join(CWD, 'md_full', 'run_5reps_1us.sh')
    with open(run_sh, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# Extend md_full seed-{1..5}_sample-0 from checkpoint to 1us total.\n')
        f.write('# Each run appends to the DCD. Arg: "parallel" (all at once) or "serial" (default).\n')
        f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\n')
        f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
        for rep in CHOSEN_REPS:
            f.write(f'd="$HERE/{rep}/"\n')
            f.write('[ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; exit 1; }\n')
            f.write('echo "EXTEND $d"\n')
            f.write('if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi\n')
        f.write('[ "$MODE" = parallel ] && wait\necho "md_full 5-rep 1us extension done."\n')
    os.chmod(run_sh, 0o755)
    print(f"Launcher: {run_sh}  ->  bash {run_sh} parallel")


if __name__ == '__main__':
    main()
