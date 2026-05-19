#!/usr/bin/env python3
"""
Extend PARP14 CALVADOS simulations from checkpoint.

Targets (round 2 — match all at 20 ns):
  - FL:           already at 20 ns, no extension needed
  - md/core/mka:  currently 10 ns -> 20 ns = +10 ns = 1,000,000 additional steps

Directory layout (nested):
  {set}/seed-{N}_sample-{M}/

CALVADOS checkpoint restart runs 'steps' ADDITIONAL steps and appends to DCD.

Usage:
    python extend_sims.py                    # update all configs
    python extend_sims.py --dry-run          # preview only
    python extend_sims.py --set fl           # update FL only
    python extend_sims.py --write-launcher   # also write run_extend.sh
"""

import os
import yaml
from argparse import ArgumentParser

CWD = os.path.dirname(os.path.abspath(__file__))

# Timestep = 0.01 ps = 10 fs -> 100,000 steps = 1 ns
# FL:  already at 20 ns (2000 frames) — skip
# MD/Core/MKA: currently at 10 ns (1000 frames) -> +10 ns -> 20 ns
#   +10 ns = 1,000,000 additional steps at wfreq=1000 -> +1000 frames -> 2000 total

EXTEND_CONFIG = {
    'fl':           {'additional_steps':           0, 'wfreq': 1_000, 'target_ns': 20},
    'md':           {'additional_steps': 1_000_000,   'wfreq': 1_000, 'target_ns': 20},
    'core':         {'additional_steps': 1_000_000,   'wfreq': 1_000, 'target_ns': 20},
    'mka':          {'additional_steps': 1_000_000,   'wfreq': 1_000, 'target_ns': 20},
    'fl_optimized': {'additional_steps':           0, 'wfreq': 1_000, 'target_ns': 20},
}

SEEDS = range(1, 6)
SAMPLES = range(0, 5)


def update_config(config_path, new_steps, new_wfreq, dry_run=False):
    """Update steps and wfreq in a config.yaml file."""
    with open(config_path) as f:
        config = yaml.safe_load(f)

    old_steps = config.get('steps', 0)
    old_wfreq = config.get('wfreq', 0)

    config['steps'] = new_steps
    config['wfreq'] = new_wfreq

    if not dry_run:
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    return old_steps, old_wfreq


def write_launcher(active_sets):
    """Write a bash launcher script for the extension runs."""
    launcher_path = os.path.join(CWD, 'run_extend.sh')
    with open(launcher_path, 'w') as f:
        f.write("""#!/bin/bash
# Extend PARP14 CALVADOS simulations from checkpoint.
# Round 2: md/core/mka 10 ns -> 20 ns (+10 ns). FL already at 20 ns (skipped).
#
# Usage:
#   bash run_extend.sh              # all sets, sequential
#   bash run_extend.sh parallel     # all sets, parallel
#   bash run_extend.sh fl           # FL only
#   bash run_extend.sh md           # md only
#   bash run_extend.sh core         # core only
#   bash run_extend.sh mka          # mka only

set -e
cd "$(dirname "$0")"

MODE="${1:-sequential}"
PIDS=()

run_sim() {
    local dir="$1"
    if [ ! -f "$dir/config.yaml" ]; then
        echo "SKIP: $dir (no config.yaml)"
        return
    fi
    if [ ! -f "$dir/restart.chk" ]; then
        echo "SKIP: $dir (no checkpoint)"
        return
    fi
    echo "EXTEND: $dir"

    if [ "$MODE" = "parallel" ]; then
        (cd "$dir" && python run.py) &
        PIDS+=($!)
    else
        (cd "$dir" && python run.py)
    fi
}

run_set() {
    local prefix="$1"
    local label="$2"
    echo "=================================================="
    echo "Extending: $label"
    echo "=================================================="
    for seed in 1 2 3 4 5; do
        for sample in 0 1 2 3 4; do
            run_sim "${prefix}/seed-${seed}_sample-${sample}"
        done
    done
}

case "$MODE" in
    fl)       run_set "fl" "Full-Length (-> 20 ns)" ;;
    md)       run_set "md" "Macrodomains (-> 20 ns)" ;;
    core)     run_set "core" "Core (-> 20 ns)" ;;
    mka)      run_set "mka" "MKA (-> 20 ns)" ;;
    parallel)
""")
        for s in active_sets:
            cfg = EXTEND_CONFIG[s]
            f.write(f'        run_set "{s}" "{s.upper()} (-> {cfg["target_ns"]} ns)"\n')
        f.write("""        ;;
    *)
""")
        for s in active_sets:
            cfg = EXTEND_CONFIG[s]
            f.write(f'        run_set "{s}" "{s.upper()} (-> {cfg["target_ns"]} ns)"\n')
        f.write("""        ;;
esac

if [ "$MODE" = "parallel" ] && [ ${#PIDS[@]} -gt 0 ]; then
    echo ""
    echo "Waiting for ${#PIDS[@]} parallel jobs..."
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
fi

echo ""
echo "=================================================="
echo "EXTENSION COMPLETE"
echo "=================================================="
""")

    os.chmod(launcher_path, 0o755)
    print(f"\n  Wrote: {launcher_path}")
    print(f"  Usage:")
    print(f"    bash run_extend.sh              # all sets, sequential")
    print(f"    bash run_extend.sh parallel     # all sets, parallel")
    for s in active_sets:
        print(f"    bash run_extend.sh {s:<15} # {s} only")


def main():
    parser = ArgumentParser(description='Extend PARP14 simulations from checkpoint')
    parser.add_argument('--dry-run', action='store_true', help='Preview changes without writing')
    parser.add_argument('--set', nargs='+', default=None,
                        choices=['fl', 'md', 'core', 'mka', 'fl_optimized'],
                        help='Specific sets to extend (default: all)')
    parser.add_argument('--write-launcher', action='store_true',
                        help='Write run_extend.sh launcher script')
    args = parser.parse_args()

    active_sets = args.set if args.set else ['fl', 'md', 'core', 'mka',
                                              'fl_optimized']

    # Filter out sets that need no extension
    need_extension = [s for s in active_sets if EXTEND_CONFIG[s]['additional_steps'] > 0]
    already_done = [s for s in active_sets if EXTEND_CONFIG[s]['additional_steps'] == 0]

    print("=" * 70)
    print("PARP14 Simulation Extension (round 2: all to 20 ns)")
    print("=" * 70)
    for s in already_done:
        cfg = EXTEND_CONFIG[s]
        print(f"  {s:<6}: ALREADY at {cfg['target_ns']} ns — skipping")
    for s in need_extension:
        cfg = EXTEND_CONFIG[s]
        add_ns = cfg['additional_steps'] * 0.01 / 1000
        print(f"  {s:<6}: +{add_ns:.0f} ns ({cfg['additional_steps']:,} steps) -> {cfg['target_ns']} ns total")

    updated = 0
    skipped = 0

    for set_key in need_extension:
        cfg = EXTEND_CONFIG[set_key]
        print(f"\n--- {set_key.upper()} (-> {cfg['target_ns']} ns) ---")

        for seed in SEEDS:
            for sample in SAMPLES:
                rep_name = f'seed-{seed}_sample-{sample}'
                sim_dir = os.path.join(CWD, set_key, rep_name)
                config_path = os.path.join(sim_dir, 'config.yaml')
                chk_path = os.path.join(sim_dir, 'restart.chk')

                if not os.path.isfile(config_path):
                    skipped += 1
                    continue

                has_chk = os.path.isfile(chk_path)
                if not has_chk:
                    print(f"  WARN: {set_key}/{rep_name} — no checkpoint, will run from scratch")

                old_steps, old_wfreq = update_config(
                    config_path, cfg['additional_steps'], cfg['wfreq'],
                    dry_run=args.dry_run
                )
                tag = "DRY RUN" if args.dry_run else "UPDATED"
                print(f"  {tag}: {set_key}/{rep_name}  steps: {old_steps:,} -> {cfg['additional_steps']:,}")
                updated += 1

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  Updated: {updated}  |  Skipped (no config): {skipped}  |  Already done: {len(already_done)*25}")
    print()
    for s in active_sets:
        cfg = EXTEND_CONFIG[s]
        if cfg['additional_steps'] > 0:
            add_ns = cfg['additional_steps'] * 0.01 / 1000
            n_frames = cfg['additional_steps'] // cfg['wfreq']
            print(f"  {s:<6}: {cfg['additional_steps']:>10,} additional steps = +{add_ns:.0f} ns = +{n_frames} frames")
        else:
            print(f"  {s:<6}: no extension needed")
    print()
    print(f"  After extension (all sets):")
    for s in active_sets:
        cfg = EXTEND_CONFIG[s]
        total_frames = cfg['target_ns'] * 100
        print(f"    {s:<6}: {cfg['target_ns']} ns total, {total_frames} frames/rep, "
              f"{total_frames * 25:,} frames cumulative (25 reps)")
    total_all = sum(EXTEND_CONFIG[s]['target_ns'] * 100 * 25 for s in active_sets)
    print(f"\n  TOTAL: {total_all:,} frames across all {len(active_sets)} sets "
          f"({total_all * 0.01:.0f} ns cumulative)")

    if args.dry_run:
        print(f"\n  ** DRY RUN — no files were modified **")
    elif need_extension:
        print(f"\n  Configs updated. Run simulations with:")
        print(f"    bash run_extend.sh parallel     # md/core/mka in parallel")
        for s in need_extension:
            print(f"    bash run_extend.sh {s:<15} # {s} only")
    else:
        print(f"\n  Nothing to extend — all sets already at target.")

    if args.write_launcher:
        write_launcher(need_extension if need_extension else active_sets)

    print(f"{'='*70}")


if __name__ == '__main__':
    main()
