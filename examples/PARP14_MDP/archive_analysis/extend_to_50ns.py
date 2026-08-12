#!/usr/bin/env python3
"""
Set up PARP14 CALVADOS simulations for 50 ns total across 25 replicates.

Target: 25 replicates x 2 ns each = 50 ns total.
Equilibration: discard first 0.5 ns of each replicate.
Effective sampling: 25 x 1.5 ns = 37.5 ns compiled.

Updates config.yaml in each simulation directory. If a checkpoint exists,
CALVADOS will run 'steps' additional steps from the checkpoint and append
to the existing DCD.

Covers:
  - 25 construct sims:   parp14_seed-{1-5}_sample-{0-4}/
  - 1 full-length sim:   parp14/

Usage:
    python extend_to_50ns.py           # update configs
    python extend_to_50ns.py --dry-run # preview changes without writing
"""

import os
import yaml
from argparse import ArgumentParser

CWD = os.path.dirname(os.path.abspath(__file__))

# Timestep is 0.01 ps = 10 fs
# 1 ns = 100,000 steps
# Target per replicate: 2 ns = 200,000 steps
TARGET_STEPS = 200_000        # 2 ns per replicate
WFREQ = 1_000                 # save every 1000 steps -> 200 frames per sim

# Full-length sim: 2 ns target as well
FL_TARGET_STEPS = 200_000
FL_WFREQ = 1_000

# Equilibration: 0.5 ns = 50,000 steps = 50 frames to discard
EQUIL_NS = 0.5
EQUIL_FRAMES = 50  # at wfreq=1000


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


def main():
    parser = ArgumentParser(description='Set up PARP14 simulations for 50 ns total')
    parser.add_argument('--dry-run', action='store_true', help='Preview changes without writing')
    args = parser.parse_args()

    print("=" * 70)
    print("PARP14 Simulations: 50 ns Total (25 x 2 ns)")
    print("=" * 70)

    updated = 0
    skipped = 0

    # --- Construct simulations (25 x seed/sample) ---
    print("\n--- Construct simulations (25 replicates) ---")
    for seed in range(1, 6):
        for sample in range(0, 5):
            sim_name = f'parp14_seed-{seed}_sample-{sample}'
            config_path = os.path.join(CWD, sim_name, 'config.yaml')
            chk_path = os.path.join(CWD, sim_name, 'restart.chk')

            if not os.path.isfile(config_path):
                print(f"  SKIP (no config): {sim_name}")
                skipped += 1
                continue

            has_chk = os.path.isfile(chk_path)
            old_steps, old_wfreq = update_config(
                config_path, TARGET_STEPS, WFREQ, dry_run=args.dry_run
            )
            status = "DRY RUN" if args.dry_run else "UPDATED"
            chk_status = "checkpoint found" if has_chk else "NO CHECKPOINT"
            print(f"  {status}: {sim_name}  steps: {old_steps:,} -> {TARGET_STEPS:,}  ({chk_status})")
            updated += 1

    # --- Full-length simulation ---
    print("\n--- Full-length simulation ---")
    fl_config = os.path.join(CWD, 'parp14', 'config.yaml')
    fl_chk = os.path.join(CWD, 'parp14', 'restart.chk')
    if os.path.isfile(fl_config):
        has_chk = os.path.isfile(fl_chk)
        old_steps, old_wfreq = update_config(
            fl_config, FL_TARGET_STEPS, FL_WFREQ, dry_run=args.dry_run
        )
        status = "DRY RUN" if args.dry_run else "UPDATED"
        chk_status = "checkpoint found" if has_chk else "NO CHECKPOINT"
        print(f"  {status}: parp14 (full-length)  steps: {old_steps:,} -> {FL_TARGET_STEPS:,}  ({chk_status})")
        updated += 1
    else:
        print(f"  SKIP: parp14/ config.yaml not found")
        skipped += 1

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Updated:  {updated}")
    print(f"  Skipped:  {skipped}")
    print(f"\n  Per replicate:  {TARGET_STEPS:,} steps = {TARGET_STEPS * 0.01 / 1000:.1f} ns")
    print(f"  Write freq:     {WFREQ} steps -> {TARGET_STEPS // WFREQ} frames per sim")
    print(f"  Equilibration:  discard first {EQUIL_NS} ns = {EQUIL_FRAMES} frames")
    print(f"\n  Protocol:")
    print(f"    25 replicates x 2.0 ns = 50.0 ns total")
    print(f"    Discard first 0.5 ns of each -> 25 x 1.5 ns = 37.5 ns effective")
    if args.dry_run:
        print(f"\n  ** DRY RUN -- no files were modified **")
    else:
        print(f"\n  To run a simulation (appends to DCD if checkpoint exists):")
        print(f"    cd parp14_seed-1_sample-0 && python run.py")
        print(f"\n  To run all 25 in parallel:")
        print(f"    for d in parp14_seed-*_sample-*; do")
        print(f"      (cd $d && python run.py) &")
        print(f"    done")
    print("=" * 70)


if __name__ == '__main__':
    main()
