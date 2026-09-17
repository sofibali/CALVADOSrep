#!/usr/bin/env python
"""
Launch the extension rounds already staged in config.yaml (steps= set to the
additional amount, restart='checkpoint') for fl, md, core, mka, norrm, noart,
md3art -- 25 replicates each, 175 jobs. (fl_optimized's own extension already
ran to completion earlier; not re-launched here.)

Bounded concurrency (default 32 concurrent replicates x 4 threads each = 128
threads) so this doesn't oversubscribe the machine alongside the dashboard
backfill batch already running. Continues past any single replicate's
failure (note it, don't retry) rather than aborting the whole run.

Usage:
    python launch_staged_extensions.py [--workers 32]
"""
import argparse
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

CWD = os.path.dirname(os.path.abspath(__file__))
PY = '/home/sbali/miniconda3/envs/CALVADOS/bin/python'
SETS = ['fl', 'md', 'core', 'mka', 'norrm', 'noart', 'md3art']
SEEDS = range(1, 6)
SAMPLES = range(0, 5)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_one(set_key, seed, sample):
    rep_dir = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
    if not os.path.isfile(os.path.join(rep_dir, 'restart.chk')):
        return set_key, seed, sample, False, 'no restart.chk'
    result = subprocess.run([PY, 'run.py'], cwd=rep_dir,
                             capture_output=True, text=True)
    ok = result.returncode == 0
    err = '' if ok else result.stderr.strip()[-400:]
    return set_key, seed, sample, ok, err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=32)
    args = ap.parse_args()

    jobs = [(s, seed, sample) for s in SETS for seed in SEEDS for sample in SAMPLES]
    log(f"Launching {len(jobs)} replicate extensions across {len(SETS)} sets, "
        f"{args.workers} concurrent (~{args.workers * 4} threads)")

    n_ok = n_fail = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, *j): j for j in jobs}
        for fut in as_completed(futures):
            set_key, seed, sample, ok, err = fut.result()
            if ok:
                n_ok += 1
            else:
                n_fail += 1
                log(f"  FAILED {set_key}/seed-{seed}_sample-{sample}: {err}")
            if (n_ok + n_fail) % 10 == 0:
                log(f"  progress: {n_ok + n_fail}/{len(jobs)} done ({n_fail} failed)")

    log(f"=== DONE: {n_ok}/{len(jobs)} OK, {n_fail} failed ===")


if __name__ == '__main__':
    main()
