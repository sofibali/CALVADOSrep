#!/usr/bin/env python
"""
Production-step bookkeeping for one slab run directory.

WHY
---
CALVADOS restarting from `restart.chk` runs `steps` ADDITIONAL steps and
appends to the DCD (sim.py, `append=True`). That is the right behaviour for
picking up a killed job, but it means re-running a queue over a tree where
some constructs already finished would silently push those past the 2e8-step
target -- and re-running a half-finished one overshoots by whatever it had
already done. Nothing in the tree guards against that.

This reads how many production steps a run has actually completed and
rewrites `config.yaml`'s `steps` to the remainder, so a queue can be
re-run over the same tree as many times as needed and every construct
still lands on exactly its target.

`slab_meta.yaml` holds the authoritative target; `config.yaml`'s `steps`
is the working value this script edits.

The step count comes from the StateDataReporter log, which is attached
*after* slab equilibration -- so its Step column is production steps only
and the 5e6 equilibration steps are correctly excluded.

USAGE
-----
    python slab_steps.py status  <run_dir>              # done / target / remaining
    python slab_steps.py prepare <run_dir>              # patch config.yaml; rc=3 if complete
    python slab_steps.py prepare <run_dir> --max-leg N  # cap this leg at N steps

`--max-leg` is what makes a run preemptible. sim.py checkpoints once per
batch and splits a leg into 10 batches, so the checkpoint interval is
leg/10 -- that is the most work a SIGTERM can destroy. A full 2e8-step leg
checkpoints only every 2e7 steps (~1 h), which is far too coarse to yield a
GPU politely; a 2e7-step leg checkpoints every 2e6 (~7 min), which is fine.
"""
import os
import sys

import yaml

# sim.py splits production into 10 checkpointed batches, so keep the value a
# multiple of 10 to avoid losing steps to integer division.
NBATCHES = 10


def read_meta(run_dir):
    fmeta = os.path.join(run_dir, 'slab_meta.yaml')
    if not os.path.isfile(fmeta):
        raise SystemExit(f'no slab_meta.yaml in {run_dir}')
    return yaml.safe_load(open(fmeta))


def steps_done(run_dir, sysname):
    """Production steps completed, from the reporter log (0 if never started)."""
    flog = os.path.join(run_dir, f'{sysname}.log')
    if not os.path.isfile(flog):
        return 0
    last = 0
    with open(flog, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                last = max(last, int(float(line.split('\t')[0])))
            except (ValueError, IndexError):
                continue
    return last


def main():
    argv = sys.argv[1:]
    max_leg = None
    if '--max-leg' in argv:
        i = argv.index('--max-leg')
        try:
            max_leg = int(float(argv[i + 1]))
        except (IndexError, ValueError):
            raise SystemExit('--max-leg needs a step count')
        del argv[i:i + 2]
    if len(argv) != 2 or argv[0] not in ('status', 'prepare'):
        raise SystemExit(__doc__.strip())
    cmd, run_dir = argv[0], os.path.abspath(argv[1])

    meta = read_meta(run_dir)
    sysname = meta['sysname']
    target = int(float(meta['steps']))
    done = steps_done(run_dir, sysname)
    remaining = max(0, target - done)

    if cmd == 'status':
        print(f'{done} {target} {remaining}')
        return

    if remaining == 0:
        print(f'  [steps] {sysname}: {done:,}/{target:,} production steps -- COMPLETE, skipping')
        sys.exit(3)

    leg = remaining if max_leg is None else min(remaining, max_leg)
    leg = max(NBATCHES, (leg // NBATCHES) * NBATCHES)

    fcfg = os.path.join(run_dir, 'config.yaml')
    cfg = yaml.safe_load(open(fcfg))
    if int(float(cfg.get('steps', 0))) != leg:
        cfg['steps'] = leg
        with open(fcfg, 'w') as fh:
            yaml.safe_dump(cfg, fh)

    tail = f' (checkpoint every {leg // NBATCHES:,})' if max_leg else ''
    if done:
        print(f'  [steps] {sysname}: resuming at {done:,}/{target:,}; '
              f'this leg runs {leg:,} steps{tail}')
    else:
        print(f'  [steps] {sysname}: fresh start, {leg:,} of {target:,} production '
              f'steps{tail} (after {int(float(cfg.get("steps_eq", 0))):,} equilibration steps)')


if __name__ == '__main__':
    main()
