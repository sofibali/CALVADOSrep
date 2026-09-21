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
    # sim.py runs `nbatches` batches of `int(steps/nbatches)` steps, so a leg
    # must be a whole multiple of NBATCHES or steps are silently dropped --
    # and a leg BELOW NBATCHES gives batch=0, i.e. the run does nothing at all,
    # `done` never advances, and the opportunistic runner relaunches forever.
    #
    # Neither rounding direction alone is safe: rounding up always would
    # overshoot the target on the last leg, rounding down always strands a
    # 1..NBATCHES-1 remainder that can never be consumed. So: round down while
    # that leaves a runnable leg, and only for a final sub-batch remainder
    # round up to NBATCHES, overshooting by <10 steps out of 2e8 rather than
    # livelocking.
    if leg >= NBATCHES:
        leg = (leg // NBATCHES) * NBATCHES
    else:
        print(f'  [steps] {sysname}: {leg} step(s) left is below sim.py\'s '
              f'{NBATCHES}-batch floor; running {NBATCHES} to finish '
              f'(overshoots target by {NBATCHES - leg})')
        leg = NBATCHES

    fcfg = os.path.join(run_dir, 'config.yaml')
    cfg = yaml.safe_load(open(fcfg))
    changed = int(float(cfg.get('steps', 0))) != leg

    # For a PREEMPTIBLE leg (--max-leg, i.e. the opportunistic runner), pin
    # logfreq to the checkpoint interval. sim.py checkpoints after each of its
    # 10 batches, so with logfreq == leg/10 the last line in {sysname}.log is
    # written at exactly the step the checkpoint holds. trim_dcd.py relies on
    # that to know which trajectory frames the checkpoint actually backs, and
    # therefore which are re-simulated leftovers to drop.
    #
    # Not done for a full-target leg: there logfreq would become 2e7 and the
    # run would log only 10 times in 2e8 steps, which is useless to monitor --
    # and a run that never yields has nothing to trim anyway.
    if max_leg is not None:
        want_logfreq = max(1, leg // NBATCHES)
        if int(float(cfg.get('logfreq', 0))) != want_logfreq:
            cfg['logfreq'] = want_logfreq
            changed = True

    if changed:
        cfg['steps'] = leg
        # Write atomically. The opportunistic runner SIGKILLs stragglers and
        # these jobs get killed on logout; a kill partway through an in-place
        # `open(fcfg,'w')` would leave config.yaml empty or half-written and
        # permanently brick the construct.
        tmp = fcfg + '.tmp'
        with open(tmp, 'w') as fh:
            yaml.safe_dump(cfg, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, fcfg)

    tail = f' (checkpoint every {leg // NBATCHES:,})' if max_leg else ''
    if done:
        print(f'  [steps] {sysname}: resuming at {done:,}/{target:,}; '
              f'this leg runs {leg:,} steps{tail}')
    else:
        print(f'  [steps] {sysname}: fresh start, {leg:,} of {target:,} production '
              f'steps{tail} (after {int(float(cfg.get("steps_eq", 0))):,} equilibration steps)')


if __name__ == '__main__':
    main()
