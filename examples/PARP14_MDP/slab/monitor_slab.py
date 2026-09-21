#!/usr/bin/env python
"""
Monitor the PARP14 slab campaign on lyra.

Reads each run's StateDataReporter log (`{sysname}.log`: Step / Speed /
Elapsed) and finds live runs in /proc, then prints progress, throughput and
an ETA per run, alongside what the shared GPUs are doing.

Works with whatever started the run -- `run_slab_queue.sh`, a bare
`python run.py`, or a resumed leg. There is no state file to go stale.

    python monitor_slab.py                 # one snapshot
    python monitor_slab.py --watch         # refresh every 60 s
    python monitor_slab.py --watch -n 15   # refresh every 15 s
    python monitor_slab.py --arm both

Phases you will see:

  equilibrating   inside the 5e6-step slab_eq pull. No reporter is attached
                  yet (sim.py builds a fresh Simulation once the centering
                  force is removed), so there is NO step output at all
                  during this phase -- a run can sit here for many minutes
                  looking idle while being perfectly healthy.
  running         production; steps, ns/day and ETA are live.
  stalled         process alive but the log has not advanced in 15 min.
  done            production steps >= the slab_meta target.
"""
import argparse

import os
import subprocess
import time
from datetime import datetime, timedelta

import yaml

SLAB_ROOT = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(SLAB_ROOT, 'logs')

ARMS = ('homotypic', 'rna')

DT_PS = 0.01  # integrator timestep, ps -- config.yaml is fixed at this


def _parse_rows(path):
    """
    [(step, speed, elapsed), ...] from a StateDataReporter log.

    Column positions are taken from the `#"Step"...` header rather than
    assumed. With `report_potential_energy: true` the reporter emits potential
    energy BEFORE speed, so fixed indices would silently return the energy as
    the speed and the speed as the elapsed time -- every rate and ETA would be
    nonsense with no error. The header is rewritten whenever append=False, so
    the last one in the file describes the current columns.
    """
    if not os.path.isfile(path):
        return []
    idx = {'step': 0, 'speed': 1, 'elapsed': 2}
    rows = []
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#'):
                cols = [c.strip().strip('"').lower() for c in line.lstrip('#').split('\t')]
                found = {}
                for i, c in enumerate(cols):
                    if c.startswith('step'):
                        found['step'] = i
                    elif c.startswith('speed'):
                        found['speed'] = i
                    elif c.startswith('elapsed'):
                        found['elapsed'] = i
                if len(found) == 3:
                    idx = found
                continue
            parts = line.split('\t')
            if len(parts) <= max(idx.values()):
                continue
            try:
                rows.append((float(parts[idx['step']]),
                             float(parts[idx['speed']]),
                             float(parts[idx['elapsed']])))
            except ValueError:
                continue
    return rows


def read_log(path):
    """(last_step, last_speed_ns_day, last_elapsed_s, n_rows) from a reporter log."""
    rows = _parse_rows(path)
    if not rows:
        return 0, 0.0, 0.0, 0
    return int(max(r[0] for r in rows)), rows[-1][1], rows[-1][2], len(rows)


def log_interval(path, window=5):
    """Median wall-clock seconds between log lines, for the stall threshold."""
    rows = _parse_rows(path)
    if len(rows) < 2:
        return None
    deltas = [b[2] - a[2] for a, b in zip(rows[-(window + 1):], rows[-window:])
              if b[2] > a[2]]
    if not deltas:
        return None
    deltas.sort()
    return deltas[len(deltas) // 2]


def recent_rate(path, window=5):
    """
    ns/day over the last `window` log intervals.

    Two traps here:

    * The Speed column is a running average since the start of the LEG, so it
      lags badly when a GPU becomes contended mid-run.
    * The Elapsed column restarts at ~0 on every leg, and an opportunistic run
      is many short legs. Differencing the first and last row of a window that
      straddles a leg boundary divides a real step count by a near-zero time
      and reports a rate several times the true one.

    So: accumulate interval by interval and drop any interval whose elapsed
    time did not advance (i.e. a leg boundary).
    """
    rows = [(r[0], r[2]) for r in _parse_rows(path)]
    if len(rows) < 2:
        return 0.0
    dsteps = dt = 0.0
    for (s0, t0), (s1, t1) in list(zip(rows, rows[1:]))[-window:]:
        if t1 <= t0 or s1 <= s0:      # leg boundary, or no progress
            continue
        dsteps += s1 - s0
        dt += t1 - t0
    if dt <= 0 or dsteps <= 0:
        return 0.0
    return (dsteps * DT_PS * 1e-3) / dt * 86400.0  # ns/day


def gpu_table():
    try:
        out = subprocess.check_output(
            ['nvidia-smi', '--query-gpu=index,name,utilization.gpu,memory.used,memory.total',
             '--format=csv,noheader,nounits'], text=True)
    except (OSError, subprocess.CalledProcessError):
        return []
    rows = []
    for line in out.strip().splitlines():
        idx, name, util, used, total = [x.strip() for x in line.split(',')]
        rows.append((int(idx), name, int(util), int(used), int(total)))
    return rows


def running_runs():
    """
    Map run directory -> {pid, gpu} for every live `python run.py`.

    Read from /proc rather than a state file, so this works no matter how the
    run was started -- run_slab_queue.sh, a bare `python run.py`, or a
    resubmitted leg -- and never goes stale.
    """
    found = {}
    for entry in os.listdir('/proc'):
        if not entry.isdigit():
            continue
        try:
            with open(f'/proc/{entry}/cmdline', 'rb') as fh:
                cmd = fh.read().decode(errors='replace').split('\0')
            if not any(c.endswith('run.py') for c in cmd):
                continue
            cwd = os.path.realpath(f'/proc/{entry}/cwd')
            gpu = '-'
            with open(f'/proc/{entry}/environ', 'rb') as fh:
                for kv in fh.read().decode(errors='replace').split('\0'):
                    if kv.startswith('CUDA_VISIBLE_DEVICES='):
                        gpu = kv.split('=', 1)[1] or '-'
                        break
            found[cwd] = {'pid': int(entry), 'gpu': gpu}
        except (OSError, PermissionError, ValueError):
            continue
    return found


def collect(arms):
    live = running_runs()
    rows = []
    for arm in arms:
        arm_dir = os.path.join(SLAB_ROOT, arm)
        if not os.path.isdir(arm_dir):
            continue
        for name in sorted(os.listdir(arm_dir)):
            d = os.path.join(arm_dir, name)
            fmeta = os.path.join(d, 'slab_meta.yaml')
            if not os.path.isfile(fmeta):
                continue
            meta = yaml.safe_load(open(fmeta))
            sysname = meta['sysname']
            target = int(float(meta['steps']))
            flog = os.path.join(d, f'{sysname}.log')
            step, avg_speed, elapsed, nrows = read_log(flog)
            rate = recent_rate(flog) or avg_speed

            job = live.get(os.path.realpath(d), {})
            alive = bool(job)

            # run_slab_queue.sh writes logs/<arm>_<construct>.log
            fout = os.path.join(LOG_DIR, f'{arm}_{name}.log')
            in_eq = False
            if alive and os.path.isfile(fout):
                try:
                    txt = open(fout, errors='replace').read()
                    # only the current leg matters -- earlier legs already
                    # passed through equilibration
                    leg = txt.rsplit('Starting slab equilibration', 1)
                    in_eq = len(leg) > 1 and 'STARTING SIMULATION' not in leg[-1]
                except OSError:
                    pass

            if step >= target:
                status = 'done'
            elif alive and in_eq:
                status = 'equilibrating'
            elif alive:
                # The stall threshold must scale with how often this run
                # actually logs. logfreq is 1e6 steps; on a contended card
                # (~830 steps/s, per the README's own table) that is ~20 min
                # between lines, so a fixed 900 s cutoff would label every
                # healthy large construct 'stalled'. Use 3x the observed
                # interval, with 900 s only as a floor.
                mtime = os.path.getmtime(flog) if os.path.isfile(flog) else 0
                iv = log_interval(flog)
                limit = max(900.0, 3.0 * iv) if iv else 900.0
                status = 'running' if time.time() - mtime < limit else 'stalled'
            elif step > 0:
                status = 'partial'
            else:
                status = 'queued'

            eta = ''
            if status == 'running' and rate > 0:
                ns_left = (target - step) * DT_PS * 1e-3
                days = ns_left / rate
                eta = str(timedelta(days=days)).split('.')[0]
            rows.append({
                'arm': arm, 'construct': name, 'beads': int(meta['beads']),
                'step': step, 'target': target, 'rate': rate,
                'status': status, 'gpu': job.get('gpu', '-'),
                'pid': job.get('pid', '-') if alive else '-', 'eta': eta,
            })
    return rows


def render(rows, show_gpu=True):
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    out = [f'PARP14 slab campaign -- {now}', '']
    out.append(f"{'ARM':<11}{'CONSTRUCT':<16}{'GPU':>4}{'STEPS':>14}{'/TARGET':>14}"
               f"{'%':>7}{'ns/day':>9}{'ETA':>14}  STATUS")
    out.append('-' * 104)
    tot_rate = 0.0
    counts = {}
    for r in rows:
        pct = 100.0 * r['step'] / r['target'] if r['target'] else 0.0
        counts[r['status']] = counts.get(r['status'], 0) + 1
        if r['status'] == 'running':
            tot_rate += r['rate']
        rate_s = f"{r['rate']:.0f}" if r['rate'] else '-'
        out.append(f"{r['arm']:<11}{r['construct']:<16}{str(r['gpu']):>4}"
                   f"{r['step']:>14,}{r['target']:>14,}{pct:>6.1f}%{rate_s:>9}"
                   f"{r['eta']:>14}  {r['status']}")
    out.append('-' * 104)
    summary = '  '.join(f'{k}={v}' for k, v in sorted(counts.items()))
    out.append(f'{summary}    aggregate {tot_rate:.0f} ns/day')

    done_bead_steps = sum(r['step'] * r['beads'] for r in rows)
    todo_bead_steps = sum((r['target'] - r['step']) * r['beads'] for r in rows)
    out.append(f'bead-steps done {done_bead_steps:.3e}   remaining {todo_bead_steps:.3e}')

    if show_gpu:
        out += ['', 'GPUs:']
        for idx, name, util, used, total in gpu_table():
            bar = '#' * (util // 10) + '.' * (10 - util // 10)
            out.append(f'  [{idx}] {name:<14} {bar} {util:>3}%  {used:>6}/{total} MiB')
        out.append('  (memory beyond our ~7 GB/run is another user -- lyra is shared)')
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--arm', choices=['homotypic', 'rna', 'both'], default='both')
    ap.add_argument('--watch', action='store_true')
    ap.add_argument('-n', '--interval', type=int, default=60)
    ap.add_argument('--no-gpu', action='store_true')
    args = ap.parse_args()
    arms = ARMS if args.arm == 'both' else (args.arm,)

    if not args.watch:
        print(render(collect(arms), show_gpu=not args.no_gpu))
        return
    try:
        while True:
            txt = render(collect(arms), show_gpu=not args.no_gpu)
            print('\033[2J\033[H' + txt, flush=True)
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print()


if __name__ == '__main__':
    main()
