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


def read_log(path):
    """(last_step, last_speed_ns_day, last_elapsed_s, n_rows) from a reporter log."""
    if not os.path.isfile(path):
        return 0, 0.0, 0.0, 0
    step = elapsed = 0.0
    speed = 0.0
    rows = 0
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                s, sp, el = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            step, speed, elapsed, rows = max(step, s), sp, el, rows + 1
    return int(step), speed, elapsed, rows


def recent_rate(path, window=5):
    """
    ns/day over the last `window` log intervals.

    The Speed column is a running average since the start of the leg, so it
    lags badly when a GPU becomes contended mid-run. This differences the
    raw Step/Elapsed columns instead, which tracks the current rate.
    """
    rows = []
    if not os.path.isfile(path):
        return 0.0
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                rows.append((float(parts[0]), float(parts[2])))
            except ValueError:
                continue
    if len(rows) < 2:
        return 0.0
    sub = rows[-(window + 1):]
    dstep = sub[-1][0] - sub[0][0]
    dt = sub[-1][1] - sub[0][1]
    if dt <= 0 or dstep <= 0:
        return 0.0
    return (dstep * DT_PS * 1e-3) / dt * 86400.0  # ns/day


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
                mtime = os.path.getmtime(flog) if os.path.isfile(flog) else 0
                status = 'running' if time.time() - mtime < 900 else 'stalled'
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
