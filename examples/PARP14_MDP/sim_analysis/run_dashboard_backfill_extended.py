#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Re-run the full dashboard-backfill recipe for the 7 sets whose trajectories
were just extended by launch_staged_extensions.py (fl, md, core, mka, norrm,
noart, md3art -- 175/175 replicates completed with 0 failures). Their
data/*.npz + figures were computed from the pre-extension trajectories, so
everything needs a refresh now that the extra steps have landed.

Same recipe/steps as run_dashboard_backfill.py, just scoped to these 7 sets.
"""
import subprocess
import sys
import time

HERE_ENV = ['conda', 'run', '-n', 'calvados']
SETS = ['fl', 'md', 'core', 'mka', 'norrm', 'noart', 'md3art']

STEPS = [
    ('run_analysis', ['bash', 'run_analysis.sh', '--set', '{sim}']),
    ('md_distances', ['python', 'figure_md_distances.py', '--set', '{sim}']),
    ('active_sites', ['python', 'analyze_active_sites.py', '--set', '{sim}']),
    ('accessibility', ['python', 'analyze_accessibility.py', '--set', '{sim}']),
    ('lys_persistence', ['python', 'figure_lys_exposed_persistence.py', '--sets', '{sim}']),
    ('dashboard', ['python', 'build_sim_dashboard.py', '--sim', '{sim}']),
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_step(cmd_template, sim):
    cmd = [c.format(sim=sim) for c in cmd_template]
    env_cmd = HERE_ENV + ['env', 'WORKERS=16'] + cmd if cmd[0] == 'bash' else HERE_ENV + cmd
    result = subprocess.run(env_cmd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def main():
    results = []
    for sim in SETS:
        log(f"=== {sim} ===")
        sim_ok = True
        for step_name, cmd_template in STEPS:
            t0 = time.time()
            rc, out, err = run_step(cmd_template, sim)
            dt = time.time() - t0
            ok = rc == 0
            sim_ok = sim_ok and ok
            log(f"  [{step_name}] {'OK' if ok else 'FAILED'} in {dt:.0f}s")
            if not ok:
                log(f"    stderr tail: {err.strip()[-600:]}")
            results.append({'sim': sim, 'step': step_name, 'ok': ok, 'seconds': round(dt)})
        log(f"=== {sim} {'DONE' if sim_ok else 'DONE WITH FAILURES'} ===\n")

    log("=== BACKFILL SUMMARY ===")
    n_ok = sum(1 for r in results if r['ok'])
    log(f"{n_ok}/{len(results)} steps completed without error")
    for r in results:
        if not r['ok']:
            log(f"  [FAILED] {r['sim']} / {r['step']}")

    return 0 if n_ok == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())
