#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full dashboard-backfill recipe for core_full_go (missing its own by_sim
dashboard -- only ever analyzed as part of multi-set comparison runs) and
md1_md2/md2_md3 (brand new 2-domain isolation-test constructs, never
analyzed). Restores the comprehensive 13-way comparison cache at the end."""
import subprocess, sys, time

HERE_ENV = ['conda', 'run', '-n', 'calvados']
SETS = ['core_full_go', 'md1_md2', 'md2_md3']

STEPS = [
    ('run_analysis', ['bash', 'run_analysis.sh', '--set', '{sim}']),
    ('md_distances', ['python', 'figure_md_distances.py', '--set', '{sim}']),
    ('active_sites', ['python', 'analyze_active_sites.py', '--set', '{sim}']),
    ('accessibility', ['python', 'analyze_accessibility.py', '--set', '{sim}']),
    ('lys_persistence', ['python', 'figure_lys_exposed_persistence.py', '--sets', '{sim}']),
    ('dashboard', ['python', 'build_sim_dashboard.py', '--sim', '{sim}']),
]

RESTORE_CACHE_CMD = [
    'python', 'figure_md_distances.py',
    '--set', 'md_full', 'mka_full', 'core_full_go',
    'kh1_art_full', 'kh1_wwe_full', 'md2_art_full', 'md2_wwe_full',
    'md3_art_full', 'md3_wwe_full', 'core_wwe_full_go', 'mka_wwe_full', 'fl_wwe_full_go', 'md1_md2', 'md2_md3',
    '--sim-folder-as', 'fl_go_5rep1us:/home/sbali/CALVADOS/examples/PARP14_MDP/fl_go',
    '--max-reps', 'fl_go_5rep1us:5', '--max-ns', 'fl_go_5rep1us:1000',
    '--only-sample', 'md_full:0',
    '--workers', '24',
]

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def run_cmd(cmd):
    env_cmd = HERE_ENV + ['env', 'WORKERS=16'] + cmd if cmd[0] == 'bash' else HERE_ENV + cmd
    result = subprocess.run(env_cmd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr

def run_step(cmd_template, sim):
    cmd = [c.format(sim=sim) for c in cmd_template]
    return run_cmd(cmd)

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
            results.append({'sim': sim, 'step': step_name, 'ok': ok})
        log(f"=== {sim} {'DONE' if sim_ok else 'DONE WITH FAILURES'} ===\n")

    log("=== Restoring comprehensive 13-way comparison cache ===")
    t0 = time.time()
    rc, out, err = run_cmd(RESTORE_CACHE_CMD)
    ok = rc == 0
    log(f"  [restore_cache] {'OK' if ok else 'FAILED'} in {time.time()-t0:.0f}s")
    if not ok:
        log(f"    stderr tail: {err.strip()[-600:]}")
    results.append({'sim': 'ALL', 'step': 'restore_cache', 'ok': ok})

    log("=== BACKFILL SUMMARY ===")
    n_ok = sum(1 for r in results if r['ok'])
    log(f"{n_ok}/{len(results)} steps completed without error")
    for r in results:
        if not r['ok']:
            log(f"  [FAILED] {r['sim']} / {r['step']}")
    return 0 if n_ok == len(results) else 1

if __name__ == '__main__':
    sys.exit(main())
