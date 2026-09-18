#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full dashboard-backfill recipe for the 9 new contiguous constructs launched
this session (all 45 replicates confirmed complete at 1us):
kh1_art_full, kh1_wwe_full, md2_art_full, md2_wwe_full, md3_art_full,
md3_wwe_full, core_wwe_full_go, mka_wwe_full, fl_wwe_full_go.

Same recipe as run_dashboard_backfill_extended.py. IMPORTANT: the per-sim
md_distances step below (figure_md_distances.py --set <sim>) OVERWRITES the
shared data/md_distances_com.npz cache with ONLY that one set each time --
necessary to populate each sim's own by_sim/<sim>/04_md_distances/ figures,
but it destroys the comprehensive 13-way comparison cache already built
(md_full+mka_full+core_full_go+9 new sets+fl_go_5rep1us). So the LAST step
here restores that full comparison cache in one final combined call.
"""
import subprocess
import sys
import time

HERE_ENV = ['conda', 'run', '-n', 'calvados']
SETS = ['kh1_art_full', 'kh1_wwe_full', 'md2_art_full', 'md2_wwe_full',
        'md3_art_full', 'md3_wwe_full', 'core_wwe_full_go', 'mka_wwe_full',
        'fl_wwe_full_go']

STEPS = [
    ('run_analysis', ['bash', 'run_analysis.sh', '--set', '{sim}']),
    ('md_distances', ['python', 'figure_md_distances.py', '--set', '{sim}']),
    ('active_sites', ['python', 'analyze_active_sites.py', '--set', '{sim}']),
    ('accessibility', ['python', 'analyze_accessibility.py', '--set', '{sim}']),
    ('lys_persistence', ['python', 'figure_lys_exposed_persistence.py', '--sets', '{sim}']),
    ('dashboard', ['python', 'build_sim_dashboard.py', '--sim', '{sim}']),
]

# Final step: restore the comprehensive cross-set comparison cache (see
# module docstring) after all per-sim md_distances calls above have
# clobbered it down to one set at a time.
RESTORE_CACHE_CMD = [
    'python', 'figure_md_distances.py',
    '--set', 'md_full', 'mka_full', 'core_full_go',
    'kh1_art_full', 'kh1_wwe_full', 'md2_art_full', 'md2_wwe_full',
    'md3_art_full', 'md3_wwe_full', 'core_wwe_full_go', 'mka_wwe_full', 'fl_wwe_full_go',
    '--sim-folder-as', 'fl_go_5rep1us:/home/sbali/CALVADOS/examples/PARP14_MDP/fl_go',
    '--max-reps', 'fl_go_5rep1us:5', '--max-ns', 'fl_go_5rep1us:1000',
    # md_full has 25 replicates total, but only seed-{1..5}_sample-0 were
    # extended to the full ~1000 ns target (sample-{1,2,3,4} are still at the
    # original ~100 ns) -- restrict to those 5 so it's length-matched to
    # every other _full/_full_go set (all genuinely 5 reps x 1000 ns).
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
            results.append({'sim': sim, 'step': step_name, 'ok': ok, 'seconds': round(dt)})
        log(f"=== {sim} {'DONE' if sim_ok else 'DONE WITH FAILURES'} ===\n")

    log("=== Restoring comprehensive 13-way comparison cache ===")
    t0 = time.time()
    rc, out, err = run_cmd(RESTORE_CACHE_CMD)
    dt = time.time() - t0
    ok = rc == 0
    log(f"  [restore_cache] {'OK' if ok else 'FAILED'} in {dt:.0f}s")
    if not ok:
        log(f"    stderr tail: {err.strip()[-600:]}")
    results.append({'sim': 'ALL', 'step': 'restore_cache', 'ok': ok, 'seconds': round(dt)})

    log("=== BACKFILL SUMMARY ===")
    n_ok = sum(1 for r in results if r['ok'])
    log(f"{n_ok}/{len(results)} steps completed without error")
    for r in results:
        if not r['ok']:
            log(f"  [FAILED] {r['sim']} / {r['step']}")

    return 0 if n_ok == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())
