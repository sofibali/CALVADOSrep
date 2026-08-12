#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch driver for make_surface_gallery.py across every requested
representative-frame state set.

Runs states SEQUENTIALLY (each state is itself CPU-heavy: APBS + several
ray-traced surfaces), continuing to the next state on failure rather than
retrying -- if a state fails, it's noted in the summary and skipped.

Usage:
    conda run -n pymol-render python run_surface_gallery_batch.py
"""
import subprocess
import sys
import time
from pathlib import Path

CWD = Path(__file__).resolve().parent
ROOT = CWD.parent
sys.path.insert(0, str(ROOT))
from _fig_layout import get_fig_dir  # noqa: E402

# (run_tag, set_key, state_pdb_dir)
RUNS = [
    ('fl_optimized_ca25_tica', 'fl_optimized',
     ROOT / 'representative_frames/2026-06-25/fl_optimized_ca25_tica_allatom'),
    ('fl_optimized_iface10_tica', 'fl_optimized',
     ROOT / 'representative_frames/2026-06-25/fl_optimized_iface10_tica_allatom'),
    ('md_full_pose', 'md',
     ROOT / 'states/md_full_pose_allatom'),
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    fig_root = get_fig_dir('09_surface_gallery')
    results = []

    for run_tag, set_key, pdb_dir in RUNS:
        state_pdbs = sorted(pdb_dir.glob('state_*_allatom.pdb'),
                             key=lambda p: int(p.stem.split('_')[1]))
        if not state_pdbs:
            log(f"WARNING: no state PDBs found in {pdb_dir}, skipping {run_tag}")
            continue
        log(f"=== {run_tag} ({set_key}): {len(state_pdbs)} states ===")

        for pdb_path in state_pdbs:
            label = pdb_path.stem.replace('_allatom', '')
            out_dir = fig_root / run_tag / label
            pse_out = out_dir / f'{label}.pse'
            t0 = time.time()
            log(f"--- {run_tag}/{label} -> {out_dir} ---")

            proc = subprocess.run(
                ['conda', 'run', '-n', 'pymol-render', 'python',
                 str(CWD / 'make_surface_gallery.py'),
                 '--set-key', set_key, '--pdb', str(pdb_path),
                 '--label', label, '--out-dir', str(out_dir),
                 '--pse-out', str(pse_out)],
                capture_output=True, text=True)

            dt = time.time() - t0
            ok = proc.returncode == 0
            log(f"    {'OK' if ok else 'FAILED'} in {dt:.0f}s "
                f"(returncode={proc.returncode})")
            if not ok:
                log(f"    stderr tail: {proc.stderr.strip()[-800:]}")
            # Always surface any FAILED lines from the child's own per-step
            # try/except logging, even on overall success -- partial
            # failures (e.g. one pocket) don't fail the returncode.
            for line in proc.stdout.splitlines():
                if 'FAILED' in line or 'SKIPPED' in line:
                    log(f"    {line.strip()}")

            results.append({
                'run_tag': run_tag, 'label': label, 'ok': ok,
                'seconds': round(dt), 'out_dir': str(out_dir),
            })

    log("=== BATCH SUMMARY ===")
    n_ok = sum(1 for r in results if r['ok'])
    log(f"{n_ok}/{len(results)} states completed without a fatal error")
    for r in results:
        status = 'OK' if r['ok'] else 'FAILED'
        log(f"  [{status}] {r['run_tag']}/{r['label']} ({r['seconds']}s) -> {r['out_dir']}")

    return 0 if n_ok == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())
