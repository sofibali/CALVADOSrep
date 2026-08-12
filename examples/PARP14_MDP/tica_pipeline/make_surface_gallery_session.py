#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build one combined PyMOL session per run_tag with every state's
domain-colored cartoon loaded as a separate, aligned object, plus named
scenes to flip between states (mirrors tica_pipeline/make_cluster_pse.py's
conventions). Each state's own full surface-property session (from
make_surface_gallery.py's --pse-out) still exists separately as
figures/09_surface_gallery/<date>/<run_tag>/state_N/state_N.pse for anyone
who wants to open one state's electrostatic/hydrophobic/face/pocket
surfaces interactively -- this combined session is for comparing
conformations side by side, not for re-hosting every heavy surface object
from every state in one file.

Usage:
    conda run -n pymol-render python make_surface_gallery_session.py --auto
    conda run -n pymol-render python make_surface_gallery_session.py \\
        --run-dir figures/09_surface_gallery/<date>/<run_tag> --set-key fl_optimized
"""
import argparse
import sys
from pathlib import Path

import pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-cq'])

CWD = Path(__file__).resolve().parent
ROOT = CWD.parent
sys.path.insert(0, str(CWD))
sys.path.insert(0, str(ROOT))
import surface_gallery_lib as sgl   # noqa: E402
import sim_registry as sr           # noqa: E402


def find_state_dirs(run_dir):
    def key(p):
        try:
            return int(p.name.split('_')[1])
        except (IndexError, ValueError):
            return p.name
    return sorted((d for d in run_dir.iterdir() if d.is_dir() and d.name.startswith('state_')), key=key)


def build_session(run_dir, set_key, out_pse):
    run_dir = Path(run_dir)
    state_dirs = find_state_dirs(run_dir)
    if not state_dirs:
        print(f"  no state_* dirs under {run_dir}")
        return False

    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)

    domains = sgl.get_domain_ranges_for_set(set_key)
    all_cat = []
    for site in sr.SITE_NAMES:
        cat, _ = sgl.get_pocket_for_set(set_key, site)
        all_cat += cat

    state_objs = []
    for sd in state_dirs:
        # The centered PDB written by prep_apbs.py has the same residue
        # numbering as the original state PDB, just recentered -- use it
        # if present (matches the surface renders' frame), else fall back.
        apbs_pdb = sd / '_apbs_work' / 'centered.pdb'
        # Fall back to the original representative-frame PDB via the label.
        obj = sd.name
        pdb_path = apbs_pdb if apbs_pdb.exists() else None
        if pdb_path is None:
            print(f"  WARNING: no centered PDB for {obj}, skipping")
            continue
        cmd.load(str(pdb_path), obj)
        cmd.hide('everything', obj)
        cmd.show('cartoon', obj)
        cmd.color('gray80', obj)
        for dname, (color_name, rgb, (cs, ce)) in domains.items():
            cmd.set_color(color_name, list(rgb))
            cmd.color(color_name, f'{obj} and resi {cs}-{ce}')
        if all_cat:
            sel = f'{obj} and resi {"+".join(map(str, all_cat))} and name CA'
            cmd.show('spheres', sel)
            cmd.color('yellow', sel)
            cmd.set('sphere_scale', 0.7, sel)
        state_objs.append(obj)

    if not state_objs:
        return False

    ref_obj = state_objs[0]
    align_sel = 'name CA'
    for obj in state_objs[1:]:
        try:
            cmd.align(f'{obj} and {align_sel}', f'{ref_obj} and {align_sel}')
        except Exception as e:
            print(f"  align {obj}->{ref_obj} failed (non-fatal): {e}")

    run_tag = run_dir.name
    cmd.group(f'{run_tag}_states', ' '.join(state_objs))
    cmd.orient(ref_obj)

    for i, obj in enumerate(state_objs, start=1):
        cmd.disable('all')
        cmd.enable(obj)
        cmd.scene(f'S{i}_{obj}', 'store')

    cmd.disable('all')
    for obj in state_objs:
        cmd.enable(obj)
        cmd.set('cartoon_transparency', 0.4, obj)
    cmd.scene('All_states_overlay', 'store')
    for obj in state_objs:
        cmd.set('cartoon_transparency', 0.0, obj)

    cmd.scene(f'S1_{state_objs[0]}', 'recall')
    out_pse = Path(out_pse)
    out_pse.parent.mkdir(parents=True, exist_ok=True)
    cmd.save(str(out_pse))
    print(f"  saved {out_pse}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', action='append', default=[])
    ap.add_argument('--set-key', action='append', default=[],
                     help='registry set key matching each --run-dir, in order')
    ap.add_argument('--auto', action='store_true')
    args = ap.parse_args()

    jobs = list(zip(args.run_dir, args.set_key))
    if args.auto:
        gallery_root = ROOT / 'figures' / '09_surface_gallery'
        latest = gallery_root / 'latest'
        # run_tag -> set_key, matching run_surface_gallery_batch.py's RUNS
        auto_map = {
            'fl_optimized_ca25_tica': 'fl_optimized',
            'fl_optimized_iface10_tica': 'fl_optimized',
            'md_full_pose': 'md',
        }
        if latest.is_dir():
            for run_dir in latest.iterdir():
                if run_dir.is_dir() and run_dir.name in auto_map:
                    jobs.append((str(run_dir), auto_map[run_dir.name]))

    if not jobs:
        print("Nothing to do: pass --run-dir/--set-key pairs or --auto")
        return 1

    out_root = ROOT / 'figures' / '09_surface_gallery' / 'latest' / 'sessions'
    ok = True
    for run_dir, set_key in jobs:
        run_dir = Path(run_dir)
        print(f"Building session for {run_dir.name} ({set_key})...")
        cmd.reinitialize()
        ok = build_session(run_dir, set_key, out_root / f'{run_dir.name}.pse') and ok

    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
