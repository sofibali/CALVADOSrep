#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process an existing per-state surface-gallery .pse (from
make_surface_gallery.py --pse-out): add named PyMOL selections for every
domain's face-classification group (active/back for MD1L1/MD2/MD3/ART,
faceA/faceB for the rest) on the mol_face object, plus a 6-angle rotation
scene set (front/back/left/right/top/bottom) so the face-classification
surface can be spun through preset views instead of just the one saved
orientation. Re-saves in place (or to --pse-out if given).

Selections are recomputed from mol_face's own coordinates with the same
classify_active_face()/classify_generic_face() logic make_surface_gallery.py
used to color it -- deterministic given the same seed, so this reproduces
the exact same groups, just also exposed as named selections.

Usage:
    conda run -n pymol-render python enhance_face_classification_pse.py \\
        --pse figures/09_surface_gallery/.../state_1/state_1.pse --set-key fl_optimized
"""
import argparse
import sys
from pathlib import Path

import pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-cq'])

CWD = Path(__file__).resolve().parent
sys.path.insert(0, str(CWD))
sys.path.insert(0, str(CWD.parent))
import surface_gallery_lib as sgl                                              # noqa: E402
import sim_registry as sr                                                      # noqa: E402
from make_surface_gallery import classify_active_face, classify_generic_face, DOMAIN_TO_SITE  # noqa: E402


def add_selections_and_scenes(obj, set_key):
    domains = sgl.get_domain_ranges_for_set(set_key)
    made = []
    for i, (dname, (color_name, base_rgb, (cs, ce))) in enumerate(domains.items()):
        site = DOMAIN_TO_SITE.get(dname)
        faces = {}
        if site:
            cat, _ = sgl.get_pocket_for_set(set_key, site)
            if cat:
                faces = classify_active_face(obj, (cs, ce), cat)
        if not faces:
            faces = classify_generic_face(obj, (cs, ce), seed=42 + i)
        if not faces:
            continue

        groups = {}
        for r, f in faces.items():
            groups.setdefault(f, []).append(r)

        safe_dname = dname.replace('-', '_')
        for fname, resids in groups.items():
            sel_name = f'{safe_dname}_{fname}'
            resi_str = "+".join(str(r) for r in sorted(resids))
            cmd.select(sel_name, f'{obj} and resi {resi_str}')
            made.append(sel_name)
    cmd.deselect()
    return made


def add_rotation_scenes(obj, prefix):
    """front/back/left/right/top/bottom scenes for `obj`, relative to its
    current cmd.orient()-ed pose (front)."""
    for name in cmd.get_names('objects'):
        cmd.disable(name)
    cmd.enable(obj)
    cmd.orient(obj)

    views = [
        ('front', None),
        ('back', ('y', 180)),
        ('left', ('y', 90)),
        ('right', ('y', -90)),
        ('top', ('x', 90)),
        ('bottom', ('x', -90)),
    ]
    base_view = cmd.get_view()
    for vname, turn in views:
        cmd.set_view(base_view)
        if turn:
            axis, angle = turn
            cmd.turn(axis, angle)
        cmd.scene(f'{prefix}_{vname}', 'store')
    cmd.set_view(base_view)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pse', required=True)
    ap.add_argument('--set-key', required=True)
    ap.add_argument('--pse-out', default=None, help='default: overwrite --pse in place')
    ap.add_argument('--obj', default='mol_face')
    args = ap.parse_args()

    cmd.load(args.pse)
    if args.obj not in cmd.get_names('objects'):
        print(f"ERROR: object '{args.obj}' not found in {args.pse}. "
              f"Objects present: {cmd.get_names('objects')}")
        return 1

    made = add_selections_and_scenes(args.obj, args.set_key)
    print(f"Created {len(made)} named selections: {', '.join(made)}")

    add_rotation_scenes(args.obj, 'face')
    print("Added rotation scenes: face_front, face_back, face_left, face_right, face_top, face_bottom")

    out_path = args.pse_out or args.pse
    cmd.scene('face_front', 'recall')
    cmd.save(out_path)
    print(f"Saved: {out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
