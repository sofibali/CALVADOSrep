#!/usr/bin/env python3
"""
Headless ray-traced rendering of PARP14 accessibility scenes.

Loads PDB frames directly (not from .pse session), applies styling,
and renders in two orientations:
  1. "MD1 view": aligned on MD1 domain, MD1 active site facing up
  2. "ART view": aligned on ART domain, ART active site facing down

Usage:
    conda run -n pymol-render python render_scenes.py
    conda run -n pymol-render python render_scenes.py --width 3600 --height 2700
    conda run -n pymol-render python render_scenes.py --objects fl_MD1_buried core_MD1_exposed
"""

import os
import sys
import numpy as np
from argparse import ArgumentParser

CWD = os.path.dirname(os.path.abspath(__file__))
FRAMES_DIR = os.path.join(CWD, 'pymol_frames')
RENDER_DIR = os.path.join(CWD, 'renders')

# ============================================================
# Domain definitions
# ============================================================

DOMAIN_UNITS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789),
    'md1l1': (790, 1004), 'md2': (1004, 1193), 'md3': (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe': (1534, 1602), 'art': (1603, 1801),
}

CONSTRUCT_UNITS = {
    'fl':    list(DOMAIN_UNITS.keys()),
    'md':    ['md1l1', 'md2', 'md3'],
    'core':  ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'mka':   ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'norrm': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'noart': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'],
}

CATALYTIC_FL = {
    'MD1': [831, 923, 962],
    'MD2': [1035, 1046, 1134, 1171],
    'MD3': [1248, 1259, 1330, 1371],
    'ART': [1684, 1705, 1706, 1722],
}

FL_DOMAIN_RANGES = {
    'MD1': (791, 978), 'MD2': (1003, 1190), 'MD3': (1216, 1387), 'ART': (1603, 1801),
}


def compute_fl_blocks(unit_names):
    ranges = sorted([DOMAIN_UNITS[u] for u in unit_names])
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def build_fl_to_construct_map(unit_names):
    fl_blocks = compute_fl_blocks(unit_names)
    segments = []
    pos = 1
    for fl_s, fl_e in fl_blocks:
        segments.append((fl_s, fl_e, pos - fl_s))
        pos += fl_e - fl_s + 1
    def map_resid(r):
        for fl_s, fl_e, off in segments:
            if fl_s <= r <= fl_e:
                return r + off
        return None
    return map_resid


def get_domain_ranges(sk):
    if sk == 'fl':
        return dict(FL_DOMAIN_RANGES)
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[sk])
    result = {}
    for dname, (fl_s, fl_e) in FL_DOMAIN_RANGES.items():
        c_s, c_e = fl_to_c(fl_s), fl_to_c(fl_e)
        if c_s is not None and c_e is not None:
            result[dname] = (c_s, c_e)
    return result


def get_catalytic(sk):
    if sk == 'fl':
        return dict(CATALYTIC_FL)
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[sk])
    result = {}
    for sname, resids in CATALYTIC_FL.items():
        mapped = [fl_to_c(r) for r in resids if fl_to_c(r) is not None]
        if mapped:
            result[sname] = mapped
    return result


def read_picks():
    picks = []
    with open(os.path.join(FRAMES_DIR, 'picks.txt')) as f:
        for line in f:
            parts = line.strip().split('\t')
            picks.append((parts[0], parts[3], parts[4]))  # sk, label, desc
    return picks


DOMAIN_COLORS = {
    'MD1': [0.40, 0.10, 0.50],
    'MD2': [0.58, 0.30, 0.68],
    'MD3': [0.72, 0.50, 0.80],
    'ART': [1.00, 0.00, 0.00],
}
CAT_GRAY = [0.55, 0.55, 0.55]


def style_object(cmd, label, sk):
    """Apply coloring and surface to a loaded object."""
    cmd.show('surface', label)
    cmd.hide('cartoon', label)
    cmd.hide('lines', label)
    cmd.color('white', label)

    dom_ranges = get_domain_ranges(sk)
    catalytic = get_catalytic(sk)

    for dname, rgb in DOMAIN_COLORS.items():
        if dname in dom_ranges:
            cname = f'{dname}_color'
            cmd.set_color(cname, rgb)
            s, e = dom_ranges[dname]
            cmd.color(cname, f'{label} and resi {s}-{e}')

    cmd.set_color('cat_gray', CAT_GRAY)
    for sname, resids in catalytic.items():
        resi_str = '+'.join(str(r) for r in resids)
        cmd.color('cat_gray', f'{label} and resi {resi_str}')


def load_and_style(cmd, label, sk):
    """Load a PDB and apply full styling."""
    pdb_path = os.path.join(FRAMES_DIR, f'{label}.pdb')
    if not os.path.isfile(pdb_path):
        return False
    cmd.load(pdb_path, label)
    style_object(cmd, label, sk)
    return True


def render_one(cmd, label, out_path, width, height):
    """Enable one object, ray trace, save."""
    cmd.disable('all')
    cmd.enable(label)
    cmd.ray(width, height)
    cmd.png(out_path, dpi=300)


def main():
    parser = ArgumentParser()
    parser.add_argument('--outdir', default=RENDER_DIR)
    parser.add_argument('--width', type=int, default=2400)
    parser.add_argument('--height', type=int, default=1800)
    parser.add_argument('--objects', nargs='*', default=None)
    args = parser.parse_args()

    import pymol
    from pymol import cmd

    pymol.finish_launching(['pymol', '-cq'])
    os.makedirs(args.outdir, exist_ok=True)

    # Global settings
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadow', 0)
    cmd.set('ray_trace_mode', 0)
    cmd.set('antialias', 2)
    cmd.set('surface_quality', 1)
    cmd.bg_color('white')

    picks = read_picks()
    sk_for = {label: sk for sk, label, _ in picks}

    # Filter
    if args.objects:
        picks = [(sk, label, desc) for sk, label, desc in picks if label in args.objects]

    # Load all objects from PDB files (not session)
    loaded = []
    for sk, label, desc in picks:
        if load_and_style(cmd, label, sk):
            loaded.append((sk, label))
            print(f"  Loaded: {label}")

    if not loaded:
        print("No objects loaded!")
        sys.exit(1)

    # Find FL reference for alignment
    ref_label = None
    ref_sk = None
    for sk, label in loaded:
        if sk == 'fl':
            ref_label = label
            ref_sk = sk
            break
    if ref_label is None:
        ref_label = loaded[0][1]
        ref_sk = loaded[0][0]

    # ============================================================
    # Orientation 1: Align on MD1, site facing up
    # ============================================================

    print(f"\n--- MD1 view ({len(loaded)} objects) ---")
    md1_dir = os.path.join(args.outdir, 'MD1_view')
    os.makedirs(md1_dir, exist_ok=True)

    ref_dom = get_domain_ranges(ref_sk)
    if 'MD1' in ref_dom:
        rs, re = ref_dom['MD1']
        ref_md1_sel = f'{ref_label} and resi {rs}-{re}'
    else:
        ref_md1_sel = ref_label

    # Align all to reference MD1
    for sk, label in loaded:
        if label == ref_label:
            continue
        tgt_dom = get_domain_ranges(sk)
        if 'MD1' in tgt_dom:
            ts, te = tgt_dom['MD1']
            try:
                cmd.align(f'{label} and resi {ts}-{te}', ref_md1_sel)
            except:
                cmd.align(label, ref_label)
        else:
            cmd.align(label, ref_label)

    # Set view: zoom to fit the largest object with generous buffer
    cmd.disable('all')
    # Enable all to find the widest extent, then zoom
    for _, label in loaded:
        cmd.enable(label)
    cmd.zoom('all', buffer=15)
    # Disable again before per-object renders
    cmd.disable('all')

    # Save the view matrix so all renders use same orientation
    view = cmd.get_view()

    for sk, label in loaded:
        cmd.set_view(view)
        out_path = os.path.join(md1_dir, f'{label}.png')
        render_one(cmd, label, out_path, args.width, args.height)
        print(f"  {label} -> MD1_view/")

    # ============================================================
    # Orientation 2: Align on ART, site facing down
    # ============================================================

    art_loaded = [(sk, label) for sk, label in loaded if 'ART' in get_domain_ranges(sk)]

    if art_loaded:
        print(f"\n--- ART view ({len(art_loaded)} objects) ---")
        art_dir = os.path.join(args.outdir, 'ART_view')
        os.makedirs(art_dir, exist_ok=True)

        # Find ART reference
        art_ref = None
        art_ref_sk = None
        for sk, label in art_loaded:
            if sk == 'fl':
                art_ref = label
                art_ref_sk = sk
                break
        if art_ref is None:
            art_ref = art_loaded[0][1]
            art_ref_sk = art_loaded[0][0]

        art_ref_dom = get_domain_ranges(art_ref_sk)
        ars, are = art_ref_dom['ART']
        art_ref_sel = f'{art_ref} and resi {ars}-{are}'

        # Re-align all ART objects on ART domain
        for sk, label in art_loaded:
            if label == art_ref:
                continue
            tgt_dom = get_domain_ranges(sk)
            if 'ART' in tgt_dom:
                ts, te = tgt_dom['ART']
                try:
                    cmd.align(f'{label} and resi {ts}-{te}', art_ref_sel)
                except:
                    cmd.align(label, art_ref)

        # Orient on ART, then flip 180 so ART faces down
        cmd.disable('all')
        for _, label in art_loaded:
            cmd.enable(label)
        cmd.zoom('all', buffer=15)
        cmd.rotate('x', 180)
        cmd.disable('all')

        art_view = cmd.get_view()

        for sk, label in art_loaded:
            cmd.set_view(art_view)
            out_path = os.path.join(art_dir, f'{label}.png')
            render_one(cmd, label, out_path, args.width, args.height)
            print(f"  {label} -> ART_view/")

    print(f"\nAll renders saved to: {args.outdir}/")
    cmd.quit()


if __name__ == '__main__':
    main()
