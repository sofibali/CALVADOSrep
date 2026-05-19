#!/usr/bin/env python3
"""
Extract representative frames and build a PyMOL session showing
active site accessibility across PARP14 constructs.

Step 1: Extract PDB frames from DCDs (needs calvados env)
Step 2: Build PyMOL session (needs pymol-render env)

Per construct: picks the replicate with the most buried and most
exposed conformation for each active site (MD1, MD2, MD3, ART).

Usage:
    conda run -n calvados python make_pymol_session.py --extract
    conda run -n pymol-render python make_pymol_session.py --build-session
"""

import os
import sys
import numpy as np
from argparse import ArgumentParser

CWD = os.path.dirname(os.path.abspath(__file__))
FRAMES_DIR = os.path.join(CWD, 'pymol_frames')

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
    'fl':           list(DOMAIN_UNITS.keys()),
    'md':           ['md1l1', 'md2', 'md3'],
    'core':         ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'mka':          ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'norrm':        ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'noart':        ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    'fl_optimized': list(DOMAIN_UNITS.keys()),
}

SETS_INFO = {
    'fl': 'parp14', 'md': 'parp14_macrodomains', 'core': 'parp14_core',
    'mka': 'parp14_mka', 'norrm': 'parp14_norrm', 'noart': 'parp14_noart',
    'fl_optimized': 'parp14',
}

CONSTRUCT_SITES = {
    'fl':           ['MD1', 'MD2', 'MD3', 'ART'],
    'md':           ['MD1', 'MD2', 'MD3'],
    'core':         ['MD1', 'MD2', 'MD3', 'ART'],
    'mka':          ['MD1', 'MD2', 'MD3', 'ART'],
    'norrm':        ['MD1', 'MD2', 'MD3', 'ART'],
    'noart':        ['MD1', 'MD2', 'MD3'],
    'fl_optimized': ['MD1', 'MD2', 'MD3', 'ART'],
}

# Catalytic residues only (FL numbering)
CATALYTIC_FL = {
    'MD1': [831, 923, 962],
    'MD2': [1035, 1046, 1134, 1171],
    'MD3': [1248, 1259, 1330, 1371],
    'ART': [1684, 1705, 1706, 1722],
}

# FL domain ranges for coloring (full domain, not trimmed)
FL_DOMAIN_RANGES = {
    'MD1': (791, 978), 'MD2': (1003, 1190), 'MD3': (1216, 1387), 'ART': (1603, 1801),
}

# ============================================================
# Residue mapping
# ============================================================

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


def get_construct_domain_ranges(sk):
    if sk == 'fl':
        return dict(FL_DOMAIN_RANGES)
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[sk])
    result = {}
    for dname, (fl_s, fl_e) in FL_DOMAIN_RANGES.items():
        c_s, c_e = fl_to_c(fl_s), fl_to_c(fl_e)
        if c_s is not None and c_e is not None:
            result[dname] = (c_s, c_e)
    return result


def get_construct_catalytic(sk):
    if sk == 'fl':
        return dict(CATALYTIC_FL)
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[sk])
    result = {}
    for sname, resids in CATALYTIC_FL.items():
        mapped = [fl_to_c(r) for r in resids]
        mapped = [r for r in mapped if r is not None]
        if mapped:
            result[sname] = mapped
    return result


# ============================================================
# Frame selection: per-site extremes from accessibility data
# ============================================================

def pick_frames():
    """Pick replicates showing extreme accessibility for each site."""
    data = np.load(os.path.join(CWD, 'data', 'accessibility_stats.npz'))
    picks = []  # (sk, seed, sample, label, desc)
    seen = set()

    for sk in ['fl', 'core', 'mka', 'norrm', 'noart', 'md', 'fl_optimized']:
        if sk not in CONSTRUCT_SITES:
            continue
        sites = CONSTRUCT_SITES[sk]
        for site in sites:
            key = f'{sk}_{site}_saa'
            if key not in data:
                continue
            vals = data[key]
            for extremum, tag in [(np.argmin, 'buried'), (np.argmax, 'exposed')]:
                idx = extremum(vals)
                seed = idx // 5 + 1
                sample = idx % 5
                label = f'{sk}_{site}_{tag}'
                if label in seen:
                    continue
                seen.add(label)
                saa = vals[idx]
                desc = f'{sk.upper()}: {site} {tag} (SAA={saa:.3f})'
                picks.append((sk, seed, sample, label, desc))

    return picks


# ============================================================
# Step 1: Extract frames
# ============================================================

def extract_frames():
    import MDAnalysis as mda

    os.makedirs(FRAMES_DIR, exist_ok=True)

    picks = pick_frames()
    print(f"  {len(picks)} frames to extract\n")

    for sk, seed, sample, label, desc in picks:
        sysname = SETS_INFO[sk]
        sim_dir = os.path.join(CWD, sk, f'seed-{seed}_sample-{sample}')
        pdb = os.path.join(sim_dir, 'top.pdb')
        dcd = os.path.join(sim_dir, f'{sysname}.dcd')

        if not os.path.isfile(dcd):
            print(f"  SKIP: {label} — no DCD")
            continue

        u = mda.Universe(pdb, dcd)
        n = len(u.trajectory)
        frame_idx = int(n * 0.75)
        u.trajectory[frame_idx]

        out_pdb = os.path.join(FRAMES_DIR, f'{label}.pdb')
        u.atoms.write(out_pdb)
        print(f"  {label}: frame {frame_idx}/{n}  {desc}")

    # Write picks list for the session builder
    picks_file = os.path.join(FRAMES_DIR, 'picks.txt')
    with open(picks_file, 'w') as f:
        for sk, seed, sample, label, desc in picks:
            f.write(f'{sk}\t{seed}\t{sample}\t{label}\t{desc}\n')
    print(f"\n  Saved {len(picks)} frames to {FRAMES_DIR}/")


# ============================================================
# Step 2: Build PyMOL session
# ============================================================

def build_session():
    import pymol
    from pymol import cmd

    pymol.finish_launching(['pymol', '-cq'])

    # Read picks
    picks_file = os.path.join(FRAMES_DIR, 'picks.txt')
    picks = []
    with open(picks_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            picks.append((parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]))

    # Custom colors
    cmd.set_color('dark_purple', [0.40, 0.10, 0.50])
    cmd.set_color('medium_purple', [0.58, 0.30, 0.68])
    cmd.set_color('light_purple', [0.72, 0.50, 0.80])
    cmd.set_color('cat_gray', [0.55, 0.55, 0.55])

    # Load all
    loaded = []
    for sk, seed, sample, label, desc in picks:
        pdb_path = os.path.join(FRAMES_DIR, f'{label}.pdb')
        if not os.path.isfile(pdb_path):
            print(f"  SKIP: {pdb_path}")
            continue
        cmd.load(pdb_path, label)
        loaded.append((sk, label))

    print(f"  Loaded {len(loaded)} objects")

    # Style each object
    domain_colors = {
        'MD1': 'dark_purple', 'MD2': 'medium_purple',
        'MD3': 'light_purple', 'ART': 'red',
    }

    for sk, label in loaded:
        cmd.color('white', label)
        cmd.show('surface', label)
        cmd.hide('cartoon', label)
        cmd.hide('lines', label)

        dom_ranges = get_construct_domain_ranges(sk)
        catalytic = get_construct_catalytic(sk)

        # Color domains
        for dname, color in domain_colors.items():
            if dname in dom_ranges:
                s, e = dom_ranges[dname]
                cmd.color(color, f'{label} and resi {s}-{e}')

        # Color catalytic residues gray
        for sname, resids in catalytic.items():
            resi_str = '+'.join(str(r) for r in resids)
            cmd.color('cat_gray', f'{label} and resi {resi_str}')

    # --- Align all to first FL object on MD1 domain ---
    ref_label = None
    ref_sk = None
    for sk, label in loaded:
        if sk == 'fl':
            ref_label = label
            ref_sk = sk
            break
    if ref_label is None and loaded:
        ref_label = loaded[0][1]
        ref_sk = loaded[0][0]

    ref_dom = get_construct_domain_ranges(ref_sk)
    if 'MD1' in ref_dom:
        rs, re = ref_dom['MD1']
        ref_sel = f'{ref_label} and resi {rs}-{re}'
    else:
        ref_sel = ref_label

    for sk, label in loaded:
        if label == ref_label:
            continue
        tgt_dom = get_construct_domain_ranges(sk)
        if 'MD1' in tgt_dom:
            ts, te = tgt_dom['MD1']
            tgt_sel = f'{label} and resi {ts}-{te}'
        else:
            tgt_sel = label
        try:
            cmd.align(tgt_sel, ref_sel)
        except:
            try:
                cmd.align(label, ref_label)
            except:
                pass

    # Global settings
    cmd.set('surface_quality', 1)
    cmd.set('transparency', 0.0)
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadow', 0)
    cmd.bg_color('white')

    # Group by construct
    for sk in ['fl', 'core', 'mka', 'norrm', 'noart', 'md', 'fl_optimized']:
        members = [label for s, label in loaded if s == sk]
        if members:
            cmd.group(sk, ' '.join(members))

    # Disable all, enable FL
    cmd.disable('all')
    cmd.enable('fl')
    cmd.zoom('fl')

    session_path = os.path.join(CWD, 'parp14_accessibility.pse')
    cmd.save(session_path)
    print(f"\n  Session saved: {session_path}")
    print(f"  Open in PyMOL, orient, then run render_scenes.py on cluster")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--extract', action='store_true')
    parser.add_argument('--build-session', action='store_true')
    args = parser.parse_args()

    if not args.extract and not args.build_session:
        parser.print_help()
        print("\n  Run --extract first (calvados env), then --build-session (pymol-render env)")
        sys.exit(1)

    if args.extract:
        print("Extracting frames...")
        extract_frames()
    if args.build_session:
        print("Building PyMOL session...")
        build_session()
