#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a PyMOL session (.pse) from the representative state PDBs produced
by cluster_states.py.

For each clustering run (e.g. md_ca25_tica), loads the K representative
state PDBs from representative_frames/<date>/<run_tag>/, aligns them all
to state_1, colors each state distinctly, highlights catalytic residues
in red, and saves a session.

Outputs:
  figures/05_clustering/<date>/<run_tag>/cluster_states.pse
  figures/05_clustering/<date>/<run_tag>/cluster_states.pml  (script)

Usage:
    python make_cluster_pse.py --auto --set md
    python make_cluster_pse.py --auto --set fl_optimized
    python make_cluster_pse.py --auto --set md --run-tag md_ca25_tica
    python make_cluster_pse.py --states representative_frames/2026-05-19/md_ca25_tica/state_*.pdb
"""
import os
import sys
import argparse
import json
from pathlib import Path

CWD = Path(__file__).resolve().parent


# ─────────────────────────────────────────────────────────────────
# Domain layouts (FL numbering for fl/fl_optimized; construct-local for others)
# ─────────────────────────────────────────────────────────────────

# PyMOL color names and FL residue ranges per domain.
# Uses the dark/medium/light purple gradient for MD1L1/MD2/MD3 + red for
# ART (consistent with make_pymol_session.py). Other domains use the
# same FL palette as figure_sasa_faces.py / figure_md_distances.py.
FL_DOMAIN_COLORS = {
    'RRM1':    ('rrm1_col',   (1, 145),    (0.12, 0.46, 0.71)),
    'RRM2':    ('rrm2_col',   (146, 224),  (0.17, 0.63, 0.17)),
    'RRM3':    ('rrm3_col',   (225, 314),  (0.68, 0.78, 0.91)),
    'KH1-KH6': ('kh16_col',   (315, 737),  (1.00, 0.50, 0.05)),
    'KH7a':    ('kh7a_col',   (738, 789),  (0.84, 0.15, 0.16)),
    # Pink → blue gradient for MD1L1/MD2/MD3, picked for max mutual
    # contrast AND distinctness from RRM1/3 (blues), KHb-KH8 (soft pink).
    'MD1L1':   ('md1_hotpink',   (790, 1004),  (0.84, 0.13, 0.43)),
    'MD2':     ('md2_navy',      (1005, 1193), (0.00, 0.20, 0.63)),
    'MD3':     ('md3_cyan',      (1207, 1388), (0.00, 0.78, 0.78)),
    'KHb-KH8': ('khb_col',    (1389, 1533), (0.89, 0.47, 0.76)),
    'WWE':     ('wwe_col',    (1534, 1602), (0.50, 0.50, 0.50)),
    'ART':     ('red',        (1603, 1801), (1.00, 0.00, 0.00)),
}

# Catalytic residues in FL numbering
CATALYTIC_FL = {
    'MD1':  [831, 923, 962],
    'MD2':  [1035, 1046, 1134, 1171],
    'MD3':  [1248, 1259, 1330, 1371],
    'ART':  [1684, 1705, 1706, 1722],
}

DOMAIN_UNITS = {
    'rrm1':    (1, 145),    'rrm2':    (146, 224), 'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),  'kh7a':    (738, 789),
    'md1l1':   (790, 1004), 'md2':     (1004, 1193), 'md3':   (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe':    (1534, 1602), 'art':   (1603, 1801),
}

CONSTRUCT_UNITS = {
    'fl':           list(DOMAIN_UNITS.keys()),
    'fl_optimized': list(DOMAIN_UNITS.keys()),
    'md':           ['md1l1', 'md2', 'md3'],
    'mka':          ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'core':         ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'norrm':        ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8',
                     'wwe', 'art'],
    'noart':        ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8',
                     'wwe'],
    'md3art':       ['md3', 'khb-kh8', 'wwe', 'art'],
}


# State colors (consistent with cluster plots — tab10)
STATE_COLORS = [
    'tv_red', 'tv_orange', 'yellow', 'limegreen', 'cyan',
    'slate', 'magenta', 'salmon', 'palegreen', 'lightblue',
]


def compute_fl_to_construct(set_key):
    """Build FL→construct residue map for a named set or fragment."""
    if set_key in ('fl', 'fl_optimized'):
        return {r: r for r in range(1, 1802)}  # identity for FL

    if set_key.startswith('frag_'):
        frag = set_key[5:]
        meta_path = CWD / 'fragments' / frag / 'metadata.json'
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
            units = meta.get('units', [])
        else:
            return {}
    else:
        units = CONSTRUCT_UNITS.get(set_key, [])

    # Compute continuous FL blocks (merging overlaps)
    ranges = sorted([DOMAIN_UNITS[u] for u in units if u in DOMAIN_UNITS])
    if not ranges:
        return {}
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    mapping = {}
    pos = 1
    for fl_s, fl_e in merged:
        for fl_r in range(fl_s, fl_e + 1):
            mapping[fl_r] = pos
            pos += 1
    return mapping


def auto_find_run_tag(set_key, date_str=None):
    """Locate latest representative_frames/<date>/<tag>/ for this set."""
    rep_root = CWD / 'representative_frames'
    if not rep_root.is_dir():
        return None

    if date_str:
        date_dirs = [rep_root / date_str]
    else:
        date_dirs = sorted(rep_root.glob('20*'), reverse=True)

    for d in date_dirs:
        if not d.is_dir():
            continue
        # Find run_tag dirs containing this set
        candidates = sorted(d.iterdir(), reverse=True)
        for c in candidates:
            if not c.is_dir():
                continue
            # Match by prefix: e.g. 'md_ca25_tica' starts with 'md_'
            if c.name.startswith(f'{set_key}_'):
                pdbs = sorted(c.glob('state_*.pdb'))
                if pdbs:
                    return c, pdbs
    return None, []


def write_pml_script(states, set_key, fig_dir, run_tag):
    """Generate a .pml script that loads, aligns, colors, and saves."""
    fl_to_c = compute_fl_to_construct(set_key)
    catalytic_construct = {}
    for site, resids in CATALYTIC_FL.items():
        catalytic_construct[site] = [
            fl_to_c[r] for r in resids if r in fl_to_c]

    # Map each domain's FL range → construct range
    domains_construct = {}
    for dname, (color_name, (fs, fe), rgb) in FL_DOMAIN_COLORS.items():
        in_construct = [fl_to_c[r] for r in range(fs, fe + 1)
                        if r in fl_to_c]
        if in_construct:
            domains_construct[dname] = (color_name, rgb,
                                         (min(in_construct),
                                          max(in_construct)))

    # Decide alignment selection: prefer ART if present, else MD1L1
    if 'ART' in domains_construct:
        align_dom = 'ART'
    elif 'MD1L1' in domains_construct:
        align_dom = 'MD1L1'
    elif 'MD3' in domains_construct:
        align_dom = 'MD3'
    elif domains_construct:
        align_dom = next(iter(domains_construct))
    else:
        align_dom = None

    if align_dom:
        cs, ce = domains_construct[align_dom][2]
        align_sel = f'resi {cs}-{ce} and name CA'
    else:
        align_sel = 'name CA'

    pml_path = fig_dir / 'cluster_states.pml'
    pse_path = fig_dir / 'cluster_states.pse'

    state_objs = [s.stem for s in states]
    ref_obj = state_objs[0] if state_objs else 'state_1'

    with open(pml_path, 'w') as f:
        f.write(f"# Cluster states PyMOL session\n")
        f.write(f"# Run tag: {run_tag}\n")
        f.write(f"# Set: {set_key}\n")
        f.write(f"# {len(states)} states\n")
        f.write(f"# Aligned on: {align_dom} (sel: {align_sel})\n\n")

        # PyMOL settings
        f.write("bg_color white\n")
        f.write("set cartoon_transparency, 0\n")
        f.write("set ray_opaque_background, 1\n")
        f.write("set ray_shadow, 0\n")
        f.write("set cartoon_trace_atoms, 1\n")
        f.write("set cartoon_tube_radius, 0.25\n")
        f.write("\n")

        # Define custom colors
        f.write("# Custom colors\n")
        defined = set()
        for dname, (color_name, rgb, _) in domains_construct.items():
            # Skip built-in PyMOL colors
            if color_name in ('red', 'green', 'blue', 'yellow', 'cyan',
                               'magenta', 'orange', 'purple'):
                continue
            if color_name in defined:
                continue
            f.write(f"set_color {color_name}, "
                    f"[{rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}]\n")
            defined.add(color_name)
        f.write("\n")

        # Load each state PDB
        f.write("# Load state PDBs\n")
        for state_pdb in states:
            obj = state_pdb.stem
            f.write(f"load {state_pdb}, {obj}\n")
        f.write("\n")

        # Align all to ref using only the chosen alignment domain
        f.write(f"# Align all to {ref_obj} using {align_dom} CAs\n")
        for obj in state_objs[1:]:
            f.write(f"align {obj} and {align_sel}, "
                    f"{ref_obj} and {align_sel}\n")
        f.write("\n")

        # Apply DOMAIN coloring to all state objects (default style)
        f.write("# Domain coloring: applied to every state object\n")
        for obj in state_objs:
            f.write(f"hide everything, {obj}\n")
            f.write(f"show cartoon, {obj}\n")
            f.write(f"color gray80, {obj}\n")
            for dname, (color_name, _, (cs, ce)) in domains_construct.items():
                f.write(f"color {color_name}, "
                        f"{obj} and resi {cs}-{ce}\n")
        f.write("\n")

        # Highlight catalytic residues with bright spheres on each state
        f.write("# Catalytic residues (yellow spheres)\n")
        all_cat = sorted({r for resids in catalytic_construct.values()
                          for r in resids})
        if all_cat:
            cat_str = '+'.join(map(str, all_cat))
            for obj in state_objs:
                f.write(f"select cat_{obj}, {obj} and resi {cat_str} "
                        f"and name CA\n")
                f.write(f"show spheres, cat_{obj}\n")
                f.write(f"color yellow, cat_{obj}\n")
                f.write(f"set sphere_scale, 0.8, cat_{obj}\n")
            f.write("\n")

        # Group state objects
        f.write("# Group\n")
        f.write(f"group {run_tag}_states, {' '.join(state_objs)}\n")
        f.write("\n")

        # ── Scenes ──
        # Start: orient on ref and the chosen domain
        f.write("# Scenes\n")
        f.write(f"orient {ref_obj} and {align_sel}\n\n")

        # Scene 1: state_1 only (domain colored — base reference view)
        f.write(f"# Scene 1: domain-colored state_1\n")
        f.write(f"disable all\n")
        f.write(f"enable {ref_obj}\n")
        f.write(f"scene S1_alone, store\n\n")

        # Overlay scenes: state_1 + state_N, with state_1 dimmed
        # to make state_N pop
        for i, obj in enumerate(state_objs[1:], start=2):
            f.write(f"# Scene: state_1 + state_{i} overlay\n")
            f.write(f"disable all\n")
            f.write(f"enable {ref_obj}\n")
            f.write(f"enable {obj}\n")
            # Dim state_1 slightly so state_N is more visible
            f.write(f"set cartoon_transparency, 0.55, {ref_obj}\n")
            f.write(f"set cartoon_transparency, 0.0, {obj}\n")
            f.write(f"scene S1_vs_S{i}, store\n")
            # Reset transparency for next scene
            f.write(f"set cartoon_transparency, 0.0, {ref_obj}\n")
            f.write("\n")

        # All states overlayed, each transparent
        f.write(f"# Scene: all states overlayed (each semi-transparent)\n")
        f.write(f"disable all\n")
        for obj in state_objs:
            f.write(f"enable {obj}\n")
            f.write(f"set cartoon_transparency, 0.4, {obj}\n")
        f.write(f"scene All_states_overlay, store\n")
        # Restore
        for obj in state_objs:
            f.write(f"set cartoon_transparency, 0.0, {obj}\n")
        f.write("\n")

        # Per-state alone (useful for screenshots)
        for i, obj in enumerate(state_objs, start=1):
            f.write(f"disable all\n")
            f.write(f"enable {obj}\n")
            f.write(f"scene S{i}_alone, store\n")
        f.write("\n")

        # Default view: S1_alone
        f.write(f"scene S1_alone, recall\n")
        f.write(f"save {pse_path}\n")

    return pml_path, pse_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--auto', action='store_true',
                        help='Auto-locate latest cluster output for --set')
    parser.add_argument('--set', dest='set_key',
                        help='Set name (md, fl_optimized, etc.) or fragment name')
    parser.add_argument('--run-tag', default=None,
                        help='Specific run_tag (e.g. md_ca25_tica). '
                             'Default: pick first matching directory.')
    parser.add_argument('--date', default=None,
                        help='Specific date dir (default: latest)')
    parser.add_argument('--states', nargs='+', type=Path,
                        help='State PDB files (overrides --auto)')
    parser.add_argument('--out-pml-only', action='store_true',
                        help='Only write the .pml script; run pymol yourself')
    args = parser.parse_args()

    # Resolve state PDBs
    if args.states:
        states = [Path(s).resolve() for s in args.states]
        run_tag = states[0].parent.name
        set_key = args.set_key or run_tag.split('_')[0]
        fig_dir = states[0].parent.parent.parent.parent / 'figures' / \
            '05_clustering' / states[0].parent.parent.name / run_tag
    elif args.auto:
        if not args.set_key:
            parser.error('--auto requires --set')
        set_key = args.set_key
        if args.run_tag:
            # User specified exact tag — find it under any date
            rep_root = CWD / 'representative_frames'
            date_dirs = sorted(rep_root.glob('20*'), reverse=True) \
                if not args.date else [rep_root / args.date]
            run_dir = None
            for d in date_dirs:
                cand = d / args.run_tag
                if cand.is_dir():
                    run_dir = cand
                    break
            if run_dir is None:
                print(f"ERROR: no run dir found for tag={args.run_tag}")
                return 1
            states = sorted(run_dir.glob('state_*.pdb'))
            run_tag = args.run_tag
        else:
            result = auto_find_run_tag(set_key, args.date)
            if result is None or not result[1]:
                print(f"ERROR: no clustering output found for set={set_key}")
                return 1
            run_dir, states = result
            run_tag = run_dir.name
        fig_dir = (CWD / 'figures' / '05_clustering' /
                   run_dir.parent.name / run_tag)
    else:
        parser.error('Specify --auto or --states')

    fig_dir.mkdir(parents=True, exist_ok=True)
    print(f"Set:     {set_key}")
    print(f"Run tag: {run_tag}")
    print(f"States:  {len(states)}")
    for s in states:
        print(f"  {s.name}")
    print(f"Output:  {fig_dir}/")

    # Translate set_key for FL-mapping logic
    if set_key not in ('fl', 'fl_optimized', 'md', 'core', 'mka',
                       'norrm', 'noart', 'md3art'):
        # Treat as fragment
        set_key_internal = f'frag_{set_key}'
    else:
        set_key_internal = set_key

    pml_path, pse_path = write_pml_script(states, set_key_internal,
                                            fig_dir, run_tag)
    print(f"\nWrote: {pml_path}")

    if args.out_pml_only:
        print(f"To generate session: pymol {pml_path}")
        return 0

    # Try to run pymol -cq directly
    import subprocess
    pymol_exe = subprocess.run(['which', 'pymol'],
                                capture_output=True, text=True).stdout.strip()
    if not pymol_exe:
        # Try the pymol-render conda env (mentioned in CLAUDE.md)
        candidates = [
            '/home/sbali/miniconda3/envs/pymol-render/bin/pymol',
            '/home/sbali/miniconda3/envs/pymol/bin/pymol',
        ]
        for c in candidates:
            if Path(c).exists():
                pymol_exe = c
                break

    if pymol_exe:
        print(f"\nRunning: {pymol_exe} -cq {pml_path}")
        result = subprocess.run([pymol_exe, '-cq', str(pml_path)],
                                 capture_output=True, text=True,
                                 timeout=300)
        if result.returncode == 0:
            print(f"Saved PSE: {pse_path}")
        else:
            print(f"PyMOL failed: {result.stderr[:500]}")
            print(f"You can run it yourself: pymol -cq {pml_path}")
    else:
        print(f"PyMOL not found in PATH. Run yourself:")
        print(f"  conda activate pymol-render")
        print(f"  pymol -cq {pml_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
