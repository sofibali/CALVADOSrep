#!/usr/bin/env python3
"""
Build PyMOL session visualizing the PARP14 CALVADOS analysis.

For each construct (fl, md, core, mka): loads a representative replicate
trajectory (seed-1_sample-0, subsampled to ~45 states), colors by domain,
highlights active site catalytic + pocket residues, draws inter-site distance
measurements, and groups everything for easy toggling.

Usage:
    conda run -n pymol-render python make_analysis_pse.py

Objects hierarchy in PyMOL:
  FL/
    FL_traj          — trajectory (cartoon, domain-colored)
    FL_MD1_cat       — MD1 catalytic residues (spheres)
    FL_MD1_pocket    — MD1 pocket residues (sticks)
    FL_MD2_cat, ...
    FL_dist_MD1_MD2  — distance dash between site COMs
    FL_dist_MD1_MD3  — ...
    FL_dist_MD1_ART
    FL_dist_MD2_MD3
    FL_dist_MD2_ART
    FL_dist_MD3_ART
  MD/
    ...
  CORE/
    ...
  MKA/
    ...
"""

import os
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ============================================================
# Domain definitions for coloring (FL numbering)
# ============================================================

DOMAINS_FL = {
    'RRM1':   (1, 145),
    'RRM2':   (146, 224),
    'RRM3':   (225, 314),
    'KH1-6':  (315, 737),
    'KH7a':   (738, 789),
    'MD1L1':  (790, 1004),
    'MD2':    (1004, 1193),
    'MD3':    (1207, 1388),
    'KHb-8':  (1389, 1533),
    'WWE':    (1534, 1602),
    'ART':    (1603, 1801),
}

DOMAIN_COLORS_RGB = {
    'RRM1':   [0.12, 0.47, 0.71],
    'RRM2':   [0.17, 0.63, 0.17],
    'RRM3':   [0.20, 0.73, 0.20],
    'KH1-6':  [1.00, 0.50, 0.05],
    'KH7a':   [0.84, 0.15, 0.16],
    'MD1L1':  [0.68, 0.78, 0.91],
    'MD2':    [0.58, 0.40, 0.74],
    'MD3':    [0.55, 0.34, 0.29],
    'KHb-8':  [0.89, 0.47, 0.76],
    'WWE':    [0.50, 0.50, 0.50],
    'ART':    [0.74, 0.74, 0.13],
}

# Active sites (FL numbering)
ACTIVE_SITES_FL = {
    'MD1': {
        'catalytic': [831, 923, 962],
        'cat_labels': ['D831', 'N923', 'D962'],
        'pocket': [822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833,
                   834, 835, 836, 919, 920, 921, 922, 923, 924, 925, 926, 927,
                   961, 962, 966],
    },
    'MD2': {
        'catalytic': [1035, 1046, 1134, 1171],
        'cat_labels': ['G1035', 'I1046', 'G1134', 'D1171'],
        'pocket': [1021, 1022, 1023, 1024, 1034, 1035, 1036, 1037, 1038, 1039,
                   1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1130, 1131,
                   1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141,
                   1170, 1171, 1175, 1178],
    },
    'MD3': {
        'catalytic': [1248, 1259, 1330, 1371],
        'cat_labels': ['G1248', 'V1259', 'G1330', 'N1371'],
        'pocket': [1235, 1236, 1237, 1247, 1248, 1249, 1250, 1251, 1252, 1253,
                   1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1302, 1303,
                   1304, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332,
                   1333, 1334, 1335, 1336, 1337, 1369, 1370, 1371, 1375],
    },
    'ART': {
        'catalytic': [1684, 1705, 1706, 1722],
        'cat_labels': ['H1684', 'Y1705', 'E1706', 'I1722'],
        'pocket': [1681, 1682, 1683, 1684, 1685, 1688, 1701, 1704, 1705, 1706,
                   1707, 1708, 1709, 1714, 1715, 1716, 1721, 1722, 1726, 1727,
                   1781],
    },
}

SITE_COLORS = {
    'MD1': [0.90, 0.10, 0.29],
    'MD2': [0.24, 0.71, 0.29],
    'MD3': [0.26, 0.39, 0.85],
    'ART': [0.96, 0.51, 0.19],
}

# ============================================================
# Construct definitions
# ============================================================

DOMAIN_UNITS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789), 'md1l1': (790, 1004),
    'md2': (1004, 1193), 'md3': (1207, 1388), 'khb-kh8': (1389, 1533),
    'wwe': (1534, 1602), 'art': (1603, 1801),
}

SETS = {
    'fl': {
        'sysname': 'parp14',
        'label': 'Full-length',
        'units': list(DOMAIN_UNITS.keys()),
        'sites': ['MD1', 'MD2', 'MD3', 'ART'],
    },
    'md': {
        'sysname': 'parp14_macrodomains',
        'label': 'Macrodomains',
        'units': ['md1l1', 'md2', 'md3'],
        'sites': ['MD1', 'MD2', 'MD3'],
    },
    'core': {
        'sysname': 'parp14_core',
        'label': 'Core',
        'units': ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
        'sites': ['MD1', 'MD2', 'MD3', 'ART'],
    },
    'mka': {
        'sysname': 'parp14_mka',
        'label': 'MKA',
        'units': ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
        'sites': ['MD1', 'MD2', 'MD3', 'ART'],
    },
    'fl_optimized': {
        'sysname': 'parp14',
        'label': 'FL (optimized restraints)',
        'units': list(DOMAIN_UNITS.keys()),
        'sites': ['MD1', 'MD2', 'MD3', 'ART'],
    },
}


def compute_fl_blocks(unit_names):
    ranges = sorted([DOMAIN_UNITS[n] for n in unit_names])
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
    construct_pos = 1
    for fl_start, fl_end in fl_blocks:
        offset = construct_pos - fl_start
        segments.append((fl_start, fl_end, offset))
        construct_pos += (fl_end - fl_start + 1)

    def map_resid(fl_resid):
        for fl_s, fl_e, off in segments:
            if fl_s <= fl_resid <= fl_e:
                return fl_resid + off
        return None
    return map_resid


def get_construct_domains(set_key):
    """Map FL domain ranges to construct numbering for coloring."""
    if set_key == 'fl':
        return DOMAINS_FL

    fl_to_c = build_fl_to_construct_map(SETS[set_key]['units'])
    fl_blocks = compute_fl_blocks(SETS[set_key]['units'])
    construct_fl_residues = set()
    for s, e in fl_blocks:
        construct_fl_residues.update(range(s, e + 1))

    mapped = {}
    for dname, (fl_s, fl_e) in DOMAINS_FL.items():
        # Find overlap with construct
        cs = fl_s
        while cs <= fl_e and cs not in construct_fl_residues:
            cs += 1
        ce = fl_e
        while ce >= cs and ce not in construct_fl_residues:
            ce -= 1
        if cs > ce:
            continue
        c_start = fl_to_c(cs)
        c_end = fl_to_c(ce)
        if c_start is not None and c_end is not None:
            mapped[dname] = (c_start, c_end)
    return mapped


def get_construct_active_sites(set_key):
    """Map FL active site residues to construct numbering."""
    if set_key == 'fl':
        return ACTIVE_SITES_FL

    fl_to_c = build_fl_to_construct_map(SETS[set_key]['units'])
    sites = {}
    for sname in SETS[set_key]['sites']:
        fl_data = ACTIVE_SITES_FL[sname]
        cat = []
        cat_labels = []
        for r, lab in zip(fl_data['catalytic'], fl_data['cat_labels']):
            cr = fl_to_c(r)
            if cr is not None:
                cat.append(cr)
                cat_labels.append(lab)
        pocket = [fl_to_c(r) for r in fl_data['pocket'] if fl_to_c(r) is not None]
        sites[sname] = {'catalytic': cat, 'cat_labels': cat_labels, 'pocket': pocket}
    return sites


# ============================================================
# Build session
# ============================================================

def setup_colors():
    """Register custom colors."""
    for dname, rgb in DOMAIN_COLORS_RGB.items():
        cname = f'dom_{dname.replace("-", "_")}'
        cmd.set_color(cname, rgb)
    for sname, rgb in SITE_COLORS.items():
        cmd.set_color(f'site_{sname}', rgb)


def build_set(set_key):
    """Load trajectory and build all visual objects for one construct."""
    info = SETS[set_key]
    prefix = set_key.upper()
    sim_dir = os.path.join(CWD, f'{set_key}_seed-1_sample-0')
    pdb_path = os.path.join(sim_dir, 'top.pdb')
    dcd_path = os.path.join(sim_dir, f'{info["sysname"]}.dcd')

    if not os.path.isfile(pdb_path) or not os.path.isfile(dcd_path):
        print(f"  SKIP {set_key}: files not found")
        return

    traj_name = f'{prefix}_traj'
    group_members = [traj_name]

    # --- Load trajectory (subsample: skip first 50 frames, every 10th) ---
    print(f"  Loading {set_key} trajectory...")
    cmd.load(pdb_path, traj_name)
    cmd.load_traj(dcd_path, traj_name, start=51, interval=10)
    n_states = cmd.count_states(traj_name)
    print(f"    {n_states} states loaded")

    # --- Style: cartoon backbone ---
    cmd.hide('everything', traj_name)
    cmd.show('cartoon', traj_name)
    cmd.set('cartoon_trace_atoms', 1, traj_name)
    cmd.set('cartoon_tube_radius', 0.4, traj_name)
    cmd.color('gray70', traj_name)

    # --- Color by domain ---
    domains = get_construct_domains(set_key)
    for dname, (s, e) in domains.items():
        cname = f'dom_{dname.replace("-", "_")}'
        cmd.color(cname, f'{traj_name} and resi {s}-{e}')

    # --- Active sites ---
    active_sites = get_construct_active_sites(set_key)
    available_sites = SETS[set_key]['sites']

    for sname in available_sites:
        site_data = active_sites[sname]
        if not site_data['catalytic']:
            continue

        # Catalytic residues: spheres with labels
        cat_resi = '+'.join(str(r) for r in site_data['catalytic'])
        cat_obj = f'{prefix}_{sname}_cat'
        cmd.select(cat_obj, f'{traj_name} and resi {cat_resi} and name CA')
        cmd.show('spheres', cat_obj)
        cmd.set('sphere_scale', 1.5, cat_obj)
        cmd.color(f'site_{sname}', cat_obj)
        # Label catalytic residues
        for r, lab in zip(site_data['catalytic'], site_data['cat_labels']):
            cmd.label(f'{traj_name} and resi {r} and name CA', f'"{lab}"')
        group_members.append(cat_obj)

        # Pocket residues: sticks (slightly transparent)
        pocket_resi = '+'.join(str(r) for r in site_data['pocket'])
        pocket_obj = f'{prefix}_{sname}_pocket'
        cmd.select(pocket_obj, f'{traj_name} and resi {pocket_resi} and name CA')
        cmd.show('sticks', pocket_obj)
        cmd.set('stick_radius', 0.3, pocket_obj)
        cmd.color(f'site_{sname}', pocket_obj)
        cmd.set('stick_transparency', 0.4, pocket_obj)
        group_members.append(pocket_obj)

    # --- Distance measurements between active site pairs ---
    # Use distance on catalytic residue selections (averages over CA atoms)
    site_sels = {}
    for sname in available_sites:
        site_data = active_sites[sname]
        if site_data['catalytic']:
            cat_resi = '+'.join(str(r) for r in site_data['catalytic'])
            sel_name = f'_tmp_{prefix}_{sname}'
            cmd.select(sel_name, f'{traj_name} and resi {cat_resi} and name CA')
            site_sels[sname] = sel_name

    # Pairwise distances
    avail = [s for s in available_sites if s in site_sels]
    for i in range(len(avail)):
        for j in range(i + 1, len(avail)):
            si, sj = avail[i], avail[j]
            dist_name = f'{prefix}_dist_{si}_{sj}'
            cmd.distance(dist_name, site_sels[si], site_sels[sj])
            cmd.set('dash_color', SITE_COLORS[si], dist_name)
            cmd.set('dash_width', 3.0, dist_name)
            cmd.set('dash_gap', 0.3, dist_name)
            cmd.set('dash_length', 0.5, dist_name)
            cmd.set('label_size', 14, dist_name)
            cmd.set('label_color', 'black', dist_name)
            group_members.append(dist_name)

    # --- Radial distance: protein COM pseudoatom to each site ---
    # Create a pseudoatom at the protein COM (state-dependent)
    com_name = f'{prefix}_COM'
    cmd.pseudoatom(com_name, f'{traj_name}', name='COM', label='COM')
    cmd.hide('everything', com_name)
    cmd.show('spheres', com_name)
    cmd.set('sphere_scale', 2.0, com_name)
    cmd.color('white', com_name)
    group_members.append(com_name)

    for sname in avail:
        radial_name = f'{prefix}_radial_{sname}'
        cmd.distance(radial_name, com_name, site_sels[sname])
        cmd.set('dash_color', 'gray50', radial_name)
        cmd.set('dash_width', 1.5, radial_name)
        cmd.set('dash_gap', 0.2, radial_name)
        cmd.set('dash_length', 0.3, radial_name)
        cmd.set('label_size', 12, radial_name)
        cmd.disable(radial_name)  # off by default (cluttered)
        group_members.append(radial_name)

    # Clean up temp selections
    for sel_name in site_sels.values():
        cmd.delete(sel_name)

    # --- Group everything ---
    cmd.group(prefix, ' '.join(group_members))

    # Sub-groups within construct
    site_members = []
    for sname in available_sites:
        cat_obj = f'{prefix}_{sname}_cat'
        pocket_obj = f'{prefix}_{sname}_pocket'
        site_members.extend([cat_obj, pocket_obj])
    if site_members:
        cmd.group(f'{prefix}_active_sites', ' '.join(site_members))

    dist_members = [m for m in group_members if '_dist_' in m]
    if dist_members:
        cmd.group(f'{prefix}_distances', ' '.join(dist_members))

    radial_members = [m for m in group_members if '_radial_' in m]
    if radial_members:
        radial_members.append(com_name)
        cmd.group(f'{prefix}_radial', ' '.join(radial_members))

    return prefix


# ============================================================
# Main
# ============================================================

print("=" * 60)
print("Building PARP14 Analysis PyMOL Session")
print("=" * 60)

setup_colors()

built = []
for set_key in ['fl', 'md', 'core', 'mka', 'fl_optimized']:
    name = build_set(set_key)
    if name:
        built.append(name)

# --- Global settings ---
cmd.set('label_color', 'black')
cmd.set('label_size', 14)
cmd.set('label_font_id', 7)
cmd.set('label_position', [0, 0, 4])
cmd.set('cartoon_fancy_helices', 0)
cmd.set('cartoon_smooth_loops', 1)
cmd.set('ray_opaque_background', 1)
cmd.bg_color('white')
cmd.set('antialias', 2)
cmd.set('dash_round_ends', 0)
cmd.set('state', 1)

# Enable only FL by default; others off
for name in built:
    if name == 'FL':
        cmd.enable(name)
    else:
        cmd.disable(name)

# Orient to FL
cmd.orient('FL_traj')
cmd.zoom('FL_traj', buffer=10)
cmd.deselect()

# Save
output_path = os.path.join(CWD, 'data', 'PARP14_analysis.pse')
cmd.save(output_path)
print(f"\nSaved: {output_path}")
print()
print("Session structure:")
print("  FL/               — Full-length (enabled by default)")
print("    FL_traj           trajectory (cartoon, domain-colored)")
print("    FL_active_sites/  catalytic spheres + pocket sticks per site")
print("    FL_distances/     inter-site distance dashes")
print("    FL_radial/        COM-to-site distances (disabled by default)")
print("  MD/               — Macrodomains only")
print("  CORE/             — Core construct")
print("  MKA/              — MKA construct")
print()
print("Usage:")
print("  Toggle groups on/off in the PyMOL object panel")
print("  Use the state slider to scrub through trajectory frames")
print("  Distances update per-frame as you scrub")

cmd.quit()
