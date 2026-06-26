#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Active site steric accessibility analysis for PARP14 CALVADOS simulations.

Question: Can a second PARP14's ART domain (~3 nm) physically approach each
active site pocket, or is it sterically blocked by the protein's own chain?

Method — Solid-Angle Accessibility (SAA):
  For each active site pocket COM, cast N_RAYS uniformly distributed rays
  outward. A ray is "blocked" if any non-pocket protein bead falls within
  PROBE_RADIUS of the ray path (cylinder test) out to MAX_DIST. The fraction
  of unblocked rays = accessibility (0 = fully buried, 1 = fully exposed).

  PROBE_RADIUS = 0.5 nm  (CG bead clearance for approach path)
  MAX_DIST     = 5.0 nm  (how far out to check for obstruction)
  N_RAYS       = 200     (Fibonacci sphere sampling)

Also computes:
  - Half-space density: fraction of a sphere shell (2-5 nm) occupied by
    protein beads around each site. Low density = accessible.
  - Approach cone angle: the widest cone (from site COM) that contains
    no protein beads, giving the effective "docking angle" available.

Usage:
    python analyze_accessibility.py                  # all sets
    python analyze_accessibility.py --set fl md      # specific sets
    python analyze_accessibility.py --nrays 500      # more rays (slower)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
import os
import warnings
from argparse import ArgumentParser
from collections import defaultdict

import sim_registry as reg

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(CWD, 'data')
# Dated category subdir for figures: figures/03_accessibility/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_PATH = str(_get_fig_dir('03_accessibility'))
os.makedirs(DATA_PATH, exist_ok=True)

SEEDS = range(1, 6)
SAMPLES = range(0, 5)
SKIP_FRAMES = 50

# Accessibility parameters
PROBE_RADIUS = 0.5    # nm — CG bead radius; tests if approach path is sterically clear
MAX_DIST = 5.0        # nm — how far to check for obstruction
N_RAYS = 200          # rays per site per frame (Fibonacci sphere)
SEQ_SEP = 10          # ignore beads within ±10 residues of pocket (own fold)

# Set definitions
SETS = {
    'fl':   {'sysname': 'parp14',              'label': 'Full-length (1801)'},
    'md':   {'sysname': 'parp14_macrodomains', 'label': 'Macrodomains (586)'},
    'core': {'sysname': 'parp14_core',         'label': 'Core (1051)'},
    'mka':  {'sysname': 'parp14_mka',          'label': 'MKA (999)'},
    'fl_optimized': {'sysname': 'parp14',      'label': 'FL (optimized restraints)'},
}

# Active sites (FL numbering)
ACTIVE_SITES_FL = {
    'MD1': {'catalytic': [831, 923, 962],
            'pocket': [822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833,
                       834, 835, 836, 919, 920, 921, 922, 923, 924, 925, 926, 927,
                       961, 962, 966]},
    'MD2': {'catalytic': [1035, 1046, 1134, 1171],
            'pocket': [1021, 1022, 1023, 1024, 1034, 1035, 1036, 1037, 1038, 1039,
                       1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1130, 1131,
                       1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141,
                       1170, 1171, 1175, 1178]},
    'MD3': {'catalytic': [1248, 1259, 1330, 1371],
            'pocket': [1235, 1236, 1237, 1247, 1248, 1249, 1250, 1251, 1252, 1253,
                       1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1302, 1303,
                       1304, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332,
                       1333, 1334, 1335, 1336, 1337, 1369, 1370, 1371, 1375]},
    'ART': {'catalytic': [1684, 1705, 1706, 1722],
            'pocket': [1681, 1682, 1683, 1684, 1685, 1688, 1701, 1704, 1705, 1706,
                       1707, 1708, 1709, 1714, 1715, 1716, 1721, 1722, 1726, 1727,
                       1781]},
}

SITE_NAMES = ['MD1', 'MD2', 'MD3', 'ART']
SITE_COLORS = {'MD1': '#e6194b', 'MD2': '#3cb44b', 'MD3': '#4363d8', 'ART': '#f58231'}

CONSTRUCT_SITES = {
    'fl': ['MD1', 'MD2', 'MD3', 'ART'],
    'md': ['MD1', 'MD2', 'MD3'],
    'core': ['MD1', 'MD2', 'MD3', 'ART'],
    'mka': ['MD1', 'MD2', 'MD3', 'ART'],
}

# Domain unit mapping (for construct numbering)
DOMAIN_UNITS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789), 'md1l1': (790, 1004),
    'md2': (1004, 1193), 'md3': (1207, 1388), 'khb-kh8': (1389, 1533),
    'wwe': (1534, 1602), 'art': (1603, 1801),
}

CONSTRUCT_UNITS = {
    'fl': list(DOMAIN_UNITS.keys()),
    'md': ['md1l1', 'md2', 'md3'],
    'core': ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'mka': ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
}

# Directories registered via --sim-folder: set_key -> absolute folder path.
EXTERNAL_DIRS = {}


def register_sim_folder(path, units=None):
    """Register an arbitrary simulation folder (from --sim-folder) as a set so
    accessibility analysis can process it. Domain units are read from the folder's
    metadata.json, the `units` arg, or (if the basename is a known set) that set.
    Returns the set_key under which it is registered.
    """
    folder = os.path.abspath(path)
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"sim folder not found: {folder}")
    resolved_units = reg._resolve_units(folder, units)
    meta = reg.read_metadata(folder) or {}
    sysname = meta.get('sysname') or reg.detect_sysname(folder)
    set_key = os.path.basename(os.path.normpath(folder))
    SETS[set_key] = {'sysname': sysname,
                     'label': meta.get('label') or set_key,
                     'color': '#444444'}
    CONSTRUCT_UNITS[set_key] = resolved_units
    CONSTRUCT_SITES[set_key] = [reg.UNIT_TO_SITE[u] for u in resolved_units
                                if u in reg.UNIT_TO_SITE]
    EXTERNAL_DIRS[set_key] = folder
    return set_key


# ============================================================
# FL-to-construct mapping (reused from analyze_all.py)
# ============================================================

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


def get_active_sites_for_set(set_key):
    if set_key == 'fl':
        return ACTIVE_SITES_FL
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[set_key])
    sites = {}
    for sname in CONSTRUCT_SITES[set_key]:
        data = ACTIVE_SITES_FL[sname]
        sites[sname] = {
            'catalytic': [r for r in (fl_to_c(x) for x in data['catalytic']) if r is not None],
            'pocket': [r for r in (fl_to_c(x) for x in data['pocket']) if r is not None],
        }
    return sites


# ============================================================
# Ray-based accessibility
# ============================================================

def fibonacci_sphere(n):
    """Generate n approximately uniform directions on unit sphere."""
    golden = (1 + np.sqrt(5)) / 2
    indices = np.arange(n)
    theta = 2 * np.pi * indices / golden
    phi = np.arccos(1 - 2 * (indices + 0.5) / n)
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi) * np.cos(theta)  # intentional: gives better coverage
    # Actually use proper formula
    dirs = np.column_stack([
        np.sin(phi) * np.cos(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(phi),
    ])
    return dirs


def compute_ray_accessibility(site_com, blocker_pos, ray_dirs, probe_radius, max_dist):
    """
    For each ray from site_com, check if any blocker bead falls within
    a cylinder of radius probe_radius along the ray, up to max_dist.

    Fully vectorized: (N_rays, M_blockers) operations.

    Returns: fraction of unblocked rays (0-1).
    """
    n_rays = len(ray_dirs)
    if len(blocker_pos) == 0:
        return 1.0

    vecs = blocker_pos - site_com  # (M, 3)
    dists_sq = np.sum(vecs**2, axis=1)

    cutoff_sq = (max_dist + probe_radius) ** 2
    mask = dists_sq < cutoff_sq
    if not np.any(mask):
        return 1.0
    vecs = vecs[mask]  # (M, 3)

    # proj[i,j] = dot(ray_i, vec_j) = projection of blocker j onto ray i
    proj = ray_dirs @ vecs.T  # (N_rays, M)

    # Ahead mask: 0 < proj < max_dist
    ahead = (proj > 0) & (proj < max_dist)

    # Perpendicular distance squared: |vec|^2 - proj^2
    vecs_sq = np.sum(vecs**2, axis=1)  # (M,)
    perp_sq = vecs_sq[None, :] - proj**2  # (N_rays, M)

    # A ray is blocked if any ahead blocker has perp_sq < probe_radius^2
    blocked = np.any(ahead & (perp_sq < probe_radius**2), axis=1)  # (N_rays,)

    return 1.0 - np.sum(blocked) / n_rays


def compute_max_cone_angle(site_com, blocker_pos, max_dist):
    """
    Find the widest unobstructed cone from site_com.
    Returns the half-angle (degrees) of the largest empty cone.

    Method: for each blocker within max_dist, compute its angular position.
    The max cone angle = largest angular gap between consecutive blockers
    on the unit sphere. Approximated by finding the direction with the
    largest minimum angle to any blocker.
    """
    vecs = blocker_pos - site_com
    dists = np.linalg.norm(vecs, axis=1)
    mask = (dists > 0) & (dists < max_dist)
    if not np.any(mask):
        return 180.0

    vecs = vecs[mask]
    dists = dists[mask]
    # Normalize to unit vectors
    unit_vecs = vecs / dists[:, None]

    # Sample many test directions and find the one with largest min angle to blockers
    test_dirs = fibonacci_sphere(500)
    # Cosine of angle between each test dir and each blocker
    cos_angles = test_dirs @ unit_vecs.T  # (500, M)
    # For each test dir, find the closest blocker (highest cosine = smallest angle)
    max_cos_per_dir = np.max(cos_angles, axis=1)  # (500,)
    # The best direction is the one where the closest blocker is farthest away
    best_idx = np.argmin(max_cos_per_dir)
    best_max_cos = max_cos_per_dir[best_idx]
    # Convert to angle
    half_angle = np.degrees(np.arccos(np.clip(best_max_cos, -1, 1)))
    return half_angle


def compute_shell_density(site_com, blocker_pos, r_inner=2.0, r_outer=5.0):
    """
    Fraction of a spherical shell (r_inner to r_outer nm) around site_com
    that is occupied by protein beads. Each bead occupies a ~0.38 nm radius
    sphere (CA bead). Returns occupancy fraction.
    """
    vecs = blocker_pos - site_com
    dists = np.linalg.norm(vecs, axis=1)
    n_in_shell = np.sum((dists >= r_inner) & (dists < r_outer))

    # Shell volume
    shell_vol = (4/3) * np.pi * (r_outer**3 - r_inner**3)
    # Each CG bead ~ sphere of radius 0.19 nm -> volume ~0.029 nm^3
    bead_vol = (4/3) * np.pi * 0.19**3
    occupied_vol = n_in_shell * bead_vol

    return occupied_vol / shell_vol


# ============================================================
# Per-frame analysis
# ============================================================

def analyze_site_frame(site_com, pocket_resids, all_pos, all_resids,
                       ray_dirs, probe_radius, max_dist, seq_sep):
    """Compute all accessibility metrics for one site in one frame."""
    # Exclude pocket residues and sequence neighbors from blockers
    pocket_set = set(pocket_resids)
    blocker_mask = np.array([
        r not in pocket_set and all(abs(r - pr) > seq_sep for pr in pocket_resids)
        for r in all_resids
    ])
    blocker_pos = all_pos[blocker_mask]

    saa = compute_ray_accessibility(site_com, blocker_pos, ray_dirs,
                                    probe_radius, max_dist)
    cone = compute_max_cone_angle(site_com, blocker_pos, max_dist)
    density = compute_shell_density(site_com, blocker_pos, r_inner=2.0, r_outer=5.0)

    return saa, cone, density


# ============================================================
# Main analysis
# ============================================================

def run_accessibility(active_sets, n_rays):
    print("\n" + "=" * 70)
    print("ACTIVE SITE STERIC ACCESSIBILITY ANALYSIS")
    print("=" * 70)
    print(f"  Probe radius:  {PROBE_RADIUS} nm (half-width of approaching domain)")
    print(f"  Max distance:  {MAX_DIST} nm")
    print(f"  N rays:        {n_rays}")
    print(f"  Seq separation: ±{SEQ_SEP} residues excluded from blockers")

    ray_dirs = fibonacci_sphere(n_rays)

    results = {}

    for set_key in active_sets:
        info = SETS[set_key]
        sites = get_active_sites_for_set(set_key)
        available_sites = CONSTRUCT_SITES[set_key]

        print(f"\n--- {info['label']} ---")

        # Per-replicate mean values
        rep_saa = {s: [] for s in available_sites}
        rep_cone = {s: [] for s in available_sites}
        rep_density = {s: [] for s in available_sites}

        n_analyzed = 0
        for seed in SEEDS:
            for sample in SAMPLES:
                sysname = info['sysname']
                # Handle fragment sets (prefix frag_), --sim-folder dirs, and nested layout
                if set_key in EXTERNAL_DIRS:
                    sim_dir = os.path.join(EXTERNAL_DIRS[set_key],
                                            f'seed-{seed}_sample-{sample}')
                elif set_key.startswith('frag_'):
                    frag_name = set_key[5:]
                    sim_dir = os.path.join(CWD, 'fragments', frag_name,
                                            f'seed-{seed}_sample-{sample}')
                else:
                    sim_dir = os.path.join(CWD, set_key,
                                            f'seed-{seed}_sample-{sample}')
                pdb = os.path.join(sim_dir, 'top.pdb')
                dcd = os.path.join(sim_dir, f'{sysname}.dcd')

                if not os.path.isfile(pdb) or not os.path.isfile(dcd):
                    continue

                u = mda.Universe(pdb, dcd)
                all_ag = u.select_atoms('all')
                all_resids = all_ag.resids

                # Precompute blocker masks and pocket atom groups (resid-based, frame-independent)
                site_info_cache = {}
                for sname in available_sites:
                    pocket_resids = sites[sname]['pocket']
                    pocket_sel = ' or '.join([f'resid {r}' for r in pocket_resids])
                    pocket_ag = u.select_atoms(pocket_sel)
                    if len(pocket_ag) == 0:
                        continue
                    pocket_set = set(pocket_resids)
                    blocker_mask = np.array([
                        r not in pocket_set and all(abs(r - pr) > SEQ_SEP for pr in pocket_resids)
                        for r in all_resids
                    ])
                    # Indices into all_ag for pocket atoms
                    pocket_idx = np.array([i for i, r in enumerate(all_resids) if r in pocket_set])
                    site_info_cache[sname] = (pocket_idx, blocker_mask)

                # Per-site accumulators for this replicate
                site_saa = {s: [] for s in available_sites}
                site_cone = {s: [] for s in available_sites}
                site_density = {s: [] for s in available_sites}

                for t, ts in enumerate(u.trajectory):
                    if t < SKIP_FRAMES:
                        continue

                    all_pos = all_ag.positions / 10.0  # nm

                    for sname in available_sites:
                        if sname not in site_info_cache:
                            continue
                        pocket_idx, blocker_mask = site_info_cache[sname]

                        site_com = np.mean(all_pos[pocket_idx], axis=0)
                        blocker_pos = all_pos[blocker_mask]

                        saa = compute_ray_accessibility(site_com, blocker_pos, ray_dirs,
                                                        PROBE_RADIUS, MAX_DIST)
                        cone = compute_max_cone_angle(site_com, blocker_pos, MAX_DIST)
                        density = compute_shell_density(site_com, blocker_pos, 2.0, 5.0)

                        site_saa[sname].append(saa)
                        site_cone[sname].append(cone)
                        site_density[sname].append(density)

                for sname in available_sites:
                    if site_saa[sname]:
                        rep_saa[sname].append(np.mean(site_saa[sname]))
                        rep_cone[sname].append(np.mean(site_cone[sname]))
                        rep_density[sname].append(np.mean(site_density[sname]))

                n_analyzed += 1
                if n_analyzed % 5 == 0:
                    print(f"    {n_analyzed}/25 replicates done")

        results[set_key] = {
            'n_analyzed': n_analyzed,
            'saa': rep_saa,
            'cone': rep_cone,
            'density': rep_density,
            'available_sites': available_sites,
        }

        print(f"  Analyzed: {n_analyzed} replicates")
        for sname in available_sites:
            saa_vals = rep_saa[sname]
            cone_vals = rep_cone[sname]
            dens_vals = rep_density[sname]
            if saa_vals:
                print(f"  {sname}: SAA={np.mean(saa_vals):.3f}±{np.std(saa_vals):.3f}"
                      f"  Cone={np.mean(cone_vals):.1f}±{np.std(cone_vals):.1f}°"
                      f"  ShellDens={np.mean(dens_vals):.4f}±{np.std(dens_vals):.4f}")

    # --- Save ---
    save_dict = {}
    for set_key, res in results.items():
        for sname in res['available_sites']:
            if res['saa'][sname]:
                save_dict[f'{set_key}_{sname}_saa'] = np.array(res['saa'][sname])
                save_dict[f'{set_key}_{sname}_cone'] = np.array(res['cone'][sname])
                save_dict[f'{set_key}_{sname}_density'] = np.array(res['density'][sname])
    np.savez(os.path.join(DATA_PATH, 'accessibility_stats.npz'), **save_dict)
    print(f"\n  Saved: accessibility_stats.npz")

    # ============================================================
    # Plots
    # ============================================================

    # --- Plot 1: SAA comparison (main result) ---
    fig, axes = plt.subplots(1, len(SITE_NAMES), figsize=(4 * len(SITE_NAMES), 5), squeeze=False)
    axes = axes[0]
    for si, sname in enumerate(SITE_NAMES):
        ax = axes[si]
        set_data = []
        set_labels = []
        set_colors = []
        for set_key in active_sets:
            res = results[set_key]
            if sname in res['available_sites'] and res['saa'][sname]:
                set_data.append(res['saa'][sname])
                set_labels.append(SETS[set_key]['label'].split('(')[0].strip())
                set_colors.append(SITE_COLORS.get(sname, '#333333'))

        if set_data:
            bp = ax.boxplot(set_data, patch_artist=True, widths=0.6)
            for i, (patch, c) in enumerate(zip(bp['boxes'], set_colors)):
                patch.set_facecolor(c)
                patch.set_alpha(0.5)
            for i, (d, c) in enumerate(zip(set_data, set_colors)):
                jitter = np.random.default_rng(42).normal(0, 0.06, len(d))
                ax.scatter(i + 1 + jitter, d, c=c, s=20, alpha=0.7,
                           edgecolors='black', linewidths=0.3)
            ax.set_xticklabels(set_labels, fontsize=8, rotation=30, ha='right')

        ax.set_ylim(0, 1)
        ax.set_title(f'{sname} Active Site', fontsize=11)
        ax.set_ylabel('Solid-Angle Accessibility')

    fig.suptitle(f'Steric Accessibility to Active Sites\n'
                 f'(fraction of approach directions unblocked by own chain, '
                 f'probe={PROBE_RADIUS} nm, range={MAX_DIST} nm)',
                 fontsize=12, y=1.06)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'accessibility_saa.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: accessibility_saa.png")

    # --- Plot 2: Max cone angle ---
    fig, axes = plt.subplots(1, len(SITE_NAMES), figsize=(4 * len(SITE_NAMES), 5), squeeze=False)
    axes = axes[0]
    for si, sname in enumerate(SITE_NAMES):
        ax = axes[si]
        set_data = []
        set_labels = []
        set_colors = []
        for set_key in active_sets:
            res = results[set_key]
            if sname in res['available_sites'] and res['cone'][sname]:
                set_data.append(res['cone'][sname])
                set_labels.append(SETS[set_key]['label'].split('(')[0].strip())
                set_colors.append(SITE_COLORS.get(sname, '#333333'))

        if set_data:
            bp = ax.boxplot(set_data, patch_artist=True, widths=0.6)
            for i, (patch, c) in enumerate(zip(bp['boxes'], set_colors)):
                patch.set_facecolor(c)
                patch.set_alpha(0.5)
            for i, (d, c) in enumerate(zip(set_data, set_colors)):
                jitter = np.random.default_rng(42).normal(0, 0.06, len(d))
                ax.scatter(i + 1 + jitter, d, c=c, s=20, alpha=0.7,
                           edgecolors='black', linewidths=0.3)
            ax.set_xticklabels(set_labels, fontsize=8, rotation=30, ha='right')

        ax.set_ylim(0, 180)
        ax.axhline(90, color='gray', ls=':', lw=1, label='hemisphere')
        ax.set_title(f'{sname} Active Site', fontsize=11)
        ax.set_ylabel('Max approach cone (degrees)')
        ax.legend(fontsize=7)

    fig.suptitle('Widest Unobstructed Approach Cone\n'
                 '(half-angle of largest empty cone from pocket COM)',
                 fontsize=12, y=1.06)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'accessibility_cone.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: accessibility_cone.png")

    # --- Plot 3: Summary heatmap ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    for ax_idx, (metric, metric_label, cmap, vrange) in enumerate([
        ('saa', 'Solid-Angle Accessibility\n(0=buried, 1=open)', 'RdYlGn', (0, 1)),
        ('cone', 'Max Cone Angle (°)\n(wider=more accessible)', 'RdYlGn', (0, 180)),
        ('density', 'Shell Density (2-5 nm)\n(lower=more accessible)', 'RdYlGn_r', (0, None)),
    ]):
        ax = axes[ax_idx]
        matrix = np.full((len(active_sets), len(SITE_NAMES)), np.nan)
        for si, set_key in enumerate(active_sets):
            res = results[set_key]
            for sj, sname in enumerate(SITE_NAMES):
                if sname in res['available_sites'] and res[metric][sname]:
                    matrix[si, sj] = np.mean(res[metric][sname])

        vmin, vmax = vrange
        if vmax is None:
            vmax = np.nanmax(matrix) * 1.1
        im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
        ax.set_yticks(range(len(active_sets)))
        ax.set_yticklabels([SETS[s]['label'] for s in active_sets], fontsize=9)
        ax.set_xticks(range(len(SITE_NAMES)))
        ax.set_xticklabels(SITE_NAMES, fontsize=11)
        for i in range(len(active_sets)):
            for j in range(len(SITE_NAMES)):
                v = matrix[i, j]
                fmt = f'{v:.3f}' if metric == 'density' else (f'{v:.2f}' if metric == 'saa' else f'{v:.0f}°')
                if not np.isnan(v):
                    ax.text(j, i, fmt, ha='center', va='center', fontsize=10,
                            fontweight='bold')
                else:
                    ax.text(j, i, '—', ha='center', va='center', fontsize=10, color='gray')
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title(metric_label, fontsize=11)

    fig.suptitle('Active Site Accessibility Summary', fontsize=14, y=1.04)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'accessibility_summary.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: accessibility_summary.png")

    # --- Final summary table ---
    print(f"\n{'='*90}")
    print("ACCESSIBILITY SUMMARY")
    print(f"{'='*90}")
    print(f"{'Set':<20} {'Site':<6} {'SAA (mean±std)':<20} {'Cone° (mean±std)':<20} {'ShellDens (mean±std)':<22}")
    print("-" * 88)
    for set_key in active_sets:
        res = results[set_key]
        for sname in SITE_NAMES:
            if sname not in res['available_sites']:
                continue
            saa = res['saa'][sname]
            cone = res['cone'][sname]
            dens = res['density'][sname]
            s_str = f"{np.mean(saa):.3f}±{np.std(saa):.3f}" if saa else "—"
            c_str = f"{np.mean(cone):.1f}±{np.std(cone):.1f}" if cone else "—"
            d_str = f"{np.mean(dens):.4f}±{np.std(dens):.4f}" if dens else "—"
            print(f"{SETS[set_key]['label']:<20} {sname:<6} {s_str:<20} {c_str:<20} {d_str:<22}")
    print(f"{'='*90}")
    print()
    print("Interpretation:")
    print("  SAA > 0.5  = majority of approach directions open — site readily accessible")
    print("  SAA < 0.3  = most directions blocked — site sterically occluded")
    print("  Cone > 90° = at least a hemisphere is free — easy docking")
    print("  Cone < 45° = only a narrow approach window — docking constrained")


# ============================================================
# CLI
# ============================================================

def discover_fragments():
    """Auto-discover fragment simulations and register them in SETS."""
    import json as _json
    frag_dir = os.path.join(CWD, 'fragments')
    if not os.path.isdir(frag_dir):
        return []
    discovered = []
    for name in sorted(os.listdir(frag_dir)):
        fdir = os.path.join(frag_dir, name)
        meta_file = os.path.join(fdir, 'metadata.json')
        if not os.path.isfile(meta_file):
            continue
        # Check at least one replicate has a trajectory
        has_traj = any(
            os.path.isfile(os.path.join(fdir, entry, 'parp14.dcd'))
            for entry in os.listdir(fdir) if entry.startswith('seed-'))
        if not has_traj:
            continue
        with open(meta_file) as f:
            meta = _json.load(f)
        set_key = f'frag_{name}'
        SETS[set_key] = {'sysname': 'parp14', 'label': name}
        # Active sites present in this fragment
        unit_to_site = {'md1l1': 'MD1', 'md2': 'MD2', 'md3': 'MD3', 'art': 'ART'}
        sites = [unit_to_site[u] for u in meta.get('units', [])
                 if u in unit_to_site]
        CONSTRUCT_SITES[set_key] = sites
        discovered.append(set_key)
    return discovered


def main():
    parser = ArgumentParser(description='Active site steric accessibility analysis')
    parser.add_argument('--set', nargs='+', default=None,
                        help="Sets to analyze (e.g. fl, fl_optimized, "
                             "md1l1_md2, or 'all')")
    parser.add_argument('--include-fragments', action='store_true',
                        help='Auto-discover and include fragments/')
    parser.add_argument('--sim-folder', nargs='+', default=None, metavar='PATH',
                        help='One or more simulation-set folders to analyze '
                             '(each containing seed-*_sample-*/ replicates). '
                             "Domain units are read from the folder's metadata.json, "
                             'or pass --units. Use this for new simulations '
                             'without editing the script.')
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units in the --sim-folder construct '
                             '(e.g. md1l1 md2 md3). Needed only when the folder '
                             'has no metadata.json.')
    parser.add_argument('--nrays', type=int, default=N_RAYS,
                        help=f'Number of probe rays (default: {N_RAYS})')
    args = parser.parse_args()

    # Register any folders passed via --sim-folder before running analysis.
    folder_keys = [register_sim_folder(f, units=args.units)
                   for f in (args.sim_folder or [])]
    for k in folder_keys:
        print(f"  Registered sim folder -> set '{k}' "
              f"(units: {', '.join(CONSTRUCT_UNITS[k])})")

    if args.include_fragments or (args.set and 'all' in args.set):
        discovered = discover_fragments()
        print(f"  Discovered {len(discovered)} fragment simulations")

    if args.set:
        if 'all' in args.set:
            active_sets = list(SETS.keys())
        else:
            active_sets = []
            for s in args.set:
                if s in SETS:
                    active_sets.append(s)
                elif f'frag_{s}' in SETS:
                    active_sets.append(f'frag_{s}')
                else:
                    print(f"  WARN: unknown set '{s}'")
        # --sim-folder always adds its folders on top of --set selection
        active_sets += [k for k in folder_keys if k not in active_sets]
    elif folder_keys:
        active_sets = folder_keys
    else:
        active_sets = ['fl', 'md', 'core', 'mka']

    run_accessibility(active_sets, args.nrays)


if __name__ == '__main__':
    main()
