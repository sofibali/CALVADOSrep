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

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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

# Cache: set_key -> replicate dirs discovered by _flat_replicate_dirs() below.
_FLAT_DIRS_CACHE = {}


def _flat_replicate_dirs(folder, sysname):
    """Every immediate subdirectory of `folder` containing a `<sysname>.dcd`.

    Fallback for --sim-folder layouts that aren't a seed-{1-5}_sample-{0-4}
    grid -- e.g. a handful of independent long runs named arbitrarily
    (fl_go's state-*_tica_seed-*_sample-*_fr*/ dirs)."""
    dirs = []
    for entry in sorted(os.listdir(folder)):
        d = os.path.join(folder, entry)
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, f'{sysname}.dcd')):
            dirs.append(d)
    return dirs


def _resolve_sim_paths(set_key, sysname, seed, sample):
    """(pdb, dcd) for one replicate, handling fragments, --sim-folder dirs
    (with a flat-directory fallback for non-grid layouts), and named sets."""
    if set_key in EXTERNAL_DIRS:
        base = EXTERNAL_DIRS[set_key]
        sim_dir = os.path.join(base, f'seed-{seed}_sample-{sample}')
        if not os.path.isdir(sim_dir):
            if set_key not in _FLAT_DIRS_CACHE:
                _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base, sysname)
            flat = _FLAT_DIRS_CACHE[set_key]
            idx = (seed - 1) * len(SAMPLES) + sample
            if idx < len(flat):
                sim_dir = flat[idx]
    elif set_key.startswith('frag_'):
        sim_dir = os.path.join(CWD, 'fragments', set_key[5:],
                                f'seed-{seed}_sample-{sample}')
    else:
        sim_dir = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
    dcd = os.path.join(sim_dir, f'{sysname}.dcd')
    for name in ('top.pdb', 'restart.pdb', 'checkpoint.pdb'):
        pdb = os.path.join(sim_dir, name)
        if os.path.isfile(pdb):
            return pdb, dcd
    return os.path.join(sim_dir, 'top.pdb'), dcd

# Accessibility parameters
PROBE_RADIUS = 0.5    # nm — CG bead radius; tests if approach path is sterically clear
# How far along each ray obstruction is tested. A SELECTED cutoff, not a derived
# one -- chosen at 3.0 nm to match the DSS CA-CA crosslink ceiling used in
# analyze_lys_contacts.py, so "reachable" means the same distance in both
# analyses. A sweep over 1-10 nm (sweep_maxdist.sh) shows the site ordering
# MD1 < MD2 < ART < MD3 < WWE is identical at every distance from 2 nm up, and
# that SAA is within a few percent of its asymptote by ~5 nm, so no conclusion
# here depends on the choice. Below 2 nm the metric saturates (MD3 and WWE both
# hit 1.000) and is unusable.
MAX_DIST = 3.0        # nm — see above; matches the DSS crosslink ceiling
N_RAYS = 200          # rays per site per frame (Fibonacci sphere)
SEQ_SEP = 10          # ignore beads within ±10 residues of pocket (own fold)

# Set definitions -- derived from sim_registry.SETS (the shared source of
# truth) rather than a hardcoded literal, which had gone stale and silently
# KeyError'd on any set added there since (norrm/noart/md3art/md_full/mka_full).
# Still a plain mutable dict: --sim-folder registration further down adds to
# this SAME object at runtime (see register_external_folder-equivalent below).
SETS = {k: {'sysname': v['sysname'], 'label': v['label']} for k, v in reg.SETS.items()}

# Active sites (FL numbering)
# Was a stale hardcoded duplicate of sim_registry.ACTIVE_SITES_FL (missing
# the later-added WWE site until this file was pointed at reg's copy
# directly) -- same stale-duplicate-dict bug class already fixed elsewhere.
ACTIVE_SITES_FL = reg.ACTIVE_SITES_FL
SITE_NAMES = reg.SITE_NAMES
SITE_COLORS = reg.SITE_COLORS

# Domain unit mapping (for construct numbering)
DOMAIN_UNITS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789), 'md1l1': (790, 1004),
    'md2': (1004, 1193), 'md3': (1207, 1388), 'khb-kh8': (1389, 1533),
    'wwe': (1534, 1602), 'art': (1603, 1801),
}

# Derived from sim_registry.SETS rather than hardcoded -- the literal version
# here had gone stale (missing norrm/noart/md3art/md_full/mka_full, which then
# KeyError'd at every use below).
CONSTRUCT_UNITS = {k: v['units'] for k, v in reg.SETS.items()}
CONSTRUCT_SITES = {k: [reg.UNIT_TO_SITE[u] for u in units if u in reg.UNIT_TO_SITE]
                    for k, units in CONSTRUCT_UNITS.items()}

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


def build_fl_to_construct_map(unit_names, set_key=None):
    # Full-length (all 11 units): identity map -- see sim_registry.py's
    # build_fl_to_construct_map for the full rationale. Without this,
    # 'fl_optimized' (all 11 units in CONSTRUCT_UNITS, not special-cased by
    # the literal 'fl' check below) silently mis-maps MD3/KHb-KH8/WWE/ART
    # active sites by the un-excised 13-residue MD2-MD3 linker gap.
    if set(unit_names) == set(DOMAIN_UNITS.keys()):
        return lambda fl_resid: fl_resid if 1 <= fl_resid <= 1801 else None
    # md_full/mka_full: same gap, but a genuine sub-range (not all 11 units)
    # that was deliberately built to keep every linker within its span -- see
    # sim_registry.CONTIGUOUS_FL_RANGE for why this can't be inferred from
    # unit_names alone.
    if set_key in reg.CONTIGUOUS_FL_RANGE:
        fl_start, fl_end = reg.CONTIGUOUS_FL_RANGE[set_key]
        offset = 1 - fl_start
        return lambda fl_resid: (fl_resid + offset
                                  if fl_start <= fl_resid <= fl_end else None)
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
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[set_key], set_key=set_key)
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

def run_accessibility(active_sets, n_rays, stride=1, target_frames=None, tag=''):
    print("\n" + "=" * 70)
    print("ACTIVE SITE STERIC ACCESSIBILITY ANALYSIS")
    print("=" * 70)
    print(f"  Probe radius:  {PROBE_RADIUS} nm (half-width of approaching domain)")
    print(f"  Max distance:  {MAX_DIST} nm")
    print(f"  N rays:        {n_rays}")
    print(f"  Seq separation: ±{SEQ_SEP} residues excluded from blockers")
    if target_frames:
        print(f"  Frame budget:  ~{target_frames} frames/replicate "
              f"(stride chosen per trajectory)")
    else:
        print(f"  Frame stride:  {stride}"
              + ("" if stride == 1 else
                 f" (2 us / 100k-frame runs -> 1 frame per {stride * 0.02:.2f} ns)"))

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
                pdb, dcd = _resolve_sim_paths(set_key, sysname, seed, sample)

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

                # Subsample post-equilibration frames. The 2 us extension runs
                # hold 100k frames at 20 ps/frame, but accessibility follows
                # domain reorientation, which decorrelates on tens of ns --
                # consecutive frames are near-identical, so analyzing all of
                # them costs ~50x more for no extra independent sampling.
                #
                # With --target-frames the stride is picked PER TRAJECTORY: the
                # sets span 3-4k-frame originals and 100k-frame extensions, and
                # one fixed stride cannot serve both (stride 50 leaves a
                # 3514-frame run with only ~69 frames).
                rep_stride = stride
                if target_frames:
                    usable = max(0, len(u.trajectory) - SKIP_FRAMES)
                    rep_stride = max(1, usable // target_frames)

                for t, ts in enumerate(u.trajectory):
                    if t < SKIP_FRAMES:
                        continue
                    if rep_stride > 1 and (t - SKIP_FRAMES) % rep_stride:
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
    # MERGE with any previously cached sets rather than overwriting the file.
    # This run only holds `active_sets`; a plain np.savez here wiped every set
    # computed by earlier invocations, so an incremental backfill (one set per
    # run -- the only tractable way to cover 200+ replicates) destroyed its own
    # results and left the npz holding just the set that happened to run last.
    npz_path = os.path.join(
        DATA_PATH,
        'accessibility_stats.npz' if not tag
        else f'accessibility_sweep{tag}.npz')
    save_dict = {}
    if os.path.isfile(npz_path):
        with np.load(npz_path, allow_pickle=True) as _old:
            save_dict = {k: _old[k] for k in _old.files}
        n_prev = len({k.rsplit('_', 2)[0] for k in save_dict})
        print(f"\n  Merging into {n_prev} previously cached set(s) in accessibility_stats.npz")
    for set_key, res in results.items():
        for sname in res['available_sites']:
            if res['saa'][sname]:
                save_dict[f'{set_key}_{sname}_saa'] = np.array(res['saa'][sname])
                save_dict[f'{set_key}_{sname}_cone'] = np.array(res['cone'][sname])
                save_dict[f'{set_key}_{sname}_density'] = np.array(res['density'][sname])
    np.savez(npz_path, **save_dict)
    print(f"  Saved: accessibility_stats.npz "
          f"({len({k.rsplit('_', 2)[0] for k in save_dict})} sets total)")

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
    # Hoisted: MAX_DIST/PROBE_RADIUS are read below as argparse defaults and
    # rebound after parsing, and Python requires the declaration first.
    global MAX_DIST, PROBE_RADIUS, FIG_PATH
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
    parser.add_argument('--max-dist', type=float, default=MAX_DIST, metavar='NM',
                        help=f'How far along each ray to test for obstruction '
                             f'(default {MAX_DIST} nm). This sets what counts as '
                             f'"blocking": it is the size of the approaching '
                             f'species you care about. Sweep it to check whether '
                             f'a conclusion depends on the choice.')
    parser.add_argument('--probe-radius', type=float, default=PROBE_RADIUS,
                        metavar='NM',
                        help=f'Radius of the cylinder swept along each ray '
                             f'(default {PROBE_RADIUS} nm).')
    parser.add_argument('--tag', default='', metavar='STR',
                        help='Write to data/accessibility_sweep<TAG>.npz instead '
                             'of the production cache. Use for parameter sweeps: '
                             'tagged KEYS in the shared file would be read '
                             'downstream as extra constructs.')
    parser.add_argument('--stride', type=int, default=1,
                        help='Analyze every Nth post-equilibration frame '
                             '(default: 1 = every frame). The 2 us runs hold '
                             '100k frames at 20 ps each; --stride 50 gives '
                             '1 ns spacing (~2000 frames/replicate), which is '
                             'still far finer than the ~10-100 ns decorrelation '
                             'time of domain motion, at ~50x lower cost.')
    parser.add_argument('--target-frames', type=int, default=None, metavar='N',
                        help='Per-replicate frame budget; the stride is chosen '
                             'PER TRAJECTORY as (n_frames - skip) // N. Prefer '
                             'this over --stride when a run spans both the '
                             '100k-frame extensions and the 3-4k-frame original '
                             'sets, where one fixed stride would either waste '
                             'compute on the long runs or leave the short ones '
                             'with only a few dozen frames. Overrides --stride.')
    args = parser.parse_args()
    if args.stride < 1:
        parser.error('--stride must be >= 1')
    if args.target_frames is not None and args.target_frames < 1:
        parser.error('--target-frames must be >= 1')

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

    FIG_PATH = str(_get_fig_dir('03_accessibility', sims=active_sets))

    MAX_DIST = args.max_dist
    PROBE_RADIUS = args.probe_radius
    run_accessibility(active_sets, args.nrays, stride=args.stride,
                      target_frames=args.target_frames, tag=args.tag)


if __name__ == '__main__':
    main()
