#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inter-domain distance violin plot for fl_optimized PARP14 simulations.

For each frame in each replicate, compute:
  - MD1L1 - MD3 center-of-mass distance
  - MD1L1 - ART center-of-mass distance
  - (Optional) minimum CA-CA distance between domain pairs

Then plot distributions as violins. The "frequency" axis comes from the
density of frames at each distance.

This is a focused alternative to the grid contact-frequency heatmap when you
want to see the full distribution shape (bimodal, narrow, etc.) for specific
inter-domain pairs.

Usage:
    python figure_md_distances.py
    python figure_md_distances.py --set fl fl_optimized   # compare multiple sets
    python figure_md_distances.py --metric min            # use min CA-CA instead of COM
"""

import os
import argparse
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent.parent
DATA_PATH = CWD / 'data'
DATA_PATH.mkdir(exist_ok=True)

# Arbitrary simulation folders registered via --sim-folder.
# Populated in the main process before any worker pool is created so that
# forked workers (Linux fork start method) inherit them.
EXTERNAL_DIRS = {}      # set_key -> absolute Path to the sim folder
EXTERNAL_SYSNAME = {}   # set_key -> dcd basename (sysname)
CONSTRUCT_UNITS = {}    # set_key -> FL domain units (drives FL->construct remap)

# Dated category subdir: figures/04_md_distances/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, str(CWD))
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_PATH = _get_fig_dir('04_md_distances')

# Domain definitions in FL numbering
DOMAINS_FL = {
    'MD1L1': (790, 1004),
    'MD2':   (1005, 1193),
    'MD3':   (1207, 1388),
    'ART':   (1603, 1801),
}

# Pairs to plot
PAIRS = [
    ('MD1L1', 'MD3'),
    ('MD1L1', 'ART'),
    # Additional pairs to add to violin row for context:
    ('MD1L1', 'MD2'),
    ('MD2',   'MD3'),
    ('MD3',   'ART'),
]

# Contact distance threshold (from cmap_traj soft tanh cutoff)
CONTACT_CUTOFF_NM = 1.0
SKIP_FRAMES = 50  # equilibration

SEEDS = range(1, 6)
SAMPLES = range(0, 5)


# ============================================================
# External simulation-folder registration (--sim-folder)
# ============================================================

def register_sim_folder(path, units=None):
    """Register an arbitrary simulation folder (from --sim-folder) as a set so
    this figure can be generated for it without editing the script. Domain units
    are read from the folder's metadata.json, the `units` arg, or (if the
    basename is a known set) that set. Returns the set_key it is registered under.

    Call this in the main process before launching the worker pool so that
    forked workers inherit EXTERNAL_DIRS / EXTERNAL_SYSNAME / CONSTRUCT_UNITS.
    """
    folder = Path(path).resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"sim folder not found: {folder}")
    resolved_units = reg._resolve_units(str(folder), units)
    meta = reg.read_metadata(str(folder)) or {}
    sysname = meta.get('sysname') or reg.detect_sysname(str(folder))
    set_key = folder.name
    EXTERNAL_DIRS[set_key] = folder
    EXTERNAL_SYSNAME[set_key] = sysname
    CONSTRUCT_UNITS[set_key] = resolved_units
    return set_key


# Cache: set_key -> replicate dirs discovered by _flat_replicate_dirs() below.
_FLAT_DIRS_CACHE = {}


def _flat_replicate_dirs(folder, sysname):
    """Every immediate subdirectory of `folder` containing a `<sysname>.dcd`.

    Fallback for --sim-folder layouts that aren't a seed-{1-5}_sample-{0-4}
    grid -- e.g. a handful of independent long runs named arbitrarily
    (fl_go's state-*_tica_seed-*_sample-*_fr*/ dirs)."""
    return [d for d in sorted(folder.iterdir())
            if d.is_dir() and (d / f'{sysname}.dcd').is_file()]


def _resolve_sim_paths(set_key, seed, sample):
    if set_key in EXTERNAL_DIRS:
        base = EXTERNAL_DIRS[set_key]
        sysname = EXTERNAL_SYSNAME[set_key]
        sim_dir = base / f'seed-{seed}_sample-{sample}'
        if not sim_dir.is_dir():
            if set_key not in _FLAT_DIRS_CACHE:
                _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base, sysname)
            flat = _FLAT_DIRS_CACHE[set_key]
            idx = (seed - 1) * len(SAMPLES) + sample
            if idx < len(flat):
                sim_dir = flat[idx]
    elif set_key.startswith('frag_'):
        sysname = 'parp14'
        sim_dir = CWD / 'fragments' / set_key[5:] / f'seed-{seed}_sample-{sample}'
    else:
        sysname = 'parp14'
        sim_dir = CWD / set_key / f'seed-{seed}_sample-{sample}'
    dcd = sim_dir / f'{sysname}.dcd'
    for name in ('top.pdb', 'restart.pdb', 'checkpoint.pdb'):
        pdb = sim_dir / name
        if pdb.is_file():
            return pdb, dcd
    return sim_dir / 'top.pdb', dcd


def domain_ranges_for_set(set_key):
    """Domain residue ranges (construct numbering) for the inter-domain pairs.

    For named/fragment sets the historical FL-numbering ranges (DOMAINS_FL) are
    used unchanged. For folders registered via --sim-folder the FL ranges are
    remapped into the construct's numbering using the resolved domain units, so
    sub-construct trajectories select the correct residues.
    """
    units = CONSTRUCT_UNITS.get(set_key)
    if not units:
        return DOMAINS_FL
    fl_to_c = reg.build_fl_to_construct_map(units)
    mapped = {}
    for dname, (fl_s, fl_e) in DOMAINS_FL.items():
        c_s = fl_to_c(fl_s)
        c_e = fl_to_c(fl_e)
        if c_s is not None and c_e is not None:
            mapped[dname] = (c_s, c_e)
    return mapped


# ============================================================
# Worker: compute distances for one replicate
# ============================================================

def compute_distances_one_replicate(args):
    """Compute per-frame inter-domain distances for one replicate.

    Returns dict {pair_label: array of distances (nm)} or None on failure.
    """
    set_key, seed, sample, metric = args
    import MDAnalysis as mda

    pdb, dcd = _resolve_sim_paths(set_key, seed, sample)
    if not pdb.is_file() or not dcd.is_file():
        return None

    u = mda.Universe(str(pdb), str(dcd))

    # Pre-select domain atom groups
    domains = domain_ranges_for_set(set_key)
    ags = {}
    for dname, (ds, de) in domains.items():
        ag = u.select_atoms(f'resid {ds}:{de}')
        if len(ag) > 0:
            ags[dname] = ag

    if not ags:
        return None

    # Initialize result arrays
    results = {f'{a}_{b}': [] for a, b in PAIRS if a in ags and b in ags}

    for ts in u.trajectory[SKIP_FRAMES:]:
        for a, b in PAIRS:
            key = f'{a}_{b}'
            if key not in results:
                continue
            ag_a = ags[a]
            ag_b = ags[b]

            if metric == 'com':
                # COM-COM distance (nm)
                # MDAnalysis positions are in Angstroms, convert to nm
                com_a = ag_a.center_of_geometry() / 10.0
                com_b = ag_b.center_of_geometry() / 10.0
                d = float(np.linalg.norm(com_a - com_b))
            else:  # metric == 'min'
                # Minimum CA-CA distance (nm)
                pos_a = ag_a.positions / 10.0  # Å → nm
                pos_b = ag_b.positions / 10.0
                diff = pos_a[:, np.newaxis, :] - pos_b[np.newaxis, :, :]
                d = float(np.sqrt(np.sum(diff ** 2, axis=2)).min())

            results[key].append(d)

    return {
        'set_key': set_key,
        'seed': seed,
        'sample': sample,
        'distances': {k: np.array(v) for k, v in results.items()},
    }


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--set', nargs='+', default=['fl_optimized'],
                        help='Sets to include (e.g. fl_optimized, fl, '
                             'or fragment names)')
    parser.add_argument('--metric', choices=['com', 'min'], default='com',
                        help='Distance metric: COM-COM (default) or min CA-CA')
    parser.add_argument('--workers', type=int, default=8,
                        help='Parallel workers (default 8)')
    parser.add_argument('--pairs', nargs='+', default=None,
                        help='Pairs to plot (e.g. MD1L1_MD3 MD1L1_ART). '
                             'Default: all 5 pairs')
    parser.add_argument('--sim-folder', nargs='+', default=None, metavar='PATH',
                        help='One or more NEW simulation-set folders to analyze '
                             '(each containing seed-*_sample-*/ replicates with '
                             'top.pdb + <sysname>.dcd). Domain units are read from '
                             "the folder's metadata.json, or pass --units.")
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units in the --sim-folder construct '
                             '(e.g. md1l1 md2 md3). Required only if the folder '
                             'has no metadata.json.')
    args = parser.parse_args()

    # Register any --sim-folder folders before building jobs / launching workers
    # so forked workers inherit EXTERNAL_DIRS / EXTERNAL_SYSNAME / CONSTRUCT_UNITS.
    external_keys = []
    if args.sim_folder:
        for folder in args.sim_folder:
            k = register_sim_folder(folder, units=args.units)
            external_keys.append(k)
            print(f"Registered --sim-folder '{folder}' as set '{k}' "
                  f"(sysname: {EXTERNAL_SYSNAME[k]}, "
                  f"units: {', '.join(CONSTRUCT_UNITS[k])})")
        # Replace the default set with the folders unless --set was given explicitly.
        if args.set == ['fl_optimized']:
            args.set = external_keys
        else:
            args.set = args.set + [k for k in external_keys if k not in args.set]

    sets = args.set
    metric = args.metric

    # Filter pairs
    if args.pairs:
        wanted = set(args.pairs)
        pair_keys = [f'{a}_{b}' for a, b in PAIRS if f'{a}_{b}' in wanted]
    else:
        pair_keys = [f'{a}_{b}' for a, b in PAIRS]

    print(f"Sets:       {sets}")
    print(f"Metric:     {metric} ({'COM-COM' if metric=='com' else 'min CA-CA'})")
    print(f"Pairs:      {pair_keys}")
    print(f"Workers:    {args.workers}")

    # Build job list
    jobs = []
    for sk in sets:
        for seed in SEEDS:
            for sample in SAMPLES:
                jobs.append((sk, seed, sample, metric))

    # Run in parallel
    print(f"\nProcessing {len(jobs)} replicates...")
    pooled = {sk: {pk: [] for pk in pair_keys} for sk in sets}

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(compute_distances_one_replicate, j): j
                   for j in jobs}
        for fut in as_completed(futures):
            res = fut.result()
            if res is None:
                continue
            sk = res['set_key']
            for k, arr in res['distances'].items():
                if k in pooled[sk]:
                    pooled[sk][k].append(arr)

    # Pool all replicate arrays into one big array per pair per set
    pooled_arrays = {}
    for sk in sets:
        pooled_arrays[sk] = {}
        for pk in pair_keys:
            if pooled[sk][pk]:
                pooled_arrays[sk][pk] = np.concatenate(pooled[sk][pk])
                n_frames = len(pooled_arrays[sk][pk])
                print(f"  {sk} {pk}: {n_frames} frames, "
                      f"mean={pooled_arrays[sk][pk].mean():.2f} nm, "
                      f"std={pooled_arrays[sk][pk].std():.2f}")
            else:
                pooled_arrays[sk][pk] = np.array([])

    # Save raw data
    out_npz = DATA_PATH / f'md_distances_{metric}.npz'
    save_dict = {}
    for sk in sets:
        for pk in pair_keys:
            save_dict[f'{sk}__{pk}'] = pooled_arrays[sk][pk]
    np.savez(out_npz, **save_dict)
    print(f"\n  Saved: {out_npz}")

    # ── Violin plot ──
    plot_violin(pooled_arrays, sets, pair_keys, metric)


# ============================================================
# Plotting
# ============================================================

def plot_violin(pooled_arrays, sets, pair_keys, metric):
    """Massive violin plot: one violin per (set, pair) combination.

    Y-axis: inter-domain distance (nm)
    X-axis: pair, grouped by set with side-by-side comparison
    Violin shape: frequency density across all frames
    Overlay: scatter of all data points (low alpha)
    """
    # Color per set
    set_colors = {
        'fl':            '#1f77b4',
        'fl_optimized':  '#aec7e8',
        'md':            '#ff7f0e',
        'core':          '#2ca02c',
        'mka':           '#d62728',
    }
    default_palette = plt.cm.tab10(np.linspace(0, 1, len(sets)))

    n_pairs = len(pair_keys)
    n_sets = len(sets)
    width_per_pair = 1.5
    fig_width = max(8, n_pairs * width_per_pair * max(n_sets, 1) + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    positions = []
    labels_minor = []
    labels_major = []
    label_positions = []

    pos = 0
    pair_centers = []
    for pk in pair_keys:
        a, b = pk.split('_')
        pair_start = pos
        for i, sk in enumerate(sets):
            arr = pooled_arrays[sk][pk]
            if len(arr) == 0:
                pos += 1
                continue
            color = set_colors.get(sk, default_palette[i])
            parts = ax.violinplot([arr], positions=[pos], widths=0.8,
                                  showmeans=True, showmedians=True,
                                  showextrema=False)
            for body in parts['bodies']:
                body.set_facecolor(color)
                body.set_alpha(0.55)
                body.set_edgecolor('black')
                body.set_linewidth(0.5)
            for key in ('cmeans', 'cmedians'):
                if key in parts:
                    parts[key].set_color('black')
                    parts[key].set_linewidth(1.2)

            # Scatter overlay (subsampled for performance)
            rng = np.random.default_rng(42)
            sub = arr if len(arr) < 5000 else rng.choice(arr, 5000, replace=False)
            jitter = rng.uniform(-0.18, 0.18, len(sub))
            ax.scatter(pos + jitter, sub, s=1.5, c=color, alpha=0.15,
                       edgecolors='none', zorder=2)

            positions.append(pos)
            labels_minor.append(sk if n_sets > 1 else '')
            pos += 1
        # Center label for the pair
        pair_centers.append((pair_start + pos - 1) / 2)
        labels_major.append(f'{a}–{b}')
        pos += 0.8  # gap between pairs

    # Contact cutoff reference
    ax.axhline(CONTACT_CUTOFF_NM, color='red', linestyle='--', alpha=0.6,
               linewidth=1.5,
               label=f'Contact cutoff ({CONTACT_CUTOFF_NM:.1f} nm)')

    # X-axis: minor ticks for sets (if multi-set), major ticks for pairs
    if n_sets > 1:
        ax.set_xticks(positions, minor=True)
        ax.set_xticklabels(labels_minor, minor=True, rotation=45, ha='right',
                           fontsize=7)
        ax.set_xticks(pair_centers, minor=False)
        ax.set_xticklabels(labels_major, minor=False, fontsize=12,
                           fontweight='bold')
        ax.tick_params(axis='x', which='major', pad=22)
    else:
        ax.set_xticks(pair_centers)
        ax.set_xticklabels(labels_major, fontsize=12, fontweight='bold')

    metric_label = ('COM–COM distance' if metric == 'com'
                    else 'min CA–CA distance')
    ax.set_ylabel(f'{metric_label} (nm)', fontsize=13)
    n_frames_per = '\n'.join(
        f'{sk}: {sum(len(pooled_arrays[sk][pk]) for pk in pair_keys)} pooled frames'
        for sk in sets)
    ax.set_title(
        f'Inter-Domain {metric_label} Distributions\n'
        f'{n_frames_per}',
        fontsize=12)

    # Legend with set colors
    legend_elements = []
    for sk in sets:
        color = set_colors.get(sk, '#666')
        legend_elements.append(Line2D([0], [0], marker='s', color='w',
                                       markerfacecolor=color, markersize=12,
                                       label=sk))
    legend_elements.append(Line2D([0], [0], color='red', linestyle='--',
                                   label=f'Contact ({CONTACT_CUTOFF_NM:.1f} nm)'))
    ax.legend(handles=legend_elements, fontsize=10, loc='upper right')

    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    sets_str = '_'.join(sets)
    out_png = FIG_PATH / f'md_distance_violin_{sets_str}_{metric}.png'
    out_svg = FIG_PATH / f'md_distance_violin_{sets_str}_{metric}.svg'
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    fig.savefig(out_svg, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_png}")
    print(f"  Saved: {out_svg}")

    # ── Per-pair detailed histogram + cumulative ──
    plot_distance_histograms(pooled_arrays, sets, pair_keys, metric)


def plot_distance_histograms(pooled_arrays, sets, pair_keys, metric):
    """Per-pair histogram + cumulative distribution."""
    set_colors = {
        'fl':            '#1f77b4',
        'fl_optimized':  '#aec7e8',
        'md':            '#ff7f0e',
        'core':          '#2ca02c',
        'mka':           '#d62728',
    }
    default_palette = plt.cm.tab10(np.linspace(0, 1, len(sets)))

    n_pairs = len(pair_keys)
    fig, axes = plt.subplots(2, n_pairs, figsize=(5 * n_pairs, 7),
                              squeeze=False)

    # Find global x range
    all_dists = []
    for sk in sets:
        for pk in pair_keys:
            all_dists.append(pooled_arrays[sk][pk])
    all_combined = np.concatenate([a for a in all_dists if len(a) > 0])
    xmax = np.percentile(all_combined, 99.5) * 1.05 if len(all_combined) else 10

    for col, pk in enumerate(pair_keys):
        a, b = pk.split('_')
        ax_hist = axes[0, col]
        ax_cdf = axes[1, col]

        for i, sk in enumerate(sets):
            arr = pooled_arrays[sk][pk]
            if len(arr) == 0:
                continue
            color = set_colors.get(sk, default_palette[i])
            ax_hist.hist(arr, bins=80, density=True, alpha=0.5,
                         color=color, label=sk, range=(0, xmax))
            # CDF
            sorted_arr = np.sort(arr)
            cdf = np.arange(1, len(sorted_arr) + 1) / len(sorted_arr)
            ax_cdf.plot(sorted_arr, cdf, color=color, label=sk, linewidth=1.5)

            # Fraction in contact
            frac_contact = np.mean(arr < CONTACT_CUTOFF_NM)
            ax_cdf.axvline(CONTACT_CUTOFF_NM, color='red', linestyle='--',
                           alpha=0.4, linewidth=1)
            ax_hist.axvline(CONTACT_CUTOFF_NM, color='red', linestyle='--',
                            alpha=0.4, linewidth=1)
            ax_cdf.text(0.98, 0.05 + i * 0.06,
                        f'{sk}: {frac_contact*100:.1f}% contact',
                        transform=ax_cdf.transAxes, fontsize=8,
                        ha='right', color=color, fontweight='bold')

        ax_hist.set_title(f'{a}–{b}', fontsize=12, fontweight='bold')
        ax_hist.set_xlabel('Distance (nm)', fontsize=10)
        ax_hist.set_ylabel('Density', fontsize=10)
        ax_hist.legend(fontsize=8, loc='upper right')
        ax_hist.set_xlim(0, xmax)

        ax_cdf.set_xlabel('Distance (nm)', fontsize=10)
        ax_cdf.set_ylabel('Cumulative fraction', fontsize=10)
        ax_cdf.set_xlim(0, xmax)
        ax_cdf.set_ylim(0, 1.02)
        ax_cdf.grid(alpha=0.3)

    metric_label = ('COM–COM' if metric == 'com'
                    else 'min CA–CA')
    fig.suptitle(f'Inter-Domain {metric_label} Distance Distributions',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    sets_str = '_'.join(sets)
    out_png = FIG_PATH / f'md_distance_hist_{sets_str}_{metric}.png'
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_png}")


if __name__ == '__main__':
    main()
