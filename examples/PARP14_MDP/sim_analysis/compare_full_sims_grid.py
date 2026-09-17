#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Domain x domain grid comparison of inter-domain COM-COM distances across the
four "full/contiguous" PARP14 constructs, matched to the same simulated
length (1000 ns/replicate):

  md_full        MD1L1-MD2-MD3 contiguous (790-1388, 599 res), 5 reps x 1000 ns
  mka_full       MD1L1-ART contiguous (790-1801, 1012 res), 5 reps x 1000 ns
  core_full_go   KH7a-ART contiguous (738-1801, 1064 res), 5 reps x 1000 ns
                 + Go-model KH7a-KHb custom restraints
  fl_go_5rep1us  Full length (1-1801, 1801 res), FIRST 5 of fl_go's 11
                 replicates, truncated to the first 1000 ns of each 2000 ns
                 run (matched length/replicate-count to the other three) +
                 the same Go-model KH7a-KHb custom restraints

Two domain x domain grids side by side (histograms | cumulative fraction),
domains MD1L1/MD2/MD3/ART on both axes: cell (row=i, col=j) for i<j holds the
domain_i-domain_j distribution (upper triangle only); diagonal cells show the
domain name; lower triangle is blank. Reads the SAME data/md_distances_com.npz
cache figure_md_distances.py already writes (pairs + per-domain mean Rg) --
run that script first with all 4 sets so the cache has everything:

    python figure_md_distances.py --set md_full mka_full core_full_go \\
        --sim-folder-as fl_go_5rep1us:/path/to/fl_go \\
        --max-reps fl_go_5rep1us:5 --max-ns fl_go_5rep1us:1000

Usage:
    python compare_full_sims_grid.py
    python compare_full_sims_grid.py --metric min
"""
import argparse
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sys as _sys
CWD = Path(__file__).resolve().parent.parent
_sys.path.insert(0, str(CWD))
from _fig_layout import get_fig_dir as _get_fig_dir

from figure_md_distances import detect_bimodality, AH_CUTOFF_NM

DATA_PATH = CWD / 'data'

DOMAINS = ['MD1L1', 'MD2', 'MD3', 'ART']

SETS = ['md_full', 'mka_full', 'core_full_go', 'fl_go_5rep1us']

# All 13 constructs from the "full/contiguous, matched to 1000 ns" family
# (the 4 above + the 9 new ones launched this session), for --sets ALL13-style
# invocations. Colors match sim_registry.SETS for consistency across figures.
ALL_SETS = SETS + ['kh1_art_full', 'kh1_wwe_full', 'md2_art_full', 'md2_wwe_full',
                    'md3_art_full', 'md3_wwe_full', 'core_wwe_full_go',
                    'mka_wwe_full', 'fl_wwe_full_go']

# Colors match sim_registry.SETS / parp14_mdp_colors.py's CONSTRUCT_GO
# (Go-restraint sets) and CONSTRUCT_NOGO (non-Go sets) -- validated
# colorblind-safe (OKLab CVD Delta E, --pairs all) within each group; the two
# groups are never drawn in the same figure so hues don't need to be distinct
# across them. See parp14_mdp_colors.py for the derivation notes.
SET_COLORS = {
    'md_full':       '#c8831c',
    'mka_full':      '#2845bd',
    'core_full_go':  '#E69F00',
    'fl_go_5rep1us': '#CC79A7',
    'kh1_art_full':      '#56B4E9',
    'kh1_wwe_full':      '#009E73',
    'md2_art_full':      '#eb79eb',
    'md2_wwe_full':      '#a82366',
    'md3_art_full':      '#1cc895',
    'md3_wwe_full':      '#2387a8',
    'core_wwe_full_go':  '#0072B2',
    'mka_wwe_full':      '#8a62e8',
    'fl_wwe_full_go':    '#D55E00',
}

# Run parameters for the table below the title. Lengths/replicate counts are
# the ANALYZED values (post skip-frames, post any --max-reps/--max-ns cap),
# not necessarily the full config.yaml `steps:` (which for checkpoint-restart
# sims is only the last increment, not the accumulated total).
RUN_PARAMS = {
    'md_full': {
        'label': 'MD1L1-MD2-MD3 contiguous',
        'span': 'FL 790-1388 (599 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic (domain cores only)',
    },
    'mka_full': {
        'label': 'MD1L1-ART contiguous',
        'span': 'FL 790-1801 (1012 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic (domain cores only)',
    },
    'core_full_go': {
        'label': 'KH7a-ART contiguous',
        'span': 'FL 738-1801 (1064 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic + Go-model KH7a-KHb (141 pairs, k=15)',
    },
    'fl_go_5rep1us': {
        'label': 'Full length (subset)',
        'span': 'FL 1-1801 (1801 res)',
        'reps': '5 of 11', 'length': '1000 of 2000 ns',
        'restraints': 'harmonic + Go-model KH7a-KHb (141 pairs, k=15)',
    },
    'kh1_art_full': {
        'label': 'KH1-6-ART contiguous', 'span': 'FL 315-1801 (1487 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic + Go-model KH7a-KHb (141 pairs, k=15)',
    },
    'kh1_wwe_full': {
        'label': 'KH1-6-WWE contiguous (no ART)', 'span': 'FL 315-1602 (1288 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic + Go-model KH7a-KHb (141 pairs, k=15)',
    },
    'md2_art_full': {
        'label': 'MD2-ART contiguous', 'span': 'FL 1004-1801 (798 res)',
        'reps': '5', 'length': '1000 ns', 'restraints': 'harmonic (domain cores only)',
    },
    'md2_wwe_full': {
        'label': 'MD2-WWE contiguous (no ART)', 'span': 'FL 1004-1602 (599 res)',
        'reps': '5', 'length': '1000 ns', 'restraints': 'harmonic (domain cores only)',
    },
    'md3_art_full': {
        'label': 'MD3-ART contiguous', 'span': 'FL 1207-1801 (595 res)',
        'reps': '5', 'length': '1000 ns', 'restraints': 'harmonic (domain cores only)',
    },
    'md3_wwe_full': {
        'label': 'MD3-WWE contiguous (no ART)', 'span': 'FL 1207-1602 (396 res)',
        'reps': '5', 'length': '1000 ns', 'restraints': 'harmonic (domain cores only)',
    },
    'core_wwe_full_go': {
        'label': 'KH7a-WWE contiguous (no ART)', 'span': 'FL 738-1602 (865 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic + Go-model KH7a-KHb (141 pairs, k=15)',
    },
    'mka_wwe_full': {
        'label': 'MD1L1-WWE contiguous (no ART)', 'span': 'FL 790-1602 (813 res)',
        'reps': '5', 'length': '1000 ns', 'restraints': 'harmonic (domain cores only)',
    },
    'fl_wwe_full_go': {
        'label': 'FL-WWE contiguous (no ART)', 'span': 'FL 1-1602 (1602 res)',
        'reps': '5', 'length': '1000 ns',
        'restraints': 'harmonic + Go-model KH7a-KHb (141 pairs, k=15)',
    },
}


def load_cache(metric, sets, pair_keys):
    """Reconstruct pooled_arrays[sk][pk] and mean_rg[sk][domain] from the
    shared md_distances_{metric}.npz cache figure_md_distances.py writes."""
    npz_path = DATA_PATH / f'md_distances_{metric}.npz'
    if not npz_path.is_file():
        raise FileNotFoundError(
            f"{npz_path} not found -- run figure_md_distances.py --metric {metric} "
            f"--set {' '.join(sets)} first (see this script's module docstring).")
    data = np.load(npz_path)
    pooled_arrays = {sk: {} for sk in sets}
    mean_rg = {sk: {} for sk in sets}
    missing = []
    for sk in sets:
        for pk in pair_keys:
            key = f'{sk}__{pk}'
            if key in data:
                pooled_arrays[sk][pk] = data[key]
            else:
                pooled_arrays[sk][pk] = np.array([])
                missing.append(key)
        for dname in DOMAINS:
            rgkey = f'{sk}__rg__{dname}'
            if rgkey in data:
                mean_rg[sk][dname] = float(data[rgkey][0])
    if missing:
        print(f"  WARNING: cache missing {len(missing)} set/pair combos "
              f"(will be blank in the grid): {missing}")
    return pooled_arrays, mean_rg


def _draw_slash_percentages(ax, renderer, x0, y, items, fontsize=8):
    """Draw items = [(text, color), ...] left to right as 'A/B/C/...', each
    segment colored, separators in gray -- a compact one-line legend instead
    of one stacked line per set (which overlaps once there are 3-4 sets).
    x0, y in axes fraction coords. Returns nothing; draws directly on ax.

    `renderer` must be a single renderer obtained via ONE fig.canvas.draw()
    up front (see build_grid_figure) -- calling canvas.draw() (a full figure
    re-render) per text segment instead made this scale as O(n_cells *
    n_sets) full-figure redraws and was minutes slow on the real grid.
    """
    x = x0
    for i, (text, color) in enumerate(items):
        if i > 0:
            sep = ax.text(x, y, '/', transform=ax.transAxes, fontsize=fontsize,
                          color='dimgray', ha='left', va='center')
            bbox = sep.get_window_extent(renderer=renderer)
            x = ax.transAxes.inverted().transform((bbox.x1, 0))[0]
        t = ax.text(x, y, text, transform=ax.transAxes, fontsize=fontsize,
                    color=color, fontweight='bold', ha='left', va='center')
        bbox = t.get_window_extent(renderer=renderer)
        x = ax.transAxes.inverted().transform((bbox.x1, 0))[0]


def _bimodality_for_cell(pair_key, pooled_arrays, sets, bimodal_cache):
    """Compute (or fetch from bimodal_cache) the bimodality check for every
    set's array for this pair. Cached because build_grid_figure calls this
    once per (kind='hist'|'cdf'), and detect_bimodality's gaussian_kde
    evaluation is O(N) and genuinely slow (~5-6s for a 500k-point array,
    much more for fl_go_full's 2.2M) -- computing it twice per pair (once
    for hist, once for cdf) was doubling that cost for nothing.
    """
    out = {}
    for sk in sets:
        key = (sk, pair_key)
        if key in bimodal_cache:
            out[sk] = bimodal_cache[key]
            continue
        arr = pooled_arrays[sk][pair_key]
        if len(arr) == 0:
            out[sk] = None
            continue
        # Subsample for the KDE shape check only (hist/CDF below still use
        # the FULL array) -- bimodality is a shape property that doesn't
        # need every frame, and gaussian_kde's evaluate cost scales with N.
        sample = arr if len(arr) <= 100_000 else np.random.RandomState(0).choice(arr, 100_000, replace=False)
        bm = detect_bimodality(sample)
        bimodal_cache[key] = bm
        out[sk] = bm
    return out


def _domain_pair_cell(ax, renderer, pair_key, pooled_arrays, sets, mean_rg, metric, xmax, kind, bimodal_cache):
    """Populate one grid cell (kind='hist' or 'cdf') for one domain pair."""
    a, b = pair_key.split('_')
    rg_floor = None
    if metric == 'com':
        sums = [mean_rg[sk][a] + mean_rg[sk][b] for sk in sets
                if a in mean_rg.get(sk, {}) and b in mean_rg.get(sk, {})]
        if sums:
            rg_floor = float(np.mean(sums))

    bimodal = _bimodality_for_cell(pair_key, pooled_arrays, sets, bimodal_cache)

    floor_items, energetic_items = [], []
    for sk in sets:
        arr = pooled_arrays[sk][pair_key]
        if len(arr) == 0:
            continue
        color = SET_COLORS.get(sk, 'black')

        if kind == 'hist':
            ax.hist(arr, bins=70, density=True, alpha=0.5, color=color,
                    range=(0, xmax))
        else:
            sorted_arr = np.sort(arr)
            cdf = np.arange(1, len(sorted_arr) + 1) / len(sorted_arr)
            ax.plot(sorted_arr, cdf, color=color, linewidth=1.5)

        # Bimodal peak positions -- colored per-set (was black in the linear
        # per-pair figure) so multiple sets don't collide visually, drawn on
        # BOTH hist and cdf (was hist-only before).
        bm = bimodal.get(sk)
        if bm and bm['is_bimodal']:
            for p in bm['peaks']:
                ax.axvline(p, color=color, linestyle=':', linewidth=1.2, alpha=0.7)

        if rg_floor is not None:
            energetic = rg_floor + AH_CUTOFF_NM
            frac_floor = np.mean(arr < rg_floor)
            frac_energetic = np.mean(arr < energetic)
            floor_items.append((f'{frac_floor*100:.0f}%', color))
            energetic_items.append((f'{frac_energetic*100:.0f}%', color))

    if rg_floor is not None:
        ax.axvline(rg_floor, color='darkorange', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(rg_floor + AH_CUTOFF_NM, color='mediumpurple', linestyle='--', alpha=0.5, linewidth=1)
        if kind == 'cdf' and floor_items:
            _draw_slash_percentages(ax, renderer, 0.02, 0.92, floor_items, fontsize=7)
            _draw_slash_percentages(ax, renderer, 0.02, 0.84, energetic_items, fontsize=7)

    ax.set_xlim(0, xmax)
    if kind == 'cdf':
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.25)


def _draw_run_param_table(fig, sets, y0=0.955):
    """Small run-parameter table below the suptitle, one column per sim."""
    cell_text = []
    for field, key in (('Construct', 'label'), ('Span', 'span'), ('Reps', 'reps'),
                        ('Length', 'length'), ('Restraints', 'restraints')):
        cell_text.append([field] + [RUN_PARAMS[sk][key] for sk in sets])
    ax = fig.add_axes([0.06, y0 - 0.14, 0.88, 0.13])
    ax.axis('off')
    table = ax.table(cellText=cell_text, loc='center', cellLoc='center',
                     colWidths=[0.10] + [0.90 / len(sets)] * len(sets))
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor('#cccccc')
        if c == 0:
            cell.set_text_props(fontweight='bold', ha='left')
            cell.set_facecolor('#f5f5f5')
        elif r == -1:
            pass
        else:
            sk = sets[c - 1]
            cell.set_facecolor(SET_COLORS.get(sk, 'white'))
            cell.set_alpha(0.15)
    # Header row: set names, colored
    for c, sk in enumerate(sets, start=1):
        table.add_cell(-1, c, width=0.90 / len(sets), height=0.16,
                       text=sk, loc='center',
                       facecolor=SET_COLORS.get(sk, 'white'))
        table[(-1, c)].set_text_props(fontweight='bold', color='black')
    table.add_cell(-1, 0, width=0.10, height=0.16, text='', loc='center',
                   facecolor='#f5f5f5')


def build_grid_figure(pooled_arrays, mean_rg, sets, metric, out_path):
    n = len(DOMAINS)
    all_dists = [pooled_arrays[sk][pk] for sk in sets
                 for pk in pooled_arrays[sk] if len(pooled_arrays[sk][pk])]
    all_combined = np.concatenate(all_dists) if all_dists else np.array([0, 10])
    xmax = np.percentile(all_combined, 99.5) * 1.05

    fig = plt.figure(figsize=(24, 13))
    subfigs = fig.subfigures(1, 2, wspace=0.04)
    subfigs[0].suptitle('Distance Histograms', fontsize=13, fontweight='bold', y=1.01)
    subfigs[1].suptitle('Cumulative Fraction', fontsize=13, fontweight='bold', y=1.01)

    # ONE renderer for the whole figure, reused by every cell's slash-percentage
    # labels (see _draw_slash_percentages) instead of a fresh full canvas.draw()
    # per label -- that was the O(n_cells * n_sets) perf bug found on first run.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bimodal_cache = {}

    xlabel = 'COM-COM distance (nm)' if metric == 'com' else 'min CA-CA distance (nm)'

    for subfig, kind in ((subfigs[0], 'hist'), (subfigs[1], 'cdf')):
        axes = subfig.subplots(n, n)
        for i in range(n):
            for j in range(n):
                ax = axes[i, j]
                if i == j:
                    ax.axis('off')
                    ax.text(0.5, 0.5, DOMAINS[i], transform=ax.transAxes,
                            fontsize=14, fontweight='bold', ha='center', va='center')
                elif i > j:
                    ax.axis('off')
                else:
                    pk = f'{DOMAINS[i]}_{DOMAINS[j]}'
                    if pk not in pooled_arrays[sets[0]]:
                        ax.axis('off')
                        continue
                    _domain_pair_cell(ax, renderer, pk, pooled_arrays, sets, mean_rg,
                                      metric, xmax, kind, bimodal_cache)
                if i == n - 1:
                    ax.set_xlabel(xlabel, fontsize=8)
                if j == 0 and i != j:
                    ax.set_ylabel('Density' if kind == 'hist' else 'Cum. fraction',
                                  fontsize=8)
                ax.tick_params(labelsize=7)

    # Legend (set -> color) as a shared strip under both grids
    handles = [plt.Line2D([0], [0], color=SET_COLORS[sk], lw=3, label=sk) for sk in sets]
    fig.legend(handles=handles, loc='lower center', ncol=len(sets), fontsize=10,
               bbox_to_anchor=(0.5, -0.01), frameon=False)

    metric_label = 'COM-COM' if metric == 'com' else 'minimum CA-CA'
    fig.suptitle(f'Inter-Domain {metric_label} Distance Comparison — Full/Contiguous Constructs, Matched to 1000 ns',
                 fontsize=15, fontweight='bold', y=1.06)
    _draw_run_param_table(fig, sets, y0=1.0)

    fig.savefig(out_path.with_suffix('.png'), dpi=180, bbox_inches='tight')
    fig.savefig(out_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_path.with_suffix('.png')}")
    print(f"  Saved: {out_path.with_suffix('.svg')}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--metric', choices=['com', 'min'], default='com')
    parser.add_argument('--sets', nargs='+', default=SETS)
    args = parser.parse_args()

    sets = args.sets
    pair_keys = [f'{DOMAINS[i]}_{DOMAINS[j]}'
                 for i in range(len(DOMAINS)) for j in range(i + 1, len(DOMAINS))]

    print(f"Sets:  {sets}")
    print(f"Pairs: {pair_keys}")
    pooled_arrays, mean_rg = load_cache(args.metric, sets, pair_keys)

    out_dir = _get_fig_dir('04_md_distances', sims=sets)
    out_path = out_dir / f'md_distance_grid_{args.metric}'
    build_grid_figure(pooled_arrays, mean_rg, sets, args.metric, out_path)


if __name__ == '__main__':
    main()
