#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Companion to compare_full_sims_contact_bars.py, for a DIFFERENT correction
scope: instead of correcting the ~6 sim-pair comparisons WITHIN one panel
(one domain-pair/threshold combination) together, this corrects each SIM-PAIR
(e.g. md_full vs mka_full) across the 12 panels it's repeated in (6 domain
pairs x 2 thresholds: floor/energetic) -- because those 12 tests are NOT
independent: all 12 draw on the SAME 5 replicate trajectories per sim
(replicate i's MD1L1-MD2 value and replicate i's MD1L1-MD3 value come from
the same simulation run).

Benjamini-Hochberg FDR (not Holm-Bonferroni) is used for this cross-panel
family: with only N=5 replicates/sim, a Mann-Whitney U test's smallest
possible two-sided p-value is fixed at 2/C(10,5) = 0.0079 (perfect
separation, the strongest possible signal) -- Holm-correcting that across 12
repeated tests (0.0079 x 12 = 0.095) would make EVERY comparison "ns"
regardless of true effect size. FDR can still recover a real, consistent
signal that shows up across several of the 12 correlated panels (expected
here, since several domain pairs share a moving domain and should co-vary).

Output: ONE wide combined figure (all sim-pairs' brackets overlaid on the
same all-4-sims-per-panel bar layout as compare_full_sims_contact_bars.py),
rather than one separate figure per sim-pair.

Usage:
    python compare_full_sims_contact_bars_pairwise_fdr.py
    python compare_full_sims_contact_bars_pairwise_fdr.py --sets md_full mka_full
"""
import argparse
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from itertools import combinations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys as _sys
CWD = Path(__file__).resolve().parent.parent
_sys.path.insert(0, str(CWD))
from _fig_layout import get_fig_dir as _get_fig_dir

from figure_md_distances import AH_CUTOFF_NM
from compare_full_sims_grid import DOMAINS, SETS, SET_COLORS
from compare_full_sims_contact_bars import (
    load_per_replicate, per_replicate_fractions, sig_stars, _draw_sig_brackets)

THRESHOLDS = [('floor', 'Rg_A+Rg_B\n("surfaces touching")'),
              ('energetic', f'Rg_A+Rg_B+{AH_CUTOFF_NM:.0f}\n("AH could still reach")')]


def compute_sim_pair_qvalues(replicate_arrays, mean_rg, sk1, sk2, pair_keys):
    """All 12 (domain-pair, threshold) Mann-Whitney U tests for ONE sim-pair,
    Benjamini-Hochberg FDR corrected together (see module docstring for why).
    Returns list of dicts (one per valid panel): pair, thresh_key, p_raw, p_fdr.
    """
    panel_data = []
    for pk in pair_keys:
        a, b = pk.split('_')
        sums = [mean_rg[sk][a] + mean_rg[sk][b] for sk in (sk1, sk2)
                if a in mean_rg.get(sk, {}) and b in mean_rg.get(sk, {})]
        rg_floor = float(np.mean(sums)) if len(sums) == 2 else None
        for thresh_key, _ in THRESHOLDS:
            if rg_floor is None:
                continue
            threshold_val = rg_floor if thresh_key == 'floor' else rg_floor + AH_CUTOFF_NM
            f1 = per_replicate_fractions(replicate_arrays[sk1].get(pk, []), threshold_val)
            f2 = per_replicate_fractions(replicate_arrays[sk2].get(pk, []), threshold_val)
            if len(f1) < 2 or len(f2) < 2:
                continue
            try:
                _, p = mannwhitneyu(f1, f2, alternative='two-sided')
            except ValueError:
                p = 1.0
            panel_data.append({'pair': pk, 'thresh_key': thresh_key, 'p_raw': p})

    if not panel_data:
        return []
    raw_p = [e['p_raw'] for e in panel_data]
    _, p_fdr, _, _ = multipletests(raw_p, alpha=0.05, method='fdr_bh')
    for e, q in zip(panel_data, p_fdr):
        e['p_fdr'] = q
    return panel_data


def build_combined_figure(replicate_arrays, mean_rg, sets, pair_keys, metric, out_path, reference=None):
    """ONE wide figure: 2 rows (floor/energetic) x len(pair_keys) columns, all
    `sets` as bars per panel, with brackets for every sim-pair -- each
    bracket's significance is that specific sim-pair's own FDR q-value
    (computed across ITS 12 panels, not the other sim-pairs').

    reference: if given, only (reference, other) sim-pairs are computed/drawn
    -- C(n,2) pairwise brackets is unreadable/expensive past a handful of
    sets (e.g. 78 brackets for 13 sets); comparing everything against one
    reference instead gives n-1 brackets per group, still fully FDR-corrected
    per sim-pair across its own 12 repeated panels.
    """
    # qvals[(sk1, sk2)][(pair, thresh_key)] = q
    qvals = {}
    sim_pairs = ([(reference, sk) for sk in sets if sk != reference] if reference
                 else list(combinations(sets, 2)))
    for sk1, sk2 in sim_pairs:
        panel_data = compute_sim_pair_qvalues(replicate_arrays, mean_rg, sk1, sk2, pair_keys)
        qvals[(sk1, sk2)] = {(e['pair'], e['thresh_key']): e['p_fdr'] for e in panel_data}
        for e in panel_data:
            print(f"    [{sk1} vs {sk2}] {e['pair']} {e['thresh_key']}: "
                  f"raw p={e['p_raw']:.4g}, FDR q={e['p_fdr']:.4g} {sig_stars(e['p_fdr'])}")

    n_pairs = len(pair_keys)
    gap = 1.5  # empty units between domain-pair groups
    fig, axes = plt.subplots(2, 1, figsize=(3.2 * n_pairs, 12), squeeze=True)

    for row, (thresh_key, thresh_label) in enumerate(THRESHOLDS):
        ax = axes[row]
        cursor = 0.0
        group_centers, group_names = [], []
        overall_y_top = 100.0

        for pk in pair_keys:
            a, b = pk.split('_')
            rg_sums = [mean_rg[sk][a] + mean_rg[sk][b] for sk in sets
                       if a in mean_rg.get(sk, {}) and b in mean_rg.get(sk, {})]
            rg_floor = float(np.mean(rg_sums)) if rg_sums else None
            if rg_floor is None:
                continue
            threshold_val = rg_floor if thresh_key == 'floor' else rg_floor + AH_CUTOFF_NM

            means, sems, colors, labels, xpos = [], [], [], [], []
            for sk in sets:
                fracs = per_replicate_fractions(replicate_arrays[sk].get(pk, []), threshold_val)
                if len(fracs) == 0:
                    continue
                means.append(fracs.mean())
                sems.append(fracs.std(ddof=1) / np.sqrt(len(fracs)) if len(fracs) > 1 else 0.0)
                colors.append(SET_COLORS.get(sk, 'gray'))
                labels.append(sk)
                xpos.append(cursor + len(labels) - 1)

            ax.bar(xpos, means, yerr=sems, color=colors, alpha=0.85, capsize=4,
                   edgecolor='black', linewidth=0.5, width=0.8)
            group_centers.append(cursor + (len(labels) - 1) / 2)
            group_names.append(f'{a}-{b}')

            # Look up each bracket's significance from ITS OWN sim-pair's FDR
            # family (computed above across that sim-pair's 12 panels), not a
            # fresh per-panel correction. i/j here are ABSOLUTE x positions
            # (this group's bar slots), not local 0..n indices, so brackets
            # land correctly on a single shared axis with many groups.
            idx_pairs = (combinations(range(len(labels)), 2) if reference is None else
                        [(i, j) for i, j in combinations(range(len(labels)), 2)
                         if labels[i] == reference or labels[j] == reference])
            tests = []
            for (i_local, j_local) in idx_pairs:
                sk1, sk2 = labels[i_local], labels[j_local]
                key = (sk1, sk2) if (sk1, sk2) in qvals else (sk2, sk1)
                q = qvals.get(key, {}).get((pk, thresh_key))
                if q is not None:
                    tests.append({'i': xpos[i_local], 'j': xpos[j_local], 'p_corr': q})

            bar_top = max((m + s for m, s in zip(means, sems)), default=0)
            y_base = max(100, bar_top + 8)
            y_step = max(9, y_base * 0.09)
            y_top = _draw_sig_brackets(ax, tests, y_base, y_step) if tests else y_base
            overall_y_top = max(overall_y_top, y_top + y_step * 0.5)

            cursor += len(labels) + gap

        ax.set_xticks(group_centers)
        ax.set_xticklabels(group_names, fontsize=10, fontweight='bold')
        ax.set_ylabel(f'% frames in contact\n({thresh_label})', fontsize=9)
        ax.set_ylim(0, overall_y_top)
        ax.set_xlim(-1, cursor - gap + 1)

    handles = [plt.Rectangle((0, 0), 1, 1, color=SET_COLORS.get(sk, 'gray')) for sk in sets]
    fig.legend(handles, sets, loc='lower center', ncol=len(sets), fontsize=10,
               bbox_to_anchor=(0.5, -0.02), frameon=False)

    metric_label = 'COM-COM' if metric == 'com' else 'minimum CA-CA'
    ref_note = f' -- brackets vs. reference "{reference}" only' if reference else ''
    fig.suptitle(f'% Frames in Contact — {metric_label} Distance, Full/Contiguous Constructs{ref_note}\n'
                 f'(error bars: SEM across replicates; brackets: pairwise Mann-Whitney U, '
                 f'Benjamini-Hochberg FDR corrected PER SIM-PAIR across its 12 repeated panels '
                 f'[NOT per-panel Holm -- see module docstring]; ns/*/**/*** = q>=0.05/<0.05/<0.01/<0.001)',
                 fontsize=13, fontweight='bold', y=1.03)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix('.png'), dpi=180, bbox_inches='tight')
    fig.savefig(out_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_path.with_suffix('.png')}")
    print(f"  Saved: {out_path.with_suffix('.svg')}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--metric', choices=['com', 'min'], default='com')
    parser.add_argument('--sets', nargs='+', default=SETS)
    parser.add_argument('--reference', default=None,
                        help='Only compute/draw brackets between this set and each other '
                             '(n-1 brackets/group instead of C(n,2)) -- required past a '
                             'handful of --sets or the bracket count becomes unreadable.')
    args = parser.parse_args()

    sets = args.sets
    pair_keys = [f'{DOMAINS[i]}_{DOMAINS[j]}'
                 for i in range(len(DOMAINS)) for j in range(i + 1, len(DOMAINS))]

    print(f"Sets:      {sets}")
    print(f"Pairs:     {pair_keys}")
    print(f"Reference: {args.reference}")
    replicate_arrays, mean_rg = load_per_replicate(args.metric, sets, pair_keys)

    out_dir = _get_fig_dir('04_md_distances', sims=sets)
    suffix = f'_vs_{args.reference}' if args.reference else '_combined'
    out_path = out_dir / f'md_distance_contact_bars_fdr{suffix}_{args.metric}'
    print("\n  Statistical test report (per sim-pair, FDR-corrected across its own 12 panels):")
    build_combined_figure(replicate_arrays, mean_rg, sets, pair_keys, args.metric, out_path,
                          reference=args.reference)


if __name__ == '__main__':
    main()
