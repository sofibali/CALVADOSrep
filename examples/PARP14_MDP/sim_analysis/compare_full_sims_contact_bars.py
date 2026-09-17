#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bar-graph companion to compare_full_sims_grid.py: for each inter-domain pair,
% of frames "in contact" by the SAME two thresholds already annotated on the
grid's CDF panels --

  floor:      COM-COM < Rg_A + Rg_B            ("surfaces touching")
  energetic:  COM-COM < Rg_A + Rg_B + AH_CUTOFF_NM   ("AH could still reach")

-- computed PER REPLICATE (not on the pooled/concatenated array), so the
across-replicate mean +/- SEM is a real error estimate, plus a Kruskal-Wallis
test across the 4 sims per (pair, threshold) and pairwise Mann-Whitney U
tests (small N=5 per group -- nonparametric, no normality assumption).

Requires figure_md_distances.py to have been run AFTER the per-replicate
caching was added (saves `{set}__{pair}__rep{i}` + `{set}__{pair}__nreps`
keys) -- re-run it if this script errors on a missing nreps key:

    python figure_md_distances.py --set md_full mka_full core_full_go \\
        --sim-folder-as fl_go_5rep1us:/path/to/fl_go \\
        --max-reps fl_go_5rep1us:5 --max-ns fl_go_5rep1us:1000

Usage:
    python compare_full_sims_contact_bars.py
"""
import argparse
import numpy as np
from pathlib import Path
from scipy.stats import kruskal, mannwhitneyu
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
from compare_full_sims_grid import DOMAINS, SETS, SET_COLORS, RUN_PARAMS

DATA_PATH = CWD / 'data'


def load_per_replicate(metric, sets, pair_keys):
    """Returns replicate_arrays[sk][pk] = list of per-replicate distance
    arrays, and mean_rg[sk][domain], from the md_distances_{metric}.npz cache
    (see figure_md_distances.py's per-replicate caching)."""
    npz_path = DATA_PATH / f'md_distances_{metric}.npz'
    if not npz_path.is_file():
        raise FileNotFoundError(f"{npz_path} not found -- run figure_md_distances.py first.")
    data = np.load(npz_path)
    replicate_arrays = {sk: {} for sk in sets}
    mean_rg = {sk: {} for sk in sets}
    for sk in sets:
        for pk in pair_keys:
            nkey = f'{sk}__{pk}__nreps'
            if nkey not in data:
                replicate_arrays[sk][pk] = []
                continue
            n = int(data[nkey][0])
            replicate_arrays[sk][pk] = [data[f'{sk}__{pk}__rep{i}'] for i in range(n)]
        for dname in DOMAINS:
            rgkey = f'{sk}__rg__{dname}'
            if rgkey in data:
                mean_rg[sk][dname] = float(data[rgkey][0])
    return replicate_arrays, mean_rg


def per_replicate_fractions(arr_list, threshold):
    """Fraction of frames < threshold, computed separately per replicate."""
    return np.array([np.mean(a < threshold) * 100 for a in arr_list if len(a) > 0])


def sig_stars(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


def _pairwise_tests(labels, per_sim_fracs):
    """All-pairs Mann-Whitney U between every sim with >=2 replicate values,
    Holm-Bonferroni corrected across this panel's own family of comparisons
    (one panel = one (domain-pair, threshold) combination; comparisons in
    OTHER panels test different hypotheses so aren't pooled into the same
    correction family). Returns list of dicts: i, j (indices into `labels`),
    p_raw, p_corr, reject.
    """
    idx_pairs = [(i, j) for i, j in combinations(range(len(labels)), 2)
                 if len(per_sim_fracs.get(labels[i], [])) > 1
                 and len(per_sim_fracs.get(labels[j], [])) > 1]
    if not idx_pairs:
        return []
    raw_p = []
    for i, j in idx_pairs:
        try:
            _, p = mannwhitneyu(per_sim_fracs[labels[i]], per_sim_fracs[labels[j]],
                                alternative='two-sided')
        except ValueError:
            p = 1.0
        raw_p.append(p)
    reject, p_corr, _, _ = multipletests(raw_p, alpha=0.05, method='holm')
    return [{'i': i, 'j': j, 'p_raw': pr, 'p_corr': pc, 'reject': rj}
            for (i, j), pr, pc, rj in zip(idx_pairs, raw_p, p_corr, reject)]


def _draw_sig_brackets(ax, tests, y_base, y_step):
    """Draw one horizontal bracket per pairwise test, stacked bottom to top,
    narrowest span first so brackets nest visually instead of crossing.
    Returns the y coordinate just above the topmost bracket."""
    order = sorted(range(len(tests)), key=lambda k: tests[k]['j'] - tests[k]['i'])
    y_top = y_base
    for rank, idx in enumerate(order):
        t = tests[idx]
        y = y_base + rank * y_step
        tick = y_step * 0.12
        ax.plot([t['i'], t['i'], t['j'], t['j']],
                [y, y + tick, y + tick, y], color='black', linewidth=0.8)
        label = sig_stars(t['p_corr'])
        ax.text((t['i'] + t['j']) / 2, y + tick, label, ha='center', va='bottom',
                fontsize=7, fontweight='bold' if label != 'ns' else 'normal',
                color='black' if label != 'ns' else 'gray')
        y_top = max(y_top, y + tick + y_step * 0.3)
    return y_top


def build_bar_figure(replicate_arrays, mean_rg, sets, pair_keys, metric, out_path):
    n_pairs = len(pair_keys)
    fig, axes = plt.subplots(2, n_pairs, figsize=(4.8 * n_pairs, 11), squeeze=False)
    thresholds = [('floor', 'Rg_A+Rg_B\n("surfaces touching")'),
                  ('energetic', f'Rg_A+Rg_B+{AH_CUTOFF_NM:.0f}\n("AH could still reach")')]

    stats_report = []

    for col, pk in enumerate(pair_keys):
        a, b = pk.split('_')
        sums = [mean_rg[sk][a] + mean_rg[sk][b] for sk in sets
                if a in mean_rg.get(sk, {}) and b in mean_rg.get(sk, {})]
        rg_floor = float(np.mean(sums)) if sums else None

        for row, (thresh_key, thresh_label) in enumerate(thresholds):
            ax = axes[row, col]
            if rg_floor is None:
                ax.axis('off')
                continue
            threshold_val = rg_floor if thresh_key == 'floor' else rg_floor + AH_CUTOFF_NM

            means, sems, colors, labels, per_sim_fracs = [], [], [], [], {}
            for sk in sets:
                arrs = replicate_arrays[sk].get(pk, [])
                if not arrs:
                    continue
                fracs = per_replicate_fractions(arrs, threshold_val)
                if len(fracs) == 0:
                    continue
                per_sim_fracs[sk] = fracs
                means.append(fracs.mean())
                sems.append(fracs.std(ddof=1) / np.sqrt(len(fracs)) if len(fracs) > 1 else 0.0)
                colors.append(SET_COLORS.get(sk, 'gray'))
                labels.append(sk)

            x = np.arange(len(labels))
            ax.bar(x, means, yerr=sems, color=colors, alpha=0.85, capsize=4,
                   edgecolor='black', linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=7)
            if col == 0:
                ax.set_ylabel(f'% frames in contact\n({thresh_label})', fontsize=8)
            if row == 0:
                ax.set_title(f'{a}-{b}', fontsize=11, fontweight='bold')

            # Every sim compared to every other sim (Mann-Whitney U, small-N
            # nonparametric -- no normality assumption; N=5 replicates/group),
            # Holm-Bonferroni corrected within this panel's family of
            # comparisons, drawn as stacked brackets above the bars.
            tests = _pairwise_tests(labels, per_sim_fracs)
            bar_top = max((m + s for m, s in zip(means, sems)), default=0)
            y_base = max(100, bar_top + 8)
            y_step = max(9, y_base * 0.09)
            if tests:
                y_top = _draw_sig_brackets(ax, tests, y_base, y_step)
            else:
                y_top = y_base
            ax.set_ylim(0, y_top + y_step * 0.5)

            for t in tests:
                stats_report.append({'pair': pk, 'threshold': thresh_key,
                                     'test': f'mann-whitney {labels[t["i"]]} vs {labels[t["j"]]} '
                                             f'(Holm-corrected, n={len(tests)} comparisons)',
                                     'p': t['p_corr'], 'p_raw': t['p_raw']})

            # Omnibus check (Kruskal-Wallis across all sims at once) kept in
            # the console report only, for cross-reference against the
            # pairwise Holm-corrected results above.
            groups = [v for v in per_sim_fracs.values() if len(v) > 1]
            if len(groups) >= 2:
                try:
                    kw_stat, kw_p = kruskal(*groups)
                except ValueError:
                    kw_p = float('nan')
                stats_report.append({'pair': pk, 'threshold': thresh_key,
                                     'test': 'kruskal-wallis (omnibus, uncorrected)', 'p': kw_p,
                                     'p_raw': kw_p})

    handles = [plt.Rectangle((0, 0), 1, 1, color=SET_COLORS.get(sk, 'gray')) for sk in sets]
    fig.legend(handles, sets, loc='lower center', ncol=len(sets), fontsize=10,
               bbox_to_anchor=(0.5, -0.02), frameon=False)

    metric_label = 'COM-COM' if metric == 'com' else 'minimum CA-CA'
    fig.suptitle(f'% Frames in Contact — {metric_label} Distance, Full/Contiguous Constructs\n'
                 f'(error bars: SEM across replicates; brackets: pairwise Mann-Whitney U, '
                 f'Holm-Bonferroni corrected per panel; ns/*/**/*** = p>=0.05/<0.05/<0.01/<0.001)',
                 fontsize=13, fontweight='bold', y=1.03)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix('.png'), dpi=180, bbox_inches='tight')
    fig.savefig(out_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_path.with_suffix('.png')}")
    print(f"  Saved: {out_path.with_suffix('.svg')}")

    print("\n  Statistical test report:")
    for r in stats_report:
        raw_note = f" (raw p={r['p_raw']:.4g})" if r.get('p_raw') is not None and r['p_raw'] != r['p'] else ""
        print(f"    [{r['pair']}] {r['threshold']}: {r['test']}: "
              f"p={r['p']:.4g}{raw_note} {sig_stars(r['p'])}")


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
    replicate_arrays, mean_rg = load_per_replicate(args.metric, sets, pair_keys)

    out_dir = _get_fig_dir('04_md_distances', sims=sets)
    out_path = out_dir / f'md_distance_contact_bars_{args.metric}'
    build_bar_figure(replicate_arrays, mean_rg, sets, pair_keys, args.metric, out_path)


if __name__ == '__main__':
    main()
