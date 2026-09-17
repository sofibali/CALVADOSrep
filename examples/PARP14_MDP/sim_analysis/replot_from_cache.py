#!/usr/bin/env python3
"""Re-render figure_md_distances.py's per-construct plots (violin,
histogram/CDF, block-averaging) straight from the existing
data/md_distances_{metric}.npz cache -- no DCD recompute.

Only exists because re-running figure_md_distances.py's main() to pick up a
color-only change would re-walk every replicate's DCD again (hours), and
would re-save the npz cache with only the sets passed in THIS invocation,
silently narrowing it for every other script that reads it (see
run_dashboard_backfill_final3.py's RESTORE_CACHE_CMD bug). This script never
touches the npz file -- read-only.

Usage:
    python replot_from_cache.py --set md_full mka_full core_full_go fl_go_5rep1us
    python replot_from_cache.py --set core_full_go kh1_art_full kh1_wwe_full \\
        core_wwe_full_go fl_wwe_full_go fl_go_5rep1us --label go_group
"""
import argparse
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

import figure_md_distances as fmd

CWD = Path(__file__).resolve().parent.parent
DATA_PATH = CWD / 'data'


def load_all(metric, sets, pair_keys):
    npz_path = DATA_PATH / f'md_distances_{metric}.npz'
    data = np.load(npz_path)
    pooled_arrays = {sk: {} for sk in sets}
    mean_rg = {sk: {} for sk in sets}
    pooled = {sk: {pk: [] for pk in pair_keys} for sk in sets}
    for sk in sets:
        for pk in pair_keys:
            key = f'{sk}__{pk}'
            pooled_arrays[sk][pk] = data[key] if key in data else np.array([])
            nreps_key = f'{sk}__{pk}__nreps'
            nreps = int(data[nreps_key][0]) if nreps_key in data else 0
            for i in range(nreps):
                rep_key = f'{sk}__{pk}__rep{i}'
                if rep_key in data:
                    pooled[sk][pk].append(data[rep_key])
        for dname in ('MD1L1', 'MD2', 'MD3', 'ART'):
            rgkey = f'{sk}__rg__{dname}'
            if rgkey in data:
                mean_rg[sk][dname] = float(data[rgkey][0])
    return pooled_arrays, mean_rg, pooled


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', nargs='+', required=True)
    ap.add_argument('--metric', choices=['com', 'min'], default='com')
    ap.add_argument('--pairs', nargs='+', default=None)
    args = ap.parse_args()

    sets = args.set
    metric = args.metric
    if args.pairs:
        wanted = set(args.pairs)
        pair_keys = [f'{a}_{b}' for a, b in fmd.PAIRS if f'{a}_{b}' in wanted]
    else:
        pair_keys = [f'{a}_{b}' for a, b in fmd.PAIRS]

    fmd.FIG_PATH = fmd._get_fig_dir('04_md_distances', sims=sets)

    pooled_arrays, mean_rg, pooled = load_all(metric, sets, pair_keys)

    print(f"Re-rendering {sets} ({metric}) from cache -> {fmd.FIG_PATH}")
    fmd.plot_violin(pooled_arrays, sets, pair_keys, metric, mean_rg)
    fmd.plot_block_averaging(pooled, sets, pair_keys, metric)
    fmd.plot_distance_histograms(pooled_arrays, sets, pair_keys, metric, mean_rg)
    print("Done.")


if __name__ == '__main__':
    main()
