#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Re-plot ONE fresh comparison figure across simulations, reading straight from
the cached data/*.npz arrays each analysis script already writes -- no
recomputation, no re-reading trajectories.

This is the "flags" counterpart to the by-simulation dashboard
(build_sim_dashboard.py shows everything about ONE sim; this shows ONE
metric across SEVERAL sims side by side).

Usage:
    python compare_sims.py --metric rg --sims fl_go fl_optimized md
    python compare_sims.py --metric ree --sims fl fl_optimized
    python compare_sims.py --metric list          # show available metrics

Output: figures/comparisons/<metric>/<sims_joined>.png (+ .svg)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CWD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CWD))
from _fig_layout import get_fig_dir  # noqa: E402
import sim_registry as reg  # noqa: E402

DATA_PATH = CWD / 'data'


def load_conf_prop(sims, prop):
    """{sim: [per-replicate mean values]} for prop in ('rg', 'ree'), from the
    conf_prop.npz cache analyze_all.py's --conf-prop module writes (merge-
    safe across separate per-sim runs -- see merge_savez in analyze_all.py)."""
    f = DATA_PATH / 'conf_prop.npz'
    if not f.is_file():
        raise FileNotFoundError(
            f"{f} not found -- run analyze_all.py --conf-prop on these sims first")
    d = np.load(f, allow_pickle=True)
    out = {}
    for sk in sims:
        prefix = f'{sk}_{prop}_rep'
        keys = sorted(k for k in d.files if k.startswith(prefix))
        if not keys:
            print(f"  WARN: no '{prop}' data for '{sk}' in {f.name} -- "
                  f"run analyze_all.py --sim-folder <path to {sk}> --conf-prop")
            continue
        out[sk] = [float(np.mean(d[k])) for k in keys]
    return out


def load_md_distance(sims, pair, metric='com'):
    """{sim: array of pooled per-frame distances} for one domain pair, from
    figure_md_distances.py's cache. NOTE: that cache file (md_distances_
    <metric>.npz) is overwritten each run with only the sets from that one
    invocation -- if a requested sim is missing, re-run figure_md_distances.py
    with ALL the sims you want to compare in one --set/--sim-folder call."""
    f = DATA_PATH / f'md_distances_{metric}.npz'
    if not f.is_file():
        raise FileNotFoundError(
            f"{f} not found -- run figure_md_distances.py --metric {metric} "
            f"--set {' '.join(sims)} first")
    d = np.load(f, allow_pickle=True)
    out = {}
    for sk in sims:
        key = f'{sk}__{pair}'
        if key not in d.files:
            print(f"  WARN: '{pair}' not cached for '{sk}' in {f.name} -- "
                  f"re-run figure_md_distances.py --metric {metric} "
                  f"--set {' '.join(sims)} (together, in one call) to refresh it")
            continue
        out[sk] = d[key]
    return out


METRICS = {
    'rg': dict(kind='conf_prop', prop='rg', ylabel='Radius of gyration (nm)'),
    'ree': dict(kind='conf_prop', prop='ree', ylabel='End-to-end distance (nm)'),
}


def sim_label(sk):
    return reg.SETS.get(sk, {}).get('label', sk)


def sim_color(sk, i, n):
    c = reg.SETS.get(sk, {}).get('color')
    return c if c else plt.cm.tab10(i / max(n - 1, 1))


def plot_conf_prop_comparison(metric_key, sims, data):
    spec = METRICS[metric_key]
    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(sims) + 2), 6))
    positions, labels, colors = [], [], []
    for i, sk in enumerate(sims):
        vals = data.get(sk)
        if not vals:
            continue
        pos = len(positions)
        color = sim_color(sk, i, len(sims))
        parts = ax.violinplot([vals], positions=[pos], showmeans=True, showextrema=True)
        for body in parts['bodies']:
            body.set_facecolor(color)
            body.set_alpha(0.6)
        rng = np.random.default_rng(42)
        jitter = rng.normal(0, 0.05, len(vals))
        ax.scatter(pos + jitter, vals, c=[color], s=25, alpha=0.8,
                  edgecolors='black', linewidths=0.5, zorder=3)
        positions.append(pos)
        labels.append(sim_label(sk))
        colors.append(color)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=10, rotation=20, ha='right')
    ax.set_ylabel(spec['ylabel'], fontsize=12)
    ax.set_title(f'{spec["ylabel"]} — {len(positions)} sets compared\n'
                 f'(each point = one replicate mean)', fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--metric', required=True,
                    help=f"One of {list(METRICS)}, a 'DOMAIN_DOMAIN' inter-domain pair "
                         f"(e.g. MD1L1_MD3), or 'list' to show available metrics")
    ap.add_argument('--sims', nargs='+', default=None,
                    help='Sim/set keys to compare (e.g. fl_go fl_optimized md)')
    ap.add_argument('--distance-metric', choices=['com', 'min'], default='com',
                    help="For a DOMAIN_DOMAIN pair metric: which figure_md_distances.py "
                         "cache to read (default: com)")
    args = ap.parse_args()

    if args.metric == 'list':
        print("Built-in metrics:", list(METRICS))
        print("Or pass a domain pair (e.g. MD1L1_MD3) to compare an "
              "inter-domain distance from figure_md_distances.py's cache.")
        return

    if not args.sims or len(args.sims) < 2:
        raise SystemExit("--sims needs at least 2 sim/set keys to compare")

    if args.metric in METRICS:
        spec = METRICS[args.metric]
        data = load_conf_prop(args.sims, spec['prop'])
        fig = plot_conf_prop_comparison(args.metric, args.sims, data)
    elif '_' in args.metric:
        data = load_md_distance(args.sims, args.metric, args.distance_metric)
        fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(args.sims) + 2), 6))
        positions, labels = [], []
        for i, sk in enumerate(args.sims):
            arr = data.get(sk)
            if arr is None or len(arr) == 0:
                continue
            pos = len(positions)
            color = sim_color(sk, i, len(args.sims))
            parts = ax.violinplot([arr], positions=[pos], showmeans=True, showextrema=False)
            for body in parts['bodies']:
                body.set_facecolor(color); body.set_alpha(0.6)
            positions.append(pos); labels.append(sim_label(sk))
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=10, rotation=20, ha='right')
        metric_label = 'COM-COM' if args.distance_metric == 'com' else 'min CA-CA'
        ax.set_ylabel(f'{args.metric.replace("_", "-")} {metric_label} distance (nm)', fontsize=12)
        ax.set_title(f'{args.metric.replace("_", "-")} distance — {len(positions)} sets compared', fontsize=12)
        ax.grid(axis='y', alpha=0.3)
        fig.tight_layout()
    else:
        raise SystemExit(f"unknown --metric '{args.metric}'. Use --metric list to see options, "
                         f"or pass a DOMAIN_DOMAIN pair (e.g. MD1L1_MD3).")

    # get_fig_dir already routes to figures/comparisons/<metric>/<sims_joined>/
    # for a multi-sim call, so the filename itself just needs to identify the
    # metric (the directory already encodes which sims).
    out_dir = get_fig_dir(args.metric, sims=args.sims)
    out_png = out_dir / f'{args.metric}.png'
    out_svg = out_dir / f'{args.metric}.svg'
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    fig.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_png}")
    print(f"  Saved: {out_svg}")


if __name__ == '__main__':
    main()
