#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Featurization sweep for MSM/TICA analysis of CALVADOS CG-MD.

Runs the same diagnostic across multiple featurization strategies in
parallel, producing ITS validation plots and (optionally) VAMP-2 scoring.
Used to pick the best featurization BEFORE doing the full clustering
analysis in cluster_states.py.

Featurization configs tested by default (in CONFIGS):
  linker_ca     positions of unrestrained linker CAs (slow dofs)
  orient        inter-domain rigid-body rotation matrices (hinge motion)
  interface_ca  CA-CA distances at domain edges (hinge contacts)
  ca_stride25   pairwise CA-CA distances within domains (intra-domain detail)
  com           inter-domain COM-COM distances (baseline)

Default lag times (linear x-axis, evenly labeled):
  ITS:  20, 40, 80, 140, 200, 400 frames = 0.2, 0.4, 0.8, 1.4, 2.0, 4.0 ns
  VAMP: 40, 80, 110 frames           = 0.4, 0.8, 1.1 ns

For each featurization:
  1. Extract or load cached features
  2. Run TICA at multiple lags with bootstrap CI
  3. (Optional --vamp) compute VAMP-2 score, either block-bootstrapped
     (biased, training-set score) or --cv k-fold cross-validated
     (unbiased, held-out-test score; recommended).
  4. Produce ITS / VAMP plots and a cross-featurization comparison plot.

Output: figures/05_clustering/<date>/sweep_its_<set>/
  ├── <name>.png            (per-featurization ITS plot)
  ├── <name>_vamp.png       (per-featurization VAMP plot, if --vamp)
  ├── comparison_IC1.png    (ITS IC1 across all featurizations)
  ├── comparison_VAMP2_by_lag.png    (grouped bars per lag)
  └── comparison_VAMP2_panels.png    (per-lag panels)

Plus cached arrays in data/its_results_*.npz and data/vamp_*_*.npz.

Usage:
    # Standard ITS sweep
    python sweep_its.py --set fl_optimized

    # With CV-VAMP-2 (recommended for picking the best featurization)
    python sweep_its.py --set fl_optimized --vamp-only --cv --cv-folds 5

    # Custom subset of featurizations
    python sweep_its.py --set fl_optimized --vamp --cv \\
        --only linker_ca orient interface_ca

    # Custom lags or skip slow configs
    python sweep_its.py --set md --its-lags 50 100 200 400 800 \\
        --skip ca_stride25

    # Replot from cache with more ICs (no recompute)
    python sweep_its.py --set fl_optimized --vamp --n-show 10 \\
        --replot-only
"""
import os
# Cap BLAS threads BEFORE numpy imports
os.environ.setdefault('OPENBLAS_NUM_THREADS', '8')
os.environ.setdefault('OMP_NUM_THREADS', '8')
os.environ.setdefault('MKL_NUM_THREADS', '8')

import sys
import argparse
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CWD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CWD))

# Reuse functions from cluster_states.py
import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
from cluster_states import (
    collect_features, compute_tica, implied_timescales,
    its_sweep, its_sweep_bootstrap, plot_its,
    vamp_sweep, vamp_sweep_bootstrap, vamp_sweep_cv, plot_vamp,
    get_construct_domain_ranges, DOMAINS_FL,
    DATA_PATH, SEEDS, SAMPLES,
)
from _fig_layout import get_fig_dir


# Featurization configs to try (display_name, args.features, stride, add_rg)
# Featurizations targeting SLOW degrees of freedom in restrained CG-MD
# (where intra-domain modes are frozen by the restraints):
#
#   - linker_ca:    positions of unrestrained linker CAs (centered at
#                   chain COM). Where the real slow motion lives.
#                   stride=5 means every 5th linker residue.
#   - orient:       inter-domain relative rotation matrices, 9 elements
#                   per pair. Captures rigid-body domain rotations.
#   - interface_ca: CA-CA distances between edge residues of adjacent
#                   domains — residues that close/open at hinges.
#                   stride=50 means use ~10 edge residues per domain.
#   - ca_stride25:  intra-domain CA-CA distances (intermediate detail).
#   - com:          COM-COM domain distances (simple baseline).
CONFIGS = [
    ('linker_ca',     'linker_ca',     5,  False),
    ('orient',        'orient',        10, False),
    ('interface_ca',  'interface_ca',  50, False),
    ('inter_exposed', 'inter_exposed', 5,  False),  # inter-domain, exposed-CA only
    ('ca_stride25',   'ca',            25, False),
    ('com',           'com',           10, False),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--set', required=True,
                        help='Set name (fl, fl_optimized, md, etc.)')
    parser.add_argument('--domains', nargs='+',
                        default=['RRM1', 'RRM2', 'RRM3',
                                  'KH1-KH6', 'KH7a',
                                  'MD1L1', 'MD2', 'MD3',
                                  'KHb-KH8', 'WWE', 'ART'],
                        help='All 11 grouped domain units by default. '
                             'Auto-filtered to those present in the set.')
    parser.add_argument('--its-lags', nargs='+', type=int,
                        default=[20, 40, 80, 140, 200, 400],
                        help='Lag times in frames (default = 0.2, 0.4, '
                             '0.8, 1.4, 2.0, 4.0 ns at 10 ps/frame)')
    parser.add_argument('--vamp-lags', nargs='+', type=int,
                        default=[40, 80, 110],
                        help='Lag times for VAMP-2 scoring (default = '
                             '0.4, 0.8, 1.1 ns at 10 ps/frame). Fewer lags '
                             '= faster.')
    parser.add_argument('--cv', action='store_true',
                        help='Use k-fold cross-validation for VAMP-2 '
                             'instead of bootstrap (recommended — '
                             'correctly handles overfitting bias).')
    parser.add_argument('--cv-folds', type=int, default=5,
                        help='k for k-fold CV (default 5). With 25 reps, '
                             'this gives 5 test reps per fold.')
    parser.add_argument('--bootstrap', type=int, default=30,
                        help='Bootstrap resamples for 95%% CI (default 30; '
                             '0 = no CI)')
    parser.add_argument('--n-components', type=int, default=10,
                        help='Top N TICA components to compute')
    parser.add_argument('--n-show', type=int, default=5,
                        help='Number of ICs to show on plot (default 5; '
                             'n_components=10 are stored, so can be up '
                             'to 10 without recomputing)')
    parser.add_argument('--replot-only', action='store_true',
                        help='Load cached ITS/VAMP results and just '
                             'regenerate plots. Skips all bootstrap.')
    parser.add_argument('--force-its', action='store_true',
                        help='Force recomputation of ITS results '
                             '(ignore cached bootstrap arrays).')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--skip', nargs='+', default=[],
                        help='Skip these configs (by name)')
    parser.add_argument('--only', nargs='+', default=None,
                        help='Only run these configs')
    parser.add_argument('--x-scale', choices=['linear', 'log'],
                        default='linear',
                        help='X-axis scale (default linear, starting at 0)')
    parser.add_argument('--vamp', action='store_true',
                        help='Also compute VAMP-2 score sweep across lags. '
                             'Use to rank featurizations objectively.')
    parser.add_argument('--vamp-only', action='store_true',
                        help='Only VAMP-2 sweep (skip ITS). Faster.')
    args = parser.parse_args()

    # Set name auto-prefix
    set_key = args.set
    if set_key not in ('fl', 'fl_optimized', 'md', 'core', 'mka',
                       'norrm', 'noart', 'md3art'):
        if (CWD / 'fragments' / set_key).is_dir():
            set_key = f'frag_{set_key}'

    # Filter requested configs
    configs = list(CONFIGS)
    if args.only:
        configs = [c for c in configs if c[0] in args.only]
    configs = [c for c in configs if c[0] not in args.skip]
    print(f"Set: {set_key}")
    print(f"Configs: {[c[0] for c in configs]}")
    print(f"Lags (frames): {args.its_lags}")
    print(f"Bootstrap: {args.bootstrap}")
    print(f"X-axis: {args.x_scale}\n")

    # Resolve effective domain list (drop domains not in this construct)
    available_ranges = get_construct_domain_ranges(set_key, args.domains)
    effective_domains = list(available_ranges.keys())
    if len(effective_domains) < 2:
        print(f"ERROR: too few domains in {set_key}: {effective_domains}")
        return 1
    print(f"Domains used: {effective_domains}\n")

    out_dir = get_fig_dir('05_clustering', subname=f'sweep_its_{set_key.replace("frag_","")}')
    print(f"Output dir: {out_dir}\n")

    # Track results for comparison plot
    all_results = {}  # config_name -> (ts_median, ts_lower, ts_upper)

    for name, features, ca_stride, add_rg in configs:
        print(f"━━━ {name} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

        # Cache key encodes features
        # Cache key encodes features + stride
        if features == 'ca':
            feat_tag_cache = f'_ca{ca_stride}'
        elif features == 'com':
            feat_tag_cache = ''
        else:
            feat_tag_cache = f'_{features}{ca_stride}'
        rg_tag = '_rg' if add_rg else ''
        cache_npz = DATA_PATH / (
            f'cluster_features_{set_key.replace("frag_","")}'
            f'{feat_tag_cache}{rg_tag}.npz')

        # Load or compute features
        if cache_npz.exists():
            print(f"  Loading cached features: {cache_npz.name}")
            cached = np.load(cache_npz, allow_pickle=True)
            X_raw = cached['features']
            metadata = list(cached['metadata'])
            if 'rep_boundaries' in cached:
                rep_boundaries = [tuple(b) for b in cached['rep_boundaries']]
            else:
                # Reconstruct
                rep_boundaries = []
                cur = None
                start = 0
                for i, m in enumerate(metadata):
                    k = (int(m[0]), int(m[1]))
                    if cur is None:
                        cur = k
                    elif k != cur:
                        rep_boundaries.append((start, i))
                        start = i
                        cur = k
                rep_boundaries.append((start, len(metadata)))
        else:
            print(f"  Extracting features (cache miss)...")
            data = collect_features(set_key, effective_domains, add_rg,
                                     args.workers, features, ca_stride)
            if data is None:
                print(f"  FAILED: no replicates loaded")
                continue
            X_raw = data['features']
            metadata = data['metadata']
            rep_boundaries = data['rep_boundaries']
            np.savez(cache_npz, features=X_raw,
                     metadata=np.array(metadata, dtype=object),
                     feature_names=np.array(data['feature_names']),
                     rep_boundaries=np.array(rep_boundaries))
            print(f"  Saved cache: {cache_npz.name}")

        # Standardize
        X_std = (X_raw - X_raw.mean(axis=0)) / (X_raw.std(axis=0) + 1e-9)
        print(f"  Feature matrix: {X_std.shape}")

        n_comp = min(args.n_components, X_std.shape[1])

        # VAMP-2 sweep (optional). --cv uses cross-validation; otherwise
        # block-bootstrap on training-set scores (overestimates real score).
        if args.vamp or args.vamp_only:
            vamp_lags = args.vamp_lags
            method = 'cv' if args.cv else 'boot'
            n_iters = (args.cv_folds if args.cv
                       else max(args.bootstrap, 20))
            vamp_cache = DATA_PATH / (
                f'vamp_{method}_{set_key.replace("frag_","")}_'
                f'{name}_b{n_iters}_n{n_comp}.npz')
            vlags_arr = np.asarray(vamp_lags)
            vamp_med = None
            if vamp_cache.exists() and not args.force_its:
                try:
                    d = np.load(vamp_cache)
                    if np.array_equal(d['lags'], vlags_arr):
                        print(f"  Loaded cached VAMP ({method}): "
                              f"{vamp_cache.name}")
                        vamp_med = d['vamp_median']
                        vamp_low = d['vamp_lower']
                        vamp_up = d['vamp_upper']
                except (KeyError, ValueError, OSError):
                    vamp_cache.unlink(missing_ok=True)
                    vamp_med = None
            if vamp_med is None:
                if args.replot_only:
                    print(f"  No VAMP cache for {name}; skipping")
                    if args.vamp_only:
                        continue
                elif args.cv:
                    print(f"  Cross-validated VAMP-2 at lags "
                          f"{[f'{l*0.01:.1f} ns' for l in vamp_lags]} "
                          f"({args.cv_folds}-fold CV)...")
                    vamp_med, vamp_low, vamp_up, _ = vamp_sweep_cv(
                        X_std, rep_boundaries, vamp_lags,
                        n_components=n_comp, n_folds=args.cv_folds)
                else:
                    print(f"  Bootstrap VAMP-2 at lags "
                          f"{[f'{l*0.01:.1f} ns' for l in vamp_lags]} "
                          f"({n_iters} resamples)...")
                    vamp_med, vamp_low, vamp_up = vamp_sweep_bootstrap(
                        X_std, rep_boundaries, vamp_lags,
                        n_components=n_comp, n_bootstrap=n_iters)
                if vamp_med is not None:
                    np.savez(vamp_cache,
                             lags=vlags_arr,
                             vamp_median=vamp_med,
                             vamp_lower=vamp_low,
                             vamp_upper=vamp_up,
                             n_components=n_comp,
                             method=method,
                             n_iters=n_iters)
                    print(f"  Cached VAMP: {vamp_cache.name}")
            if vamp_med is not None:
                plot_vamp(vamp_lags, vamp_med, vamp_low, vamp_up,
                          frame_dt_ns=0.01,
                          outpath=out_dir / f'{name}_vamp',
                          x_scale=args.x_scale, label=name)
                print(f"  VAMP-2 by lag (ns):  " + '  '.join(
                    f'{lag*0.01:.1f}={v:.3f}'
                    for lag, v in zip(vamp_lags, vamp_med)))

        if args.vamp_only:
            # Stash for comparison, skip ITS
            all_results[name] = ('vamp', vamp_med, vamp_low, vamp_up,
                                  vamp_lags)
            continue

        # ITS sweep — cache results to disk so replots are instant
        its_cache = DATA_PATH / (
            f'its_results_{set_key.replace("frag_","")}_'
            f'{name}_b{args.bootstrap}_n{n_comp}.npz')
        # Cache hit requires same lags
        lags_arr = np.asarray(args.its_lags)
        if (its_cache.exists() and not args.force_its
                and not args.force_its):
            try:
                d = np.load(its_cache)
                if (np.array_equal(d['lags'], lags_arr)
                        and d['ts_matrix'].shape[1] >= args.n_show):
                    print(f"  Loaded cached ITS: {its_cache.name}")
                    ts_med = d['ts_matrix']
                    ts_low = d['ts_lower']
                    ts_up = d['ts_upper']
                else:
                    raise ValueError('cache mismatch')
            except (KeyError, ValueError, OSError):
                # Cache invalid — recompute
                its_cache.unlink(missing_ok=True)
                ts_med = None

            if ts_med is None:
                pass  # fall through to compute
        else:
            ts_med = None

        if ts_med is None:
            if args.replot_only:
                print(f"  No cache for {name}; skipping (--replot-only)")
                continue
            if args.bootstrap > 0:
                print(f"  Bootstrap ITS ({args.bootstrap} resamples, "
                      f"{len(args.its_lags)} lags)...")
                ts_med, ts_low, ts_up = its_sweep_bootstrap(
                    X_std, rep_boundaries, args.its_lags,
                    n_components=n_comp, n_bootstrap=args.bootstrap)
            else:
                print(f"  Point estimate ITS (no bootstrap)...")
                ts_med, _ = its_sweep(X_std, rep_boundaries, args.its_lags,
                                       n_components=n_comp)
                ts_low = ts_med
                ts_up = ts_med

            # Save cache
            np.savez(its_cache,
                     lags=lags_arr,
                     ts_matrix=ts_med, ts_lower=ts_low, ts_upper=ts_up,
                     n_components=n_comp, n_bootstrap=args.bootstrap)
            print(f"  Cached ITS: {its_cache.name}")

        # Print top IC1-3 at each lag
        print(f"  {'Lag(ns)':<10} " + ' '.join(
            f'IC{i+1}(ns)         ' for i in range(min(3, n_comp))))
        for i, lag in enumerate(args.its_lags):
            parts = [f"{lag * 0.01:<10.2f}"]
            for m in range(min(3, n_comp)):
                med = ts_med[i, m]
                lo = ts_low[i, m]
                hi = ts_up[i, m]
                if np.isnan(med):
                    parts.append('    nan          ')
                else:
                    parts.append(f"{med:6.2f} [{lo:5.2f},{hi:5.2f}]")
            print('  ' + ' '.join(parts))

        # Plot
        plot_its(args.its_lags, ts_med, frame_dt_ns=0.01,
                 outpath=out_dir / name, n_show=args.n_show,
                 units='ps',
                 ts_lower=ts_low if args.bootstrap > 0 else None,
                 ts_upper=ts_up if args.bootstrap > 0 else None,
                 x_scale=args.x_scale,
                 title=f'ITS — {name} ({set_key})')

        all_results[name] = ('its', ts_med, ts_low, ts_up)
        print(f"  Saved: {out_dir / (name + '.png')}\n")

    # Comparison plot — handles ITS (IC1) and VAMP-2 modes
    if len(all_results) >= 2 and not args.vamp_only:
        print("━━━ Comparison plot (IC1 across featurizations) ━━━━━━━━━━")
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = plt.cm.tab10(np.linspace(0, 1, len(all_results)))
        lag_axis_ps = np.array(args.its_lags) * 0.01 * 1000  # frames → ps
        for i, (name, val) in enumerate(all_results.items()):
            if val[0] != 'its':
                continue
            _, ts_med, ts_low, ts_up = val
            ic1_med = ts_med[:, 0] * 1000  # ns → ps
            ic1_low = ts_low[:, 0] * 1000
            ic1_up = ts_up[:, 0] * 1000
            valid = ~np.isnan(ic1_med)
            if not valid.any():
                continue
            if args.bootstrap > 0:
                ax.fill_between(lag_axis_ps[valid], ic1_low[valid],
                                ic1_up[valid], color=colors[i], alpha=0.2)
            ax.plot(lag_axis_ps[valid], ic1_med[valid], 'o-',
                    color=colors[i], linewidth=2, markersize=7, label=name)

        # Trust region
        x_extend = np.linspace(0, lag_axis_ps.max() * 1.05, 50)
        ax.fill_between(x_extend, 0, x_extend, color='gray', alpha=0.1)
        ax.plot(x_extend, x_extend, 'k--', linewidth=1, alpha=0.4)

        ax.set_ylabel('IC1 timescale / ps', fontsize=12)
        ax.set_title(f'ITS comparison — IC1 across featurizations\n'
                     f'({set_key}, 95% CI shaded)', fontsize=13)
        if args.x_scale == 'linear':
            ax.set_xlim(0, lag_axis_ps.max() * 1.05)
            ax.set_xticks(lag_axis_ps)
            if lag_axis_ps.max() >= 1000:  # ≥ 1 ns
                ax.set_xticklabels([f'{v/1000:.1f}' for v in lag_axis_ps])
                ax.set_xlabel('Lag time (ns)', fontsize=12)
            else:
                ax.set_xticklabels([f'{v:.0f}' for v in lag_axis_ps])
                ax.set_xlabel('Lag time / ps', fontsize=12)
        else:
            ax.set_xscale('log')
            ax.set_xlabel('Lag time / ps', fontsize=12)
        ax.set_yscale('log')
        ax.legend(fontsize=10, loc='best')
        ax.grid(alpha=0.3, which='both')
        plt.tight_layout()
        fig.savefig(out_dir / 'comparison_IC1.png', dpi=200,
                     bbox_inches='tight')
        fig.savefig(out_dir / 'comparison_IC1.svg', bbox_inches='tight')
        plt.close()
        print(f"  Saved: {out_dir / 'comparison_IC1.png'}")

    # VAMP-2 comparison plot — bars grouped by lag, colored by featurization
    if (args.vamp or args.vamp_only) and len(all_results) >= 2:
        print("━━━ VAMP-2 comparison: features at each lag ━━━━━━━━━━━━━━━")

        # Collect per-featurization VAMP arrays at the same lags
        vamp_results = {}  # name -> (lags, med, lo, hi)
        vamp_lags = args.vamp_lags
        for name, val in all_results.items():
            if val[0] == 'vamp':
                _, vmed, vlow, vup, vlags = val
                vamp_results[name] = (vlags, vmed, vlow, vup)
            else:
                # Re-load and compute VAMP for this featurization
                feat = next((c for c in configs if c[0] == name), None)
                if feat is None:
                    continue
                _, features, ca_stride, add_rg = feat
                # Cache key encodes features + stride
                if features == 'ca':
                    feat_tag_cache = f'_ca{ca_stride}'
                elif features == 'com':
                    feat_tag_cache = ''
                else:
                    feat_tag_cache = f'_{features}{ca_stride}'
                rg_tag = '_rg' if add_rg else ''
                cache_npz = DATA_PATH / (
                    f'cluster_features_{set_key.replace("frag_","")}'
                    f'{feat_tag_cache}{rg_tag}.npz')
                if not cache_npz.exists():
                    continue
                cached = np.load(cache_npz, allow_pickle=True)
                Xr = cached['features']
                meta = list(cached['metadata'])
                if 'rep_boundaries' in cached:
                    rb = [tuple(b) for b in cached['rep_boundaries']]
                else:
                    rb = [(0, len(meta))]
                X_s = (Xr - Xr.mean(axis=0)) / (Xr.std(axis=0) + 1e-9)
                n_c = min(args.n_components, X_s.shape[1])
                vmed, vlow, vup = vamp_sweep_bootstrap(
                    X_s, rb, vamp_lags,
                    n_components=n_c,
                    n_bootstrap=max(args.bootstrap, 20))
                vamp_results[name] = (vamp_lags, vmed, vlow, vup)

        if vamp_results:
            # Grouped bar chart: x = lag, bars within group = featurizations
            feat_names = list(vamp_results.keys())
            n_feat = len(feat_names)
            n_lag = len(vamp_lags)
            colors = plt.cm.tab10(np.linspace(0, 1, n_feat))
            bar_width = 0.8 / n_feat

            fig, ax = plt.subplots(figsize=(max(8, 1.5 * n_lag * n_feat / 5), 6))
            x_center = np.arange(n_lag)

            for fi, name in enumerate(feat_names):
                _, vmed, vlow, vup = vamp_results[name]
                # Position offset for grouped bars
                offset = (fi - n_feat / 2 + 0.5) * bar_width
                xs = x_center + offset
                # Asymmetric error bars (CI bounds → distance from median)
                yerr_low = vmed - vlow
                yerr_up = vup - vmed
                bars = ax.bar(xs, vmed, bar_width,
                               yerr=[yerr_low, yerr_up],
                               color=colors[fi], alpha=0.85,
                               edgecolor='black', linewidth=0.5,
                               label=name,
                               error_kw=dict(elinewidth=1, capsize=3))
                # Value labels on bars
                for x, v in zip(xs, vmed):
                    if not np.isnan(v):
                        ax.text(x, v, f'{v:.2f}', ha='center', va='bottom',
                                fontsize=7)

            ax.set_xticks(x_center)
            ax.set_xticklabels([f'{lag * 0.01:.1f} ns' for lag in vamp_lags],
                                fontsize=11)
            ax.set_xlabel('Lag time', fontsize=12)
            ax.set_ylabel('VAMP-2 score', fontsize=12)
            ax.set_title(f'VAMP-2 score per featurization (grouped by lag)\n'
                         f'{set_key} — 95% CI as error bars '
                         '(higher = better)', fontsize=13)
            ax.legend(fontsize=10, loc='best',
                       title='Featurization')
            ax.grid(axis='y', alpha=0.3)
            ax.set_axisbelow(True)

            plt.tight_layout()
            fig.savefig(out_dir / 'comparison_VAMP2_by_lag.png', dpi=200,
                         bbox_inches='tight')
            fig.savefig(out_dir / 'comparison_VAMP2_by_lag.svg',
                         bbox_inches='tight')
            plt.close()
            print(f"  Saved: {out_dir / 'comparison_VAMP2_by_lag.png'}")

            # Also: per-lag panels (one subplot per lag, features as bars)
            fig, axes = plt.subplots(1, n_lag,
                                       figsize=(4 * n_lag, 5), sharey=True,
                                       squeeze=False)
            for li, lag in enumerate(vamp_lags):
                ax = axes[0, li]
                meds = [vamp_results[n][1][li] for n in feat_names]
                lows = [vamp_results[n][2][li] for n in feat_names]
                ups = [vamp_results[n][3][li] for n in feat_names]
                yerr_low = [m - l for m, l in zip(meds, lows)]
                yerr_up = [u - m for u, m in zip(ups, meds)]
                xs = np.arange(n_feat)
                bars = ax.bar(xs, meds, color=colors[:n_feat],
                               yerr=[yerr_low, yerr_up],
                               edgecolor='black', linewidth=0.5,
                               error_kw=dict(elinewidth=1.2, capsize=4))
                for x, v in zip(xs, meds):
                    if not np.isnan(v):
                        ax.text(x, v, f'{v:.2f}', ha='center', va='bottom',
                                fontsize=8, fontweight='bold')
                ax.set_xticks(xs)
                ax.set_xticklabels(feat_names, rotation=30, ha='right',
                                    fontsize=9)
                ax.set_title(f'Lag {lag * 0.01:.1f} ns', fontsize=12)
                ax.grid(axis='y', alpha=0.3)
                ax.set_axisbelow(True)
                if li == 0:
                    ax.set_ylabel('VAMP-2 score', fontsize=11)

            fig.suptitle(f'VAMP-2 across featurizations ({set_key})\n'
                         '95% CI as error bars — higher = better',
                         fontsize=13)
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            fig.savefig(out_dir / 'comparison_VAMP2_panels.png', dpi=200,
                         bbox_inches='tight')
            fig.savefig(out_dir / 'comparison_VAMP2_panels.svg',
                         bbox_inches='tight')
            plt.close()
            print(f"  Saved: {out_dir / 'comparison_VAMP2_panels.png'}")

    print(f"\nDone. All plots in: {out_dir}/")
    return 0


if __name__ == '__main__':
    sys.exit(main())
