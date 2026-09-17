#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Writer-vs-eraser dominance prediction for PARP14 constructs.

QUESTION
--------
For a given PARP14 construct, is the WRITER (ART, ADP-ribosyltransferase) or the
ERASER (MD1, ADP-ribose hydrolase) catalytic site more able to engage substrate?

WHAT THIS IS -- AND IS NOT
--------------------------
CALVADOS is a coarse-grained, implicit-solvent, one-bead-per-residue model. It
has NO chemistry: no NAD+, no ADP-ribose, no transition states, no rate
constants. So this CANNOT predict catalytic rates.

What it CAN do is score STERIC OPPORTUNITY: across the conformational ensemble,
how much of the approach space around each catalytic pocket is left open by the
protein's own chain. A pocket that is buried by its own domains cannot be
engaged by a substrate no matter how good its chemistry; a pocket that is
solvent-exposed at least has the opportunity. Treat the output as a
hypothesis-generating RANKING over constructs, not a calibrated activity
prediction. See --help and the README notes for the calibration caveat.

METRIC
------
Primary: SAA (solid-angle accessibility) from analyze_accessibility.py -- the
fraction of 200 outward rays from the pocket COM that escape to 5 nm without
passing within 0.5 nm of a non-pocket bead. 0 = fully occluded, 1 = fully open.

Dominance score, for constructs carrying BOTH sites:

    D = (SAA_ART - SAA_MD1) / (SAA_ART + SAA_MD1)

    D > 0  -> writer site more sterically available  (writer-leaning)
    D < 0  -> eraser site more sterically available  (eraser-leaning)
    D ~ 0  -> balanced

D is a normalized contrast rather than a raw difference so constructs whose
sites are both open and both closed stay comparable.

PAIRED +/-ART ANALYSIS
----------------------
The simulation set contains matched construct pairs differing only by whether
ART is present (e.g. kh1_art_full vs kh1_wwe_full). Comparing MD1/MD2/MD3
accessibility across such a pair isolates how much the writer domain occludes
the reader/eraser sites -- a steric-coupling measurement that single constructs
cannot give.

Usage:
    python predict_enzyme_dominance.py                  # all cached sets
    python predict_enzyme_dominance.py --metric cone    # score on cone angle
    python predict_enzyme_dominance.py --set fl norrm core
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir

NPZ = os.path.join(DATA_PATH, 'accessibility_stats.npz')

WRITER_SITE = 'ART'   # ADP-ribosyltransferase  -- adds the mark
ERASER_SITE = 'MD1'   # macrodomain 1 hydrolase -- removes the mark

# Construct pairs differing ONLY by presence of the ART (writer) domain.
# (with_ART, without_ART, protocol_matched).
#
# All eight pairs are a clean ART-only difference in DOMAIN CONTENT (verified
# against sim_registry units). They are NOT all comparable, though: a pair is
# only interpretable as "the steric effect of ART" when both members were run
# under the SAME simulation protocol. The two fl pairs are not -- `fl` and
# `fl_optimized` are ~3.5k-frame runs under the original/optimized restraint
# schemes, while `fl_wwe_full_go` is a 100k-frame contiguous Go-model run, so
# any difference there confounds ART removal with restraint scheme and run
# length. They are kept (the contrast is still worth eyeballing) but marked
# unmatched, excluded from the headline figure, and flagged in the report.
ART_PAIRS = [
    ('fl',               'fl_wwe_full_go',    False),
    ('fl_optimized',     'fl_wwe_full_go',    False),
    ('norrm',            'noart',             True),
    ('core_full_go',     'core_wwe_full_go',  True),
    ('kh1_art_full',     'kh1_wwe_full',      True),
    ('mka_full',         'mka_wwe_full',      True),
    ('md2_art_full',     'md2_wwe_full',      True),
    ('md3_art_full',     'md3_wwe_full',      True),
]

METRIC_LABEL = {
    'saa':     'Solid-angle accessibility (fraction of open approach directions)',
    'cone':    'Widest unobstructed approach cone (degrees)',
    'density': 'Shell density of blocking beads (lower = more accessible)',
}


# ============================================================
# Data access
# ============================================================

def load_cache():
    """set_key -> site -> metric -> per-replicate array, from the cached npz."""
    if not os.path.isfile(NPZ):
        sys.exit(f"ERROR: {NPZ} not found.\n"
                 f"Run the accessibility backfill first:\n"
                 f"    ./backfill_accessibility.sh\n"
                 f"or for one set:\n"
                 f"    python analyze_accessibility.py --set <name> --stride 50")
    out = {}
    with np.load(NPZ, allow_pickle=True) as d:
        for key in d.files:
            # keys look like  <set_key>_<SITE>_<metric>
            set_key, site, metric = key.rsplit('_', 2)
            out.setdefault(set_key, {}).setdefault(site, {})[metric] = np.asarray(d[key])
    return out


def label_for(set_key):
    """Short, plot-friendly label. sim_registry labels carry long parentheticals.

    Used in the printed report, where the descriptive text is helpful. Figures
    use the set key instead (see `axis_label_for`): the trimmed labels are not
    unique -- `fl` -> "Full Length" and `fl_optimized` -> "FL" are impossible to
    tell apart on an axis, and the set key is what every command takes anyway.
    """
    if set_key in reg.SETS:
        lab = reg.SETS[set_key]['label']
        return lab.split('(')[0].strip().rstrip(',')
    return set_key


def axis_label_for(set_key):
    """Unambiguous figure-axis label: the set key, as used on the CLI."""
    return set_key


def bootstrap_ci(values, statfn=np.mean, n_boot=10000, alpha=0.05, seed=0):
    """Percentile bootstrap CI over replicates. Returns (lo, hi)."""
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(values), size=(n_boot, len(values)))
    stats = statfn(values[idx], axis=1)
    return (float(np.percentile(stats, 100 * alpha / 2)),
            float(np.percentile(stats, 100 * (1 - alpha / 2))))


def dominance(writer_vals, eraser_vals):
    """Normalized writer-vs-eraser contrast, computed replicate-wise.

    Replicates are paired (same trajectory), so pairing before taking the
    contrast preserves the within-replicate correlation instead of dividing
    two independently-averaged means.
    """
    w = np.asarray(writer_vals, dtype=float)
    e = np.asarray(eraser_vals, dtype=float)
    n = min(len(w), len(e))
    w, e = w[:n], e[:n]
    denom = w + e
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.where(denom > 0, (w - e) / denom, np.nan)
    return d


# ============================================================
# Analysis
# ============================================================

def build_table(cache, metric, sets_filter=None):
    rows = []
    for set_key in sorted(cache):
        if sets_filter and set_key not in sets_filter:
            continue
        sites = cache[set_key]
        row = {'set': set_key, 'label': label_for(set_key)}
        for site in reg.SITE_NAMES:
            vals = sites.get(site, {}).get(metric)
            row[f'{site}_mean'] = float(np.mean(vals)) if vals is not None else np.nan
            row[f'{site}_sd'] = float(np.std(vals)) if vals is not None else np.nan
            row[f'{site}_n'] = int(len(vals)) if vals is not None else 0

        has_w = WRITER_SITE in sites and metric in sites[WRITER_SITE]
        has_e = ERASER_SITE in sites and metric in sites[ERASER_SITE]
        if has_w and has_e:
            d = dominance(sites[WRITER_SITE][metric], sites[ERASER_SITE][metric])
            d = d[~np.isnan(d)]
            row['D'] = float(np.mean(d)) if len(d) else np.nan
            row['D_sd'] = float(np.std(d)) if len(d) else np.nan
            lo, hi = bootstrap_ci(d)
            row['D_lo'], row['D_hi'] = lo, hi
            row['n_rep'] = len(d)
            row['call'] = ('writer-leaning' if lo > 0 else
                           'eraser-leaning' if hi < 0 else 'balanced')
        else:
            row['D'] = row['D_sd'] = row['D_lo'] = row['D_hi'] = np.nan
            row['n_rep'] = 0
            missing = [s for s, present in
                       ((WRITER_SITE, has_w), (ERASER_SITE, has_e)) if not present]
            row['call'] = f'n/a (no {"/".join(missing)})'
        rows.append(row)
    return rows


def paired_art_effect(cache, metric):
    """For each matched +/-ART pair, the change in each non-ART site's metric."""
    out = []
    for with_art, without_art, matched in ART_PAIRS:
        if with_art not in cache or without_art not in cache:
            continue
        for site in ('MD1', 'MD2', 'MD3', 'WWE'):
            a = cache[with_art].get(site, {}).get(metric)
            b = cache[without_art].get(site, {}).get(metric)
            if a is None or b is None:
                continue
            ma, mb = float(np.mean(a)), float(np.mean(b))
            # delta > 0  => site MORE accessible when ART is present
            # delta < 0  => ART occludes this site
            delta = ma - mb
            # NOTE: pooled_sd is the spread across replicates of the SAME
            # construct, and every replicate starts from the same AF3/AF2 model,
            # so that spread is tiny (~1e-3 SAA). Dividing by it inflates
            # `effect_size` into the tens -- it is a precision ratio, NOT a
            # Cohen's d, and must not be read as one. `delta` is the
            # interpretable number; effect_size only ranks confidence.
            pooled_sd = float(np.sqrt((np.var(a) + np.var(b)) / 2))
            out.append({
                'with_art': with_art, 'without_art': without_art, 'site': site,
                'protocol_matched': matched,
                'mean_with': ma, 'mean_without': mb, 'delta': delta,
                'pooled_sd': pooled_sd,
                'effect_size': delta / pooled_sd if pooled_sd > 0 else np.nan,
            })
    return out


# ============================================================
# Figures
# ============================================================

def fig_dominance_ranking(rows, metric, figdir):
    scored = [r for r in rows if not np.isnan(r['D'])]
    if not scored:
        print("  (no construct carries both ART and MD1 -- skipping ranking figure)")
        return
    scored.sort(key=lambda r: r['D'])
    labels = [axis_label_for(r['set']) for r in scored]
    vals = np.array([r['D'] for r in scored])
    lo = np.array([r['D_lo'] for r in scored])
    hi = np.array([r['D_hi'] for r in scored])
    colors = ['#f58231' if v > 0 else '#e6194b' for v in vals]

    fig, ax = plt.subplots(figsize=(9, 0.45 * len(scored) + 2.6))
    y = np.arange(len(scored))
    ax.barh(y, vals, color=colors, edgecolor='black', linewidth=0.6, zorder=3)
    ax.errorbar(vals, y, xerr=[vals - lo, hi - vals], fmt='none',
                ecolor='black', capsize=3, lw=1.1, zorder=4)
    ax.axvline(0, color='black', lw=1.2, zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(f'D  =  (ART - MD1) / (ART + MD1)     [{metric}]', fontsize=10)
    ax.set_title('Writer (ART) vs eraser (MD1) steric availability\n'
                 'positive = writer site more open;  negative = eraser site more open',
                 fontsize=11)
    ax.grid(axis='x', alpha=0.3, zorder=0)
    lim = max(0.05, float(np.nanmax(np.abs(np.concatenate([lo, hi])))) * 1.25)
    ax.set_xlim(-lim, lim)
    ax.text(0.99, 0.01, 'bars = mean over replicates, whiskers = 95% bootstrap CI',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=7.5, color='#555')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'writer_eraser_dominance.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: writer_eraser_dominance.png/.svg")


def fig_site_heatmap(rows, metric, figdir):
    sets = [r for r in rows if any(r[f'{s}_n'] for s in reg.SITE_NAMES)]
    if not sets:
        return
    sets.sort(key=lambda r: r['set'])
    mat = np.array([[r[f'{s}_mean'] for s in reg.SITE_NAMES] for r in sets], dtype=float)

    fig, ax = plt.subplots(figsize=(1.15 * len(reg.SITE_NAMES) + 4.5,
                                    0.42 * len(sets) + 2.4))
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#dddddd')
    im = ax.imshow(np.ma.masked_invalid(mat), aspect='auto', cmap=cmap)
    ax.set_xticks(range(len(reg.SITE_NAMES)))
    ax.set_xticklabels(reg.SITE_NAMES, fontsize=11)
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels([axis_label_for(r['set']) for r in sets], fontsize=9)
    for i in range(len(sets)):
        for j in range(len(reg.SITE_NAMES)):
            v = mat[i, j]
            if np.isnan(v):
                ax.text(j, i, 'absent', ha='center', va='center', fontsize=7, color='#777')
            else:
                # contrast-aware label colour
                norm = (v - np.nanmin(mat)) / max(1e-9, np.nanmax(mat) - np.nanmin(mat))
                ax.text(j, i, f'{v:.3f}' if metric != 'cone' else f'{v:.0f}',
                        ha='center', va='center', fontsize=8,
                        color='white' if norm < 0.55 else 'black')
    ax.set_title(f'Per-site accessibility across constructs\n{METRIC_LABEL[metric]}',
                 fontsize=11)
    plt.colorbar(im, ax=ax, shrink=0.75)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'site_accessibility_heatmap.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: site_accessibility_heatmap.png/.svg")


def fig_art_effect(pairs, metric, figdir):
    # Only protocol-matched pairs belong in this figure: for the unmatched ones
    # the delta mixes ART removal with a different restraint scheme/run length.
    pairs = [p for p in pairs if p['protocol_matched']]
    if not pairs:
        print("  (no protocol-matched +/-ART construct pairs cached "
              "-- skipping ART-effect figure)")
        return
    labels = [f"{axis_label_for(p['with_art'])}\nvs {axis_label_for(p['without_art'])}"
              f"  [{p['site']}]" for p in pairs]
    deltas = np.array([p['delta'] for p in pairs])
    colors = ['#2ca02c' if d > 0 else '#d62728' for d in deltas]

    fig, ax = plt.subplots(figsize=(9.5, 0.42 * len(pairs) + 2.8))
    y = np.arange(len(pairs))
    ax.barh(y, deltas, color=colors, edgecolor='black', linewidth=0.6, zorder=3)
    ax.axvline(0, color='black', lw=1.2, zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel(f'Δ {metric}  (with ART  −  without ART)', fontsize=10)
    ax.set_title('Steric effect of the writer (ART) domain on the other sites\n'
                 'negative = ART occludes that site;  positive = ART opens it up\n'
                 '(protocol-matched construct pairs only)',
                 fontsize=11)
    ax.grid(axis='x', alpha=0.3, zorder=0)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'art_domain_steric_effect.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: art_domain_steric_effect.png/.svg")


# ============================================================
# Output
# ============================================================

def write_csvs(rows, pairs, metric):
    os.makedirs(DATA_PATH, exist_ok=True)
    main = os.path.join(DATA_PATH, f'enzyme_dominance_{metric}.csv')
    cols = (['set', 'label', 'n_rep', 'D', 'D_sd', 'D_lo', 'D_hi', 'call'] +
            [f'{s}_{k}' for s in reg.SITE_NAMES for k in ('mean', 'sd', 'n')])
    with open(main, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in sorted(rows, key=lambda r: (np.isnan(r['D']), r['D'])):
            w.writerow(r)
    print(f"  Saved: {os.path.relpath(main, CWD)}")

    if pairs:
        pth = os.path.join(DATA_PATH, f'art_steric_effect_{metric}.csv')
        with open(pth, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()))
            w.writeheader()
            w.writerows(pairs)
        print(f"  Saved: {os.path.relpath(pth, CWD)}")


def print_report(rows, pairs, metric):
    scored = [r for r in rows if not np.isnan(r['D'])]
    unscored = [r for r in rows if np.isnan(r['D'])]

    print(f"\n{'=' * 94}")
    print(f"WRITER (ART) vs ERASER (MD1) STERIC DOMINANCE   [metric: {metric}]")
    print('=' * 94)
    if scored:
        print(f"{'Construct':<34}{'n':>3}  {'ART':>7} {'MD1':>7}  "
              f"{'D':>7}  {'95% CI':>17}  {'call':<16}")
        print('-' * 94)
        for r in sorted(scored, key=lambda r: r['D']):
            print(f"{r['label'][:33]:<34}{r['n_rep']:>3}  "
                  f"{r[f'{WRITER_SITE}_mean']:>7.3f} {r[f'{ERASER_SITE}_mean']:>7.3f}  "
                  f"{r['D']:>+7.3f}  [{r['D_lo']:>+6.3f},{r['D_hi']:>+6.3f}]  {r['call']:<16}")
    else:
        print("  No cached construct carries BOTH the ART and MD1 sites.")
    if unscored:
        print(f"\n  Not scoreable ({len(unscored)}): "
              + ', '.join(f"{r['label'][:26]} [{r['call']}]" for r in unscored[:8])
              + (' ...' if len(unscored) > 8 else ''))

    if pairs:
        print(f"\n{'=' * 94}")
        print("STERIC EFFECT OF THE ART DOMAIN  (matched +/-ART construct pairs)")
        print('=' * 94)
        print(f"{'pair':<46}{'site':<6}{'with':>8}{'without':>9}{'delta':>9}{'d/SD':>9}")
        print('-' * 94)
        print("  (delta is the interpretable quantity. d/SD divides by the "
              "replicate spread of\n   identical starting structures, so it runs "
              "to the tens -- a precision ratio,\n   not a Cohen's d.)")
        matched = [p for p in pairs if p['protocol_matched']]
        unmatched = [p for p in pairs if not p['protocol_matched']]
        for p in sorted(matched, key=lambda p: p['delta']):
            pair_lbl = f"{label_for(p['with_art'])[:20]} vs {label_for(p['without_art'])[:20]}"
            eff = p['effect_size']
            print(f"{pair_lbl:<46}{p['site']:<6}{p['mean_with']:>8.3f}"
                  f"{p['mean_without']:>9.3f}{p['delta']:>+9.3f}"
                  f"{(f'{eff:>+9.1f}' if np.isfinite(eff) else '        -')}")
        if unmatched:
            print(f"\n  NOT protocol-matched -- ART removal is confounded with a different")
            print(f"  restraint scheme / run length. Reported for reference only, and")
            print(f"  excluded from art_domain_steric_effect.png:")
            for p in sorted(unmatched, key=lambda p: p['delta']):
                pair_lbl = (f"{label_for(p['with_art'])[:20]} vs "
                            f"{label_for(p['without_art'])[:20]}")
                print(f"    {pair_lbl:<44}{p['site']:<6}{p['mean_with']:>8.3f}"
                      f"{p['mean_without']:>9.3f}{p['delta']:>+9.3f}")

    print(f"\n{'=' * 94}")
    print("HOW TO READ THIS")
    print('=' * 94)
    print("  D > 0        writer (ART) pocket has more open approach space than the eraser (MD1)")
    print("  D < 0        eraser (MD1) pocket is the more available one")
    print("  CI excludes 0 -> the lean is consistent across replicates")
    print()
    print("  This ranks STERIC OPPORTUNITY, not catalytic rate. CALVADOS is coarse-grained")
    print("  with no NAD+, no ADP-ribose and no chemistry, so it cannot say how fast either")
    print("  reaction runs -- only whether the pocket is physically reachable in the ensemble.")
    print("  Replicate spread is small because all replicates start from the same AF3/AF2")
    print("  models: that is precision, not independent validation of the underlying structure.")
    print()
    print("  To turn this into a calibrated prediction you need an experimental activity")
    print("  readout (in vitro ADPr transfer / hydrolysis) for at least a few constructs.")
    print('=' * 94)


def main():
    ap = ArgumentParser(description='Predict writer-vs-eraser steric dominance '
                                    'across PARP14 constructs.')
    ap.add_argument('--metric', default='saa', choices=('saa', 'cone', 'density'),
                    help='Accessibility metric to score on (default: saa)')
    ap.add_argument('--set', nargs='+', default=None,
                    help='Restrict to these set keys (default: every cached set)')
    args = ap.parse_args()

    cache = load_cache()
    print(f"  Loaded {len(cache)} cached set(s) from "
          f"{os.path.relpath(NPZ, CWD)}: {', '.join(sorted(cache))}")

    rows = build_table(cache, args.metric, sets_filter=args.set)
    pairs = paired_art_effect(cache, args.metric)

    # This is a cross-construct synthesis, so `sims` would be every cached set.
    # Joining ~20 set keys produces a ~250-char directory name that brushes the
    # filesystem limit and breaks outright as sets are added, so fall back to
    # the shared dated layout once the comparison spans more than a few sets.
    sims = sorted(cache)
    figdir = str(_get_fig_dir('03_accessibility', subname='enzyme_dominance',
                              sims=sims if len(sims) <= 3 else None))
    fig_dominance_ranking(rows, args.metric, figdir)
    fig_site_heatmap(rows, args.metric, figdir)
    fig_art_effect(pairs, args.metric, figdir)
    write_csvs(rows, pairs, args.metric)
    print_report(rows, pairs, args.metric)
    print(f"\n  Figures -> {os.path.relpath(figdir, CWD)}")


if __name__ == '__main__':
    main()
