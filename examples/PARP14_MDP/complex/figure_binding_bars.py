#!/usr/bin/env python3
"""Bar figures: bound fraction, contact lifetime, encounter rate -- and what tau means.

Reads analysis/tau_full_resolution.csv (stride-1 lifetimes). It must be the
full-resolution table, not episode_stats.csv: tau measured at stride 10 is
pinned to the sampling interval rather than to the physics, so plotting that
version would put an artifact on an axis.

Colours: separated vs docked, #4A3AA7 / #C8481A (validated all-pairs,
CVD dE 25.4, normal-vision dE 31.7, both above 3:1 contrast).
"""
import sys
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

SEP, DOCK = '#4A3AA7', '#C8481A'
INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'muted': '#898781',
       'gridline': '#E1E0D9'}
NICE = {'parp14': 'PARP14', 'parp9': 'PARP9', 'dtx3l': 'DTX3L'}


def pretty(pair):
    a, b = pair.split('-')
    return f'{NICE.get(a, a)}–{NICE.get(b, b)}'


def load():
    f = OUT / 'tau_full_resolution.csv'
    if not f.is_file():
        sys.exit(f'{f} not found -- run tau_full_resolution.py first.')
    d = pd.read_csv(f)
    d['arm'] = np.where(d.set.str.endswith('_docked'), 'docked', 'separated')
    d['base'] = d.set.str.replace('_docked', '', regex=False)
    d['label'] = d.base + '\n' + d.pair.map(pretty)
    return d


def grouped(ax, d, col, err=None, ylab='', title='', logy=False):
    """Paired separated/docked bars, one group per (set, pair)."""
    keys = (d.groupby('label')[col].max().sort_values(ascending=False).index.tolist())
    x = np.arange(len(keys)); w = 0.38
    for k, (arm, colr) in enumerate((('separated', SEP), ('docked', DOCK))):
        sub = d[d.arm == arm].set_index('label').reindex(keys)
        vals = sub[col].values
        e = sub[err].values if err else None
        ax.bar(x + (k - 0.5) * w, np.nan_to_num(vals), w, color=colr, label=arm,
               yerr=e if err else None, capsize=2,
               error_kw=dict(lw=0.8, ecolor=INK['secondary']))
    ax.set_xticks(x); ax.set_xticklabels(keys, fontsize=6, rotation=45, ha='right')
    ax.set_ylabel(ylab, fontsize=8)
    ax.set_title(title, fontsize=9, loc='left')
    if logy:
        ax.set_yscale('log')
    ax.grid(axis='y', color=INK['gridline'], lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)
    return keys


def main():
    ap = ArgumentParser(); ap.add_argument('--demo', default=None,
                                           help='csv of the sampling-dependence demo')
    a = ap.parse_args()
    d = load()

    fig, axes = plt.subplots(3, 1, figsize=(11, 12))
    grouped(axes[0], d, 'bound_frac', 'bound_frac_ci95',
            'fraction of frames in contact',
            'Bound fraction — how much of the run the chains spend touching '
            '(min CA-CA < 1.0 nm; error bars 95% CI over replicates)')
    axes[0].legend(fontsize=7, frameon=False, ncol=2)

    grouped(axes[1], d, 'tau_survival_ns', None, 'contact lifetime τ (ns)',
            'Contact lifetime — survival-corrected mean, measured at the full '
            '0.05 ns trajectory resolution')
    axes[1].legend(fontsize=7, frameon=False, ncol=2)

    grouped(axes[2], d, 'episodes_per_us', None, 'association events per µs',
            'Encounter rate — how often a contact forms. High rate with short τ '
            'is the signature of a diffusive encounter complex')
    axes[2].legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / '11_binding_bars.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('11_binding_bars.png')

    # ---- what tau means -------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    if a.demo and Path(a.demo).is_file():
        dm = pd.read_csv(a.demo)
        ax.plot(dm.sampling_ns, dm.tau_mean_ns, 'o-', color=SEP, label='τ mean')
        ax2 = ax.twiny() if False else None
        ax.plot(dm.sampling_ns, dm.bound_frac * 10, 's--', color=DOCK,
                label='bound fraction (×10)')
        ax.set_xscale('log'); ax.set_xlabel('sampling interval (ns)')
        ax.set_ylabel('value')
        ax.legend(fontsize=7, frameon=False)
    ax.set_title('What τ means: it moves with the sampling interval,\n'
                 'bound fraction does not', fontsize=9, loc='left')
    ax.grid(color=INK['gridline'], lw=0.6); ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    ax.scatter(d.episodes_per_us, d.tau_survival_ns, s=34, color=SEP,
               edgecolor='white', linewidth=0.6)
    for _, r in d.iterrows():
        ax.annotate(r.base, (r.episodes_per_us, r.tau_survival_ns), fontsize=5,
                    xytext=(3, 3), textcoords='offset points', color=INK['secondary'])
    ax.set_xlabel('association events per µs'); ax.set_ylabel('contact lifetime τ (ns)')
    ax.set_title('Every set sits in the same corner: contacts form often\n'
                 'and break fast — encounter complexes, not complexes',
                 fontsize=9, loc='left')
    ax.grid(color=INK['gridline'], lw=0.6); ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / '12_what_tau_means.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('12_what_tau_means.png')


if __name__ == '__main__':
    main()
