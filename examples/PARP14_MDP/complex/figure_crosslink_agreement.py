#!/usr/bin/env python3
"""The campaign's headline figure: crosslink agreement across every model.

Reads analysis/crosslink_agreement_all.csv (score_all.py).

Panel A  per-frame % within 30 A, grouped by crosslink subset. The subsets are
         nested and ordered by evidential weight (ALL 36 -> reproducible 21 ->
         all-4-replicate 8), so a model that is actually right should rise to
         the right within each group.
Panel B  ensemble coverage: % of crosslinks satisfied at all, against % reaching
         1% of frames. The gap between them is the "contains it but does not
         weight it" signature -- high coverage, no occupancy.

Log x on panel A: the models span 0.00% to 88%, four orders of magnitude, and a
linear axis renders everything below 1% as an invisible stub.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'grid': '#E1E0D9'}
# validated 3-hue set, ordered by evidential weight
SUB = [('ALL (36)', '#4A3AA7'), ('reproducible (21)', '#C8481A'),
       ('all-4-replicate (8)', '#117C60')]
FLOOR = 0.004      # so exact zeros are visible on a log axis


def main():
    d = pd.read_csv(OUT / 'crosslink_agreement_all.csv')
    order = ['CALVADOS unrestrained (separated)', 'CALVADOS unrestrained (docked)',
             'CALVADOS ternary (separated)', 'CALVADOS ternary (docked)',
             'HyRes unrestrained (separated)', 'HyRes matched (separated)',
             'CALVADOS +3 XL Go k=15', 'CALVADOS +3 XL harmonic k=5',
             'CALVADOS +3 XL harmonic k=20', 'CALVADOS +3 XL harmonic k=100',
             'HyRes matched (docked)']
    order = [o for o in order if o in set(d.model)]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4),
                             gridspec_kw={'width_ratios': [1.55, 1]})

    ax = axes[0]
    y = np.arange(len(order)); h = 0.26
    for k, (sub, col) in enumerate(SUB):
        s = d[d.subset == sub].set_index('model').reindex(order)
        v = s.per_frame_30A.fillna(0).values
        ax.barh(y + (1 - k) * h, np.maximum(v, FLOOR), height=h, color=col,
                label=sub, zorder=3)
        for yy, val in zip(y + (1 - k) * h, v):
            ax.text(max(val, FLOOR) * 1.25, yy, f'{val:.2f}' if val < 10 else f'{val:.0f}',
                    va='center', fontsize=5.8, color=INK['secondary'])
    ax.set_xscale('log'); ax.set_xlim(FLOOR * 0.7, 400)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=7.5); ax.invert_yaxis()
    ax.set_xlabel('% of frames with the crosslink within 30 Å  (log scale)', fontsize=8.5)
    ax.set_title('A  Crosslink agreement, by evidential weight of the subset',
                 fontsize=9.5, loc='left', weight='bold', pad=10)
    ax.legend(fontsize=7, frameon=False, loc='upper center',
              bbox_to_anchor=(0.5, -0.13), ncol=3)
    ax.grid(axis='x', color=INK['grid'], lw=0.6); ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    s = d[d.subset == 'reproducible (21)'].set_index('model').reindex(order)
    ax.scatter(s.ever_30A, s.pct1_30A, s=70, color='#4A3AA7',
               edgecolor='white', linewidth=0.8, zorder=3)
    lim = [-4, 108]
    ax.plot(lim, lim, ls='--', lw=0.9, color=INK['secondary'], zorder=1)
    # models with zero occupancy all land on y=0 and their labels collide;
    # fan those out vertically with leader lines instead of overprinting
    zero = s[s.pct1_30A <= 0.01].sort_values('ever_30A')
    nz = s[s.pct1_30A > 0.01]
    for m, r in nz.iterrows():
        short = (m.replace('CALVADOS ', 'CAL ').replace('HyRes ', 'HyR ')
                  .replace('unrestrained', 'unrestr').replace('harmonic ', 'h'))
        ax.annotate(short, (r.ever_30A, r.pct1_30A), fontsize=6,
                    xytext=(6, 5), textcoords='offset points', color=INK['secondary'])
    for i, (m, r) in enumerate(zero.iterrows()):
        short = (m.replace('CALVADOS ', 'CAL ').replace('HyRes ', 'HyR ')
                  .replace('unrestrained', 'unrestr').replace('harmonic ', 'h'))
        ax.annotate(short, (r.ever_30A, r.pct1_30A), fontsize=6,
                    xytext=(-2, 10 + 11*i), textcoords='offset points',
                    color=INK['secondary'], ha='right',
                    arrowprops=dict(arrowstyle='-', lw=0.5, color=INK['secondary']))
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('% of crosslinks satisfied at all', fontsize=8.5)
    ax.set_ylabel('% satisfied in ≥1% of frames', fontsize=8.5)
    ax.set_title('B  Coverage vs occupancy  (reproducible subset)\n'
                 'far below the diagonal = contains the geometry, does not weight it',
                 fontsize=9.5, loc='left', weight='bold', pad=10)
    ax.grid(color=INK['grid'], lw=0.6); ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('PARP9–DTX3L simulations vs the 51 published BS3 crosslinks '
                 '(Ashok et al. 2022, Biochem J 479:289)',
                 fontsize=10, x=0.005, y=1.02, ha='left', color=INK['primary'])
    fig.tight_layout()
    fig.savefig(FIG / '14_crosslink_agreement.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('14_crosslink_agreement.png')


if __name__ == '__main__':
    main()
