#!/usr/bin/env python3
"""Bar-graph versions of the campaign's numeric tables.

Same numbers as the tables in the write-up, in the form that makes the ranking
readable at a glance.

16_campaign_bars   A sampling per run (engine-coloured)
                   B crosslink agreement, per-frame % within 30 A, by subset
17_result_bars     C ensemble coverage: satisfied ever vs in >=1% of frames
                   D the transfer result -- restraining 3 crosslinks moves 43
                   E HyRes fold stability and association
                   F crosslinks satisfied over single trajectories (faceted)

Sources: analysis/campaign_inventory.csv, analysis/crosslink_agreement_all.csv,
analysis/pymol/*_timecourse.csv. Panels D and E mirror hand-built tables in the
write-up; their numbers are literals here with the provenance noted inline.

Log x on panel B: models span 0.00% to 88%, and a linear axis renders everything
under 1% as an invisible stub.

Panel F is faceted rather than a four-series grouped chart: four categorical hues
cannot clear the all-pairs CVD check, so identity moves to the facet title and
every facet shares one hue.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)

INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'grid': '#E1E0D9'}
# every set below clears validate_palette.py --mode light --pairs all
SUBSET = [('ALL (36)', '#4A3AA7'), ('reproducible (21)', '#C8481A'),
          ('all-4-replicate (8)', '#117C60')]
ENGINE = {'CALVADOS': '#4A3AA7', 'HyRes': '#C8481A'}
PAIR2 = ['#1E6FB8', '#117C60']                   # two-series panels
EVER, PCT1 = '#8E7BD6', '#4A3AA7'                # validated ordinal 2-step
ONE = '#4A3AA7'                                  # single hue for the facets
FLOOR = 0.004                                    # so exact zeros show on a log axis


def bare(ax, xgrid=True, ygrid=False):
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(INK['grid'])
    ax.tick_params(colors=INK['secondary'], length=0, labelsize=9)
    if xgrid:
        ax.xaxis.grid(True, color=INK['grid'], lw=.8)
    if ygrid:
        ax.yaxis.grid(True, color=INK['grid'], lw=.8)
    ax.set_axisbelow(True)


def title(ax, tag, text, sub=None):
    """Title above the axes, subtitle under it -- never overlapping."""
    ax.text(0, 1.085 if sub else 1.035, f'{tag}   {text}', transform=ax.transAxes,
            ha='left', va='bottom', fontsize=11.5, color=INK['primary'],
            fontweight='semibold')
    if sub:
        ax.text(0, 1.028, sub, transform=ax.transAxes, ha='left', va='bottom',
                fontsize=9, color=INK['secondary'])


SHORT = {'CALVADOS unrestrained (separated)': 'unrestrained, separated',
         'CALVADOS unrestrained (docked)': 'unrestrained, docked',
         'CALVADOS ternary (separated)': 'ternary, separated',
         'CALVADOS ternary (docked)': 'ternary, docked',
         'CALVADOS +3 XL harmonic k=5': '+3 XL harmonic k=5',
         'CALVADOS +3 XL harmonic k=20': '+3 XL harmonic k=20',
         'CALVADOS +3 XL harmonic k=100': '+3 XL harmonic k=100',
         'CALVADOS +3 XL Go k=15': '+3 XL Gō k=15',
         'CALVADOS +3 XL Go tight r=1.5 k=50': '+3 XL Gō r=1.5 k=50',
         'CALVADOS +3 XL Go tight r=1.5 k=100': '+3 XL Gō r=1.5 k=100',
         'CALVADOS +3 XL Go tight r=1.2 k=50': '+3 XL Gō r=1.2 k=50',
         'CALVADOS +3 XL Go tight r=1.2 k=100': '+3 XL Gō r=1.2 k=100',
         'HyRes unrestrained (separated)': 'HyRes unrestrained',
         'HyRes matched (separated)': 'HyRes matched, separated',
         'HyRes matched (docked)': 'HyRes matched, docked'}


def fig16():
    inv = pd.read_csv(OUT / 'campaign_inventory.csv')
    agr = pd.read_csv(OUT / 'crosslink_agreement_all.csv')

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(17.0, 9.8), gridspec_kw={'width_ratios': [1, 1.1]})

    # --- A  sampling per run ------------------------------------------------
    inv = inv.sort_values('total_us')
    y = np.arange(len(inv))
    axA.barh(y, inv.total_us, color=[ENGINE[e] for e in inv.engine], height=.68)
    axA.set_yticks(y); axA.set_yticklabels(inv.run, fontsize=8.6)
    for yi, (us, n) in enumerate(zip(inv.total_us, inv.n_reps)):
        axA.text(us + .15, yi, f'{us:g} µs  ({n} rep{"s" if n > 1 else ""})',
                 va='center', fontsize=8.3, color=INK['secondary'])
    axA.set_xlim(0, inv.total_us.max() * 1.45)
    axA.set_xlabel('trajectory per run (µs)', fontsize=9.5, color=INK['secondary'])
    bare(axA)
    axA.legend(handles=[Patch(fc=c, ec='none', label=e) for e, c in ENGINE.items()],
               frameon=False, loc='lower right', fontsize=9.5)
    title(axA, 'A', 'Sampling per run',
          f'{inv.n_reps.sum()} replicates, {inv.total_us.sum():.1f} µs total'
          '  ·  every run 500 ns/replicate in a 40 nm box')

    # --- B  agreement, per-frame -------------------------------------------
    order = [m for m in SHORT if m in set(agr.model)][::-1]   # narrative order top-down
    h, y0 = .26, np.arange(len(order))
    for i, (name, col) in enumerate(SUBSET):
        raw = [agr[(agr.model == m) & (agr.subset == name)].per_frame_30A.iloc[0]
               for m in order]
        v = [max(x, FLOOR) for x in raw]
        pos = y0 + (i - 1) * h
        axB.barh(pos, v, height=h * .9, color=col, label=name)
        for yy, xx, rr in zip(pos, v, raw):
            axB.text(xx * 1.14, yy, f'{rr:g}', va='center', fontsize=7.4,
                     color=INK['secondary'])
    axB.set_xscale('log'); axB.set_xlim(FLOOR * .7, 400)
    axB.set_xticks([0.01, 0.1, 1, 10, 100])
    axB.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
    axB.set_yticks(y0); axB.set_yticklabels([SHORT[m] for m in order], fontsize=8.8)
    axB.set_xlabel('mean % of frames with Cα–Cα < 30 Å   (log scale)',
                   fontsize=9.5, color=INK['secondary'])
    bare(axB)
    axB.legend(frameon=False, loc='lower right', fontsize=9.3, ncol=3,
               title='crosslink subset', title_fontsize=9.3,
               bbox_to_anchor=(1.0, 1.075))
    title(axB, 'B', 'Crosslink agreement, per-frame',
          'bars pinned at the left edge are exact zeros')

    fig.suptitle('PARP9 / DTX3L campaign — sampling and crosslink agreement',
                 x=0.006, y=1.0, ha='left', fontsize=14.5,
                 color=INK['primary'], fontweight='semibold')
    fig.tight_layout(rect=[0, 0, 1, 0.95], w_pad=3.0)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'16_campaign_bars.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('wrote', FIG / '16_campaign_bars.png')


def fig17():
    agr = pd.read_csv(OUT / 'crosslink_agreement_all.csv')
    fig = plt.figure(figsize=(15.5, 11.0))
    gs = fig.add_gridspec(2, 2, hspace=.46, wspace=.26)
    axC = fig.add_subplot(gs[0, 0]); axD = fig.add_subplot(gs[0, 1])
    axE = fig.add_subplot(gs[1, 0])
    gsF = gs[1, 1].subgridspec(1, 4, wspace=.16)

    # --- C  ensemble coverage ----------------------------------------------
    order = [m for m in SHORT if m in set(agr.model)][::-1]
    sub = 'reproducible (21)'
    ev = [agr[(agr.model == m) & (agr.subset == sub)].ever_30A.iloc[0] for m in order]
    p1 = [agr[(agr.model == m) & (agr.subset == sub)].pct1_30A.iloc[0] for m in order]
    y, h = np.arange(len(order)), .38
    axC.barh(y + h / 2, ev, height=h * .92, color=EVER, label='satisfied ever')
    axC.barh(y - h / 2, p1, height=h * .92, color=PCT1, label='in ≥1% of frames')
    axC.set_yticks(y); axC.set_yticklabels([SHORT[m] for m in order], fontsize=8.4)
    axC.set_xlim(0, 108)
    axC.set_xlabel('% of the 30 reproducible crosslinks', fontsize=9.5,
                   color=INK['secondary'])
    bare(axC)
    axC.legend(frameon=False, loc='lower right', fontsize=9, ncol=2,
               bbox_to_anchor=(1.0, 1.075))
    title(axC, 'C', 'Ensemble coverage',
          'a wide gap = the geometry is visited, but never populated')

    # --- D  the transfer result --------------------------------------------
    # binding_go/p9_dtx3l_xl_h20, split by whether the pair is one of the 3 restrained
    pf, pc = [99.3, 17.6], [100.0, 83.7]
    x, w = np.arange(2), .34
    axD.bar(x - w / 2, pf, width=w * .95, color=PAIR2[0], label='per-frame < 30 Å')
    axD.bar(x + w / 2, pc, width=w * .95, color=PAIR2[1], label='in ≥1% of frames')
    for xi, (a, b) in enumerate(zip(pf, pc)):
        axD.text(xi - w / 2, a + 2, f'{a}%', ha='center', fontsize=10.5,
                 color=INK['primary'], fontweight='semibold')
        axD.text(xi + w / 2, b + 2, f'{b}%', ha='center', fontsize=10.5,
                 color=INK['primary'], fontweight='semibold')
    axD.set_xticks(x)
    axD.set_xticklabels(['restrained\n(the 3 driving pairs)',
                         'NOT restrained\n(the other 43)'], fontsize=9.8)
    axD.set_xlim(-.6, 1.6); axD.set_ylim(0, 118)
    axD.set_ylabel('% of crosslinks', fontsize=9.5, color=INK['secondary'])
    bare(axD, xgrid=False, ygrid=True)
    axD.legend(frameon=False, loc='upper center', fontsize=9, ncol=2,
               bbox_to_anchor=(.42, 1.0))
    title(axD, 'D', 'Restraining 3 crosslinks moves the other 43',
          'CALVADOS + 3 FL–FL crosslinks, harmonic k = 20')

    # --- E  HyRes fold stability and association ----------------------------
    # both series are percentages, so they share one axis -- no second y-scale
    runs = ['HyRes\nunrestrained', 'HyRes matched\nseparated', 'HyRes matched\ndocked']
    infl_n = [13, 0, 0]                                   # of 14 domains
    infl = [100 * n / 14 for n in infl_n]
    bound = [0.0, 74.2, 100.0]                            # mean over replicates
    blo, bhi = [0.0, 40.5, 100.0], [0.0, 95.8, 100.0]     # min-max over replicates
    worst = ['1.69×', '1.09×', '1.08×']
    x, w = np.arange(3), .34
    axE.bar(x - w / 2, infl, width=w * .95, color=PAIR2[0],
            label='domains inflated > 25% of native Rg')
    axE.bar(x + w / 2, bound, width=w * .95, color=PAIR2[1],
            label='bound fraction (mean of replicates)')
    axE.errorbar(x + w / 2, bound,
                 yerr=[np.array(bound) - np.array(blo), np.array(bhi) - np.array(bound)],
                 fmt='none', ecolor=INK['primary'], elinewidth=1.4, capsize=4)
    for xi, (a, b, n, wst) in enumerate(zip(infl, bound, infl_n, worst)):
        axE.text(xi - w / 2, a + 2, f'{n} of 14', ha='center', fontsize=9.5,
                 color=INK['primary'], fontweight='semibold')
        axE.text(xi + w / 2, max(b, bhi[xi]) + 2, f'{b:g}%', ha='center', fontsize=9.5,
                 color=INK['primary'], fontweight='semibold')
        axE.text(xi, -13, f'worst inflation {wst}', ha='center', fontsize=8.4,
                 color=INK['secondary'])
    axE.set_xticks(x); axE.set_xticklabels(runs, fontsize=9.3)
    axE.set_ylim(0, 112); axE.set_ylabel('%', fontsize=9.5, color=INK['secondary'])
    bare(axE, xgrid=False, ygrid=True)
    # E's two labels are long, so the legend goes under the axes rather than
    # across the title
    axE.legend(frameon=False, loc='upper center', fontsize=8.8, ncol=2,
               bbox_to_anchor=(.5, -.135))
    title(axE, 'E', 'HyRes: folds hold only with matched restraints',
          'whisker = min–max over the 3 replicates')

    # --- F  crosslinks satisfied along single trajectories (faceted) --------
    FACETS = [('p9_dtx3l_rep1_timecourse.csv',        'CALVADOS\nunrestrained'),
              ('p9_dtx3l_xl_h20_rep1_timecourse.csv', 'CALVADOS\n+3 XL k=20'),
              ('separated_rep2_timecourse.csv',       'HyRes matched\nseparated'),
              ('docked_rep1_timecourse.csv',          'HyRes matched\ndocked')]
    axF0, fax = None, []
    for i, (fn, lab) in enumerate(FACETS):
        a = fig.add_subplot(gsF[0, i], sharey=axF0)
        axF0 = axF0 or a
        d = pd.read_csv(OUT / 'pymol' / fn)
        xs = np.arange(len(d))
        a.bar(xs, d.n_sat, width=.72, color=ONE)
        for xi, v in zip(xs, d.n_sat):
            a.text(xi, v + .4, str(v), ha='center', fontsize=8,
                   color=INK['primary'] if v else INK['secondary'])
        a.set_xticks(xs); a.set_xticklabels([f'{v:g}' for v in d.ns], fontsize=7.6,
                                            rotation=90)
        a.set_ylim(0, 21)
        a.text(.5, 1.02, lab, transform=a.transAxes, ha='center', va='bottom',
               fontsize=8.8, color=INK['primary'])
        bare(a, xgrid=False, ygrid=True)
        if i:
            a.tick_params(labelleft=False)
        else:
            a.set_ylabel('crosslinks satisfied (of 36)', fontsize=9.5,
                         color=INK['secondary'])
        fax.append(a)
    box = [a.get_position() for a in fax]
    fig.text((box[0].x0 + box[-1].x1) / 2, box[0].y0 - .055, 'time (ns)',
             ha='center', fontsize=9.5, color=INK['secondary'])
    fig.text(box[0].x0, box[0].y1 + .085, 'F   Crosslinks satisfied along one trajectory',
             ha='left', fontsize=11.5, color=INK['primary'], fontweight='semibold')
    fig.text(box[0].x0, box[0].y1 + .058,
             'one replicate per run, frames at 0 / 25 / 50 / 75 / 100% of 500 ns',
             ha='left', fontsize=9, color=INK['secondary'])

    fig.suptitle('PARP9 / DTX3L campaign — the result tables as bars',
                 x=0.006, y=1.0, ha='left', fontsize=14.5,
                 color=INK['primary'], fontweight='semibold')
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'17_result_bars.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('wrote', FIG / '17_result_bars.png')


if __name__ == '__main__':
    fig16(); fig17()
