#!/usr/bin/env python3
"""The dimension-controlled comparison, and the SASD correction.

Reads analysis/decoy_pairs.csv (score_decoys.py) and analysis/sasd_rescore.csv
(score_sasd.py). Writes analysis/decoy_auc_ci.csv and figure 22.

Panel A is the one that matters. AUC of the published crosslinks against the
4,716 other inter-chain lysine pairs, inside each model's own distribution.
Rank-based, so a model that merely holds its chains compact gains nothing: the
decoys move with the crosslinks. The three pairs that binding_go restrains are
EXCLUDED, since a restrained pair sits at the top of its own ranking by
construction and would flatter those runs.

Only 18 positives survive that exclusion, so the bootstrap intervals are wide
and most models overlap. The figure is drawn to make the overlap visible rather
than to rank models that cannot be separated.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'grid': '#E1E0D9'}
RESTRAINED = {(557, 401), (632, 401), (557, 363)}      # the 3 binding_go pairs
B = 2000


def auc(score, y):
    y = np.asarray(y, bool)
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(score).rank().to_numpy()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def auc_ci(score, y, seed=0):
    rng = np.random.default_rng(seed)
    score, y = np.asarray(score), np.asarray(y, bool)
    pi, ni = np.where(y)[0], np.where(~y)[0]
    out = np.empty(B)
    for k in range(B):
        i = np.concatenate([rng.choice(pi, len(pi), True),
                            rng.choice(ni, len(ni), True)])
        out[k] = auc(score[i], y[i])
    return np.nanpercentile(out, [2.5, 97.5])


def table():
    d = pd.read_csv(OUT / 'decoy_pairs.csv')
    rows = []
    for m, g in d.groupby('model'):
        g = g.reset_index(drop=True)
        keep = ~g.apply(lambda r: (r.res_a, r.res_b) in RESTRAINED, axis=1)
        g = g[keep]
        a = auc(g.pct_within, g.is_repro)
        lo, hi = auc_ci(g.pct_within.to_numpy(), g.is_repro.to_numpy())
        xl, dec = g[g.is_repro], g[~g.is_xl]
        rows.append(dict(model=m, n_pos=int(g.is_repro.sum()),
                         n_decoy=int((~g.is_xl).sum()),
                         auc=round(a, 3), lo=round(lo, 3), hi=round(hi, 3),
                         xl_median=round(float(xl.pct_within.median()), 3),
                         decoy_median=round(float(dec.pct_within.median()), 3)))
    t = pd.DataFrame(rows).sort_values('auc', ascending=False)
    t.to_csv(OUT / 'decoy_auc_ci.csv', index=False)
    return t


SHORT = {'CALVADOS unrestrained (separated)': 'CALVADOS unrestrained, separated',
         'CALVADOS unrestrained (docked)': 'CALVADOS unrestrained, docked',
         'CALVADOS ternary (docked)': 'CALVADOS ternary, docked',
         'CALVADOS +3 XL harmonic k=20': 'CALVADOS +3 XL harmonic k=20',
         'CALVADOS +3 XL Go r=1.5 k=100': 'CALVADOS +3 XL Gō r=1.5 k=100',
         'CALVADOS +21 XL IMProv-tiered': 'CALVADOS +21 XL tiered',
         'HyRes unrestrained': 'HyRes unrestrained',
         'HyRes matched (separated)': 'HyRes matched, separated',
         'HyRes matched (docked)': 'HyRes matched, docked'}
COL = lambda m: '#C8481A' if m.startswith('HyRes') else '#4A3AA7'


def bare(ax):
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(INK['grid'])
    ax.tick_params(colors=INK['secondary'], length=0, labelsize=9)
    ax.xaxis.grid(True, color=INK['grid'], lw=.8); ax.set_axisbelow(True)


def main():
    t = table()
    sas = pd.read_csv(OUT / 'sasd_rescore.csv') if (OUT / 'sasd_rescore.csv').is_file() else None
    fig, ax = plt.subplots(1, 3, figsize=(18.5, 6.0),
                           gridspec_kw={'width_ratios': [1.25, 1, 1]})

    # --- A  AUC with bootstrap CI ------------------------------------------
    t2 = t.sort_values('auc')
    y = np.arange(len(t2))
    a = ax[0]
    a.barh(y, t2.auc, color=[COL(m) for m in t2.model], height=.62)
    a.errorbar(t2.auc, y, xerr=[t2.auc - t2.lo, t2.hi - t2.auc], fmt='none',
               ecolor=INK['primary'], elinewidth=1.3, capsize=4)
    a.axvline(.5, color=INK['primary'], lw=1.4, ls=(0, (4, 3)))
    a.text(.503, len(t2) - .4, 'chance', fontsize=8.5, color=INK['secondary'])
    a.set_yticks(y); a.set_yticklabels([SHORT.get(m, m) for m in t2.model], fontsize=9)
    a.set_xlim(.3, 1.0)
    a.set_xlabel('AUC: published crosslinks vs 4,716 lysine decoys, '
                 'within each model', fontsize=9.5, color=INK['secondary'])
    bare(a)
    a.text(0, 1.055, 'A   Dimension-controlled discrimination',
           transform=a.transAxes, fontsize=11.5, color=INK['primary'],
           fontweight='semibold')
    a.text(0, 1.012, 'rank-based, so chain compaction cancels  ·  '
           'the 3 restrained pairs excluded  ·  18 positives, 95% CI',
           transform=a.transAxes, fontsize=8.8, color=INK['secondary'])

    # --- B  absolute satisfaction: the compaction effect --------------------
    b = ax[1]
    w, yy = .38, np.arange(len(t2))
    FLOOR = 0.004
    b.barh(yy + w/2, np.maximum(t2.xl_median, FLOOR), height=w*.92,
           color='#4A3AA7', label='published crosslinks')
    b.barh(yy - w/2, np.maximum(t2.decoy_median, FLOOR), height=w*.92,
           color='#8E7BD6', label='lysine decoys')
    b.set_xscale('log'); b.set_xlim(FLOOR*.7, 200)
    b.set_yticks(yy); b.set_yticklabels([])
    b.set_xlabel('median % of frames within 30 Å  (log)', fontsize=9.5,
                 color=INK['secondary'])
    bare(b)
    b.legend(frameon=False, loc='lower right', fontsize=9, bbox_to_anchor=(1, 1.0), ncol=2)
    b.text(0, 1.055, 'B   What the raw numbers show', transform=b.transAxes,
           fontsize=11.5, color=INK['primary'], fontweight='semibold')
    b.text(0, 1.012, 'decoys rise with the crosslinks — that gap is compaction',
           transform=b.transAxes, fontsize=8.8, color=INK['secondary'])

    # --- C  SASD correction --------------------------------------------------
    c = ax[2]
    if sas is not None:
        g = (sas.groupby('model')[['pct_euclid', 'pct_sasd']].mean()
             .reindex(t2.model).dropna())
        yy = np.arange(len(g))
        c.barh(yy + w/2, np.maximum(g.pct_euclid, FLOOR), height=w*.92,
               color='#117C60', label='Euclidean')
        c.barh(yy - w/2, np.maximum(g.pct_sasd, FLOOR), height=w*.92,
               color='#7FCBAE', label='SASD')
        for i, (e, s_) in enumerate(zip(g.pct_euclid, g.pct_sasd)):
            if e > 0:
                c.text(max(e, FLOOR)*1.25, i + w/2, f'−{100*(e-s_)/e:.0f}%',
                       va='center', fontsize=7.5, color=INK['secondary'])
        c.set_xscale('log'); c.set_xlim(FLOOR*.7, 200)
        c.set_yticks(yy); c.set_yticklabels([])
        c.set_xlabel('mean % of frames within 30 Å  (log)', fontsize=9.5,
                     color=INK['secondary'])
        bare(c)
        c.legend(frameon=False, loc='lower right', fontsize=9,
                 bbox_to_anchor=(1, 1.0), ncol=2)
    c.text(0, 1.055, 'C   Solvent-accessible distance', transform=c.transAxes,
           fontsize=11.5, color=INK['primary'], fontweight='semibold')
    c.text(0, 1.012, 'paths that pass through protein removed',
           transform=c.transAxes, fontsize=8.8, color=INK['secondary'])

    fig.tight_layout(w_pad=2.0)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'22_decoy_control.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print(t.to_string(index=False))
    print('\nwrote', FIG / '22_decoy_control.png', 'and', OUT / 'decoy_auc_ci.csv')


if __name__ == '__main__':
    main()
