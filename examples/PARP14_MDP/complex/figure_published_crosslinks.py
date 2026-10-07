#!/usr/bin/env python3
"""Every published PARP9/DTX3L crosslink, mapped onto the domain architecture.

SOURCE
------
Ashok Y, Vela-Rodriguez C, Yang C, Alanen HI, Liu F, Paschal BM, Lehtio L.
"Reconstitution of the DTX3L-PARP9 complex reveals determinants for
high-affinity heterodimerization and multimeric assembly."
Biochem J. 2022;479(3):289-304. doi:10.1042/BCJ20210722. PMID 35037691.
Supplementary bcj-2021-0722_supp.xlsx (identical to the bioRxiv media-2.xlsx,
same md5). BS3, four replicates per construct.

300 crosslinks in total across three constructs:
    FL PARP9 + D3          DTX3L 230-510 only
    FL PARP9 + D3RD        DTX3L D3 + RING + DTC
    FL PARP9 + FL DTX3L    both full length

51 of those are INTER-protein. An earlier pass of mine reported only 3, because
it read the full-length sheet alone -- the truncated-DTX3L constructs yield 26
and 22 inter-protein crosslinks respectively, and they are the informative ones:
removing DTX3L's N-terminal region concentrates the crosslinking onto the D3
interface instead of spreading it over the flexible N terminus.

CAVEAT the '>=2/4 YES/No' reproducibility column is blank in the file's cell
VALUES -- the call appears to be carried by cell fill colour, which is not read
here. So no crosslink is filtered on reproducibility in this figure; treat
every point as "detected", not "detected reproducibly".

Usage:  python figure_published_crosslinks.py
"""
from pathlib import Path
import numpy as np, pandas as pd, yaml
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
SRC = HERE / 'data' / 'published_xlinks_v2.csv'
OUT = HERE / 'analysis'
FIG = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

# validated 3-colour categorical set (style guide PROTEIN slot; all-pairs PASS,
# worst CVD dE 9.2, worst normal-vision dE 24.5)
CONSTRUCT_COLOR = {'FL PARP9 + D3': '#4A3AA7',
                   'FL PARP9 + D3RD': '#C8481A',
                   'FL PARP9 + FL DTX3L': '#117C60'}
INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'muted': '#898781',
       'grid': '#E1E0D9'}
LEN = {'parp9': 854, 'dtx3l': 740}
NAMES = {'parp9': ['½KH1a', 'MD1', 'MD2', '½KH1b', 'KH2', 'ART'],
         'dtx3l': ['RRM', 'KH1', 'KH2', 'KH3', 'KH4', 'KH5', 'RING', 'DTC']}


def domains():
    d = yaml.safe_load(open(HERE / 'binding' / 'p9_dtx3l' / 'input' / 'domains.yaml'))
    return {c: list(zip(NAMES[c], d[c])) for c in ('parp9', 'dtx3l')}


def draw_axis_domains(ax, doms, axis, lo, hi):
    """Shade domain blocks and label them along one axis."""
    for nm, (a, b) in doms:
        if axis == 'x':
            ax.axvspan(a, b, color=INK['grid'], alpha=0.55, lw=0, zorder=0)
            ax.text((a+b)/2, hi*0.985, nm, ha='center', va='top',
                    fontsize=5.5, rotation=90, color=INK['secondary'], zorder=2)
        else:
            ax.axhspan(a, b, color=INK['grid'], alpha=0.55, lw=0, zorder=0)
            ax.text(hi*0.99, (a+b)/2, nm, ha='right', va='center',
                    fontsize=5.5, color=INK['secondary'], zorder=2)


def main():
    d = pd.read_csv(SRC)
    doms = domains()
    inter = d[d.block == 'P9-Dtx3L'].copy()
    intra9 = d[d.block == 'PARP9'].copy()
    intraD = d[d.block.isin(['DTX3L', 'DTX3L-DTX3L'])].copy()

    fig = plt.figure(figsize=(13.5, 6.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1, 1], height_ratios=[1, 0.05],
                          hspace=0.05, wspace=0.34)

    # ---- A: inter-protein map -------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    draw_axis_domains(ax, doms['parp9'], 'x', 0, LEN['dtx3l'])
    draw_axis_domains(ax, doms['dtx3l'], 'y', 0, LEN['parp9'])
    for c, sub in inter.groupby('construct'):
        ax.scatter(sub.res_a, sub.res_b, s=46, color=CONSTRUCT_COLOR[c],
                   edgecolor='white', linewidth=0.7, label=f'{c}  (n={len(sub)})',
                   zorder=3, alpha=0.92)
    fl = inter[inter.construct == 'FL PARP9 + FL DTX3L']
    offs = [(8, 10), (8, -14), (-52, 10)]
    for (_, r), o in zip(fl.iterrows(), offs):
        ax.annotate(f'K{r.res_a}–K{r.res_b}', (r.res_a, r.res_b), fontsize=6,
                    xytext=o, textcoords='offset points', color=INK['primary'],
                    zorder=4, arrowprops=dict(arrowstyle='-', lw=0.5,
                                              color=INK['secondary']))
    ax.set_xlim(0, LEN['parp9']); ax.set_ylim(0, LEN['dtx3l'])
    ax.set_xlabel('PARP9 residue', fontsize=8.5)
    ax.set_ylabel('DTX3L residue', fontsize=8.5)
    ax.set_title(f'A  Inter-protein crosslinks  (n={len(inter)})',
                 fontsize=9.5, loc='left', weight='bold', pad=10)
    ax.legend(fontsize=6.5, frameon=False, loc='upper left')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=7)

    # ---- B/C: intra-protein maps ----------------------------------------
    for k, (ax_i, sub, chain, title) in enumerate((
            (gs[0, 1], intra9, 'parp9', 'B  Intra-PARP9'),
            (gs[0, 2], intraD, 'dtx3l', 'C  Intra-DTX3L'))):
        ax = fig.add_subplot(ax_i)
        L = LEN[chain]
        draw_axis_domains(ax, doms[chain], 'x', 0, L)
        draw_axis_domains(ax, doms[chain], 'y', 0, L)
        ax.plot([0, L], [0, L], color=INK['muted'], lw=0.6, ls='--', zorder=1)
        for c, s2 in sub.groupby('construct'):
            lo = np.minimum(s2.res_a, s2.res_b); hi = np.maximum(s2.res_a, s2.res_b)
            ax.scatter(lo, hi, s=26, color=CONSTRUCT_COLOR[c], edgecolor='white',
                       linewidth=0.5, zorder=3, alpha=0.9)
            ax.scatter(hi, lo, s=26, color=CONSTRUCT_COLOR[c], edgecolor='white',
                       linewidth=0.5, zorder=3, alpha=0.9)
        ax.set_xlim(0, L); ax.set_ylim(0, L)
        ax.set_xlabel(f'{chain.upper()} residue', fontsize=8.5)
        ax.set_ylabel(f'{chain.upper()} residue', fontsize=8.5)
        ax.set_title(f'{title}  (n={len(sub)})', fontsize=9.5, loc='left', weight='bold', pad=10)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=7)

    fig.suptitle('Published BS3 crosslinks, PARP9 / DTX3L  —  Ashok et al. 2022, '
                 'Biochem J 479:289 (doi 10.1042/BCJ20210722)\n'
                 'All 300 crosslinks from the supplement; labelled pairs are the three '
                 'inter-protein crosslinks of the full-length + full-length condition',
                 fontsize=9, x=0.005, y=1.04, ha='left', color=INK['primary'])
    fig.savefig(FIG / '13_published_crosslinks.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('13_published_crosslinks.png')
    print(f'  inter {len(inter)} | intra-PARP9 {len(intra9)} | intra-DTX3L {len(intraD)}')


if __name__ == '__main__':
    main()
