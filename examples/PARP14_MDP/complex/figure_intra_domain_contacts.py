#!/usr/bin/env python3
"""Intra-chain inter-domain contacts: how each protein is packed against itself.

The companion to the interface maps. Same two quantities -- % of frames within
1.0 nm, and the mean of the minimum Ca-Ca distance -- but between two domains of
the SAME chain, so they describe each protein's own tertiary arrangement and
whether the force field or the partner changes it.

No conditioning here. An intra-chain pair is always part of one molecule, so
there is no "associated" subset to restrict to; every frame counts.

TWO THINGS TO DISCOUNT WHEN READING THESE MAPS

1. Sequence-adjacent domains are in contact for free. Several linkers in this
   system are only a few residues long -- PARP9 1/2KH1b->KH2 is 3 residues,
   DTX3L KH2->KH3 is 2, KH4->KH5 is 1, RING->DTC is 0 -- so those pairs cannot
   be far apart whatever the force field does. Adjacent cells are outlined in
   grey and labelled with their linker length.

2. PARP9 1/2KH1a <-> 1/2KH1b is not a tertiary contact at all: they are the two
   halves of ONE KH1 fold, split in sequence by the macrodomains. CALVADOS joins
   them with 85 custom restraints; prepare_hyres_matched.py does not build them.
   That cell is outlined in red and is the direct read-out of the gap.

Usage:  python figure_intra_domain_contacts.py
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
INK = {'primary': '#0B0B0B', 'secondary': '#52514E'}

CHAINS = {
    'parp9': (['1/2KH1a', 'MD1', 'MD2', '1/2KH1b', 'KH2', 'ART'],
              [(60, 96), (107, 296), (311, 493), (516, 539), (543, 622), (644, 823)]),
    'dtx3l': (['RRM', 'KH1', 'KH2', 'KH3', 'KH4', 'KH5', 'RING', 'DTC'],
              [(11, 87), (138, 201), (233, 302), (305, 362), (378, 447),
               (449, 507), (556, 607), (608, 737)]),
}
SPLIT_KH1 = ('1/2KH1a', '1/2KH1b')        # one fold, not a tertiary contact
# label -> the models in domain_contact_maps.csv it pools over
MODELS = [('CALVADOS\nunrestrained',      ['binding/p9_dtx3l']),
          ('CALVADOS\n+3 XL k=20',        ['binding_go/p9_dtx3l_xl_h20']),
          ('HyRes matched\nseparated',    [f'HyRes separated rep{i}' for i in (1, 2, 3)]),
          ('HyRes matched\ndocked',       [f'HyRes docked rep{i}' for i in (1, 2, 3)])]


def pool(d, models, chain, names, col):
    """Frame-weighted mean of one column over the given models, as a matrix."""
    mats, ws = [], []
    for m in models:
        g = d[(d.model == m) & (d.chain_i == chain)]
        if not len(g):
            continue
        M = np.full((len(names), len(names)), np.nan)
        for _, r in g.iterrows():
            i, j = names.index(r.domain_i), names.index(r.domain_j)
            M[max(i, j), min(i, j)] = r[col]        # lower triangle
        mats.append(M); ws.append(g.n_frames.iloc[0])
    if not mats:
        return None
    return np.average(mats, axis=0, weights=ws)


def main():
    d = pd.read_csv(OUT / 'domain_contact_maps.csv')
    d = d[d.kind == 'intra']

    nrow = 2 * len(CHAINS)
    fig, axes = plt.subplots(nrow, len(MODELS), figsize=(4.1 * len(MODELS), 19.5))
    for ci, (chain, (names, blocks)) in enumerate(CHAINS.items()):
        gaps = {(k + 1, k): blocks[k + 1][0] - blocks[k][1] - 1
                for k in range(len(names) - 1)}
        for which, (col, cmap, lo, hi, what, unit, fmt) in enumerate([
                ('pct_contact', 'Purples', 0, 100, 'contact occupancy',
                 '% of frames < 1.0 nm', '{:.0f}'),
                ('mean_min_nm', 'Greens_r', 0, 8, 'closest approach',
                 'mean min Cα–Cα (nm)', '{:.1f}')]):
            r = 2 * ci + which
            for c, (label, models) in enumerate(MODELS):
                ax = axes[r, c]
                M = pool(d, models, chain, names, col)
                if M is None:
                    ax.set_axis_off(); continue
                im = ax.imshow(M, cmap=cmap, aspect='auto', vmin=lo, vmax=hi)
                ax.set_xticks(range(len(names)))
                ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
                ax.set_yticks(range(len(names)))
                ax.set_yticklabels(names if c == 0 else [], fontsize=8)
                ax.tick_params(colors=INK['secondary'], length=0)
                for y in range(len(names)):
                    for x in range(len(names)):
                        v = M[y, x]
                        if np.isnan(v):
                            continue
                        frac = (v - lo) / (hi - lo)
                        dark = frac < .45 if cmap.endswith('_r') else frac > .55
                        ax.text(x, y, fmt.format(v), ha='center', va='center',
                                fontsize=7, color='white' if dark else INK['primary'])
                # sequence-adjacent pairs: in contact for free
                for (y, x), g in gaps.items():
                    ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fill=False,
                                           edgecolor=INK['secondary'], lw=1.1))
                    if which == 0:
                        ax.text(x + .44, y - .40, f'{g}', ha='right', va='top',
                                fontsize=5.6, color=INK['secondary'])
                # the split KH1: one fold, unrestrained in HyRes
                if chain == 'parp9':
                    y, x = names.index(SPLIT_KH1[1]), names.index(SPLIT_KH1[0])
                    ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fill=False,
                                           edgecolor='#B0342A', lw=2.2, zorder=4))
                if r % 2 == 0:
                    ax.set_title(label, fontsize=10, color=INK['primary'],
                                 fontweight='semibold', pad=6)
                if c == len(MODELS) - 1:
                    fig.colorbar(im, ax=ax, fraction=.045, pad=.03)
            axes[r, 0].text(0, 1.22 if r % 2 == 0 else 1.045,
                            f'{chain.upper()} — {what}, {unit}'
                            if r % 2 else f'{chain.upper()} — {what}, {unit}',
                            transform=axes[r, 0].transAxes, ha='left', va='bottom',
                            fontsize=11.5, color=INK['primary'], fontweight='semibold')
    fig.text(0.006, -0.004,
             'Grey outline = domains adjacent in sequence, labelled with the linker '
             'length in residues; short linkers force contact whatever the force '
             'field does. Red outline = PARP9 1/2KH1a ↔ 1/2KH1b, the two halves of '
             'one KH1 fold: CALVADOS joins them with 85 custom restraints, '
             'prepare_hyres_matched.py does not build them.',
             ha='left', va='top', fontsize=9, color=INK['secondary'])
    fig.tight_layout(rect=[0, 0.005, 1, 0.985], h_pad=4.0)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'21_intra_domain_contacts.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('wrote', FIG / '21_intra_domain_contacts.png')


if __name__ == '__main__':
    main()
