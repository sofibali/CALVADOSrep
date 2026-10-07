#!/usr/bin/env python3
"""Domain contact maps for the two restrained HyRes arms, per replicate.

These are the only HyRes runs worth reading as a force-field result: the
unrestrained run never associated at all. Both arms carry the within-block
harmonic network built by prepare_hyres_matched.py and NO inter-chain
restraint, so whatever interface appears is the force field's own.

Three columns per arm:
    occupancy        % of associated frames with min Ca-Ca < 1.0 nm, pooled
    closest approach mean of that minimum over associated frames, pooled
    reproducibility  how many of the 3 replicates put the pair above 50%

The third column is the one that matters for this campaign: the separated arm's
bound fraction ranges 40.5-95.8% across its replicates, so a pooled map could be
one replicate's interface averaged with two replicates of nothing.

CAVEAT drawn on the figure: PARP9's split KH1 is NOT held together in HyRes.
prepare_hyres_matched.py builds pairs within each domains.yaml block only, so
the 85 custom restraints that join blocks [60,96] and [516,539] into one fold
are missing; the halves reach 5.0 nm. The 1/2KH1a and 1/2KH1b rows therefore
describe two free half-domains, not a KH1 fold, and are hatched.

Usage:  python figure_hyres_domain_contacts.py
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
P9 = ['1/2KH1a', 'MD1', 'MD2', '1/2KH1b', 'KH2', 'ART']
DX = ['RRM', 'KH1', 'KH2', 'KH3', 'KH4', 'KH5', 'RING', 'DTC']
UNFOLDED = {'1/2KH1a', '1/2KH1b'}          # split KH1, unrestrained in HyRes
ARMS = [('separated', 'HyRes matched, separated'),
        ('docked', 'HyRes matched, docked')]


def grid(sub, col):
    return (sub.pivot(index='domain_i', columns='domain_j', values=col)
            .reindex(index=P9, columns=DX).to_numpy(float))


def main():
    d = pd.read_csv(OUT / 'domain_contact_maps.csv')
    d = d[(d.kind == 'inter') & (d.chain_i == 'parp9') & (d.chain_j == 'dtx3l')]

    fig, axes = plt.subplots(2, 3, figsize=(17.0, 9.0))
    for r, (arm, label) in enumerate(ARMS):
        reps = [d[d.model == f'HyRes {arm} rep{i}'] for i in (1, 2, 3)]
        reps = [g for g in reps if len(g)]
        if not reps:
            continue
        # frame-weighted pool over replicates
        w = np.array([g.n_frames.iloc[0] for g in reps], float)
        occ = np.average([grid(g, 'pct_contact_bound') for g in reps], axis=0,
                         weights=w)
        dist = np.average([grid(g, 'mean_min_bound_nm') for g in reps], axis=0,
                          weights=w)
        repro = np.sum([grid(g, 'pct_contact_bound') > 50 for g in reps], axis=0)
        assoc = [100 * g.n_bound_frames.iloc[0] / g.n_frames.iloc[0] for g in reps]

        panels = [(occ, 'Purples', 0, 100, 'contact occupancy',
                   '% of associated frames < 1.0 nm', '{:.0f}'),
                  (dist, 'Greens_r', 0, 6, 'closest approach',
                   'mean min Cα–Cα (nm)', '{:.1f}'),
                  (repro, 'Blues', 0, 3, 'reproducibility',
                   'replicates of 3 with the pair > 50%', '{:.0f}')]
        for c, (M, cmap, lo, hi, title, unit, fmt) in enumerate(panels):
            ax = axes[r, c]
            im = ax.imshow(M, cmap=cmap, aspect='auto', vmin=lo, vmax=hi)
            ax.set_xticks(range(len(DX))); ax.set_xticklabels(DX, rotation=45,
                                                              ha='right', fontsize=8.5)
            ax.set_yticks(range(len(P9))); ax.set_yticklabels(P9, fontsize=8.5)
            ax.tick_params(colors=INK['secondary'], length=0)
            for y in range(M.shape[0]):
                for x in range(M.shape[1]):
                    v = M[y, x]
                    if np.isnan(v):
                        continue
                    frac = (v - lo) / (hi - lo)
                    dark = frac < .45 if cmap.endswith('_r') else frac > .55
                    ax.text(x, y, fmt.format(v), ha='center', va='center',
                            fontsize=7.2, color='white' if dark else INK['primary'])
            # hatch the two half-domains of the unrestrained split KH1
            for y, name in enumerate(P9):
                if name in UNFOLDED:
                    ax.add_patch(Rectangle((-.5, y - .5), len(DX), 1, fill=False,
                                           hatch='///', edgecolor='#B0342A',
                                           linewidth=0, zorder=3, alpha=.55))
            ax.set_title(f'{title}\n{unit}', loc='left', fontsize=10,
                         color=INK['primary'], fontweight='semibold', pad=6)
            if c == 0:
                ax.set_ylabel('PARP9', fontsize=9.5, color=INK['secondary'])
            if r == 1:
                ax.set_xlabel('DTX3L', fontsize=9.5, color=INK['secondary'])
            fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
        axes[r, 0].text(0, 1.30, f'{label}  ·  associated in '
                        + ' / '.join(f'{a:.0f}%' for a in assoc) + ' of frames '
                        '(rep 1 / 2 / 3)', transform=axes[r, 0].transAxes,
                        ha='left', va='bottom', fontsize=11.5,
                        color=INK['primary'], fontweight='semibold')
    fig.text(0.006, -0.012,
             'Hatched rows: PARP9 1/2KH1a and 1/2KH1b are the two halves of one '
             'KH1 fold, joined in CALVADOS by 85 custom restraints that '
             'prepare_hyres_matched.py does not build. In HyRes they separate to '
             '5.0 nm, so these rows are two free half-domains, not a domain. '
             'No inter-chain restraints are present in either arm.',
             ha='left', va='top', fontsize=9, color=INK['secondary'])
    fig.tight_layout(rect=[0, 0.01, 1, 0.96], h_pad=5.0)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'20_hyres_domain_contacts.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('wrote', FIG / '20_hyres_domain_contacts.png')


if __name__ == '__main__':
    main()
