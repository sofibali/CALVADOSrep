#!/usr/bin/env python3
"""Which Lys and Asp/Glu are solvent-exposed, per residue, along each sequence.

DMTMM couples a carboxyl (Asp/Glu) to a primary amine (Lys), so those two
classes are the entire reactive inventory for the crosslinking work. A buried
residue cannot be crosslinked at any distance, so what matters is not how many
the sequence contains but how many are on the surface, and where.

One lane per residue class per chain, one vertical line per residue at its
sequence position: light where buried, saturated where exposed. Reading down
the figure gives PARP14, then PARP9, then DTX3L.

RSA is computed with freesasa on the same all-atom input structures the
simulations were built from, normalised by the Chothia Gly-X-Gly
extended-tripeptide max ASA, and thresholded exactly as
sim_analysis/figure_sasa_faces.py does -- that table and those cutoffs are
imported from it rather than restated, so the two analyses cannot drift:

    RSA < 0.05          buried
    0.05 <= RSA < 0.20  partial
    RSA >= 0.20         surface

A caveat worth carrying: this is the exposure of ONE static structure, the AF3
/ AF2 model used to seed the simulations. Side chains move, and a residue that
is partial here may be reactive in the ensemble. The CG trajectories cannot
settle it -- a one-bead-per-residue model has no side chain to expose -- so
treat the counts as a structural prior, not as a reactivity measurement.

Usage:  python figure_crosslinkable_residues.py
        -> analysis/crosslinkable_residues.csv
           analysis/figures/23_crosslinkable_residues.png
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import freesasa
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parent / 'sim_analysis'))
from figure_sasa_faces import MAX_ASA_CHOTHIA, RSA_BURIED_MAX, RSA_PARTIAL_MAX

INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'grid': '#E1E0D9'}
# one hue per residue class; within a hue, buried is pale and exposed saturated,
# so the lane reads as an ordinal ramp rather than four competing colours
LYS_BURIED, LYS_EXPOSED = '#C3B6E8', '#4A3AA7'
ACID_BURIED, ACID_EXPOSED = '#F0B79A', '#C8481A'

CHAINS = [
    ('PARP14', 'binding/p14_p9/input/parp14.pdb',
     [('RRM1', 1, 145), ('RRM2', 146, 224), ('RRM3', 225, 314),
      ('KH1-6', 315, 737), ('KH7a', 738, 789), ('MD1', 790, 1003),
      ('MD2', 1004, 1193), ('MD3', 1207, 1388), ('KHb-8', 1389, 1533),
      ('WWE', 1534, 1602), ('ART', 1603, 1801)]),
    ('PARP9', 'binding/p9_dtx3l/input/parp9.pdb',
     [('1/2KH1a', 60, 96), ('MD1', 107, 296), ('MD2', 311, 493),
      ('1/2KH1b', 516, 539), ('KH2', 543, 622), ('ART', 644, 823)]),
    ('DTX3L', 'binding/p9_dtx3l/input/dtx3l.pdb',
     [('RRM', 11, 87), ('KH1', 138, 201), ('KH2', 233, 302), ('KH3', 305, 362),
      ('KH4', 378, 447), ('KH5', 449, 507), ('RING', 556, 607), ('DTC', 608, 737)]),
]
CLASSES = [('Lys', ('LYS',), LYS_BURIED, LYS_EXPOSED),
           ('Asp/Glu', ('ASP', 'GLU'), ACID_BURIED, ACID_EXPOSED)]


def rsa_table(pdb):
    """resid -> (resname, RSA). freesasa's own classifier, Chothia denominator."""
    s = freesasa.Structure(str(pdb))
    r = freesasa.calc(s)
    acc = {}
    for i in range(s.nAtoms()):
        key = int(s.residueNumber(i))
        acc.setdefault(key, [s.residueName(i).strip(), 0.0])
        acc[key][1] += r.atomArea(i)
    out = {}
    for resid, (name, sasa) in acc.items():
        mx = MAX_ASA_CHOTHIA.get(name)
        if mx:
            out[resid] = (name, sasa / mx)
    return out


def main():
    rows = []
    for chain, pdb, _ in CHAINS:
        for resid, (name, rsa) in sorted(rsa_table(HERE / pdb).items()):
            cls = ('Lys' if name == 'LYS'
                   else 'Asp/Glu' if name in ('ASP', 'GLU') else None)
            if cls is None:
                continue
            rows.append(dict(chain=chain, resid=resid, resname=name,
                             residue_class=cls, rsa=round(rsa, 4),
                             burial=('buried' if rsa < RSA_BURIED_MAX else
                                     'partial' if rsa < RSA_PARTIAL_MAX else 'surface'),
                             exposed=rsa >= RSA_PARTIAL_MAX))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'crosslinkable_residues.csv', index=False)

    print(f"{'chain':8s} {'class':8s} {'total':>6s} {'exposed':>8s} {'partial':>8s} "
          f"{'buried':>7s} {'% exposed':>10s}")
    for chain, _, _ in CHAINS:
        for cls, _, _, _ in CLASSES:
            g = df[(df.chain == chain) & (df.residue_class == cls)]
            n = len(g); e = int(g.exposed.sum())
            p = int((g.burial == 'partial').sum()); b = int((g.burial == 'buried').sum())
            print(f'{chain:8s} {cls:8s} {n:6d} {e:8d} {p:8d} {b:7d} {100*e/n:9.1f}%')

    # ---- figure: one lane per class per chain, a line per residue -----------
    nlane = len(CHAINS) * len(CLASSES)
    fig, axes = plt.subplots(nlane, 1, figsize=(15.5, 1.55 * nlane))
    k = 0
    for chain, _, domains in CHAINS:
        sub = df[df.chain == chain]
        xmax = sub.resid.max()
        for cls, _, c_bur, c_exp in CLASSES:
            ax = axes[k]; k += 1
            g = sub[sub.residue_class == cls]
            # domain bands first, so the residue lines sit on top
            for j, (dn, lo, hi) in enumerate(domains):
                ax.axvspan(lo, hi, color=INK['grid'],
                           alpha=.75 if j % 2 == 0 else .4, lw=0, zorder=0)
                if cls == 'Lys':
                    ax.text((lo + hi) / 2, 1.14, dn, ha='center', va='bottom',
                            fontsize=6.6, color=INK['secondary'], clip_on=False)
            for _, r in g.iterrows():
                ax.vlines(r.resid, 0, 1, lw=.85, zorder=2,
                          color=c_exp if r.exposed else c_bur)
            n, e = len(g), int(g.exposed.sum())
            ax.set_xlim(0, xmax + 1); ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_ylabel(f'{chain}\n{cls}', rotation=0, ha='right', va='center',
                          fontsize=9, color=INK['primary'])
            ax.text(1.004, .5, f'{e}/{n}  ({100*e/n:.0f}% exposed)',
                    transform=ax.transAxes, va='center', fontsize=8.2,
                    color=INK['secondary'])
            for s in ('top', 'right', 'left'):
                ax.spines[s].set_visible(False)
            ax.spines['bottom'].set_color(INK['grid'])
            ax.tick_params(colors=INK['secondary'], length=2, labelsize=7.5)
            if k % 2 == 0:
                ax.set_xlabel('residue number', fontsize=8.5,
                              color=INK['secondary'], labelpad=1)

    handles = [Line2D([], [], color=LYS_EXPOSED, lw=2.4, label='Lys, exposed'),
               Line2D([], [], color=LYS_BURIED, lw=2.4, label='Lys, buried/partial'),
               Line2D([], [], color=ACID_EXPOSED, lw=2.4, label='Asp/Glu, exposed'),
               Line2D([], [], color=ACID_BURIED, lw=2.4, label='Asp/Glu, buried/partial')]
    axes[0].legend(handles=handles, frameon=False, fontsize=8.6, ncol=4,
                   loc='lower left', bbox_to_anchor=(0, 1.42))
    fig.text(0.006, 1.0, 'DMTMM-reactive residues along each sequence',
             fontsize=13.5, color=INK['primary'], fontweight='semibold', va='bottom')
    fig.text(0.006, -0.012,
             f'One line per residue at its position. Exposed = RSA ≥ '
             f'{RSA_PARTIAL_MAX:.2f} by freesasa on the all-atom input structure, '
             f'Chothia Gly-X-Gly normalised. Grey bands are restrained domains; '
             f'gaps between them are linkers. Single-structure exposure, not an '
             f'ensemble measurement.',
             fontsize=8.2, color=INK['secondary'], va='top')
    fig.tight_layout(rect=[0, 0.005, 0.93, 0.985], h_pad=1.5)
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'23_crosslinkable_residues.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('\nwrote', OUT / 'crosslinkable_residues.csv')
    print('wrote', FIG / '23_crosslinkable_residues.png')


if __name__ == '__main__':
    main()
