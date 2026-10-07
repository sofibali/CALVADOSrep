#!/usr/bin/env python3
"""Chain- and disorder-resolved render panel.

The domain-coloured renders (render_states.pml) answer "which domain is this?".
This one answers "which protein is this, and is this bit folded or disordered?" --
PARP9 in blues, DTX3L in magentas, folded cores pale and the CALVADOS-unrestrained
(disordered) regions saturated.

Source PNGs: analysis/pymol/chain_*.png (render_chains.pml).
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.image as mpimg

HERE = Path(__file__).resolve().parent
PY = HERE / 'analysis' / 'pymol'
FIG = HERE / 'analysis' / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
INK = {'primary': '#0B0B0B', 'secondary': '#52514E'}

# Validated --pairs all (worst CVD dE 14.3 protan, normal-vision 20.2).
P9_FOLD, P9_DIS = '#63A2FE', '#0052B4'
DTX_FOLD, DTX_DIS = '#CC629C', '#8F1563'
XL = '#1A1A1A'

PANELS = [
    ('chain_setup_separated.png',    'A', 'Setup: separated',
     'chains 17 nm apart, 40 nm box'),
    ('chain_setup_docked.png',       'B', 'Setup: AF3 docked',
     'same box, complex pose'),
    (None, None, None, None),                       # legend cell
    ('chain_unrestrained_best.png',  'C', 'CALVADOS unrestrained',
     'best frame  ·  6 of 36 crosslinks'),
    ('chain_restrained_best.png',    'D', 'CALVADOS + 3 XL (harmonic k=20)',
     'best frame  ·  18 of 36 crosslinks'),
    ('chain_hyres_docked_250ns.png', 'E', 'HyRes matched, docked',
     '250 ns  ·  17 of 36 crosslinks'),
]


def trim(img, pad=8):
    """Crop the white margin PyMOL leaves around the render."""
    ink = (img[..., :3] < 0.985).any(axis=2)
    if not ink.any():
        return img
    ys, xs = np.where(ink)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def main():
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 9.2))
    for ax, (fn, tag, title, sub) in zip(axes.ravel(), PANELS):
        ax.set_axis_off()
        if fn is None:
            continue
        p = PY / fn
        if not p.exists():
            ax.text(.5, .5, f'missing\n{fn}', ha='center', va='center',
                    color=INK['secondary'], fontsize=9, transform=ax.transAxes)
            continue
        ax.imshow(trim(mpimg.imread(p)))
        ax.set_anchor('N')          # top-align, so titles sit on one line per row
        ax.set_title(f'{tag}   {title}', loc='left', fontsize=11.5,
                     color=INK['primary'], fontweight='semibold', pad=2)
        ax.text(0, -0.035, sub, transform=ax.transAxes, ha='left', va='top',
                fontsize=9.5, color=INK['secondary'])

    leg = axes.ravel()[2]
    handles = [Patch(fc=P9_FOLD,  ec='none', label='PARP9  folded core'),
               Patch(fc=P9_DIS,   ec='none', label='PARP9  disordered'),
               Patch(fc=DTX_FOLD, ec='none', label='DTX3L  folded core'),
               Patch(fc=DTX_DIS,  ec='none', label='DTX3L  disordered'),
               Line2D([], [], color=XL, lw=2, ls=(0, (3, 2)),
                      label='satisfied crosslink (< 30 Å)')]
    leg.legend(handles=handles, loc='center', frameon=False, fontsize=11.5,
               labelspacing=1.05, handlelength=1.9, handleheight=1.25,
               borderpad=0)
    leg.text(.5, .955, 'Colour key', ha='center', va='top', fontsize=11.5,
             color=INK['primary'], fontweight='semibold',
             transform=leg.transAxes)
    leg.text(.5, .075,
             '"Folded" = inside a restrained block of\n'
             'binding/p9_dtx3l/input/domains.yaml.\n'
             'Everything else is unrestrained in CALVADOS.',
             ha='center', va='bottom', fontsize=8.8, color=INK['secondary'],
             linespacing=1.5, transform=leg.transAxes)

    fig.suptitle('PARP9 / DTX3L by chain and disorder', x=0.008, y=0.985,
                 ha='left', fontsize=14.5, color=INK['primary'],
                 fontweight='semibold')
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'15_chain_disorder_renders.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('wrote', FIG / '15_chain_disorder_renders.png')


if __name__ == '__main__':
    main()
