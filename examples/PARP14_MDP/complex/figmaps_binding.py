#!/usr/bin/env python3
"""Bound-window contact maps + difference maps for the binding campaign.

Replaces figmaps.py, which is hardwired to the old system names
('parp14_parp9_dtx3l', ...) and the old variant tag 'prod_dom-all_kgo4'. Those
keys match nothing in this campaign, so the panel list filtered to empty and the
script died on a 0-column GridSpec.

Three comparisons, following PROPOSAL.md section 7:
  A  every set's bound-window interface
  B  docked - separated        (is the bound fraction thermodynamic, or hysteresis?)
  C  ternary - dimer           (does the third chain reshape the interface?)

Hetero-vs-homo is shown side by side in panel A rather than as a difference: a
homodimer map is (chain x itself) so its axes do not match the heterodimer's.
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

OUT = Path(__file__).resolve().parent / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
maps = np.load(OUT / 'data' / 'bound_maps.npy', allow_pickle=True)[0]
plt.rcParams.update({'font.size': 8, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

V = 'binding'


def heat(ax, M, bi, bj, title, cmap='magma_r', norm=None,
         lab='contacts / frame while bound'):
    im = ax.imshow(M, cmap=cmap, aspect='auto', norm=norm)
    ax.set_xticks(range(len(bj))); ax.set_xticklabels(bj, rotation=90, fontsize=5)
    ax.set_yticks(range(len(bi))); ax.set_yticklabels(bi, fontsize=5)
    ax.set_title(title, fontsize=6.5, loc='left')
    cb = plt.colorbar(im, ax=ax, fraction=0.046)
    cb.ax.tick_params(labelsize=5); cb.set_label(lab, fontsize=5)


# ---- A: all bound-window interfaces ----
panels = sorted(maps.keys())
if panels:
    ncol = 3
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow), squeeze=False)
    for ax, k in zip(axes.ravel(), panels):
        v = maps[k]
        heat(ax, v['map'], v['bi'], v['bj'], f"{k[0]} · {k[2]}  (n={v['nsel']} reps)")
        ax.set_ylabel(v['ci'], fontsize=6); ax.set_xlabel(v['cj'], fontsize=6)
    for ax in axes.ravel()[len(panels):]:
        ax.axis('off')
    fig.suptitle('Inter-chain domain contacts, conditioned on the chains being bound',
                 fontsize=9, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(FIG / '04_all_bound_maps.png'); plt.close(fig)
    print(f'04_all_bound_maps.png: {len(panels)} panels')


def diff_panels(comps, fname, suptitle):
    """comps: list of (key_a, key_b, title) -> map_a - map_b, shapes permitting."""
    ok = []
    for ka, kb, t in comps:
        if ka in maps and kb in maps and maps[ka]['map'].shape == maps[kb]['map'].shape:
            ok.append((ka, kb, t))
    if not ok:
        print(f'{fname}: no comparable pairs yet, skipped')
        return
    ncol = min(2, len(ok))
    nrow = int(np.ceil(len(ok) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 4.0 * nrow), squeeze=False)
    for ax, (ka, kb, t) in zip(axes.ravel(), ok):
        D = maps[ka]['map'] - maps[kb]['map']
        m = float(np.abs(D).max()) or 1.0
        heat(ax, D, maps[ka]['bi'], maps[ka]['bj'], t, cmap='RdBu_r',
             norm=TwoSlopeNorm(vcenter=0, vmin=-m, vmax=m),
             lab='Δ contacts / frame while bound')
        ax.set_ylabel(maps[ka]['ci'], fontsize=6); ax.set_xlabel(maps[ka]['cj'], fontsize=6)
    for ax in axes.ravel()[len(ok):]:
        ax.axis('off')
    fig.suptitle(suptitle, fontsize=9, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(FIG / fname); plt.close(fig)
    print(f'{fname}: {len(ok)} panels')


# ---- B: docked - separated ----
pairs_by_set = {}
for (s, v, p) in maps:
    pairs_by_set.setdefault(s, []).append(p)
comps_b = []
for s in sorted(pairs_by_set):
    if s.endswith('_docked'):
        continue
    for p in sorted(pairs_by_set[s]):
        comps_b.append(((f'{s}_docked', V, p), (s, V, p),
                        f'{s} · {p}\ndocked start − separated start'))
diff_panels(comps_b, '05_docked_minus_separated.png',
            'Start-state dependence. Flat (white) = the two arms converge, so the bound\n'
            'fraction is thermodynamic. Structure = hysteresis, and 500 ns is too short.')

# ---- C: ternary - dimer ----
comps_c = [(('ternary', V, 'parp9-dtx3l'),  ('p9_dtx3l', V, 'parp9-dtx3l'),
            'PARP9–DTX3L: in the trimer − on its own'),
           (('ternary', V, 'parp14-parp9'), ('p14_p9', V, 'parp14-parp9'),
            'PARP14–PARP9: in the trimer − on its own'),
           (('ternary', V, 'parp14-dtx3l'), ('p14_dtx3l', V, 'parp14-dtx3l'),
            'PARP14–DTX3L: in the trimer − on its own')]
diff_panels(comps_c, '06_ternary_minus_dimer.png',
            'Does the third chain reshape each interface?')
