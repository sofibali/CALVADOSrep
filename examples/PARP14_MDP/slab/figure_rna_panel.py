#!/usr/bin/env python
"""
Figures for slab run 2 (the +/-RNA panel).

    python figure_rna_panel.py

Reads the three tables written by analyze_rna_contacts.py plus
contact_exchange.csv, and writes PNG + SVG into slab/figures/.

Figure order is deliberate: the validity control (Fig 2) comes before the
contact map (Fig 4), because the contact map is only interpretable for systems
the control shows to be liquid.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.join(HERE, 'analysis')
FIG = os.path.join(HERE, 'figures')
os.makedirs(FIG, exist_ok=True)

ORDER = ['fl', 'kh1_art_full', 'core_full_go', 'md_full']
LABEL = {'fl': 'FL\n(1-1801)', 'kh1_art_full': 'KH1-ART\n(315-1801)',
         'core_full_go': 'core\n(738-1801)', 'md_full': 'MD1-MD3\n(790-1388)'}
ARMC = {'homotypic': '#4C72B0', 'rna': '#C44E52'}
DOM_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1', 'MD2', 'MD3',
             'KHb-KH8', 'WWE', 'ART']


def save(fig, name):
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(FIG, f'{name}.{ext}'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', name)


def fig1_condensation(cond):
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    w = 0.36
    for k, arm in enumerate(['homotypic', 'rna']):
        xs, ys, es = [], [], []
        for i, c in enumerate(ORDER):
            s = cond[(cond.construct == c) & (cond.arm == arm)].slab_end_nm.dropna()
            if s.empty:
                continue
            xs.append(i + (k - 0.5) * w); ys.append(s.mean()); es.append(s.std())
            ax.scatter([xs[-1]] * len(s), s, color=ARMC[arm], s=22, zorder=3,
                       edgecolor='white', linewidth=.6)
        ax.bar(xs, ys, w, yerr=es, capsize=4, color=ARMC[arm], alpha=.55,
               label=('+RNA' if arm == 'rna' else 'protein only'))
    ax.axhline(7, color='k', ls=':', lw=1)
    ax.text(len(ORDER) - .45, 9, 'equilibrated slab (~7 nm)', fontsize=7.5, ha='right')
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LABEL[c] for c in ORDER], fontsize=8.5)
    ax.set_ylabel('slab width after 2 µs (nm)')
    ax.set_title('Only MD1–MD3 stays condensed; everything with ART disperses\n'
                 '100 mM residues, n=3 seeds', fontsize=10)
    ax.legend(frameon=False, fontsize=8.5)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    save(fig, 'rnapanel_1_condensation')


def fig2_exchange(ex):
    fig, ax = plt.subplots(figsize=(7.5, 4.3))
    w = 0.36
    for k, arm in enumerate(['homotypic', 'rna']):
        xs, ys, es = [], [], []
        for i, c in enumerate(ORDER):
            s = ex[(ex.construct == c) & (ex.arm == arm)].retention_10ns
            if s.empty:
                continue
            xs.append(i + (k - 0.5) * w); ys.append(s.mean()); es.append(s.std())
            ax.scatter([xs[-1]] * len(s), s, color=ARMC[arm], s=22, zorder=3,
                       edgecolor='white', linewidth=.6)
        ax.bar(xs, ys, w, yerr=es, capsize=4, color=ARMC[arm], alpha=.55,
               label=('+RNA' if arm == 'rna' else 'protein only'))
    ax.axhspan(0, .25, color='#2ca02c', alpha=.10)
    ax.axhspan(.25, 1, color='grey', alpha=.10)
    ax.text(3.45, .12, 'LIQUID\ncontacts exchange', fontsize=7.5, ha='right', color='#1a7a1a')
    ax.text(3.45, .80, 'FROZEN\nsame pairs persist', fontsize=7.5, ha='right', color='#555')
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LABEL[c] for c in ORDER], fontsize=8.5)
    ax.set_ylabel('inter-chain contacts retained after ~10 ns')
    ax.set_ylim(0, 1)
    ax.set_title('Validity control: only MD1–MD3 is a liquid\n'
                 'the dispersed constructs keep the SAME contacts from frame 0 — '
                 'a kinetic trap, not an interaction', fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc='upper left')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    save(fig, 'rnapanel_2_contact_exchange')


def fig3_rna_domains(rna):
    r = rna[rna.arm == 'rna']
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=False)
    for ax, c in zip(axes, ORDER):
        s = r[r.construct == c]
        if s.empty:
            ax.axis('off'); continue
        g = s.groupby('domain').rna_contacts_per_residue.agg(['mean', 'std'])
        g = g.reindex([d for d in DOM_ORDER if d in g.index])
        cols = ['#C44E52' if d == 'MD1' else ('#DD8452' if d in ('KH7a', 'WWE') else '#8C8C8C')
                for d in g.index]
        ax.barh(range(len(g)), g['mean'], xerr=g['std'], color=cols, alpha=.85, capsize=2)
        ax.set_yticks(range(len(g)))
        ax.set_yticklabels(g.index, fontsize=8)
        ax.invert_yaxis()
        ax.set_title(LABEL[c].replace('\n', ' '), fontsize=9.5)
        ax.set_xlabel('RNA contacts per residue', fontsize=8.5)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    fig.suptitle('MD1 is the dominant RNA-contacting domain in every construct — '
                 'the RRMs are the weakest', fontsize=11, y=1.04)
    save(fig, 'rnapanel_3_rna_contacts_by_domain')


def fig4_mdfull_map(pairs):
    s = pairs[pairs.construct == 'md_full']
    doms = ['MD1', 'MD2', 'MD3']
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    mats = {}
    for arm in ['homotypic', 'rna']:
        M = np.zeros((3, 3))
        g = s[s.arm == arm].groupby(['domain_i', 'domain_j']).norm_per_residue_pair.mean()
        for i, a in enumerate(doms):
            for j, b in enumerate(doms):
                key = (a, b) if (a, b) in g.index else (b, a)
                M[i, j] = g.get(key, np.nan) * 1e6
        mats[arm] = M
    vmax = max(np.nanmax(mats['homotypic']), np.nanmax(mats['rna']))
    for ax, arm, t in zip(axes[:2], ['homotypic', 'rna'],
                          ['protein only', '+RNA']):
        im = ax.imshow(mats[arm], cmap='viridis', vmin=0, vmax=vmax)
        ax.set_title(t, fontsize=10)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f'{mats[arm][i,j]:.0f}', ha='center', va='center',
                        color='w' if mats[arm][i, j] < vmax * .6 else 'k', fontsize=9)
        ax.set_xticks(range(3)); ax.set_xticklabels(doms)
        ax.set_yticks(range(3)); ax.set_yticklabels(doms)
        fig.colorbar(im, ax=ax, fraction=.046, label='contacts / residue-pair (×10⁶)')
    d = mats['rna'] - mats['homotypic']
    lim = np.nanmax(np.abs(d))
    im = axes[2].imshow(d, cmap='RdBu_r', vmin=-lim, vmax=lim)
    axes[2].set_title('+RNA − protein only', fontsize=10)
    for i in range(3):
        for j in range(3):
            axes[2].text(j, i, f'{d[i,j]:+.0f}', ha='center', va='center', fontsize=9)
    axes[2].set_xticks(range(3)); axes[2].set_xticklabels(doms)
    axes[2].set_yticks(range(3)); axes[2].set_yticklabels(doms)
    fig.colorbar(im, ax=axes[2], fraction=.046, label='Δ contacts')
    fig.suptitle('Inter-chain domain contacts inside the MD1–MD3 condensate\n'
                 'RNA redistributes the network toward MD2–MD3', fontsize=11, y=1.06)
    save(fig, 'rnapanel_4_mdfull_interchain_map')


def main():
    cond = pd.read_csv(os.path.join(AN, 'condensation_summary.csv'))
    rna = pd.read_csv(os.path.join(AN, 'rna_contacts_by_domain.csv'))
    pairs = pd.read_csv(os.path.join(AN, 'interchain_domain_contacts.csv'))
    ex = pd.read_csv(os.path.join(AN, 'contact_exchange.csv'))
    fig1_condensation(cond)
    fig2_exchange(ex)
    fig3_rna_domains(rna)
    fig4_mdfull_map(pairs)
    print('figures in', FIG)


if __name__ == '__main__':
    main()
