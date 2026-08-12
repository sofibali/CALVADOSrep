#!/usr/bin/env python
"""
Figures from analyze_af3_calvados_modules.py outputs.

Reads:
  parp14/analysis/structure_analysis/aggregated/{site_sasa,site_wcn,site_saa,
        site_geometry,site_radial,site_contacts,combo_topology,combo_dmap,combo_fnc}.parquet

Writes (PNG + SVG) under parp14/analysis/structure_analysis/figures/:
  AS1_site_sasa.{png,svg}              violin per active site, all AF3 instances
  AS2_site_wcn.{png,svg}               violin per active site, WCN at 8 Å
  AS3_site_saa.{png,svg}               violin per active site, accessible_frac (CG-comparable)
  AS4_site_cone_angle.{png,svg}        violin per active site, max free-cone angle
  AS5_site_pair_distances.{png,svg}    inter-site COM distances (4x4 grid of histograms)
  AS6_site_radial_position.{png,svg}   site COM-to-chain-COM, normalised by Rg
  AS7_site_contacts_heatmap.{png,svg}  site x atS_domain mean contact count
  AS8_sasa_vs_saa.{png,svg}            scatter: atomistic SASA vs CG-style SAA per active site
  T1_rg_per_construct_size.{png,svg}   Rg distribution by construct size (boxplot)
  T2_dmap_FL.{png,svg}                 11x11 mean inter-domain COM distance (FL combo)
  T3_dmap_difference.{png,svg}         change vs FL: median(dmap_construct) - dmap_FL per pair
  T4_fnc_per_size.{png,svg}            FNC distribution vs construct size
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = Path('/home/sbali/CALVADOS/parp14')
A = ROOT / 'analysis' / 'structure_analysis' / 'aggregated'
F = ROOT / 'analysis' / 'structure_analysis' / 'figures'
F.mkdir(exist_ok=True)
ATS_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH_N', 'KH_Ca',
             'MD1', 'MD2', 'MD3', 'KH_b', 'WWE', 'ART']
SITE_ORDER = ['MD1', 'MD2', 'MD3', 'ART']
SITE_COLORS = {'MD1': '#e6194b', 'MD2': '#3cb44b', 'MD3': '#4363d8', 'ART': '#f58231'}
CMAP_BLUE = mcolors.LinearSegmentedColormap.from_list('WhiteBlue', ['#ffffff', '#08519c'])


def save(fig, name):
    fig.savefig(F / f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(F / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {name}.png + .svg')


def _violin(ax, data_dict, value_col, title, ylabel, sites=SITE_ORDER):
    data = [data_dict[s] for s in sites]
    parts = ax.violinplot(data, positions=range(len(sites)), showmedians=True,
                          showextrema=False, widths=0.85)
    for body, s in zip(parts['bodies'], sites):
        body.set_facecolor(SITE_COLORS[s]); body.set_alpha(0.5)
        body.set_edgecolor(SITE_COLORS[s])
    parts['cmedians'].set_color('black')
    ax.set_xticks(range(len(sites))); ax.set_xticklabels(sites)
    ax.set_ylabel(ylabel); ax.set_title(title)
    ax.grid(axis='y', alpha=0.3)


# ============================================================================
def plot_site_sasa():
    df = pd.read_parquet(A / 'site_sasa.parquet')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, col, title in [
        (axes[0], 'mean_cat_sasa_rel', 'Catalytic-residue mean rel. SASA'),
        (axes[1], 'frac_cat_exposed', 'Fraction of catalytic residues exposed (rel.SASA>0.20)'),
    ]:
        d = {s: df.loc[df.site == s, col].dropna().values for s in SITE_ORDER}
        _violin(ax, d, col, title, col)
    fig.suptitle('Atomistic active-site exposure (freesasa)\nall AF3 models that contain each site')
    fig.tight_layout()
    save(fig, 'AS1_site_sasa')


def plot_site_wcn():
    df = pd.read_parquet(A / 'site_wcn.parquet')
    fig, ax = plt.subplots(figsize=(7, 5))
    d = {s: df.loc[df.site == s, 'wcn_8A'].dropna().values for s in SITE_ORDER}
    _violin(ax, d, 'wcn_8A',
            'Weighted coord. number at active site (8 Å, |i-j|>10 excluded)',
            'WCN_8Å (sum 1/r over residues)')
    fig.tight_layout()
    save(fig, 'AS2_site_wcn')


def plot_site_saa():
    df = pd.read_parquet(A / 'site_saa.parquet')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    d_acc = {s: df.loc[df.site == s, 'accessible_frac'].dropna().values for s in SITE_ORDER}
    d_cone = {s: df.loc[df.site == s, 'max_cone_deg'].dropna().values for s in SITE_ORDER}
    _violin(axes[0], d_acc, 'accessible_frac',
            'SAA — fraction of unblocked rays (probe 5 Å, max 50 Å)\n'
            'matches CG analyze_accessibility.py params',
            'Accessible fraction (0–1)')
    _violin(axes[1], d_cone, 'max_cone_deg',
            'Max free-cone half-angle around site COM',
            'Cone half-angle (deg)')
    fig.tight_layout()
    save(fig, 'AS3_site_saa')


def plot_site_distances():
    df = pd.read_parquet(A / 'site_geometry.parquet')
    pairs = sorted(df[['site_i', 'site_j']].drop_duplicates().itertuples(index=False))
    n = len(pairs)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.5 * cols, 3 * rows))
    axes_flat = axes.flat if rows > 1 else axes
    for ax, (si, sj) in zip(axes_flat, pairs):
        sub = df[(df.site_i == si) & (df.site_j == sj)]
        ax.hist(sub.value, bins=50, color='#1f77b4', alpha=0.7)
        ax.set_title(f'{si} ↔ {sj}  (n={len(sub)})')
        ax.set_xlabel('COM distance (Å)')
        ax.grid(alpha=0.3)
    for ax in list(axes_flat)[len(pairs):]:
        ax.axis('off')
    fig.suptitle('Inter-active-site COM distances across AF3 ensemble')
    fig.tight_layout()
    save(fig, 'AS5_site_pair_distances')


def plot_site_radial():
    df = pd.read_parquet(A / 'site_radial.parquet')
    fig, ax = plt.subplots(figsize=(7, 5))
    d = {s: df.loc[df.site_i == s, 'ratio'].dropna().values for s in SITE_ORDER}
    _violin(ax, d, 'ratio', 'Site COM distance to chain COM, normalised by Rg',
            'd(site COM, chain COM) / Rg', sites=SITE_ORDER)
    ax.axhline(1.0, color='gray', ls='--', alpha=0.5)
    fig.tight_layout()
    save(fig, 'AS6_site_radial_position')


def plot_site_contacts_heatmap():
    df = pd.read_parquet(A / 'site_contacts.parquet')
    M = np.zeros((len(SITE_ORDER), len(ATS_ORDER)))
    for i, s in enumerate(SITE_ORDER):
        sub = df[df.site_i == s]
        means = sub.groupby('site_j')['value'].mean()
        for j, d in enumerate(ATS_ORDER):
            if d in means.index:
                M[i, j] = means[d]
    fig, ax = plt.subplots(figsize=(9, 4))
    im = ax.imshow(M, cmap=CMAP_BLUE, aspect='auto')
    ax.set_xticks(range(len(ATS_ORDER))); ax.set_xticklabels(ATS_ORDER, rotation=45)
    ax.set_yticks(range(len(SITE_ORDER))); ax.set_yticklabels(SITE_ORDER)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v >= 0.5:
                ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=7,
                        color='white' if v > np.nanmax(M) * 0.6 else 'black')
    cbar = plt.colorbar(im, ax=ax, fraction=0.04)
    cbar.set_label('mean # contacting residues (Cα<10 Å, |i-j|>10)')
    ax.set_title('Per-active-site contact profile from each atS domain')
    fig.tight_layout()
    save(fig, 'AS7_site_contacts_heatmap')


def plot_sasa_vs_saa():
    sasa = pd.read_parquet(A / 'site_sasa.parquet')
    saa = pd.read_parquet(A / 'site_saa.parquet')
    merged = sasa.merge(saa, on=['combo', 'model', 'site'])
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
    for ax, s in zip(axes, SITE_ORDER):
        sub = merged[merged.site == s]
        if sub.empty: ax.axis('off'); continue
        ax.scatter(sub.mean_cat_sasa_rel, sub.accessible_frac,
                   s=4, c=SITE_COLORS[s], alpha=0.3)
        ax.set_xlabel('Atomistic mean rel. SASA (catalytic residues)')
        ax.set_title(f'{s}  (n={len(sub)})')
        ax.grid(alpha=0.3)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
    axes[0].set_ylabel('CG-style SAA (frac unblocked rays)')
    fig.suptitle('Atomistic SASA vs CG-style SAA per active site')
    fig.tight_layout()
    save(fig, 'AS8_sasa_vs_saa')


# ============================================================================
def plot_topology_rg():
    top = pd.read_parquet(A / 'combo_topology.parquet')
    rg = top[top.kind == 'rg'].copy()
    # construct size from n_res column if available else compute
    fig, ax = plt.subplots(figsize=(8, 5))
    rg['size_bin'] = pd.cut(rg.n_res, bins=[0, 300, 600, 900, 1200, 1500, 1900],
                            labels=['<300', '300-600', '600-900',
                                    '900-1200', '1200-1500', '1500+'])
    groups = [rg.loc[rg.size_bin == b, 'value'].values
              for b in rg.size_bin.cat.categories]
    ax.boxplot(groups, labels=list(rg.size_bin.cat.categories), showfliers=False,
               patch_artist=True,
               boxprops=dict(facecolor='#1f77b4', alpha=0.4),
               medianprops=dict(color='black'))
    ax.set_xlabel('Construct size (# residues)')
    ax.set_ylabel('Rg (Å)')
    ax.set_title('Construct compactness vs size')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, 'T1_rg_per_construct_size')


def plot_dmap_FL_and_diff():
    dmap = pd.read_parquet(A / 'combo_dmap.parquet')
    # FL combo dmap (largest n_res via combo with all 11 atS domains)
    # Identify FL combo by combo containing all atS domains
    combo_doms = dmap.groupby('combo')[['a', 'b']].apply(
        lambda g: set(g.a) | set(g.b)).rename('doms')
    fl_combos = combo_doms[combo_doms.apply(lambda s: set(ATS_ORDER).issubset(s))].index.tolist()
    if not fl_combos:
        print('  no FL-spanning combo found in dmap'); return
    fl_combo = fl_combos[0]
    fl_means = dmap[dmap.combo == fl_combo].groupby(['a', 'b'])['value'].mean()
    M_fl = np.full((len(ATS_ORDER), len(ATS_ORDER)), np.nan)
    for (a, b), v in fl_means.items():
        i, j = ATS_ORDER.index(a), ATS_ORDER.index(b)
        M_fl[i, j] = v; M_fl[j, i] = v

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(M_fl, cmap=CMAP_BLUE, aspect='auto')
    for i in range(len(ATS_ORDER)):
        for j in range(len(ATS_ORDER)):
            if i == j: continue
            v = M_fl[i, j]
            if not np.isnan(v) and v < 80:
                ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=6,
                        color='white' if v > np.nanmax(M_fl) * 0.6 else 'black')
    ax.set_xticks(range(len(ATS_ORDER))); ax.set_xticklabels(ATS_ORDER, rotation=45)
    ax.set_yticks(range(len(ATS_ORDER))); ax.set_yticklabels(ATS_ORDER)
    plt.colorbar(im, ax=ax, fraction=0.04, label='mean COM distance (Å)')
    ax.set_title(f'FL combo inter-atS-domain COM distance map\n{fl_combo}')
    fig.tight_layout()
    save(fig, 'T2_dmap_FL')

    # Difference: median(other constructs) - FL
    other = dmap[dmap.combo != fl_combo].groupby(['a', 'b'])['value'].median()
    M_other = np.full((len(ATS_ORDER), len(ATS_ORDER)), np.nan)
    for (a, b), v in other.items():
        if a in ATS_ORDER and b in ATS_ORDER:
            i, j = ATS_ORDER.index(a), ATS_ORDER.index(b)
            M_other[i, j] = v; M_other[j, i] = v
    M_diff = M_other - M_fl
    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.nanmax(np.abs(M_diff))
    im = ax.imshow(M_diff, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
    for i in range(len(ATS_ORDER)):
        for j in range(len(ATS_ORDER)):
            if i == j: continue
            v = M_diff[i, j]
            if not np.isnan(v):
                ax.text(j, i, f'{v:+.0f}', ha='center', va='center', fontsize=6,
                        color='black')
    ax.set_xticks(range(len(ATS_ORDER))); ax.set_xticklabels(ATS_ORDER, rotation=45)
    ax.set_yticks(range(len(ATS_ORDER))); ax.set_yticklabels(ATS_ORDER)
    plt.colorbar(im, ax=ax, fraction=0.04, label='Δ COM distance (Å, median other - FL)')
    ax.set_title('Inter-domain distance change vs FL\n(positive = farther in other constructs)')
    fig.tight_layout()
    save(fig, 'T3_dmap_difference')


def plot_fnc_per_size():
    fnc = pd.read_parquet(A / 'combo_fnc.parquet')
    fnc = fnc[fnc.kind == 'fnc'].copy()
    # bring n_res from combo_topology
    rg = pd.read_parquet(A / 'combo_topology.parquet')
    rg = rg[rg.kind == 'rg'][['combo', 'model', 'n_res']]
    fnc = fnc.merge(rg, on=['combo', 'model'], how='left')
    fnc['size_bin'] = pd.cut(fnc.n_res, bins=[0, 300, 600, 900, 1200, 1500, 1900],
                             labels=['<300', '300-600', '600-900',
                                     '900-1200', '1200-1500', '1500+'])
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = [fnc.loc[fnc.size_bin == b, 'value'].values
              for b in fnc.size_bin.cat.categories]
    ax.boxplot(groups, labels=list(fnc.size_bin.cat.categories), showfliers=False,
               patch_artist=True,
               boxprops=dict(facecolor='#2ca02c', alpha=0.4),
               medianprops=dict(color='black'))
    ax.set_xlabel('Construct size (# residues)')
    ax.set_ylabel('FNC vs FL_AF3 (Cα<10 Å)')
    ax.set_title('Fraction of FL native contacts preserved in each construct')
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, 'T4_fnc_per_size')


def main():
    print('Loading + plotting...')
    if (A / 'site_sasa.parquet').exists():
        plot_site_sasa()
    if (A / 'site_wcn.parquet').exists():
        plot_site_wcn()
    if (A / 'site_saa.parquet').exists():
        plot_site_saa()
    if (A / 'site_geometry.parquet').exists():
        plot_site_distances()
    if (A / 'site_radial.parquet').exists():
        plot_site_radial()
    if (A / 'site_contacts.parquet').exists():
        plot_site_contacts_heatmap()
    if (A / 'site_sasa.parquet').exists() and (A / 'site_saa.parquet').exists():
        plot_sasa_vs_saa()
    if (A / 'combo_topology.parquet').exists():
        plot_topology_rg()
    if (A / 'combo_dmap.parquet').exists():
        plot_dmap_FL_and_diff()
    if (A / 'combo_fnc.parquet').exists():
        plot_fnc_per_size()
    print(f'Figures in {F}')


if __name__ == '__main__':
    main()
