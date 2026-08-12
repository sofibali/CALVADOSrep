#!/usr/bin/env python
"""
Plots from the structure analysis tables.

Reads:
    parp14/analysis/structure_analysis/aggregated/{all_summary,all_dompairs,all_kpairs}.parquet
    parp14/analysis/structure_analysis/crystal/crystal_*.parquet

Writes (PNG + SVG):
    parp14/analysis/structure_analysis/figures/
        01_rmsd_per_domain.{png,svg}        violin per atS domain, crystals overlaid
        02_rmsd_vs_size.{png,svg}           construct size vs domain RMSD (per domain)
        03_dompair_contacts_mean.{png,svg}  11x11 heatmap of mean sc8 inter-domain contacts
        04_dompair_contacts_exposed.{png,svg} 11x11 mean exposed sc8 contacts
        05_KE_KD_dompair.{png,svg}          11x11 mean K-acidic pair count (exposed)
        06_K_acidic_FL_heatmap.{png,svg}    FL residue x FL residue, persistence of K-acidic exposed contacts
        07_top_KE_KD_residue_pairs.{png,svg} bar chart of top K-D/E pairs (residue level)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = Path('/home/sbali/CALVADOS/parp14')
A = ROOT / 'analysis' / 'structure_analysis' / 'aggregated'
X = ROOT / 'analysis' / 'structure_analysis' / 'crystal'
F = ROOT / 'analysis' / 'structure_analysis' / 'figures'
F.mkdir(exist_ok=True, parents=True)

ATS_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH_N', 'KH_Ca',
             'MD1', 'MD2', 'MD3', 'KH_b', 'WWE', 'ART']
CMAP_BLUE = mcolors.LinearSegmentedColormap.from_list('WhiteBlue', ['#ffffff', '#08519c'])
CRYSTAL_COLOR = '#d62728'   # red
AF3_COLOR = '#1f77b4'       # blue

XTAL_LABELS = {'1x4r': '1X4R', '3goy': '3GOY', '3q6z': '3Q6Z', '3vfq': '3VFQ'}


def save(fig, name):
    fig.savefig(F / f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(F / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {name}.png + .svg')


# ============================================================================
# 1. RMSD per atS domain (violin), crystals overlaid
# ============================================================================
def plot_rmsd_per_domain(summary):
    af3 = summary[summary['source'] == 'af3']
    xtal = summary[summary['source'] == 'crystal']
    fig, ax = plt.subplots(figsize=(11, 5))
    data = [af3.loc[af3['atS_domain'] == d, 'rmsd_to_FLref'].values for d in ATS_ORDER]
    parts = ax.violinplot(data, positions=range(len(ATS_ORDER)), showmedians=True,
                          showextrema=False, widths=0.85)
    for body in parts['bodies']:
        body.set_facecolor(AF3_COLOR); body.set_alpha(0.35)
        body.set_edgecolor(AF3_COLOR)
    parts['cmedians'].set_color(AF3_COLOR)
    # Crystal overlay
    for _, r in xtal.iterrows():
        if r['atS_domain'] not in ATS_ORDER: continue
        x = ATS_ORDER.index(r['atS_domain'])
        ax.scatter(x, r['rmsd_to_FLref'], s=80, c=CRYSTAL_COLOR, marker='D',
                   edgecolor='black', linewidth=0.7, zorder=5,
                   label='crystal' if x == 0 or 'crystal' not in [t.get_label() for t in ax.collections] else None)
        ax.text(x + 0.15, r['rmsd_to_FLref'], XTAL_LABELS.get(r['combo'], r['combo'].upper()),
                fontsize=8, va='center', color=CRYSTAL_COLOR)
    ax.set_xticks(range(len(ATS_ORDER)))
    ax.set_xticklabels(ATS_ORDER, rotation=30)
    ax.set_ylabel('RMSD to FL_AF3 reference (Å)')
    ax.set_title('Per-domain Cα RMSD across all AF3 constructs (crystals = red diamonds)')
    ax.grid(axis='y', alpha=0.3)
    save(fig, '01_rmsd_per_domain')


# ============================================================================
# 2. RMSD vs construct size (per domain)
# ============================================================================
def plot_rmsd_vs_size(summary):
    af3 = summary[summary['source'] == 'af3'].copy()
    # construct size = #atS domains in combo
    sizes = af3.groupby('combo')['atS_domain'].nunique().rename('n_domains_in_combo')
    af3 = af3.join(sizes, on='combo')
    fig, axes = plt.subplots(3, 4, figsize=(14, 9), sharex=True)
    for ax, d in zip(axes.flat, ATS_ORDER + [None]):
        if d is None:
            ax.axis('off'); continue
        sub = af3[af3['atS_domain'] == d]
        if sub.empty: ax.axis('off'); continue
        # Boxplot per n_domains
        groups = sorted(sub['n_domains_in_combo'].unique())
        data = [sub.loc[sub['n_domains_in_combo'] == g, 'rmsd_to_FLref'].values for g in groups]
        ax.boxplot(data, positions=groups, widths=0.6, showfliers=False,
                   patch_artist=True,
                   boxprops=dict(facecolor=AF3_COLOR, alpha=0.4),
                   medianprops=dict(color='black'))
        ax.set_title(d); ax.grid(alpha=0.3)
    fig.supxlabel('# atS domains in construct')
    fig.supylabel('RMSD to FL_AF3 (Å)')
    fig.suptitle('Domain RMSD vs construct size')
    fig.tight_layout()
    save(fig, '02_rmsd_vs_size')


# ============================================================================
# 3. & 4. & 5. Domain-pair heatmaps
# ============================================================================
def _dompair_matrix(dp, value_col):
    """Build symmetric NxN matrix (N=11) where M[i,j] = mean value across all instances
       containing both dom_i and dom_j."""
    N = len(ATS_ORDER)
    M = np.full((N, N), np.nan)
    g = dp.groupby(['dom_i', 'dom_j'])[value_col].mean()
    for (di, dj), v in g.items():
        if di in ATS_ORDER and dj in ATS_ORDER:
            i = ATS_ORDER.index(di); j = ATS_ORDER.index(dj)
            M[i, j] = v; M[j, i] = v
    return M


def plot_dompair_heatmaps(dompairs):
    af3 = dompairs[dompairs['source'] == 'af3']
    for col, fname, title in [
        ('sc_n_8A', '03_dompair_contacts_mean',
         'Mean inter-domain contacts (sc heavy <8 Å) across all constructs'),
        ('sc_n_8A_exposed', '04_dompair_contacts_exposed',
         'Mean exposed-only inter-domain contacts (sc <8 Å, both res rel-SASA>20%)'),
        ('sc8_KE_or_KD_exposed', '05_KE_KD_dompair',
         'Mean exposed K-acidic inter-domain contacts (sc <8 Å)'),
    ]:
        M = _dompair_matrix(af3, col)
        fig, ax = plt.subplots(figsize=(7.5, 6.2))
        im = ax.imshow(M, cmap=CMAP_BLUE, aspect='auto')
        # gray for NaN (domain pair never co-occurred)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                               facecolor='#e8e8e8', edgecolor='none'))
                elif M[i, j] >= 0.5:
                    ax.text(j, i, f'{M[i, j]:.0f}', ha='center', va='center',
                            fontsize=7, color='black' if M[i, j] < np.nanmax(M) * 0.6 else 'white')
        ax.set_xticks(range(len(ATS_ORDER))); ax.set_xticklabels(ATS_ORDER, rotation=45)
        ax.set_yticks(range(len(ATS_ORDER))); ax.set_yticklabels(ATS_ORDER)
        ax.set_title(title, fontsize=10)
        cbar = plt.colorbar(im, ax=ax, fraction=0.04)
        cbar.set_label('mean count per construct')
        save(fig, fname)


# ============================================================================
# 6. Symmetric residue × residue heatmaps for K-acidic and K-K contacts
# ============================================================================
def _build_symmetric_heatmap(sub, bin_size=25, n_fl=1801):
    """Build a symmetric residue×residue contact-count matrix in `bin_size` bins.
    `sub` must have res_i_fl and res_j_fl columns. Each contact is counted ONCE
    (canonicalised so i<=j) and the matrix is mirrored across the diagonal."""
    i = sub.res_i_fl.values; j = sub.res_j_fl.values
    lo = np.minimum(i, j); hi = np.maximum(i, j)
    ibin = (lo - 1) // bin_size
    jbin = (hi - 1) // bin_size
    nb = (n_fl // bin_size) + 1
    M = np.zeros((nb, nb))
    np.add.at(M, (ibin, jbin), 1)        # upper triangle holds counts
    # Mirror: keep diagonal as-is (within-bin self-contacts) and reflect off-diagonal
    M_sym = M + M.T
    di = np.diag_indices(nb)
    M_sym[di] = M[di]                    # remove diagonal double-count
    return M_sym, nb


def _annotate_domain_boundaries(ax, nb, bin_size):
    import json
    atsmap = json.load(open(ROOT / 'analysis' / 'structure_analysis' / 'reference' / 'atS_map.json'))
    for d in ATS_ORDER:
        s, e = atsmap['atS_ranges'][d]
        mid = (s + e) / 2
        ax.axvline(e, color='gray', alpha=0.35, lw=0.4)
        ax.axhline(e, color='gray', alpha=0.35, lw=0.4)
        ax.text(mid, -nb * bin_size * 0.04, d, ha='center', va='top',
                fontsize=7, rotation=45)
        ax.text(-nb * bin_size * 0.04, mid, d, ha='right', va='center', fontsize=7)


def plot_K_pair_residue_heatmaps(kpairs, bin_size=25):
    """Symmetric residue×residue heatmaps for both K-D/E and K-K, in linear and log."""
    af3 = kpairs[kpairs['source'] == 'af3']
    f = af3[af3.exposed_i & af3.exposed_j & (af3.sc_dist_A < 8.0)]

    panels = [
        ('K-acidic',
         f[(((f.aa_i == 'K') & f.aa_j.isin(['D', 'E'])) |
            ((f.aa_j == 'K') & f.aa_i.isin(['D', 'E'])))],
         'Exposed K–D/E contacts (sc<8Å)'),
        ('K-K',
         f[(f.aa_i == 'K') & (f.aa_j == 'K')],
         'Exposed K–K contacts (sc<8Å)'),
    ]

    for label, sub, title in panels:
        if sub.empty:
            print(f'  no {label} contacts; skipping'); continue
        M, nb = _build_symmetric_heatmap(sub, bin_size=bin_size)
        n_combos = sub.combo.nunique()

        for scale_name, transform, cbar_label in [
            ('linear', lambda x: x, 'count of contact instances (across all constructs × models)'),
            ('log', lambda x: np.log1p(x),
             'log(1 + n contacts) — log compresses heavy-tailed counts so rare\n'
             'pairs are visible alongside the most common ones'),
        ]:
            fig, ax = plt.subplots(figsize=(8, 7.2))
            im = ax.imshow(transform(M), cmap=CMAP_BLUE, origin='lower', aspect='equal',
                           extent=[0, nb * bin_size, 0, nb * bin_size])
            cbar = plt.colorbar(im, ax=ax, fraction=0.04)
            cbar.set_label(cbar_label)
            _annotate_domain_boundaries(ax, nb, bin_size)
            ax.set_xlabel('FL residue position')
            ax.set_ylabel('FL residue position')
            ax.set_title(f'{title}, symmetric residue×residue map\n'
                         f'pooled across {n_combos} constructs, {bin_size}-residue bins '
                         f'({scale_name} scale)', fontsize=10)
            tag = '06_K_acidic' if label == 'K-acidic' else '06b_K_K'
            save(fig, f'{tag}_FL_heatmap_{scale_name}')


# ============================================================================
# 7. Top K-acidic residue pairs (bar chart)
# ============================================================================
def plot_top_KE_KD_pairs(kpairs, top_n=30):
    af3 = kpairs[kpairs['source'] == 'af3']
    is_ke = (((af3.aa_i == 'K') & af3.aa_j.isin(['D', 'E'])) |
             ((af3.aa_j == 'K') & af3.aa_i.isin(['D', 'E'])))
    sub = af3[is_ke & af3.exposed_i & af3.exposed_j & (af3.sc_dist_A < 8.0)]
    if sub.empty:
        print('  no K-acidic exposed pairs; skipping 07')
        return
    Kres = np.where(sub.aa_i == 'K', sub.res_i_fl, sub.res_j_fl)
    Ares = np.where(sub.aa_i == 'K', sub.res_j_fl, sub.res_i_fl)
    Aaa = np.where(sub.aa_i == 'K', sub.aa_j, sub.aa_i)
    df = pd.DataFrame({'K': Kres, 'A': Ares, 'Aaa': Aaa,
                       'combo': sub.combo.values})
    # Count per pair across UNIQUE combos (= how many constructs see this pair)
    pair_combos = df.groupby(['K', 'A', 'Aaa'])['combo'].nunique().reset_index(name='n_combos')
    pair_total = df.groupby(['K', 'A', 'Aaa']).size().reset_index(name='n_models')
    merged = pair_combos.merge(pair_total, on=['K', 'A', 'Aaa']).sort_values('n_combos', ascending=False)
    top = merged.head(top_n).copy()
    top['label'] = top.apply(lambda r: f'K{int(r.K)}–{r.Aaa}{int(r.A)}', axis=1)

    fig, ax = plt.subplots(figsize=(10, max(5, top_n * 0.25)))
    bars = ax.barh(top['label'][::-1], top['n_combos'][::-1], color=AF3_COLOR, alpha=0.8)
    for i, (_, r) in enumerate(top[::-1].iterrows()):
        ax.text(r['n_combos'] + 5, i, f'{int(r["n_models"])} models',
                va='center', fontsize=7, color='gray')
    ax.set_xlabel('# constructs in which the contact appears (any model)')
    ax.set_title(f'Top {top_n} exposed K–D/E pairs (sc<8 Å), ranked by construct frequency')
    ax.grid(axis='x', alpha=0.3)
    save(fig, '07_top_KE_KD_residue_pairs')


def main():
    print('Loading aggregated tables...')
    summary = pd.read_parquet(A / 'all_summary.parquet')
    dompairs = pd.read_parquet(A / 'all_dompairs.parquet')
    kpairs = pd.read_parquet(A / 'all_kpairs.parquet')
    print(f'  summary={summary.shape}, dompairs={dompairs.shape}, kpairs={kpairs.shape}')

    print('Plot 1: RMSD per domain (violin)')
    plot_rmsd_per_domain(summary)
    print('Plot 2: RMSD vs construct size')
    plot_rmsd_vs_size(summary)
    print('Plots 3-5: Inter-domain contact heatmaps')
    plot_dompair_heatmaps(dompairs)
    print('Plots 6 / 6b: K-acidic and K-K residue×residue heatmaps (linear + log)')
    plot_K_pair_residue_heatmaps(kpairs)
    print('Plot 7: Top K-D/E residue pairs')
    plot_top_KE_KD_pairs(kpairs, top_n=40)
    print(f'\nFigures in {F}')


if __name__ == '__main__':
    main()
