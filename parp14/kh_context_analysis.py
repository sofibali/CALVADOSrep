#!/usr/bin/env python
"""
E (revised) — KH-domain stability conditioned on KH7a presence.

Two specific comparisons:
  (i) atS KH_b RMSD with vs without KH7a in construct.
      atS KH_b covers FL residues 1392-1500 (= KHb proper).
  (ii) "KH8 region" RMSD with vs without KH7a.
      Defined as FL residues 1501-1533 — these are inside CSV unit `khb-kh8` but
      outside the atS KH_b boundary, so KH8 doesn't appear in all_summary. We
      compute its Cα RMSD vs the FL reference on-the-fly from per_combo residues.

Also: a heatmap showing mean RMSD for ALL atS domains conditioned on each CSV
unit being present vs absent (rows = atS domains, cols = CSV units, value =
ΔRMSD = mean(without unit) - mean(with unit)). Positive ΔRMSD means the domain
is more stable when that unit is present.

Inputs:  all_summary.parquet, per_combo/<combo>__residues.parquet
Outputs: figures/E1_KHb_KH8_vs_KH7a.{png,svg}
         figures/E2_atS_RMSD_by_unit_presence.{png,svg}
         aggregated/kh8_rmsd.parquet
"""
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/home/sbali/CALVADOS/parp14')
PC = ROOT / 'analysis' / 'structure_analysis' / 'per_combo'
A = ROOT / 'analysis' / 'structure_analysis' / 'aggregated'
REF_DIR = ROOT / 'analysis' / 'structure_analysis' / 'reference'
F = ROOT / 'analysis' / 'structure_analysis' / 'figures'
F.mkdir(exist_ok=True)

ATS_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH_N', 'KH_Ca',
             'MD1', 'MD2', 'MD3', 'KH_b', 'WWE', 'ART']

# CSV unit tokens (lowercase, as they appear in combo dir names after md1l1->md1)
CSV_UNIT_ORDER = ['rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a', 'md1', 'md2', 'md3',
                  'khb-kh8', 'wwe', 'art']
CSV_UNIT_LABELS = {'rrm1': 'RRM1', 'rrm2': 'RRM2', 'rrm3': 'RRM3',
                   'kh1-kh6': 'KH1-KH6', 'kh7a': 'KH7a', 'md1': 'MD1L1',
                   'md2': 'MD2', 'md3': 'MD3', 'khb-kh8': 'KHb-KH8',
                   'wwe': 'WWE', 'art': 'ART'}

# KH8 region boundary (FL residues 1501-1533) — within KHb-KH8 CSV unit, outside atS KH_b
KH8_FL_START = 1501
KH8_FL_END = 1533


def kabsch_rmsd(P, Q):
    if P.shape != Q.shape or P.shape[0] < 3:
        return np.nan
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, axis=1))))


def parse_units_from_combo(name):
    """Extract CSV unit tokens from a combo dir name (lowercase)."""
    parts = name.lower().split('_')
    valid = set(CSV_UNIT_ORDER)
    return set(p for p in parts if p in valid)


def compute_kh8_ref():
    """Get FL_AF3 reference Cα for residues KH8_FL_START..KH8_FL_END from the FL combo's
    seed-1_sample-0 residues parquet."""
    am = json.load(open(REF_DIR / 'atS_map.json'))
    fl_dir_name = Path(am['fl_combo_dir']).name
    # Strip timestamp suffix to get canonical combo name (file naming uses canonical)
    canonical = re.sub(r'_2\d{7}_\d{6}$', '', fl_dir_name)
    fl_residues = PC / f'{canonical}__residues.parquet'
    if not fl_residues.exists():
        # fall back to glob
        cands = list(PC.glob('*kh1-kh6*kh7a*md1*md2*md3*khb-kh8*wwe*art*__residues.parquet'))
        if not cands:
            raise RuntimeError('FL combo residues parquet not found')
        fl_residues = cands[0]
    df = pd.read_parquet(fl_residues)
    sub = df[(df.model == 'seed-1_sample-0') &
             (df.res_fl >= KH8_FL_START) & (df.res_fl <= KH8_FL_END)].sort_values('res_fl')
    if len(sub) < 5:
        raise RuntimeError(f'FL ref has only {len(sub)} KH8 residues')
    coords = sub[['ca_x', 'ca_y', 'ca_z']].values
    res_fl = sub['res_fl'].values
    return res_fl, coords


def compute_kh8_rmsd_for_combos():
    """For each combo containing the KHb-KH8 CSV unit, compute KH8 Cα RMSD vs FL ref
    for every model. Returns DataFrame: combo, model, n_res, kh8_rmsd, has_kh7a."""
    ref_res, ref_coords = compute_kh8_ref()
    print(f'  KH8 reference: {len(ref_res)} residues, FL {ref_res.min()}-{ref_res.max()}')

    rows = []
    files = sorted(PC.glob('*__residues.parquet'))
    n_processed = 0
    for p in files:
        combo = p.name.replace('__residues.parquet', '')
        units = parse_units_from_combo(combo)
        if 'khb-kh8' not in units:
            continue
        has_kh7a = 'kh7a' in units
        df = pd.read_parquet(p, columns=['model', 'res_fl', 'ca_x', 'ca_y', 'ca_z'])
        sub = df[(df.res_fl >= KH8_FL_START) & (df.res_fl <= KH8_FL_END)]
        for mdl, g in sub.groupby('model'):
            g = g.sort_values('res_fl')
            if len(g) < 5:
                continue
            # Match residues to ref
            common = np.isin(g['res_fl'].values, ref_res)
            if common.sum() < 5:
                continue
            ref_idx = np.array([np.where(ref_res == r)[0][0] for r in g['res_fl'].values[common]])
            P = g[['ca_x', 'ca_y', 'ca_z']].values[common]
            Q = ref_coords[ref_idx]
            rows.append({
                'combo': combo, 'model': mdl,
                'n_res': int(len(P)), 'kh8_rmsd': kabsch_rmsd(P, Q),
                'has_kh7a': has_kh7a,
            })
        n_processed += 1
    df = pd.DataFrame(rows)
    out = A / 'kh8_rmsd.parquet'
    df.to_parquet(out, index=False)
    print(f'  saved {out} ({len(df)} model rows from {n_processed} KHb-KH8 combos)')
    return df


def plot_KHb_KH8_vs_KH7a(kh8_df):
    """E1: violin/box of KH_b RMSD and KH8 RMSD, split by KH7a presence."""
    summary = pd.read_parquet(A / 'all_summary.parquet')
    af3 = summary[summary['source'] == 'af3']

    # KH_b: from summary, label combos by KH7a presence
    khb = af3[af3.atS_domain == 'KH_b'].copy()
    khb['has_kh7a'] = khb['combo'].apply(lambda c: 'kh7a' in parse_units_from_combo(c))

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (label, df, val_col) in zip(axes, [
        ('atS KH_b (FL 1392–1500)', khb, 'rmsd_to_FLref'),
        ('KH8 region (FL 1501–1533)', kh8_df, 'kh8_rmsd'),
    ]):
        with_kh = df[df.has_kh7a][val_col].dropna().values
        without_kh = df[~df.has_kh7a][val_col].dropna().values
        # Mann-Whitney U test
        try:
            from scipy.stats import mannwhitneyu
            stat, pval = mannwhitneyu(with_kh, without_kh, alternative='two-sided')
        except Exception:
            pval = np.nan
        parts = ax.violinplot([with_kh, without_kh], positions=[0, 1],
                              showmedians=True, showextrema=False, widths=0.85)
        colors = ['#1f77b4', '#ff7f0e']
        for body, c in zip(parts['bodies'], colors):
            body.set_facecolor(c); body.set_alpha(0.4); body.set_edgecolor(c)
        parts['cmedians'].set_color('black')
        # Overlay sample means
        for x, vals, c in [(0, with_kh, colors[0]), (1, without_kh, colors[1])]:
            ax.scatter([x], [np.mean(vals)], marker='D', s=70, color=c,
                       edgecolor='black', linewidth=0.7, zorder=10)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['KH7a present', 'KH7a absent'])
        ax.set_ylabel('Cα RMSD to FL_AF3 (Å)')
        ax.set_title(f'{label}\nN_with={len(with_kh)}, N_without={len(without_kh)}, '
                     f'p={pval:.2e}')
        ax.grid(axis='y', alpha=0.3)
        # Annotate medians
        ax.text(0, np.median(with_kh), f' med={np.median(with_kh):.2f}',
                fontsize=8, va='center', ha='left')
        ax.text(1, np.median(without_kh), f' med={np.median(without_kh):.2f}',
                fontsize=8, va='center', ha='left')
    fig.suptitle('KH_b and KH8 stability: effect of KH7a presence', fontsize=12)
    fig.tight_layout()
    fig.savefig(F / 'E1_KHb_KH8_vs_KH7a.png', dpi=200, bbox_inches='tight')
    fig.savefig(F / 'E1_KHb_KH8_vs_KH7a.svg', bbox_inches='tight')
    plt.close(fig)
    print('  wrote E1_KHb_KH8_vs_KH7a.png + .svg')


def plot_atS_by_unit_presence():
    """E2: heatmap rows=atS domains, cols=CSV units, value=ΔRMSD (without - with)."""
    summary = pd.read_parquet(A / 'all_summary.parquet')
    af3 = summary[summary['source'] == 'af3'].copy()
    # mark presence of each CSV unit per combo
    unit_present = {u: af3.combo.apply(lambda c: u in parse_units_from_combo(c)) for u in CSV_UNIT_ORDER}
    M = np.full((len(ATS_ORDER), len(CSV_UNIT_ORDER)), np.nan)
    Pmat = np.full_like(M, np.nan)
    for i, d in enumerate(ATS_ORDER):
        sub = af3[af3.atS_domain == d]
        for j, u in enumerate(CSV_UNIT_ORDER):
            present = unit_present[u][sub.index]
            with_ = sub.loc[present, 'rmsd_to_FLref'].values
            without_ = sub.loc[~present, 'rmsd_to_FLref'].values
            if len(with_) < 30 or len(without_) < 30:
                continue
            M[i, j] = np.median(without_) - np.median(with_)
            try:
                from scipy.stats import mannwhitneyu
                _, p = mannwhitneyu(with_, without_, alternative='two-sided')
                Pmat[i, j] = p
            except Exception:
                pass
    # Heatmap, diverging colormap centered at 0
    fig, ax = plt.subplots(figsize=(9, 6.5))
    vmax = np.nanmax(np.abs(M))
    im = ax.imshow(M, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           facecolor='#e8e8e8', edgecolor='none'))
                continue
            # text cell
            sig = ' *' if (Pmat[i, j] is not None and Pmat[i, j] < 0.001) else ''
            ax.text(j, i, f'{M[i, j]:+.2f}{sig}', ha='center', va='center',
                    fontsize=7, color='black')
    ax.set_xticks(range(len(CSV_UNIT_ORDER)))
    ax.set_xticklabels([CSV_UNIT_LABELS[u] for u in CSV_UNIT_ORDER], rotation=45, ha='right')
    ax.set_yticks(range(len(ATS_ORDER))); ax.set_yticklabels(ATS_ORDER)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04)
    cbar.set_label('ΔRMSD (Å) = median(without unit) − median(with unit)\npositive → stabilised by unit')
    ax.set_xlabel('CSV unit (presence/absence in construct)')
    ax.set_ylabel('atS domain (whose RMSD is measured)')
    ax.set_title('Effect of each CSV-unit presence on each atS-domain RMSD\n(* = Mann-Whitney p<0.001)')
    fig.tight_layout()
    fig.savefig(F / 'E2_atS_RMSD_by_unit_presence.png', dpi=200, bbox_inches='tight')
    fig.savefig(F / 'E2_atS_RMSD_by_unit_presence.svg', bbox_inches='tight')
    plt.close(fig)
    print('  wrote E2_atS_RMSD_by_unit_presence.png + .svg')


if __name__ == '__main__':
    print('Computing KH8 region RMSD per combo...')
    kh8_df = compute_kh8_rmsd_for_combos()
    print('Plot E1: KH_b and KH8 vs KH7a presence...')
    plot_KHb_KH8_vs_KH7a(kh8_df)
    print('Plot E2: atS RMSD by CSV-unit presence (heatmap)...')
    plot_atS_by_unit_presence()
    print('Done.')
