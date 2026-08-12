#!/usr/bin/env python
"""
D — Per-combo intra-ensemble RMSF.

For each combo, for each atS domain present, compute the per-residue RMSF of Cα across
the (up to 25) AF3 models AFTER per-domain Kabsch alignment to the combo's mean structure.

This separates "AF3 sampling noise within a combo" from "construct-dependent rearrangement
relative to FL_AF3". Domains with low intra-combo RMSF but high RMSD-to-FLref are genuinely
adopting a different fold/orientation in that construct.

Inputs:  parp14/analysis/structure_analysis/per_combo/<combo>__residues.parquet
         parp14/analysis/structure_analysis/aggregated/all_summary.parquet

Outputs: parp14/analysis/structure_analysis/aggregated/all_intra_rmsf.parquet
            columns: combo, atS_domain, n_models, mean_rmsf_A, max_rmsf_A
         figures/D1_intra_vs_FLrmsd.{png,svg}
            scatter per atS domain: x = mean intra-combo RMSF, y = RMSD to FL ref
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/home/sbali/CALVADOS/parp14')
PC = ROOT / 'analysis' / 'structure_analysis' / 'per_combo'
A = ROOT / 'analysis' / 'structure_analysis' / 'aggregated'
F = ROOT / 'analysis' / 'structure_analysis' / 'figures'
F.mkdir(exist_ok=True)

ATS_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH_N', 'KH_Ca',
             'MD1', 'MD2', 'MD3', 'KH_b', 'WWE', 'ART']


def kabsch(P, Q):
    """Optimal rotation aligning P onto Q (centered already)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    return Vt.T @ D @ U.T


def rmsf_for_combo(parquet_path):
    """Returns list of dicts: combo, atS_domain, n_models, mean_rmsf_A, max_rmsf_A."""
    combo = parquet_path.name.replace('__residues.parquet', '')
    df = pd.read_parquet(parquet_path)
    rows = []
    for d, g in df.groupby('atS_domain'):
        if d == 'none':
            continue
        # Pivot: rows = models, cols = residue (use res_construct as id), vals = ca xyz
        pivot = g.pivot_table(index='model', columns='res_construct',
                              values=['ca_x', 'ca_y', 'ca_z'])
        # pivot has multi-index columns ('ca_x', resi); reshape to (n_models, n_res, 3)
        models = pivot.index.tolist()
        residues = pivot['ca_x'].columns.tolist()
        nM, nR = len(models), len(residues)
        if nM < 2 or nR < 3:
            continue
        coords = np.stack([pivot[c].values for c in ('ca_x', 'ca_y', 'ca_z')], axis=-1)
        # coords shape: (nM, nR, 3)
        # Kabsch align all to first model
        ref = coords[0] - coords[0].mean(0)
        aligned = np.empty_like(coords)
        aligned[0] = ref
        for m in range(1, nM):
            P = coords[m] - coords[m].mean(0)
            R = kabsch(P, ref)
            aligned[m] = P @ R.T
        # Iteratively re-align to mean (1 pass usually enough for tight ensembles)
        for _ in range(2):
            mean = aligned.mean(0)
            mean = mean - mean.mean(0)
            for m in range(nM):
                P = aligned[m]
                R = kabsch(P, mean)
                aligned[m] = P @ R.T
        mean = aligned.mean(0)
        # RMSF per residue = sqrt(mean over models of ||x_m - mean||^2)
        diff2 = np.sum((aligned - mean) ** 2, axis=2)  # (nM, nR)
        rmsf = np.sqrt(diff2.mean(0))
        rows.append({
            'combo': combo, 'atS_domain': d, 'n_models': int(nM),
            'mean_rmsf_A': float(np.mean(rmsf)),
            'max_rmsf_A': float(np.max(rmsf)),
        })
    return rows


def run_all(workers=16):
    files = sorted(PC.glob('*__residues.parquet'))
    print(f'Computing intra-combo RMSF for {len(files)} combos using {workers} workers...')
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(rmsf_for_combo, p): p for p in files}
        n = 0
        for fut in as_completed(futs):
            try:
                rows.extend(fut.result())
            except Exception as e:
                print(f'  err {futs[fut].name}: {e}')
            n += 1
            if n % 200 == 0:
                print(f'  [{n}/{len(files)}]')
    df = pd.DataFrame(rows)
    out = A / 'all_intra_rmsf.parquet'
    df.to_parquet(out, index=False)
    print(f'Saved {out} ({len(df)} rows)')
    return df


def plot_intra_vs_FL(df_rmsf):
    summary = pd.read_parquet(A / 'all_summary.parquet')
    af3 = summary[summary['source'] == 'af3']
    # mean RMSD-to-FLref per (combo, domain) across the up-to-25 models
    rmsd_mean = af3.groupby(['combo', 'atS_domain'])['rmsd_to_FLref'].mean().reset_index()
    merged = df_rmsf.merge(rmsd_mean, on=['combo', 'atS_domain'])
    fig, axes = plt.subplots(3, 4, figsize=(14, 10), sharex=True, sharey=True)
    for ax, d in zip(axes.flat, ATS_ORDER + [None]):
        if d is None:
            ax.axis('off'); continue
        sub = merged[merged.atS_domain == d]
        if sub.empty:
            ax.set_title(d); ax.axis('off'); continue
        ax.scatter(sub['mean_rmsf_A'], sub['rmsd_to_FLref'],
                   s=8, alpha=0.4, c='#1f77b4', edgecolor='none')
        ax.plot([0, 10], [0, 10], '--', color='gray', alpha=0.5, lw=0.7)
        ax.set_title(f'{d}  (n={len(sub)})')
        ax.grid(alpha=0.3)
    fig.supxlabel('Intra-combo RMSF (Å)  — AF3 sampling spread within one combo')
    fig.supylabel('RMSD to FL_AF3 (Å) — construct vs FL ref')
    fig.suptitle('Intra-combo flexibility vs construct-dependent shift')
    fig.tight_layout()
    fig.savefig(F / 'D1_intra_vs_FLrmsd.png', dpi=200, bbox_inches='tight')
    fig.savefig(F / 'D1_intra_vs_FLrmsd.svg', bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {F}/D1_intra_vs_FLrmsd.png + .svg')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--plot-only', action='store_true')
    args = ap.parse_args()
    if args.plot_only:
        df = pd.read_parquet(A / 'all_intra_rmsf.parquet')
    else:
        df = run_all(workers=args.workers)
    plot_intra_vs_FL(df)
