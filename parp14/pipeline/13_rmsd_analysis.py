#!/usr/bin/env python3
"""
13_rmsd_analysis.py - Pairwise RMSD analysis across all domain compositions

For each composition (e.g. kh1-kh6_art), loads all 25 AF3 predicted structures
(5 seeds x 5 samples), computes pairwise CA-RMSD with proper Kabsch alignment,
and produces:
  1. A CSV table of per-composition RMSD statistics
  2. A heatmap of mean pairwise RMSD across compositions

Usage:
    python 13_rmsd_analysis.py --input ../alphafold_outputs --output ../results/rmsd
    python 13_rmsd_analysis.py --input ../alphafold_outputs --output ../results/rmsd --jobs 8
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio.PDB import MMCIFParser
from multiprocessing import Pool, cpu_count
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# Style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['figure.dpi'] = 150


# ── Kabsch-aligned RMSD ─────────────────────────────────────────────────────

def kabsch_rmsd(P, Q):
    """
    Compute RMSD between coordinate sets P and Q after optimal
    superposition (Kabsch algorithm).  Both arrays must be (N, 3).
    """
    assert P.shape == Q.shape
    n = P.shape[0]

    # Centre both
    p0 = P - P.mean(axis=0)
    q0 = Q - Q.mean(axis=0)

    # Covariance matrix
    H = p0.T @ q0

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.eye(3)
    sign_matrix[2, 2] = np.sign(d)

    # Optimal rotation
    R = Vt.T @ sign_matrix @ U.T

    # Rotate P onto Q
    p_rotated = p0 @ R.T

    # RMSD
    diff = p_rotated - q0
    rmsd = np.sqrt((diff ** 2).sum() / n)
    return rmsd


# ── Structure loading ────────────────────────────────────────────────────────

def get_ca_coords(cif_path):
    """Extract CA coordinates from an mmCIF file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('s', str(cif_path))
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
    return np.array(coords) if coords else None


# ── Per-composition worker ───────────────────────────────────────────────────

def process_composition(args):
    """Process one composition directory: load structures, compute pairwise RMSD."""
    comp_dir, comp_name = args

    cif_files = sorted(Path(comp_dir).glob('seed-*/model.cif'))
    if len(cif_files) < 2:
        return None

    # Load all CA coordinate sets
    coord_sets = []
    labels = []
    for cif in cif_files:
        coords = get_ca_coords(cif)
        if coords is not None and len(coords) > 0:
            coord_sets.append(coords)
            labels.append(cif.parent.name)

    if len(coord_sets) < 2:
        return None

    # Verify all have same number of atoms
    n_atoms = len(coord_sets[0])
    valid = [(c, l) for c, l in zip(coord_sets, labels) if len(c) == n_atoms]
    if len(valid) < 2:
        return None
    coord_sets, labels = zip(*valid)
    coord_sets = list(coord_sets)
    labels = list(labels)

    # Compute pairwise RMSD matrix
    n = len(coord_sets)
    rmsd_matrix = np.zeros((n, n))
    pairwise_values = []

    for i, j in combinations(range(n), 2):
        rmsd = kabsch_rmsd(coord_sets[i], coord_sets[j])
        rmsd_matrix[i, j] = rmsd
        rmsd_matrix[j, i] = rmsd
        pairwise_values.append(rmsd)

    pairwise_values = np.array(pairwise_values)

    # Split into within-seed and between-seed
    within_seed = []
    between_seed = []
    for i, j in combinations(range(n), 2):
        seed_i = labels[i].split('_')[0]  # e.g. "seed-1"
        seed_j = labels[j].split('_')[0]
        rmsd = rmsd_matrix[i, j]
        if seed_i == seed_j:
            within_seed.append(rmsd)
        else:
            between_seed.append(rmsd)

    result = {
        'composition': comp_name,
        'n_structures': n,
        'n_atoms': n_atoms,
        'mean_rmsd': pairwise_values.mean(),
        'median_rmsd': np.median(pairwise_values),
        'std_rmsd': pairwise_values.std(),
        'min_rmsd': pairwise_values.min(),
        'max_rmsd': pairwise_values.max(),
        'within_seed_mean': np.mean(within_seed) if within_seed else np.nan,
        'between_seed_mean': np.mean(between_seed) if between_seed else np.nan,
        'rmsd_matrix': rmsd_matrix,
        'labels': labels,
    }
    return result


# ── Figures ──────────────────────────────────────────────────────────────────

def make_heatmap(df, output_dir):
    """Heatmap of mean RMSD organised by domain count and composition."""

    # Count domains for sorting
    def count_domains(name):
        parts = name.replace('-', '_').lower().split('_')
        domains = set()
        for p in parts:
            if p in ('rrm1', 'rrm2', 'rrm3', 'kh7a', 'md1', 'md1l1', 'md2', 'md3', 'wwe', 'art'):
                domains.add(p)
            elif p.startswith('kh') and p not in ('khb',):
                domains.add(p)
        return len(domains)

    df = df.copy()
    df['n_domains'] = df['composition'].apply(count_domains)
    df = df.sort_values(['n_domains', 'mean_rmsd'], ascending=[True, True]).reset_index(drop=True)

    # Truncate long labels
    max_label = 40
    df['label'] = df['composition'].apply(lambda s: s[:max_label] + '…' if len(s) > max_label else s)

    # ── Main heatmap: one row per composition (no cell text for speed) ───
    n = len(df)
    fig_height = max(8, n * 0.18)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    cols = ['mean_rmsd', 'median_rmsd', 'min_rmsd', 'max_rmsd',
            'within_seed_mean_rmsd', 'between_seed_mean_rmsd']
    col_labels = ['Mean', 'Median', 'Min', 'Max', 'Within\nSeed', 'Between\nSeed']
    matrix = df[cols].values

    im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')
    ax.set_yticks(range(n))
    ax.set_yticklabels(df['label'], fontsize=max(2, min(6, 300 // n)))
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.set_title('Pairwise CA-RMSD (Å) Across Domain Compositions', fontweight='bold', fontsize=14)

    plt.colorbar(im, ax=ax, shrink=0.4, label='RMSD (Å)')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsd_heatmap.png'), bbox_inches='tight', dpi=200)
    plt.close()
    print("  Saved: rmsd_heatmap")

    # ── Bar chart: only top 40 + bottom 20 to keep readable ──────────────
    df_top = df.sort_values('mean_rmsd', ascending=False).head(40)
    df_bot = df.sort_values('mean_rmsd', ascending=True).head(20)
    df_show = pd.concat([df_top, df_bot]).drop_duplicates().sort_values('mean_rmsd', ascending=False)
    n_show = len(df_show)

    fig, ax = plt.subplots(figsize=(12, max(8, n_show * 0.28)))
    cmap_colors = plt.cm.viridis(np.linspace(0, 1, df['n_domains'].max() + 1))

    ax.barh(range(n_show), df_show['mean_rmsd'].values, xerr=df_show['std_rmsd'].values,
            color=[cmap_colors[d] for d in df_show['n_domains']],
            capsize=2, edgecolor='grey', linewidth=0.3)

    ax.set_yticks(range(n_show))
    ax.set_yticklabels(df_show['label'].values, fontsize=7)
    ax.set_xlabel('Mean Pairwise CA-RMSD (Å)')
    ax.set_title('Top 40 Most Variable & Top 20 Most Consistent Compositions',
                 fontweight='bold')
    ax.invert_yaxis()

    from matplotlib.patches import Patch
    unique_counts = sorted(df_show['n_domains'].unique())
    legend_elements = [Patch(facecolor=cmap_colors[c], label=f'{c} domains') for c in unique_counts]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsd_barchart.png'), bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(output_dir, 'rmsd_barchart.svg'), bbox_inches='tight')
    plt.close()
    print("  Saved: rmsd_barchart")

    # ── Box plot: RMSD by domain count ───────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    domain_counts = sorted(df['n_domains'].unique())
    box_data = [df[df['n_domains'] == c]['mean_rmsd'].values for c in domain_counts]
    box_colors = plt.cm.viridis(np.linspace(0, 1, max(domain_counts) + 1))

    bp = ax.boxplot(box_data, positions=domain_counts, widths=0.6, patch_artist=True)
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(box_colors[domain_counts[i]])
        patch.set_alpha(0.7)

    for c in domain_counts:
        vals = df[df['n_domains'] == c]['mean_rmsd'].values
        jitter = np.random.normal(0, 0.08, len(vals))
        ax.scatter(c + jitter, vals, alpha=0.4, s=15, color=box_colors[c], zorder=3, edgecolors='none')

    ax.set_xlabel('Number of Domains in Structure')
    ax.set_ylabel('Mean Pairwise CA-RMSD (Å)')
    ax.set_title('Structural Variability vs Domain Count', fontweight='bold')
    ax.set_xticks(domain_counts)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsd_vs_domain_count.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'rmsd_vs_domain_count.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: rmsd_vs_domain_count")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='RMSD analysis across all AF3 compositions')
    parser.add_argument('--input', '-i', default='../alphafold_outputs',
                        help='Directory containing AF3 output folders')
    parser.add_argument('--output', '-o', default='../results/rmsd',
                        help='Output directory for results')
    parser.add_argument('--jobs', '-j', type=int, default=4,
                        help='Number of parallel workers')
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Discover compositions
    comp_dirs = sorted([
        d for d in input_dir.iterdir()
        if d.is_dir() and any(d.glob('seed-*/model.cif'))
    ])
    print(f"Found {len(comp_dirs)} compositions with structures")

    # Process in parallel
    tasks = [(str(d), d.name) for d in comp_dirs]
    print(f"Computing pairwise RMSD with {args.jobs} workers...")

    with Pool(args.jobs) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(process_composition, tasks)):
            if result is not None:
                results.append(result)
            if (i + 1) % 50 == 0 or (i + 1) == len(tasks):
                print(f"  Processed {i+1}/{len(tasks)} compositions ({len(results)} with data)")

    if not results:
        print("No results. Check that structures exist.")
        return

    # Build summary table
    rows = []
    for r in results:
        rows.append({
            'composition': r['composition'],
            'n_structures': r['n_structures'],
            'n_atoms': r['n_atoms'],
            'mean_rmsd': round(r['mean_rmsd'], 3),
            'median_rmsd': round(r['median_rmsd'], 3),
            'std_rmsd': round(r['std_rmsd'], 3),
            'min_rmsd': round(r['min_rmsd'], 3),
            'max_rmsd': round(r['max_rmsd'], 3),
            'within_seed_mean_rmsd': round(r['within_seed_mean'], 3) if not np.isnan(r['within_seed_mean']) else np.nan,
            'between_seed_mean_rmsd': round(r['between_seed_mean'], 3) if not np.isnan(r['between_seed_mean']) else np.nan,
        })

    df = pd.DataFrame(rows).sort_values('mean_rmsd', ascending=False).reset_index(drop=True)

    # Save table
    csv_path = os.path.join(output_dir, 'rmsd_summary.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nSaved RMSD table: {csv_path}  ({len(df)} compositions)")

    # Print top/bottom
    print(f"\n{'='*80}")
    print("Top 10 most variable compositions (highest mean RMSD):")
    print(df.head(10)[['composition', 'n_structures', 'mean_rmsd', 'std_rmsd']].to_string(index=False))

    print(f"\nTop 10 most consistent compositions (lowest mean RMSD):")
    print(df.tail(10)[['composition', 'n_structures', 'mean_rmsd', 'std_rmsd']].to_string(index=False))
    print(f"{'='*80}")

    # Generate figures
    print("\nGenerating figures...")
    make_heatmap(df, output_dir)

    print(f"\nAll results saved to {output_dir}")


if __name__ == '__main__':
    main()
