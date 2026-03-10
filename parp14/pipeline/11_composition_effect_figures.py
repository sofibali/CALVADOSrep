#!/usr/bin/env python
"""
11_composition_effect_figures.py - Generate publication-quality figures showing
how domain composition influences MD1, MD2, MD3, and ART active sites.

Generates:
1. Heatmap: Domain composition vs active site metrics
2. Boxplots: Active site accessibility by composition category
3. Network diagram: Inter-domain distances
4. Ridge plots: Metric distributions across compositions
5. Seed variability plots: Conformational flexibility

Usage:
    python 11_composition_effect_figures.py --input ../results/per_structure --output ../figures
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set publication style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['figure.dpi'] = 150

# Color palette for domains
DOMAIN_COLORS = {
    'MD1': '#1f77b4',
    'MD2': '#ff7f0e',
    'MD3': '#2ca02c',
    'ART': '#d62728',
    'WWE': '#9467bd',
    'KH7a': '#8c564b',
    'KH1-6': '#e377c2',
    'RRM1': '#7f7f7f',
    'RRM2': '#bcbd22',
    'RRM3': '#17becf'
}


def load_data(input_dir):
    """Load all analysis results."""
    data = {}

    files = {
        'active_sites': 'active_site_metrics.csv',
        'active_sites_summary': 'active_site_summary.csv',
        'interdomain': 'interdomain_distances.csv',
        'interdomain_stats': 'interdomain_statistics.csv',
        'key_domain_summary': 'key_domain_summary.csv',
        'fpocket': 'fpocket_results.csv',
        'confidence': 'confidence_scores.csv'
    }

    for key, filename in files.items():
        filepath = os.path.join(input_dir, filename)
        if os.path.exists(filepath):
            data[key] = pd.read_csv(filepath)
            print(f"  Loaded {key}: {len(data[key])} rows")
        else:
            print(f"  Warning: {filename} not found")

    return data


def fig1_composition_effect_heatmap(data, output_dir):
    """
    Figure 1: Heatmap showing how each metric changes with domain composition
    for each key domain (MD1, MD2, MD3, ART).
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites']

    # Metrics to show
    metrics = {
        'pocket_total_sasa': 'Pocket SASA (Å²)',
        'catalytic_mean_sasa': 'Catalytic SASA (Å²)',
        'pocket_mean_plddt': 'Pocket pLDDT',
        'fpocket_pocket_volume': 'Pocket Volume (Å³)',
        'fpocket_druggability': 'Druggability'
    }

    # Filter to available metrics
    available_metrics = {k: v for k, v in metrics.items() if k in df.columns}

    if len(available_metrics) == 0:
        print("  No metrics available for heatmap")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df[df['active_site'] == site]

        if len(site_df) == 0:
            ax.text(0.5, 0.5, f'No data for {site}', ha='center', va='center')
            ax.set_title(site)
            continue

        # Group by composition category
        grouped = site_df.groupby('composition_category')[list(available_metrics.keys())].mean()

        if len(grouped) == 0:
            continue

        # Normalize for heatmap
        normalized = (grouped - grouped.min()) / (grouped.max() - grouped.min() + 1e-10)

        # Create heatmap
        sns.heatmap(normalized, ax=ax, cmap='RdYlBu_r', center=0.5,
                   xticklabels=[available_metrics[m] for m in normalized.columns],
                   yticklabels=True, annot=grouped.round(1), fmt='.0f',
                   cbar_kws={'label': 'Normalized'})

        ax.set_title(f'{site} Active Site', fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Composition')

        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig1_composition_effect_heatmap.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig1_composition_effect_heatmap.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig1_composition_effect_heatmap")


def fig2_accessibility_boxplots(data, output_dir):
    """
    Figure 2: Boxplots showing active site accessibility (SASA) by composition.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites']

    if 'pocket_total_sasa' not in df.columns:
        print("  Skipping: pocket_total_sasa not available")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df[df['active_site'] == site].copy()

        if len(site_df) == 0:
            ax.text(0.5, 0.5, f'No data for {site}', ha='center', va='center')
            continue

        # Sort by median SASA
        order = site_df.groupby('composition_category')['pocket_total_sasa'].median().sort_values().index

        # Create boxplot
        sns.boxplot(data=site_df, x='composition_category', y='pocket_total_sasa',
                   ax=ax, order=order, palette='Set2')

        # Add individual points
        sns.stripplot(data=site_df, x='composition_category', y='pocket_total_sasa',
                     ax=ax, order=order, color='black', alpha=0.3, size=3)

        ax.set_title(f'{site} Active Site Accessibility', fontweight='bold')
        ax.set_xlabel('Domain Composition')
        ax.set_ylabel('Pocket SASA (Å²)')

        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # Add sample sizes
        for i, comp in enumerate(order):
            n = len(site_df[site_df['composition_category'] == comp])
            ax.text(i, ax.get_ylim()[1], f'n={n}', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig2_accessibility_boxplots.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig2_accessibility_boxplots.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig2_accessibility_boxplots")


def fig3_interdomain_distance_heatmap(data, output_dir):
    """
    Figure 3: Heatmap showing mean inter-domain distances.
    Creates two versions: one with MD1, one with MD1L1.
    Domains are ordered by sequence position (N-terminus to C-terminus).
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain data not available")
        return

    df = data['interdomain']

    # Domain order by sequence position (N-terminus to C-terminus)
    SEQUENCE_ORDER_MD1 = ['RRM1', 'RRM2', 'RRM3', 'KH1-6', 'KH7a', 'MD1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']
    SEQUENCE_ORDER_MD1L1 = ['RRM1', 'RRM2', 'RRM3', 'KH1-6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']
    # Versions excluding RRM domains
    SEQUENCE_ORDER_MD1_noRRM = ['KH1-6', 'KH7a', 'MD1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']
    SEQUENCE_ORDER_MD1L1_noRRM = ['KH1-6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

    # Get mean distances across all structures
    mean_distances = df.groupby(['domain1', 'domain2'])['d_com'].mean().reset_index()

    # Get all domains in data
    all_domains = set(mean_distances['domain1'].unique()) | set(mean_distances['domain2'].unique())

    def create_heatmap(domain_order, exclude_domain, suffix, title_suffix):
        # Filter to domains in this order that exist in data
        domains = [d for d in domain_order if d in all_domains and d != exclude_domain]
        n_domains = len(domains)

        if n_domains < 2:
            return

        dist_matrix = np.zeros((n_domains, n_domains))
        dist_matrix[:] = np.nan

        domain_to_idx = {d: i for i, d in enumerate(domains)}

        for _, row in mean_distances.iterrows():
            d1, d2 = row['domain1'], row['domain2']
            if d1 in domain_to_idx and d2 in domain_to_idx:
                i = domain_to_idx[d1]
                j = domain_to_idx[d2]
                dist_matrix[i, j] = row['d_com']
                dist_matrix[j, i] = row['d_com']

        # Fill diagonal with 0
        np.fill_diagonal(dist_matrix, 0)

        fig, ax = plt.subplots(figsize=(10, 8))

        # Create mask for upper triangle (to show only lower)
        mask = np.triu(np.ones_like(dist_matrix, dtype=bool), k=1)

        sns.heatmap(dist_matrix, ax=ax, mask=mask, cmap='viridis_r',
                   xticklabels=domains, yticklabels=domains,
                   annot=True, fmt='.1f', cbar_kws={'label': 'Distance (Å)'})

        ax.set_title(f'Mean Inter-Domain Center-of-Mass Distances{title_suffix}', fontweight='bold')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'fig3_interdomain_distances_{suffix}.svg'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'fig3_interdomain_distances_{suffix}.png'), bbox_inches='tight', dpi=300)
        plt.close()
        print(f"  Saved: fig3_interdomain_distances_{suffix}")

    # Create MD1 version (excluding MD1L1)
    create_heatmap(SEQUENCE_ORDER_MD1, 'MD1L1', 'MD1', ' (MD1)')

    # Create MD1L1 version (excluding MD1)
    create_heatmap(SEQUENCE_ORDER_MD1L1, 'MD1', 'MD1L1', ' (MD1L1)')

    # Create versions excluding RRM domains
    create_heatmap(SEQUENCE_ORDER_MD1_noRRM, 'MD1L1', 'MD1_noRRM', ' (MD1, no RRM)')
    create_heatmap(SEQUENCE_ORDER_MD1L1_noRRM, 'MD1', 'MD1L1_noRRM', ' (MD1L1, no RRM)')


def fig3b_interdomain_distance_std_heatmap(data, output_dir):
    """
    Figure 3b: Heatmap showing standard deviation of inter-domain distances.
    Creates two versions: one with MD1, one with MD1L1.
    Domains are ordered by sequence position (N-terminus to C-terminus).
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain std data not available")
        return

    df = data['interdomain']

    # Domain order by sequence position (N-terminus to C-terminus)
    SEQUENCE_ORDER_MD1 = ['RRM1', 'RRM2', 'RRM3', 'KH1-6', 'KH7a', 'MD1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']
    SEQUENCE_ORDER_MD1L1 = ['RRM1', 'RRM2', 'RRM3', 'KH1-6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']
    # Versions excluding RRM domains
    SEQUENCE_ORDER_MD1_noRRM = ['KH1-6', 'KH7a', 'MD1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']
    SEQUENCE_ORDER_MD1L1_noRRM = ['KH1-6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

    # Get std of distances across all structures
    std_distances = df.groupby(['domain1', 'domain2'])['d_com'].std().reset_index()

    # Get all domains in data
    all_domains = set(std_distances['domain1'].unique()) | set(std_distances['domain2'].unique())

    def create_std_heatmap(domain_order, exclude_domain, suffix, title_suffix):
        # Filter to domains in this order that exist in data
        domains = [d for d in domain_order if d in all_domains and d != exclude_domain]
        n_domains = len(domains)

        if n_domains < 2:
            return

        std_matrix = np.zeros((n_domains, n_domains))
        std_matrix[:] = np.nan

        domain_to_idx = {d: i for i, d in enumerate(domains)}

        for _, row in std_distances.iterrows():
            d1, d2 = row['domain1'], row['domain2']
            if d1 in domain_to_idx and d2 in domain_to_idx:
                i = domain_to_idx[d1]
                j = domain_to_idx[d2]
                std_matrix[i, j] = row['d_com']
                std_matrix[j, i] = row['d_com']

        # Fill diagonal with 0
        np.fill_diagonal(std_matrix, 0)

        fig, ax = plt.subplots(figsize=(10, 8))

        # Create mask for upper triangle (to show only lower)
        mask = np.triu(np.ones_like(std_matrix, dtype=bool), k=1)

        # Use a different colormap for variation (e.g., Reds or YlOrRd)
        sns.heatmap(std_matrix, ax=ax, mask=mask, cmap='YlOrRd',
                   xticklabels=domains, yticklabels=domains,
                   annot=True, fmt='.1f', cbar_kws={'label': 'Std Dev (Å)'})

        ax.set_title(f'Variation in Inter-Domain Center-of-Mass Distances{title_suffix}', fontweight='bold')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'fig3b_interdomain_distances_std_{suffix}.svg'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'fig3b_interdomain_distances_std_{suffix}.png'), bbox_inches='tight', dpi=300)
        plt.close()
        print(f"  Saved: fig3b_interdomain_distances_std_{suffix}")

    # Create MD1 version (excluding MD1L1)
    create_std_heatmap(SEQUENCE_ORDER_MD1, 'MD1L1', 'MD1', ' (MD1)')

    # Create MD1L1 version (excluding MD1)
    create_std_heatmap(SEQUENCE_ORDER_MD1L1, 'MD1', 'MD1L1', ' (MD1L1)')

    # Create versions excluding RRM domains
    create_std_heatmap(SEQUENCE_ORDER_MD1_noRRM, 'MD1L1', 'MD1_noRRM', ' (MD1, no RRM)')
    create_std_heatmap(SEQUENCE_ORDER_MD1L1_noRRM, 'MD1', 'MD1L1_noRRM', ' (MD1L1, no RRM)')


def fig4_pocket_volume_by_composition(data, output_dir):
    """
    Figure 4: Ridge plot showing pocket volume distribution across compositions.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites']

    if 'fpocket_pocket_volume' not in df.columns:
        print("  Skipping: fpocket_pocket_volume not available")
        return

    # Filter out NaN volumes
    df_vol = df[df['fpocket_pocket_volume'].notna()].copy()

    if len(df_vol) == 0:
        print("  Skipping: no pocket volume data")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df_vol[df_vol['active_site'] == site]

        if len(site_df) < 10:
            ax.text(0.5, 0.5, f'Insufficient data for {site}', ha='center', va='center')
            continue

        # Get composition categories
        categories = site_df['composition_category'].unique()

        # Create violin plot
        parts = ax.violinplot([site_df[site_df['composition_category'] == cat]['fpocket_pocket_volume'].values
                               for cat in categories], positions=range(len(categories)),
                              showmeans=True, showmedians=True)

        # Color the violins
        for pc in parts['bodies']:
            pc.set_facecolor(DOMAIN_COLORS.get(site, '#1f77b4'))
            pc.set_alpha(0.7)

        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(categories, rotation=45, ha='right')
        ax.set_title(f'{site} Pocket Volume by Composition', fontweight='bold')
        ax.set_ylabel('Pocket Volume (Å³)')
        ax.set_xlabel('Composition')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig4_pocket_volume_distribution.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig4_pocket_volume_distribution.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig4_pocket_volume_distribution")


def fig5_seed_variability(data, output_dir):
    """
    Figure 5: Seed-to-seed variability in metrics (conformational flexibility).
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites']

    # Calculate variance across seeds for each structure
    variance_df = df.groupby(['structure', 'active_site']).agg({
        'pocket_total_sasa': ['mean', 'std'],
        'pocket_mean_plddt': ['mean', 'std'],
    }).reset_index()

    variance_df.columns = ['structure', 'active_site', 'sasa_mean', 'sasa_std',
                          'plddt_mean', 'plddt_std']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # SASA variability
    ax = axes[0]
    for site in ['MD1', 'MD2', 'MD3', 'ART']:
        site_data = variance_df[variance_df['active_site'] == site]
        if len(site_data) > 0:
            ax.scatter(site_data['sasa_mean'], site_data['sasa_std'],
                      label=site, alpha=0.6, c=DOMAIN_COLORS.get(site, '#1f77b4'), s=30)

    ax.set_xlabel('Mean Pocket SASA (Å²)')
    ax.set_ylabel('SASA Std Dev (across seeds)')
    ax.set_title('SASA Variability Across Seeds', fontweight='bold')
    ax.legend()

    # pLDDT variability
    ax = axes[1]
    for site in ['MD1', 'MD2', 'MD3', 'ART']:
        site_data = variance_df[variance_df['active_site'] == site]
        if len(site_data) > 0:
            ax.scatter(site_data['plddt_mean'], site_data['plddt_std'],
                      label=site, alpha=0.6, c=DOMAIN_COLORS.get(site, '#1f77b4'), s=30)

    ax.set_xlabel('Mean Pocket pLDDT')
    ax.set_ylabel('pLDDT Std Dev (across seeds)')
    ax.set_title('Confidence Variability Across Seeds', fontweight='bold')
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig5_seed_variability.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig5_seed_variability.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig5_seed_variability")


def fig6_composition_summary(data, output_dir):
    """
    Figure 6: Summary showing key domain metrics by number of neighboring domains.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    metrics = [
        ('pocket_total_sasa', 'Pocket SASA (Å²)'),
        ('catalytic_mean_sasa', 'Catalytic SASA (Å²)'),
        ('pocket_mean_plddt', 'Pocket pLDDT'),
        ('fpocket_pocket_volume', 'Pocket Volume (Å³)')
    ]

    for ax, (metric, label) in zip(axes.flatten(), metrics):
        if metric not in df.columns:
            ax.text(0.5, 0.5, f'{label} not available', ha='center', va='center')
            continue

        # Group by active site and number of neighbors
        for site in ['MD1', 'MD2', 'MD3', 'ART']:
            site_df = df[df['active_site'] == site]
            if len(site_df) == 0:
                continue

            grouped = site_df.groupby('n_neighbors')[metric].agg(['mean', 'std'])

            ax.errorbar(grouped.index, grouped['mean'], yerr=grouped['std'],
                       label=site, marker='o', capsize=3, alpha=0.8,
                       color=DOMAIN_COLORS.get(site, '#1f77b4'))

        ax.set_xlabel('Number of Neighboring Domains')
        ax.set_ylabel(label)
        ax.set_title(f'{label} vs Domain Context', fontweight='bold')
        ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig6_composition_summary.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig6_composition_summary.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig6_composition_summary")


def fig7_contact_density_network(data, output_dir):
    """
    Figure 7: Network diagram showing domain-domain contact densities.
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain data not available")
        return

    df = data['interdomain']

    # Calculate mean contact density
    contacts = df.groupby(['domain1', 'domain2']).agg({
        'contact_density': 'mean',
        'd_com': 'mean',
        'n_contacts': 'mean'
    }).reset_index()

    fig, ax = plt.subplots(figsize=(12, 10))

    # Domain positions (approximate linear layout)
    domain_order = ['RRM1', 'RRM2', 'RRM3', 'KH1-6', 'KH7a', 'MD1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

    # Filter to domains present in data
    present_domains = set(contacts['domain1'].unique()) | set(contacts['domain2'].unique())
    domain_order = [d for d in domain_order if d in present_domains]

    # Position domains in a circle
    n = len(domain_order)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    positions = {d: (np.cos(a), np.sin(a)) for d, a in zip(domain_order, angles)}

    # Draw edges (contacts)
    max_density = contacts['contact_density'].max()

    for _, row in contacts.iterrows():
        d1, d2 = row['domain1'], row['domain2']
        if d1 not in positions or d2 not in positions:
            continue

        x1, y1 = positions[d1]
        x2, y2 = positions[d2]

        # Line width proportional to contact density
        width = 1 + 10 * (row['contact_density'] / max_density)

        # Color based on distance (blue = close, red = far)
        color_val = min(row['d_com'] / 100, 1)  # Normalize
        color = plt.cm.RdYlBu_r(color_val)

        ax.plot([x1, x2], [y1, y2], color=color, linewidth=width, alpha=0.6, zorder=1)

    # Draw nodes (domains)
    for domain, (x, y) in positions.items():
        color = DOMAIN_COLORS.get(domain, '#1f77b4')
        circle = plt.Circle((x, y), 0.12, color=color, zorder=2)
        ax.add_patch(circle)
        ax.text(x, y, domain, ha='center', va='center', fontsize=8,
               fontweight='bold', color='white', zorder=3)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Inter-Domain Contact Network\n(Line width = contact density, Color = distance)',
                fontweight='bold')

    # Add colorbar for distance
    sm = plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(0, 100))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, label='Distance (Å)')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig7_contact_network.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig7_contact_network.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig7_contact_network")


def main():
    parser = argparse.ArgumentParser(description='Generate publication figures')
    parser.add_argument('--input', '-i', default='../results/per_structure',
                       help='Input directory with analysis results')
    parser.add_argument('--output', '-o', default='../figures/domain_composition_effects',
                       help='Output directory for figures')
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    data = load_data(args.input)

    if not data:
        print("No data loaded. Run analysis scripts first.")
        return

    print("\nGenerating figures...")

    fig1_composition_effect_heatmap(data, output_dir)
    fig2_accessibility_boxplots(data, output_dir)
    fig3_interdomain_distance_heatmap(data, output_dir)
    fig3b_interdomain_distance_std_heatmap(data, output_dir)
    fig4_pocket_volume_by_composition(data, output_dir)
    fig5_seed_variability(data, output_dir)
    fig6_composition_summary(data, output_dir)
    fig7_contact_density_network(data, output_dir)

    print(f"\nAll figures saved to {output_dir}")


if __name__ == '__main__':
    main()
