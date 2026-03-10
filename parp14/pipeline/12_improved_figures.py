#!/usr/bin/env python
"""
12_improved_figures.py - Improved visualizations for PARP14 analysis

Creates:
1. Violin plots with dots colored by number of domains
2. Scatter plots showing relationships between domain count and metrics
3. Clear composition labeling
4. Relationship between domain count and interdomain contacts

Usage:
    python 12_improved_figures.py --input ../results/per_structure --output ../figures
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['figure.dpi'] = 150

# Color palette for number of domains
DOMAIN_COUNT_CMAP = plt.cm.viridis

# Domain colors
SITE_COLORS = {
    'MD1': '#1f77b4',
    'MD2': '#ff7f0e',
    'MD3': '#2ca02c',
    'ART': '#d62728'
}


def count_domains_in_structure(structure_name):
    """Count number of domain groups in a structure name."""
    name_lower = structure_name.lower()

    # Domain patterns to look for
    domain_patterns = [
        'rrm1', 'rrm2', 'rrm3',
        'kh1-kh6', 'kh1_kh2', 'kh7a', 'khb-kh8', 'khb_kh8',
        'md1l1', 'md1', 'md2', 'md3',
        'wwe', 'art'
    ]

    count = 0
    counted = set()

    for pattern in domain_patterns:
        if pattern in name_lower and pattern not in counted:
            # Avoid double counting md1 and md1l1
            if pattern == 'md1' and 'md1l1' in name_lower:
                continue
            # Avoid double counting kh variations
            if pattern == 'kh1_kh2' and 'kh1-kh6' in counted:
                continue
            if pattern == 'khb_kh8' and 'khb-kh8' in counted:
                continue
            counted.add(pattern)
            count += 1

    return count


def load_data(input_dir):
    """Load analysis results."""
    data = {}

    files = {
        'active_sites': 'active_site_metrics.csv',
        'interdomain': 'interdomain_distances.csv',
    }

    for key, filename in files.items():
        filepath = os.path.join(input_dir, filename)
        if os.path.exists(filepath):
            data[key] = pd.read_csv(filepath)
            print(f"  Loaded {key}: {len(data[key])} rows")

    return data


def fig1_violin_accessibility_by_domains(data, output_dir):
    """
    Violin plots of active site accessibility with dots colored by domain count.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites'].copy()

    # Add domain count
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df[df['active_site'] == site].copy()

        if len(site_df) == 0 or 'pocket_total_sasa' not in site_df.columns:
            ax.text(0.5, 0.5, f'No data for {site}', ha='center', va='center')
            continue

        # Get unique domain counts and sort, filter those with data
        domain_counts = sorted(site_df['n_domains'].unique())

        # Filter to domain counts with actual data
        valid_counts = []
        valid_data = []
        for n in domain_counts:
            data_n = site_df[site_df['n_domains'] == n]['pocket_total_sasa'].dropna().values
            if len(data_n) > 0:
                valid_counts.append(n)
                valid_data.append(data_n)

        if len(valid_data) == 0:
            ax.text(0.5, 0.5, f'No valid data for {site}', ha='center', va='center')
            continue

        domain_counts = valid_counts

        # Create violin plot
        parts = ax.violinplot(
            valid_data,
            positions=domain_counts,
            showmeans=True,
            showmedians=True,
            widths=0.8
        )

        # Color violins
        for pc in parts['bodies']:
            pc.set_facecolor(SITE_COLORS[site])
            pc.set_alpha(0.3)

        # Add scatter points colored by domain count
        for n in domain_counts:
            subset = site_df[site_df['n_domains'] == n]['pocket_total_sasa'].dropna()
            jitter = np.random.normal(0, 0.1, len(subset))
            color = DOMAIN_COUNT_CMAP(n / max(domain_counts))
            ax.scatter(n + jitter, subset, c=[color], alpha=0.5, s=15,
                      edgecolors='none')

        ax.set_xlabel('Number of Domains in Structure')
        ax.set_ylabel('Pocket SASA (Å²)')
        ax.set_title(f'{site} Active Site Accessibility vs Domain Count', fontweight='bold')
        ax.set_xticks(domain_counts)

        # Add sample sizes
        for n in domain_counts:
            count = len(site_df[site_df['n_domains'] == n])
            ax.text(n, ax.get_ylim()[1] * 0.98, f'n={count}',
                   ha='center', va='top', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig1_violin_accessibility_by_domains.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig1_violin_accessibility_by_domains.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig1_violin_accessibility_by_domains")


def fig2_scatter_sasa_vs_plddt(data, output_dir):
    """
    Scatter plot of SASA vs pLDDT colored by domain count.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites'].copy()
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df[df['active_site'] == site].copy()

        if len(site_df) == 0:
            continue

        if 'pocket_total_sasa' not in site_df.columns or 'pocket_mean_plddt' not in site_df.columns:
            continue

        # Remove NaN
        site_df = site_df.dropna(subset=['pocket_total_sasa', 'pocket_mean_plddt'])

        scatter = ax.scatter(
            site_df['pocket_mean_plddt'],
            site_df['pocket_total_sasa'],
            c=site_df['n_domains'],
            cmap=DOMAIN_COUNT_CMAP,
            alpha=0.6,
            s=20,
            edgecolors='none'
        )

        ax.set_xlabel('Pocket Mean pLDDT (Confidence)')
        ax.set_ylabel('Pocket SASA (Å²)')
        ax.set_title(f'{site}: Accessibility vs Confidence', fontweight='bold')

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('# Domains')

        # Add trend line
        z = np.polyfit(site_df['pocket_mean_plddt'], site_df['pocket_total_sasa'], 1)
        p = np.poly1d(z)
        x_line = np.linspace(site_df['pocket_mean_plddt'].min(), site_df['pocket_mean_plddt'].max(), 100)
        ax.plot(x_line, p(x_line), 'r--', alpha=0.8, label=f'Trend')

        # Correlation
        corr = site_df['pocket_mean_plddt'].corr(site_df['pocket_total_sasa'])
        ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
               fontsize=10, va='top')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig2_scatter_sasa_vs_plddt.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig2_scatter_sasa_vs_plddt.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig2_scatter_sasa_vs_plddt")


def fig3_domain_count_vs_metrics(data, output_dir):
    """
    Summary plot: how metrics change with increasing domain count.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites'].copy()
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    metrics = [
        ('pocket_total_sasa', 'Pocket SASA (Å²)'),
        ('pocket_mean_plddt', 'Pocket pLDDT'),
        ('catalytic_mean_sasa', 'Catalytic Residue SASA (Å²)'),
        ('pocket_fraction_exposed', 'Fraction Exposed (>40Å²)')
    ]

    for ax, (metric, label) in zip(axes.flatten(), metrics):
        if metric not in df.columns:
            ax.text(0.5, 0.5, f'{label} not available', ha='center', va='center')
            continue

        for site in ['MD1', 'MD2', 'MD3', 'ART']:
            site_df = df[df['active_site'] == site]

            # Group by domain count
            grouped = site_df.groupby('n_domains')[metric].agg(['mean', 'std', 'count'])
            grouped = grouped[grouped['count'] >= 10]  # Filter low counts

            ax.errorbar(
                grouped.index,
                grouped['mean'],
                yerr=grouped['std'] / np.sqrt(grouped['count']),  # SEM
                label=site,
                marker='o',
                capsize=3,
                color=SITE_COLORS[site],
                linewidth=2,
                markersize=8
            )

        ax.set_xlabel('Number of Domains in Structure')
        ax.set_ylabel(label)
        ax.set_title(f'{label} vs Domain Context', fontweight='bold')
        ax.legend(loc='best')
        ax.set_xticks(sorted(df['n_domains'].unique()))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig3_domain_count_vs_metrics.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig3_domain_count_vs_metrics.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig3_domain_count_vs_metrics")


def fig4_interdomain_contacts_vs_domain_count(data, output_dir):
    """
    How interdomain contacts change with total domain count.
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain data not available")
        return

    df = data['interdomain'].copy()
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Plot 1: Contact density vs domain count
    ax = axes[0]
    grouped = df.groupby('n_domains').agg({
        'contact_density': ['mean', 'std'],
        'n_contacts': ['mean', 'std']
    })

    ax.errorbar(
        grouped.index,
        grouped[('contact_density', 'mean')],
        yerr=grouped[('contact_density', 'std')],
        marker='o',
        capsize=3,
        color='#1f77b4',
        linewidth=2,
        markersize=8
    )
    ax.set_xlabel('Number of Domains in Structure')
    ax.set_ylabel('Mean Contact Density')
    ax.set_title('Inter-Domain Contact Density vs Domain Count', fontweight='bold')

    # Plot 2: Number of contacts vs domain count
    ax = axes[1]
    ax.errorbar(
        grouped.index,
        grouped[('n_contacts', 'mean')],
        yerr=grouped[('n_contacts', 'std')],
        marker='s',
        capsize=3,
        color='#ff7f0e',
        linewidth=2,
        markersize=8
    )
    ax.set_xlabel('Number of Domains in Structure')
    ax.set_ylabel('Mean Number of Contacts')
    ax.set_title('Inter-Domain Contacts vs Domain Count', fontweight='bold')

    # Plot 3: Scatter of all domain pairs colored by total domains
    ax = axes[2]
    scatter = ax.scatter(
        df['d_com'],
        df['contact_density'],
        c=df['n_domains'],
        cmap=DOMAIN_COUNT_CMAP,
        alpha=0.3,
        s=10,
        edgecolors='none'
    )
    ax.set_xlabel('Center-of-Mass Distance (Å)')
    ax.set_ylabel('Contact Density')
    ax.set_title('Distance vs Contact Density\n(colored by # domains)', fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('# Domains')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig4_interdomain_contacts_vs_domains.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig4_interdomain_contacts_vs_domains.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig4_interdomain_contacts_vs_domains")


def fig5_composition_explained(data, output_dir):
    """
    Figure explaining what composition categories mean with example structures.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites'].copy()
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df[df['active_site'] == site].copy()

        if len(site_df) == 0 or 'pocket_total_sasa' not in site_df.columns:
            continue

        # Get categories and their descriptions
        categories = site_df['composition_category'].unique()

        # Create a better display name
        cat_data = []
        for cat in categories:
            cat_df = site_df[site_df['composition_category'] == cat]
            mean_sasa = cat_df['pocket_total_sasa'].mean()
            std_sasa = cat_df['pocket_total_sasa'].std()
            mean_domains = cat_df['n_domains'].mean()
            count = len(cat_df)

            # Get an example structure
            example = cat_df['structure'].iloc[0] if len(cat_df) > 0 else ''

            cat_data.append({
                'category': cat,
                'mean_sasa': mean_sasa,
                'std_sasa': std_sasa,
                'mean_domains': mean_domains,
                'count': count,
                'example': example
            })

        cat_df = pd.DataFrame(cat_data).sort_values('mean_domains')

        # Create grouped bar chart
        x = np.arange(len(cat_df))
        width = 0.6

        bars = ax.bar(x, cat_df['mean_sasa'], width, yerr=cat_df['std_sasa'],
                     capsize=5, color=SITE_COLORS[site], alpha=0.7,
                     edgecolor='black', linewidth=1)

        # Add domain count annotation
        for i, (_, row) in enumerate(cat_df.iterrows()):
            ax.text(i, row['mean_sasa'] + row['std_sasa'] + 20,
                   f"{row['mean_domains']:.1f} domains\nn={row['count']}",
                   ha='center', va='bottom', fontsize=9)

        ax.set_xlabel('Composition Category')
        ax.set_ylabel('Pocket SASA (Å²)')
        ax.set_title(f'{site} Active Site by Composition\n(with mean domain count)', fontweight='bold')

        # Better x-tick labels
        labels = [c.replace('_', '\n').replace('with', '+') for c in cat_df['category']]
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)

        ax.set_ylim(0, ax.get_ylim()[1] * 1.2)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig5_composition_explained.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig5_composition_explained.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig5_composition_explained")


def fig6_all_sites_violin_combined(data, output_dir):
    """
    Combined violin plot showing all active sites together, colored by domain count.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites'].copy()
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    if 'pocket_total_sasa' not in df.columns:
        print("  Skipping: pocket_total_sasa not available")
        return

    fig, ax = plt.subplots(figsize=(14, 8))

    sites = ['MD1', 'MD2', 'MD3', 'ART']
    positions = []
    all_data = []
    all_colors = []
    labels = []

    pos = 0
    for site in sites:
        site_df = df[df['active_site'] == site]
        domain_counts = sorted(site_df['n_domains'].unique())

        for n in domain_counts:
            subset = site_df[site_df['n_domains'] == n]['pocket_total_sasa'].dropna()
            if len(subset) > 5:
                all_data.append(subset.values)
                positions.append(pos)
                all_colors.append(DOMAIN_COUNT_CMAP(n / df['n_domains'].max()))
                labels.append(f'{site}\n({n}d)')
                pos += 1

        pos += 0.5  # Gap between sites

    if len(all_data) == 0:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center')
        plt.savefig(os.path.join(output_dir, 'fig6_all_sites_violin.svg'), bbox_inches='tight')
        plt.close()
        return

    # Create violin plot
    parts = ax.violinplot(all_data, positions=positions, showmeans=True, showmedians=True, widths=0.8)

    # Color each violin
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(all_colors[i])
        pc.set_alpha(0.7)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha='right')
    ax.set_ylabel('Pocket SASA (Å²)')
    ax.set_title('Active Site Accessibility Across All Domains\n(grouped by domain count)', fontweight='bold')

    # Add site separators
    site_boundaries = [0]
    pos = 0
    for site in sites:
        site_df = df[df['active_site'] == site]
        pos += len(site_df['n_domains'].unique()) + 0.5
        site_boundaries.append(pos - 0.25)

    for b in site_boundaries[1:-1]:
        ax.axvline(b, color='gray', linestyle='--', alpha=0.5)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=DOMAIN_COUNT_CMAP,
                               norm=plt.Normalize(df['n_domains'].min(), df['n_domains'].max()))
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, label='Number of Domains')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig6_all_sites_violin.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig6_all_sites_violin.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig6_all_sites_violin")


def fig7_bimodal_investigation(data, output_dir):
    """
    Investigate the bimodal distributions seen in Figure 2.
    """
    if 'active_sites' not in data:
        print("  Skipping: active_sites data not available")
        return

    df = data['active_sites'].copy()
    df['n_domains'] = df['structure'].apply(count_domains_in_structure)

    # Check if structure contains the L1 loop region
    df['has_L1'] = df['structure'].str.contains('md1l1', case=False)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, site in enumerate(['MD1', 'MD2', 'MD3', 'ART']):
        ax = axes[idx]
        site_df = df[df['active_site'] == site].copy()

        if len(site_df) == 0 or 'pocket_total_sasa' not in site_df.columns:
            continue

        # Split by L1 presence
        with_l1 = site_df[site_df['has_L1']]['pocket_total_sasa'].dropna()
        without_l1 = site_df[~site_df['has_L1']]['pocket_total_sasa'].dropna()

        # Create overlapping histograms
        bins = np.linspace(
            site_df['pocket_total_sasa'].min(),
            site_df['pocket_total_sasa'].max(),
            50
        )

        ax.hist(with_l1, bins=bins, alpha=0.6, label=f'With MD1L1 (n={len(with_l1)})',
               color='#1f77b4', density=True)
        ax.hist(without_l1, bins=bins, alpha=0.6, label=f'Without MD1L1 (n={len(without_l1)})',
               color='#ff7f0e', density=True)

        ax.set_xlabel('Pocket SASA (Å²)')
        ax.set_ylabel('Density')
        ax.set_title(f'{site}: SASA Distribution by MD1L1 Presence', fontweight='bold')
        ax.legend()

        # Add statistics
        if len(with_l1) > 0 and len(without_l1) > 0:
            ax.text(0.95, 0.95,
                   f'With L1: μ={with_l1.mean():.0f}, σ={with_l1.std():.0f}\n'
                   f'Without: μ={without_l1.mean():.0f}, σ={without_l1.std():.0f}',
                   transform=ax.transAxes, ha='right', va='top', fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig7_bimodal_investigation.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig7_bimodal_investigation.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig7_bimodal_investigation")


def fig8_interdomain_distance_violin(data, output_dir):
    """
    Violin plots of interdomain distances for MD1, MD2, MD3, WWE domain pairs
    vs number of domains in structure.
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain data not available")
        return

    df = data['interdomain'].copy()

    # Define domain pairs of interest (involving MD1, MD2, MD3, WWE)
    # Using both MD1 and MD1L1
    target_domains = ['MD1', 'MD1L1', 'MD2', 'MD3', 'WWE', 'ART']

    # Define specific pairs to plot
    domain_pairs = [
        ('MD1', 'MD2'), ('MD1L1', 'MD2'),
        ('MD2', 'MD3'),
        ('MD3', 'WWE'),
        ('WWE', 'ART'),
        ('MD1', 'MD3'), ('MD1L1', 'MD3'),
    ]

    # Create a combined pair column for filtering
    df['pair'] = df.apply(lambda r: tuple(sorted([r['domain1'], r['domain2']])), axis=1)

    # Filter to only pairs of interest
    target_pairs = [tuple(sorted(p)) for p in domain_pairs]
    df_filtered = df[df['pair'].isin(target_pairs)].copy()

    if len(df_filtered) == 0:
        print("  Skipping: no matching domain pairs found")
        return

    # Create readable pair names
    def format_pair(pair):
        return f"{pair[0]}-{pair[1]}"

    df_filtered['pair_name'] = df_filtered['pair'].apply(format_pair)

    # Get unique pairs that have data
    available_pairs = df_filtered['pair_name'].unique()
    n_pairs = len(available_pairs)

    if n_pairs == 0:
        print("  Skipping: no domain pairs with data")
        return

    # Create figure with subplots
    n_cols = min(3, n_pairs)
    n_rows = (n_pairs + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))

    if n_pairs == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, pair_name in enumerate(sorted(available_pairs)):
        ax = axes[idx]
        pair_df = df_filtered[df_filtered['pair_name'] == pair_name].copy()

        if len(pair_df) == 0:
            ax.set_visible(False)
            continue

        # Get domain counts
        domain_counts = sorted(pair_df['n_domains'].unique())

        # Filter to counts with data
        valid_counts = []
        valid_data = []
        for n in domain_counts:
            data_n = pair_df[pair_df['n_domains'] == n]['d_com'].dropna().values
            if len(data_n) > 0:
                valid_counts.append(n)
                valid_data.append(data_n)

        if len(valid_counts) == 0:
            ax.set_visible(False)
            continue

        # Create violin plot
        parts = ax.violinplot(valid_data, positions=valid_counts, showmeans=True, showmedians=True)

        # Color violins by domain count
        norm = plt.Normalize(min(valid_counts), max(valid_counts))
        for i, pc in enumerate(parts['bodies']):
            color = DOMAIN_COUNT_CMAP(norm(valid_counts[i]))
            pc.set_facecolor(color)
            pc.set_alpha(0.7)

        # Style the violin parts
        for partname in ['cmeans', 'cmedians', 'cbars', 'cmins', 'cmaxes']:
            if partname in parts:
                parts[partname].set_color('black')
                parts[partname].set_linewidth(1)

        # Add scatter points
        for i, n in enumerate(valid_counts):
            data_n = pair_df[pair_df['n_domains'] == n]['d_com'].dropna().values
            color = DOMAIN_COUNT_CMAP(norm(n))
            jitter = np.random.normal(0, 0.08, len(data_n))
            ax.scatter(n + jitter, data_n, alpha=0.3, s=10, color=color, zorder=3)

        ax.set_xlabel('Number of Domains in Structure')
        ax.set_ylabel('Center-of-Mass Distance (Å)')
        ax.set_title(f'{pair_name} Distance vs Domain Count', fontweight='bold')
        ax.set_xticks(valid_counts)

        # Add sample sizes
        for n in valid_counts:
            count = len(pair_df[pair_df['n_domains'] == n])
            ax.text(n, ax.get_ylim()[1] * 0.98, f'n={count}',
                   ha='center', va='top', fontsize=8)

    # Hide unused axes
    for idx in range(n_pairs, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig8_interdomain_distance_violin.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig8_interdomain_distance_violin.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig8_interdomain_distance_violin")


def fig9_md1_md3_kh7a_khb_comparison(data, output_dir):
    """
    Violin plots comparing MD1-MD3 distances for structures with KH7a,
    split by presence/absence of KHb-KH8.
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain data not available")
        return

    df = data['interdomain'].copy()

    # Check for KH7a and KHb-KH8 in structure name
    def has_kh7a(structure):
        return 'kh7a' in structure.lower()

    def has_khb_kh8(structure):
        return 'khb-kh8' in structure.lower() or 'khb_kh8' in structure.lower()

    df['has_KH7a'] = df['structure'].apply(has_kh7a)
    df['has_KHb_KH8'] = df['structure'].apply(has_khb_kh8)

    # Filter for structures with KH7a
    df_kh7a = df[df['has_KH7a']].copy()

    if len(df_kh7a) == 0:
        print("  Skipping: no structures with KH7a found")
        return

    # Create label for KHb-KH8 presence
    df_kh7a['KHb_KH8_status'] = df_kh7a['has_KHb_KH8'].map({True: 'With KHb-KH8', False: 'Without KHb-KH8'})

    # Filter for MD1-MD3 and MD1L1-MD3 pairs
    def is_md1_md3_pair(row):
        pair = tuple(sorted([row['domain1'], row['domain2']]))
        return pair in [('MD1', 'MD3'), ('MD1L1', 'MD3')]

    df_pairs = df_kh7a[df_kh7a.apply(is_md1_md3_pair, axis=1)].copy()

    if len(df_pairs) == 0:
        print("  Skipping: no MD1-MD3 pairs found in KH7a structures")
        return

    # Create pair name
    df_pairs['pair_name'] = df_pairs.apply(
        lambda r: f"{r['domain1']}-{r['domain2']}" if r['domain1'] < r['domain2']
        else f"{r['domain2']}-{r['domain1']}", axis=1
    )

    # Get unique pairs
    pairs = sorted(df_pairs['pair_name'].unique())

    # Create figure
    fig, axes = plt.subplots(1, len(pairs), figsize=(6*len(pairs), 6))
    if len(pairs) == 1:
        axes = [axes]

    colors = {'With KHb-KH8': '#2ecc71', 'Without KHb-KH8': '#e74c3c'}

    for idx, pair_name in enumerate(pairs):
        ax = axes[idx]
        pair_df = df_pairs[df_pairs['pair_name'] == pair_name]

        # Prepare data for violin plot
        categories = ['Without KHb-KH8', 'With KHb-KH8']
        plot_data = []
        positions = []
        valid_categories = []

        for i, cat in enumerate(categories):
            cat_data = pair_df[pair_df['KHb_KH8_status'] == cat]['d_com'].dropna().values
            if len(cat_data) > 0:
                plot_data.append(cat_data)
                positions.append(i)
                valid_categories.append(cat)

        if len(plot_data) == 0:
            ax.set_visible(False)
            continue

        # Create violin plot
        parts = ax.violinplot(plot_data, positions=positions, showmeans=True, showmedians=True)

        # Color violins
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[valid_categories[i]])
            pc.set_alpha(0.7)

        # Style the violin parts
        for partname in ['cmeans', 'cmedians', 'cbars', 'cmins', 'cmaxes']:
            if partname in parts:
                parts[partname].set_color('black')
                parts[partname].set_linewidth(1)

        # Add scatter points with jitter
        for i, cat in enumerate(valid_categories):
            cat_data = pair_df[pair_df['KHb_KH8_status'] == cat]['d_com'].dropna().values
            jitter = np.random.normal(0, 0.05, len(cat_data))
            ax.scatter(positions[i] + jitter, cat_data, alpha=0.3, s=15,
                      color=colors[cat], zorder=3, edgecolors='none')

        ax.set_ylabel('Center-of-Mass Distance (Å)')
        ax.set_title(f'{pair_name} Distance\n(structures with KH7a)', fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels(valid_categories)

        # Add sample sizes and statistics
        for i, cat in enumerate(valid_categories):
            cat_data = pair_df[pair_df['KHb_KH8_status'] == cat]['d_com'].dropna()
            n = len(cat_data)
            mean = cat_data.mean()
            std = cat_data.std()
            ax.text(positions[i], ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02,
                   f'n={n}\nμ={mean:.1f}\nσ={std:.1f}',
                   ha='center', va='bottom', fontsize=9)

        # Add statistical comparison if both groups present
        if len(valid_categories) == 2:
            from scipy import stats
            group1 = pair_df[pair_df['KHb_KH8_status'] == 'Without KHb-KH8']['d_com'].dropna()
            group2 = pair_df[pair_df['KHb_KH8_status'] == 'With KHb-KH8']['d_com'].dropna()
            if len(group1) > 1 and len(group2) > 1:
                t_stat, p_val = stats.ttest_ind(group1, group2)
                significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                ax.text(0.5, 0.95, f'p={p_val:.2e} ({significance})',
                       transform=ax.transAxes, ha='center', va='top', fontsize=10,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig9_md1_md3_kh7a_khb_comparison.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig9_md1_md3_kh7a_khb_comparison.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig9_md1_md3_kh7a_khb_comparison")


def fig10_active_site_distance_vs_length(data, output_dir):
    """
    Scatter plots of MD1/MD2/MD3/ART pairwise distances vs structure length (n_residues).
    All 6 possible pairs in a 2x3 grid.
    """
    if 'interdomain' not in data:
        print("  Skipping: interdomain data not available")
        return

    df = data['interdomain'].copy()

    # Define the 4 domains and all 6 pairs
    targets = ['MD1', 'MD2', 'MD3', 'ART']
    pairs = []
    for i in range(len(targets)):
        for j in range(i + 1, len(targets)):
            pairs.append((targets[i], targets[j]))

    # Pair colors: blend the two domain colors
    pair_colors = {
        ('MD1', 'MD2'): '#8b6b5e',
        ('MD1', 'MD3'): '#1a8a6c',
        ('MD1', 'ART'): '#6a4f6e',
        ('MD2', 'MD3'): '#96931d',
        ('MD2', 'ART'): '#e8541b',
        ('MD3', 'ART'): '#507c28',
    }

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    axes = axes.flatten()

    for idx, (d1, d2) in enumerate(pairs):
        ax = axes[idx]

        # Filter for this pair (order may be swapped in data)
        mask = ((df['domain1'] == d1) & (df['domain2'] == d2)) | \
               ((df['domain1'] == d2) & (df['domain2'] == d1))
        pair_df = df[mask].dropna(subset=['d_com', 'n_residues'])

        if len(pair_df) == 0:
            ax.text(0.5, 0.5, f'No data for\n{d1}-{d2}',
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
            ax.set_title(f'{d1}-{d2}', fontweight='bold')
            continue

        color = pair_colors[(d1, d2)]

        # Scatter plot
        ax.scatter(pair_df['n_residues'], pair_df['d_com'],
                   alpha=0.15, s=8, color=color, edgecolors='none')

        # Add trend line (LOWESS-style using binned means)
        bin_edges = np.linspace(pair_df['n_residues'].min(), pair_df['n_residues'].max(), 20)
        bin_centers = []
        bin_means = []
        bin_stds = []
        for k in range(len(bin_edges) - 1):
            mask_bin = (pair_df['n_residues'] >= bin_edges[k]) & \
                       (pair_df['n_residues'] < bin_edges[k + 1])
            bin_data = pair_df.loc[mask_bin, 'd_com']
            if len(bin_data) >= 3:
                bin_centers.append((bin_edges[k] + bin_edges[k + 1]) / 2)
                bin_means.append(bin_data.mean())
                bin_stds.append(bin_data.std())

        if len(bin_centers) >= 2:
            bin_centers = np.array(bin_centers)
            bin_means = np.array(bin_means)
            bin_stds = np.array(bin_stds)
            ax.plot(bin_centers, bin_means, color='black', linewidth=2, zorder=5, label='Binned mean')
            ax.fill_between(bin_centers, bin_means - bin_stds, bin_means + bin_stds,
                           alpha=0.2, color=color, zorder=4, label='±1 SD')

        # Linear regression for overall trend (only if x values vary)
        from scipy import stats
        n = len(pair_df)
        if pair_df['n_residues'].nunique() > 1:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                pair_df['n_residues'], pair_df['d_com'])
            x_line = np.array([pair_df['n_residues'].min(), pair_df['n_residues'].max()])
            ax.plot(x_line, slope * x_line + intercept, '--', color='darkred',
                    linewidth=1.5, zorder=6, label=f'Linear fit (R²={r_value**2:.3f})')
            stat_text = f'n={n}\nslope={slope:.3f} Å/res\np={p_value:.2e}'
        else:
            stat_text = f'n={n}\n(single structure length)'

        ax.set_xlabel('Structure Length (residues)')
        ax.set_ylabel('Center-of-Mass Distance (Å)')
        ax.set_title(f'{d1}–{d2} Active Site Distance', fontweight='bold')
        ax.legend(fontsize=8, loc='upper left')

        # Annotate with stats
        ax.text(0.98, 0.02, stat_text,
                transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig10_active_site_distance_vs_length.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'fig10_active_site_distance_vs_length.png'), bbox_inches='tight', dpi=300)
    plt.close()
    print("  Saved: fig10_active_site_distance_vs_length")


def main():
    parser = argparse.ArgumentParser(description='Generate improved figures')
    parser.add_argument('--input', '-i', default='../results/per_structure',
                       help='Input directory with analysis results')
    parser.add_argument('--output', '-o', default='../figures/improved',
                       help='Output directory for figures')
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    data = load_data(args.input)

    if not data:
        print("No data loaded.")
        return

    print("\nGenerating improved figures...")

    fig1_violin_accessibility_by_domains(data, output_dir)
    fig2_scatter_sasa_vs_plddt(data, output_dir)
    fig3_domain_count_vs_metrics(data, output_dir)
    fig4_interdomain_contacts_vs_domain_count(data, output_dir)
    fig5_composition_explained(data, output_dir)
    fig6_all_sites_violin_combined(data, output_dir)
    fig7_bimodal_investigation(data, output_dir)
    fig8_interdomain_distance_violin(data, output_dir)
    fig9_md1_md3_kh7a_khb_comparison(data, output_dir)
    fig10_active_site_distance_vs_length(data, output_dir)

    print(f"\nAll figures saved to {output_dir}")


if __name__ == '__main__':
    main()
