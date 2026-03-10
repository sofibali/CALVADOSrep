#!/usr/bin/env python
"""
10_active_site_analysis_fast.py - Optimized active site accessibility analysis

Efficiently processes large SASA/pLDDT CSV files by:
1. Pre-grouping data by structure/seed
2. Using vectorized pandas operations
3. Processing in batches

Usage:
    python 10_active_site_analysis_fast.py --sasa ../sasa_per_residue.csv --plddt ./pLDDT_per_residue.csv --output ../results/per_structure/active_site_metrics.csv
"""

import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Active site definitions
ACTIVE_SITES = {
    'MD1': {
        'domain_range': (790, 981),
        'catalytic': [831, 923, 962],
        'pocket': list(range(822, 837)) + list(range(919, 928)) + [961, 966],
    },
    'MD2': {
        'domain_range': (1005, 1193),
        'catalytic': [1035, 1046, 1134, 1171],
        'pocket': list(range(1021, 1025)) + list(range(1034, 1048)) +
                  list(range(1130, 1142)) + [1170, 1171, 1175, 1178],
    },
    'MD3': {
        'domain_range': (1207, 1388),
        'catalytic': [1248, 1259, 1330, 1371],
        'pocket': list(range(1235, 1238)) + list(range(1247, 1262)) +
                  list(range(1302, 1305)) + list(range(1324, 1338)) +
                  list(range(1369, 1372)) + [1375],
    },
    'ART': {
        'domain_range': (1603, 1801),
        'catalytic': [1684, 1705, 1706, 1722],
        'pocket': list(range(1681, 1686)) + [1688, 1701] + list(range(1704, 1710)) +
                  list(range(1714, 1717)) + [1721, 1722, 1726, 1727, 1781],
    }
}


def check_domain_in_structure(structure_name):
    """Check which domains are present in structure based on name."""
    name_lower = structure_name.lower()

    domain_patterns = {
        'MD1': ['md1', 'md1l1'],
        'MD2': ['md2'],
        'MD3': ['md3'],
        'ART': ['art'],
        'WWE': ['wwe'],
        'KH7a': ['kh7a'],
    }

    present = {}
    for domain, patterns in domain_patterns.items():
        present[domain] = any(p in name_lower for p in patterns)

    return present


def get_composition_category(domains_present, target_domain):
    """Categorize composition relative to target domain."""
    if target_domain == 'MD1':
        if domains_present.get('MD2', False) and domains_present.get('MD3', False):
            return 'MD1_with_MD2_MD3'
        elif domains_present.get('MD2', False):
            return 'MD1_with_MD2'
        elif domains_present.get('KH7a', False):
            return 'MD1_with_KH7a'
        else:
            return 'MD1_alone'
    elif target_domain == 'MD2':
        if domains_present.get('MD1', False) and domains_present.get('MD3', False):
            return 'MD2_with_MD1_MD3'
        elif domains_present.get('MD1', False):
            return 'MD2_with_MD1'
        elif domains_present.get('MD3', False):
            return 'MD2_with_MD3'
        else:
            return 'MD2_alone'
    elif target_domain == 'MD3':
        if domains_present.get('MD2', False) and domains_present.get('ART', False):
            return 'MD3_with_MD2_ART'
        elif domains_present.get('MD2', False):
            return 'MD3_with_MD2'
        elif domains_present.get('ART', False):
            return 'MD3_with_ART'
        else:
            return 'MD3_alone'
    elif target_domain == 'ART':
        if domains_present.get('MD3', False) and domains_present.get('WWE', False):
            return 'ART_with_MD3_WWE'
        elif domains_present.get('WWE', False):
            return 'ART_with_WWE'
        elif domains_present.get('MD3', False):
            return 'ART_with_MD3'
        else:
            return 'ART_alone'
    return f'{target_domain}_other'


def process_structure_group(sasa_group, plddt_group, structure, seed):
    """Process a single structure/seed group for all active sites."""
    results = []

    domains_present = check_domain_in_structure(structure)
    n_neighbors = sum(1 for d, v in domains_present.items() if v)

    for site_name, site_def in ACTIVE_SITES.items():
        # Skip if this domain isn't in the structure
        if not domains_present.get(site_name, False):
            continue

        composition_category = get_composition_category(domains_present, site_name)

        result = {
            'structure': structure,
            'seed': seed,
            'active_site': site_name,
            'composition_category': composition_category,
            'n_neighbors': n_neighbors,
        }

        # Add domain presence flags
        for d, v in domains_present.items():
            result[f'has_{d}'] = v

        # Get SASA for pocket residues
        if sasa_group is not None:
            pocket_sasa = sasa_group[sasa_group['residue_number'].isin(site_def['pocket'])]
            if len(pocket_sasa) > 0:
                result['pocket_total_sasa'] = pocket_sasa['sasa'].sum()
                result['pocket_mean_sasa'] = pocket_sasa['sasa'].mean()
                result['pocket_max_sasa'] = pocket_sasa['sasa'].max()
                result['pocket_n_exposed'] = (pocket_sasa['sasa'] > 40).sum()
                result['pocket_fraction_exposed'] = (pocket_sasa['sasa'] > 40).mean()

            # Catalytic SASA
            catalytic_sasa = sasa_group[sasa_group['residue_number'].isin(site_def['catalytic'])]
            if len(catalytic_sasa) > 0:
                result['catalytic_total_sasa'] = catalytic_sasa['sasa'].sum()
                result['catalytic_mean_sasa'] = catalytic_sasa['sasa'].mean()

        # Get pLDDT for pocket residues
        if plddt_group is not None:
            pocket_plddt = plddt_group[plddt_group['residue_number'].isin(site_def['pocket'])]
            if len(pocket_plddt) > 0:
                result['pocket_mean_plddt'] = pocket_plddt['plddt'].mean()
                result['pocket_min_plddt'] = pocket_plddt['plddt'].min()
                result['pocket_std_plddt'] = pocket_plddt['plddt'].std()
                result['pocket_n_high_conf'] = (pocket_plddt['plddt'] >= 90).sum()
                result['pocket_fraction_high_conf'] = (pocket_plddt['plddt'] >= 90).mean()

            # Domain pLDDT
            domain_start, domain_end = site_def['domain_range']
            domain_plddt = plddt_group[
                (plddt_group['residue_number'] >= domain_start) &
                (plddt_group['residue_number'] <= domain_end)
            ]
            if len(domain_plddt) > 0:
                result['domain_mean_plddt'] = domain_plddt['plddt'].mean()
                result['domain_std_plddt'] = domain_plddt['plddt'].std()

        results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(description='Analyze active site accessibility (optimized)')
    parser.add_argument('--sasa', default='../sasa_per_residue.csv',
                       help='SASA per-residue CSV file')
    parser.add_argument('--plddt', default='./pLDDT_per_residue.csv',
                       help='pLDDT per-residue CSV file')
    parser.add_argument('--output', '-o', default='../results/per_structure/active_site_metrics.csv',
                       help='Output CSV file')
    parser.add_argument('--summary-output', default='../results/per_structure/active_site_summary.csv',
                       help='Output summary CSV file')
    args = parser.parse_args()

    print("Loading data files...")

    # Load data with only needed columns
    print("  Loading SASA data...")
    sasa_df = pd.read_csv(args.sasa, usecols=['structure', 'seed', 'residue_number', 'sasa'])
    print(f"    Loaded {len(sasa_df)} SASA records")

    print("  Loading pLDDT data...")
    plddt_df = pd.read_csv(args.plddt, usecols=['structure', 'seed', 'residue_number', 'plddt'])
    print(f"    Loaded {len(plddt_df)} pLDDT records")

    # Get unique structure/seed combinations
    sasa_groups = sasa_df.groupby(['structure', 'seed'])
    plddt_groups = plddt_df.groupby(['structure', 'seed'])

    # Get all unique structure/seed pairs
    sasa_keys = set(sasa_groups.groups.keys())
    plddt_keys = set(plddt_groups.groups.keys())
    all_keys = sasa_keys | plddt_keys

    print(f"\nProcessing {len(all_keys)} structure/seed combinations...")

    all_results = []

    for structure, seed in tqdm(all_keys, desc="Analyzing"):
        # Get data for this structure/seed
        sasa_group = sasa_groups.get_group((structure, seed)) if (structure, seed) in sasa_keys else None
        plddt_group = plddt_groups.get_group((structure, seed)) if (structure, seed) in plddt_keys else None

        results = process_structure_group(sasa_group, plddt_group, structure, seed)
        all_results.extend(results)

    # Save results
    if all_results:
        df = pd.DataFrame(all_results)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False)
        print(f"\nSaved {len(df)} active site records to {args.output}")

        # Calculate summary by composition category
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        summary_cols = [c for c in numeric_cols if c not in ['n_neighbors']]

        summary = df.groupby(['active_site', 'composition_category']).agg({
            **{col: ['mean', 'std', 'count'] if col == summary_cols[0] else ['mean', 'std']
               for col in summary_cols if col in df.columns}
        }).reset_index()

        # Flatten columns
        summary.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col
                          for col in summary.columns.values]
        summary.to_csv(args.summary_output, index=False)
        print(f"Saved summary to {args.summary_output}")

        # Print summary
        print("\n=== Active Site Summary ===")
        for site in ACTIVE_SITES.keys():
            site_data = df[df['active_site'] == site]
            if len(site_data) > 0:
                print(f"\n{site}:")
                print(f"  Structures: {len(site_data)}")
                print(f"  Composition categories: {site_data['composition_category'].nunique()}")
                if 'pocket_total_sasa' in site_data.columns:
                    print(f"  Mean pocket SASA: {site_data['pocket_total_sasa'].mean():.1f} Å²")
                if 'pocket_mean_plddt' in site_data.columns:
                    print(f"  Mean pocket pLDDT: {site_data['pocket_mean_plddt'].mean():.1f}")
    else:
        print("No results generated.")


if __name__ == '__main__':
    main()
