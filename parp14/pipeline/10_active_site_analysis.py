#!/usr/bin/env python
"""
10_active_site_analysis.py - Comprehensive active site accessibility analysis

Combines multiple metrics for active sites (MD1, MD2, MD3, ART):
1. SASA (from existing analysis)
2. pLDDT confidence (from existing analysis)
3. fpocket pocket metrics (from 08_fpocket_analysis.py)
4. Domain composition context

Outputs unified analysis showing how domain composition affects active site accessibility.

Usage:
    python 10_active_site_analysis.py --sasa ../sasa_per_residue.csv --plddt ../pipeline/pLDDT_per_residue.csv --fpocket ../results/per_structure/fpocket_results.csv --output ../results/per_structure/active_site_metrics.csv
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
        'all_residues': list(range(790, 982))
    },
    'MD2': {
        'domain_range': (1005, 1193),
        'catalytic': [1035, 1046, 1134, 1171],
        'pocket': list(range(1021, 1025)) + list(range(1034, 1048)) +
                  list(range(1130, 1142)) + [1170, 1171, 1175, 1178],
        'all_residues': list(range(1005, 1194))
    },
    'MD3': {
        'domain_range': (1207, 1388),
        'catalytic': [1248, 1259, 1330, 1371],
        'pocket': list(range(1235, 1238)) + list(range(1247, 1262)) +
                  list(range(1302, 1305)) + list(range(1324, 1338)) +
                  list(range(1369, 1372)) + [1375],
        'all_residues': list(range(1207, 1389))
    },
    'ART': {
        'domain_range': (1603, 1801),
        'catalytic': [1684, 1705, 1706, 1722],
        'pocket': list(range(1681, 1686)) + [1688, 1701] + list(range(1704, 1710)) +
                  list(range(1714, 1717)) + [1721, 1722, 1726, 1727, 1781],
        'all_residues': list(range(1603, 1802))
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
        'RRM1': ['rrm1'],
        'RRM2': ['rrm2'],
        'RRM3': ['rrm3'],
        'KH1-6': ['kh1', 'kh2', 'kh3', 'kh4', 'kh5', 'kh6'],
        'KH7a': ['kh7a'],
        'KHb-KH8': ['khb', 'kh8'],
        'WWE': ['wwe']
    }

    present = {}
    for domain, patterns in domain_patterns.items():
        present[domain] = any(p in name_lower for p in patterns)

    return present


def get_composition_category(domains_present, target_domain):
    """Categorize composition relative to target domain."""
    # Count how many other domains are present
    other_domains = [d for d, v in domains_present.items() if v and d != target_domain]

    # Categorize by neighbors
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
        if domains_present.get('MD2', False) and domains_present.get('WWE', False):
            return 'MD3_with_MD2_WWE'
        elif domains_present.get('MD2', False):
            return 'MD3_with_MD2'
        elif domains_present.get('WWE', False) or domains_present.get('ART', False):
            return 'MD3_with_downstream'
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


def load_sasa_for_structure(sasa_df, structure, seed, residues):
    """Extract SASA values for specific residues."""
    mask = (sasa_df['structure'] == structure) & (sasa_df['seed'] == seed)
    struct_data = sasa_df[mask]

    if len(struct_data) == 0:
        return None

    residue_sasa = struct_data[struct_data['residue_number'].isin(residues)]

    if len(residue_sasa) == 0:
        return None

    return {
        'total_sasa': residue_sasa['sasa'].sum(),
        'mean_sasa': residue_sasa['sasa'].mean(),
        'max_sasa': residue_sasa['sasa'].max(),
        'min_sasa': residue_sasa['sasa'].min(),
        'n_exposed': (residue_sasa['sasa'] > 40).sum(),  # >40 Å² = exposed
        'fraction_exposed': (residue_sasa['sasa'] > 40).mean()
    }


def load_plddt_for_structure(plddt_df, structure, seed, residues):
    """Extract pLDDT values for specific residues."""
    mask = (plddt_df['structure'] == structure) & (plddt_df['seed'] == seed)
    struct_data = plddt_df[mask]

    if len(struct_data) == 0:
        return None

    residue_plddt = struct_data[struct_data['residue_number'].isin(residues)]

    if len(residue_plddt) == 0:
        return None

    return {
        'mean_plddt': residue_plddt['plddt'].mean(),
        'min_plddt': residue_plddt['plddt'].min(),
        'std_plddt': residue_plddt['plddt'].std(),
        'n_high_conf': (residue_plddt['plddt'] >= 90).sum(),
        'n_low_conf': (residue_plddt['plddt'] < 70).sum(),
        'fraction_high_conf': (residue_plddt['plddt'] >= 90).mean()
    }


def load_fpocket_for_structure(fpocket_df, structure, seed, active_site):
    """Extract fpocket metrics for active site."""
    if fpocket_df is None or len(fpocket_df) == 0:
        return None

    mask = ((fpocket_df['structure'] == structure) &
            (fpocket_df['seed'] == seed) &
            (fpocket_df['mapped_active_site'] == active_site))
    site_data = fpocket_df[mask]

    if len(site_data) == 0:
        return None

    # Take best pocket (highest overlap)
    best = site_data.loc[site_data['overlap_fraction'].idxmax()]

    return {
        'pocket_volume': best.get('volume'),
        'druggability': best.get('druggability_score'),
        'pocket_score': best.get('score'),
        'hydrophobicity': best.get('hydrophobicity_score'),
        'overlap_fraction': best.get('overlap_fraction'),
        'catalytic_covered': best.get('catalytic_covered'),
        'n_pockets': len(site_data)
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze active site accessibility')
    parser.add_argument('--sasa', default='../sasa_per_residue.csv',
                       help='SASA per-residue CSV file')
    parser.add_argument('--plddt', default='../pipeline/pLDDT_per_residue.csv',
                       help='pLDDT per-residue CSV file')
    parser.add_argument('--fpocket', default='../results/per_structure/fpocket_results.csv',
                       help='fpocket results CSV file')
    parser.add_argument('--output', '-o', default='../results/per_structure/active_site_metrics.csv',
                       help='Output CSV file')
    parser.add_argument('--summary-output', default='../results/per_structure/active_site_summary.csv',
                       help='Output summary CSV file')
    parser.add_argument('--chunk-size', type=int, default=100000,
                       help='Chunk size for reading large CSV files')
    args = parser.parse_args()

    print("Loading data files...")

    # Load fpocket results (if exists)
    fpocket_df = None
    if os.path.exists(args.fpocket):
        fpocket_df = pd.read_csv(args.fpocket)
        print(f"  Loaded {len(fpocket_df)} fpocket records")
    else:
        print(f"  Warning: fpocket results not found at {args.fpocket}")

    # Get unique structure/seed combinations from SASA file
    print("  Scanning SASA file for structures...")
    structures = set()

    # Read in chunks to handle large file
    for chunk in pd.read_csv(args.sasa, chunksize=args.chunk_size,
                            usecols=['structure', 'seed']):
        for _, row in chunk.drop_duplicates().iterrows():
            structures.add((row['structure'], row['seed']))

    print(f"  Found {len(structures)} structure/seed combinations")

    results = []

    # Process each structure
    print("\nAnalyzing active sites...")

    # Load data in chunks for memory efficiency
    sasa_chunks = pd.read_csv(args.sasa, chunksize=args.chunk_size)
    plddt_chunks = pd.read_csv(args.plddt, chunksize=args.chunk_size)

    # Concatenate chunks (for smaller datasets) or process iteratively
    # For very large files, would need more sophisticated chunked processing
    print("  Loading SASA data...")
    sasa_df = pd.concat(pd.read_csv(args.sasa, chunksize=args.chunk_size))

    print("  Loading pLDDT data...")
    plddt_df = pd.concat(pd.read_csv(args.plddt, chunksize=args.chunk_size))

    for structure, seed in tqdm(structures, desc="Processing"):
        domains_present = check_domain_in_structure(structure)

        # Analyze each active site
        for site_name, site_def in ACTIVE_SITES.items():
            # Skip if this domain isn't in the structure
            if not domains_present.get(site_name, False):
                continue

            # Get composition category
            composition_category = get_composition_category(domains_present, site_name)

            # Count neighboring domains
            n_neighbors = sum(1 for d, v in domains_present.items()
                            if v and d != site_name and d in ['MD1', 'MD2', 'MD3', 'ART', 'WWE', 'KH7a'])

            # Base result
            result = {
                'structure': structure,
                'seed': seed,
                'active_site': site_name,
                'composition_category': composition_category,
                'n_neighbors': n_neighbors,
                **{f'has_{d}': v for d, v in domains_present.items()}
            }

            # Get SASA for pocket residues
            sasa_pocket = load_sasa_for_structure(sasa_df, structure, seed, site_def['pocket'])
            if sasa_pocket:
                result.update({f'pocket_{k}': v for k, v in sasa_pocket.items()})

            # Get SASA for catalytic residues
            sasa_catalytic = load_sasa_for_structure(sasa_df, structure, seed, site_def['catalytic'])
            if sasa_catalytic:
                result.update({f'catalytic_{k}': v for k, v in sasa_catalytic.items()})

            # Get pLDDT for pocket residues
            plddt_pocket = load_plddt_for_structure(plddt_df, structure, seed, site_def['pocket'])
            if plddt_pocket:
                result.update({f'pocket_{k}': v for k, v in plddt_pocket.items()})

            # Get pLDDT for whole domain
            plddt_domain = load_plddt_for_structure(plddt_df, structure, seed, site_def['all_residues'])
            if plddt_domain:
                result.update({f'domain_{k}': v for k, v in plddt_domain.items()})

            # Get fpocket metrics
            fpocket_metrics = load_fpocket_for_structure(fpocket_df, structure, seed, site_name)
            if fpocket_metrics:
                result.update({f'fpocket_{k}': v for k, v in fpocket_metrics.items()})

            results.append(result)

    # Save results
    if results:
        df = pd.DataFrame(results)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False)
        print(f"\nSaved {len(df)} active site records to {args.output}")

        # Calculate summary by composition category
        summary = df.groupby(['active_site', 'composition_category']).agg({
            'pocket_total_sasa': ['mean', 'std', 'count'],
            'pocket_mean_sasa': ['mean', 'std'],
            'catalytic_mean_sasa': ['mean', 'std'],
            'pocket_mean_plddt': ['mean', 'std'],
            'fpocket_pocket_volume': ['mean', 'std'],
            'fpocket_druggability': ['mean', 'std'],
            'n_neighbors': ['mean']
        }).reset_index()

        # Flatten columns
        summary.columns = ['_'.join(col).strip('_') for col in summary.columns.values]
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
                if 'fpocket_pocket_volume' in site_data.columns:
                    vol = site_data['fpocket_pocket_volume'].dropna()
                    if len(vol) > 0:
                        print(f"  Mean pocket volume: {vol.mean():.1f} Å³")
    else:
        print("No results generated.")


if __name__ == '__main__':
    main()
