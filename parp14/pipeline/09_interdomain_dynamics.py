#!/usr/bin/env python
"""
09_interdomain_dynamics.py - Analyze inter-domain distances and contacts across domain compositions

This script calculates:
1. Center-of-mass distances between all domain pairs
2. Minimum residue-residue distances between domains
3. Contact counts (residues within 8Å)
4. Contact density (normalized by domain sizes)
5. Variance across seeds (conformational flexibility proxy)

Groups results by domain composition to answer:
- How does domain X change when domain Y is present/absent?

Usage:
    python 09_interdomain_dynamics.py --input ../alphafold_outputs --output ../results/per_structure/interdomain_distances.csv
"""

import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.spatial.distance import cdist
from Bio.PDB import MMCIFParser
import warnings
warnings.filterwarnings('ignore')

# Domain definitions with residue ranges
DOMAINS = {
    'RRM1': (1, 145),
    'RRM2': (146, 224),
    'RRM3': (225, 314),
    'KH1-6': (315, 737),
    'KH7a': (738, 789),
    'MD1': (790, 981),
    'MD1L1': (790, 1004),  # Extended version
    'MD2': (1005, 1193),
    'MD3': (1207, 1388),
    'KHb-KH8': (1389, 1533),
    'WWE': (1534, 1602),
    'ART': (1603, 1801)
}

# Key domains of interest
KEY_DOMAINS = ['MD1', 'MD2', 'MD3', 'ART']

# Contact distance threshold
CONTACT_THRESHOLD = 8.0  # Angstroms


def parse_structure(cif_path):
    """Parse mmCIF file and extract CA coordinates by residue number."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_path)

    residue_coords = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    resnum = residue.id[1]
                    if 'CA' in residue:
                        residue_coords[resnum] = residue['CA'].get_coord()

    return residue_coords


def identify_domains_in_structure(residue_coords, structure_name):
    """Identify which domains are present based on residue coverage."""
    residue_nums = set(residue_coords.keys())
    domains_present = {}

    for domain_name, (start, end) in DOMAINS.items():
        domain_residues = set(range(start, end + 1))
        overlap = residue_nums & domain_residues
        coverage = len(overlap) / len(domain_residues) if domain_residues else 0

        # Consider domain present if >50% coverage
        if coverage > 0.5:
            domains_present[domain_name] = {
                'start': start,
                'end': end,
                'coverage': coverage,
                'residues': sorted(overlap)
            }

    return domains_present


def calculate_domain_metrics(residue_coords, domains_present):
    """Calculate inter-domain distance and contact metrics."""
    results = []

    domain_names = list(domains_present.keys())

    for i, domain1 in enumerate(domain_names):
        d1_info = domains_present[domain1]
        d1_residues = d1_info['residues']
        d1_coords = np.array([residue_coords[r] for r in d1_residues if r in residue_coords])

        if len(d1_coords) == 0:
            continue

        # Calculate domain 1 center of mass
        d1_com = np.mean(d1_coords, axis=0)

        # Self metrics (domain properties)
        d1_rg = np.sqrt(np.mean(np.sum((d1_coords - d1_com)**2, axis=1)))

        for j, domain2 in enumerate(domain_names):
            if i >= j:  # Skip self and duplicates
                continue

            d2_info = domains_present[domain2]
            d2_residues = d2_info['residues']
            d2_coords = np.array([residue_coords[r] for r in d2_residues if r in residue_coords])

            if len(d2_coords) == 0:
                continue

            # Calculate domain 2 center of mass
            d2_com = np.mean(d2_coords, axis=0)

            # COM distance
            com_distance = np.linalg.norm(d1_com - d2_com)

            # Pairwise distances
            dist_matrix = cdist(d1_coords, d2_coords)

            # Min distance
            min_distance = np.min(dist_matrix)

            # Mean distance
            mean_distance = np.mean(dist_matrix)

            # Contact count (pairs within threshold)
            n_contacts = np.sum(dist_matrix < CONTACT_THRESHOLD)

            # Contact density (normalized)
            contact_density = n_contacts / (len(d1_coords) * len(d2_coords))

            # Interface residue counts
            interface_d1 = np.sum(np.min(dist_matrix, axis=1) < CONTACT_THRESHOLD)
            interface_d2 = np.sum(np.min(dist_matrix, axis=0) < CONTACT_THRESHOLD)

            results.append({
                'domain1': domain1,
                'domain2': domain2,
                'd_com': com_distance,
                'd_min': min_distance,
                'd_mean': mean_distance,
                'n_contacts': n_contacts,
                'contact_density': contact_density,
                'interface_residues_d1': interface_d1,
                'interface_residues_d2': interface_d2,
                'd1_size': len(d1_coords),
                'd2_size': len(d2_coords),
                'd1_rg': d1_rg
            })

    return results


def get_composition_string(domains_present):
    """Generate a canonical composition string."""
    # Order domains by their start position
    ordered = sorted(domains_present.keys(),
                    key=lambda d: domains_present[d]['start'])
    return '_'.join(ordered)


def analyze_structure(args):
    """Analyze a single structure."""
    structure_dir, seed_dir = args

    structure_name = os.path.basename(structure_dir)
    seed_name = os.path.basename(seed_dir)

    cif_path = os.path.join(seed_dir, 'model.cif')
    if not os.path.exists(cif_path):
        return None

    try:
        # Parse structure
        residue_coords = parse_structure(cif_path)

        if len(residue_coords) == 0:
            return None

        # Identify domains
        domains_present = identify_domains_in_structure(residue_coords, structure_name)

        if len(domains_present) < 2:
            # Need at least 2 domains for inter-domain analysis
            return None

        # Calculate metrics
        domain_metrics = calculate_domain_metrics(residue_coords, domains_present)

        # Get composition
        composition = get_composition_string(domains_present)

        # Add metadata to results
        results = []
        for metric in domain_metrics:
            metric.update({
                'structure': structure_name,
                'seed': seed_name,
                'composition': composition,
                'n_domains': len(domains_present),
                'n_residues': len(residue_coords)
            })

            # Add flags for key domains
            for key_domain in KEY_DOMAINS:
                metric[f'has_{key_domain}'] = key_domain in domains_present

            results.append(metric)

        return results

    except Exception as e:
        print(f"Error processing {structure_name}/{seed_name}: {e}")
        return None


def calculate_composition_statistics(df):
    """Calculate statistics grouped by composition and domain pair."""
    # Group by composition and domain pair
    grouped = df.groupby(['composition', 'domain1', 'domain2'])

    stats = grouped.agg({
        'd_com': ['mean', 'std', 'min', 'max'],
        'd_min': ['mean', 'std', 'min', 'max'],
        'n_contacts': ['mean', 'std'],
        'contact_density': ['mean', 'std'],
        'interface_residues_d1': ['mean'],
        'interface_residues_d2': ['mean'],
    }).reset_index()

    # Flatten column names
    stats.columns = ['_'.join(col).strip('_') for col in stats.columns.values]

    return stats


def calculate_key_domain_summary(df):
    """Calculate summary statistics for key domains (MD1, MD2, MD3, ART)."""
    summaries = []

    for key_domain in KEY_DOMAINS:
        # Filter to structures containing this domain
        has_domain = df[(df['domain1'] == key_domain) | (df['domain2'] == key_domain)]

        if len(has_domain) == 0:
            continue

        # Get unique compositions
        compositions = has_domain['composition'].unique()

        for comp in compositions:
            comp_data = has_domain[has_domain['composition'] == comp]

            # Get metrics where this domain is involved
            d1_data = comp_data[comp_data['domain1'] == key_domain]
            d2_data = comp_data[comp_data['domain2'] == key_domain]

            # Combine
            domain_data = pd.concat([
                d1_data[['d_com', 'd_min', 'contact_density', 'interface_residues_d1']].rename(
                    columns={'interface_residues_d1': 'interface_residues'}),
                d2_data[['d_com', 'd_min', 'contact_density', 'interface_residues_d2']].rename(
                    columns={'interface_residues_d2': 'interface_residues'})
            ])

            # Determine which other domains are present
            other_domains = set()
            for d in comp.split('_'):
                if d != key_domain:
                    other_domains.add(d)

            summaries.append({
                'key_domain': key_domain,
                'composition': comp,
                'other_domains': ', '.join(sorted(other_domains)),
                'n_other_domains': len(other_domains),
                'd_com_mean': domain_data['d_com'].mean(),
                'd_com_std': domain_data['d_com'].std(),
                'd_min_mean': domain_data['d_min'].mean(),
                'contact_density_mean': domain_data['contact_density'].mean(),
                'interface_residues_mean': domain_data['interface_residues'].mean(),
                'n_observations': len(domain_data)
            })

    return pd.DataFrame(summaries)


def main():
    parser = argparse.ArgumentParser(description='Analyze inter-domain distances and contacts')
    parser.add_argument('--input', '-i', default='../alphafold_outputs',
                       help='Input directory with AlphaFold3 outputs')
    parser.add_argument('--output', '-o', default='../results/per_structure/interdomain_distances.csv',
                       help='Output CSV file path')
    parser.add_argument('--stats-output', default='../results/per_structure/interdomain_statistics.csv',
                       help='Output statistics CSV file path')
    parser.add_argument('--summary-output', default='../results/per_structure/key_domain_summary.csv',
                       help='Output key domain summary CSV file path')
    parser.add_argument('--workers', '-w', type=int, default=8,
                       help='Number of parallel workers')
    args = parser.parse_args()

    input_dir = Path(args.input)

    # Find all structure/seed combinations
    tasks = []
    for structure_dir in sorted(input_dir.iterdir()):
        if not structure_dir.is_dir():
            continue

        for seed_dir in sorted(structure_dir.iterdir()):
            if not seed_dir.is_dir() or not seed_dir.name.startswith('seed-'):
                continue

            if (seed_dir / 'model.cif').exists():
                tasks.append((str(structure_dir), str(seed_dir)))

    print(f"Found {len(tasks)} structures to analyze")

    # Process in parallel
    all_results = []

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(analyze_structure, task): task for task in tasks}

        for future in tqdm(as_completed(futures), total=len(futures), desc="Analyzing"):
            result = future.result()
            if result:
                all_results.extend(result)

    # Save results
    if all_results:
        df = pd.DataFrame(all_results)

        # Ensure output directories exist
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)

        # Save main results
        df.to_csv(args.output, index=False)
        print(f"\nSaved {len(df)} inter-domain records to {args.output}")

        # Calculate and save statistics
        stats = calculate_composition_statistics(df)
        stats.to_csv(args.stats_output, index=False)
        print(f"Saved composition statistics to {args.stats_output}")

        # Calculate and save key domain summary
        summary = calculate_key_domain_summary(df)
        summary.to_csv(args.summary_output, index=False)
        print(f"Saved key domain summary to {args.summary_output}")

        # Print summary
        print("\n=== Summary ===")
        print(f"Structures analyzed: {df['structure'].nunique()}")
        print(f"Unique compositions: {df['composition'].nunique()}")
        print(f"Domain pairs analyzed: {len(df['domain1'].unique())} x {len(df['domain2'].unique())}")

        print("\n=== Key Domain Presence ===")
        for domain in KEY_DOMAINS:
            count = df[f'has_{domain}'].sum() / len(df['seed'].unique())
            print(f"  {domain}: {count:.0f} structures")

    else:
        print("No results generated.")


if __name__ == '__main__':
    main()
