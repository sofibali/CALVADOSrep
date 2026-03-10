#!/usr/bin/env python3
"""
Regenerate domain_distances.csv from existing distance matrices.

This script reads the pre-computed .npz distance matrix files and calculates
domain-domain COM distances without re-processing the original structures.
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path
import csv


def read_domains(domain_file):
    """Read domain boundaries from CSV file."""
    domains = {}
    
    with open(domain_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            domain_name = row['Domain'].strip()
            start = int(row['Start'].strip())
            end = int(row['End'].strip())
            domains[domain_name] = (start, end)
    
    return domains


def map_structure_to_domains(structure_name, domains):
    """Map structure filename to actual domain boundaries."""
    tokens = structure_name.lower().replace('-', '_').split('_')
    tokens_set = set(tokens)
    
    matched_domains = []
    
    for domain_name in domains.keys():
        domain_lower = domain_name.lower()
        
        # Handle compound domains like "KH1-KH6" or "KHb-KH8"
        if '-' in domain_name:
            components = domain_lower.replace('-', '_').split('_')
            if any(comp in tokens_set for comp in components):
                matched_domains.append(domain_name)
        else:
            # Simple domain - direct match
            if domain_lower in tokens_set:
                matched_domains.append(domain_name)
    
    return matched_domains


def calculate_domain_distances_from_matrix(dist_matrix, residues, domains, matched_domains):
    """Calculate inter-domain COM distances from distance matrix."""
    # Create mapping from residue number to index in matrix
    residue_to_idx = {rn: i for i, rn in enumerate(residues)}
    
    # Extract coordinates from distance matrix using MDS-like approach
    # For COM calculation, we just need relative distances
    domain_coms = {}
    
    for domain_name in matched_domains:
        if domain_name in domains:
            start, end = domains[domain_name]
            # Find indices of residues in this domain
            domain_indices = [residue_to_idx[rn] for rn in range(start, end + 1) 
                            if rn in residue_to_idx]
            
            if domain_indices:
                # Store the indices for distance calculation
                domain_coms[domain_name] = domain_indices
    
    # Calculate pairwise domain distances
    distances = {}
    domain_list = sorted(domain_coms.keys())
    
    for i in range(len(domain_list)):
        for j in range(i + 1, len(domain_list)):
            dom1, dom2 = domain_list[i], domain_list[j]
            
            # Calculate average distance between all residue pairs
            indices1 = domain_coms[dom1]
            indices2 = domain_coms[dom2]
            
            # Extract sub-matrix of inter-domain distances
            inter_distances = []
            for idx1 in indices1:
                for idx2 in indices2:
                    inter_distances.append(dist_matrix[idx1, idx2])
            
            if inter_distances:
                # Use minimum distance as representative (closest approach)
                # Or could use mean for average separation
                avg_dist = np.mean(inter_distances)
                distances[f'{dom1}-{dom2}'] = avg_dist
    
    return distances


def main():
    # Paths
    domains_file = '/home/sbali/parp14/Domain_boundries.csv'
    distance_matrices_dir = '/home/sbali/parp14/distance_analysis/distance_matrices'
    output_file = '/home/sbali/parp14/distance_analysis/domain_distances.csv'
    
    print("Loading domain boundaries...")
    domains = read_domains(domains_file)
    print(f"  Found {len(domains)} domains\n")
    
    print("Scanning distance matrices...")
    npz_files = sorted(Path(distance_matrices_dir).glob('*.npz'))
    print(f"  Found {len(npz_files)} distance matrix files\n")
    
    all_domain_distances = []
    processed = 0
    skipped = 0
    
    print("Calculating domain distances...")
    for npz_path in npz_files:
        try:
            # Load distance matrix (allow_pickle for compatibility)
            data = np.load(npz_path, allow_pickle=True)
            dist_matrix = data['distance_matrix']
            residues = data['residue_numbers']
            
            # Parse structure and seed names from filename
            # Format: kh1_kh2_art_seed-1_sample-0_distances.npz
            filename = npz_path.stem  # Remove .npz
            filename = filename.replace('_distances', '')  # Remove _distances
            
            # Split by seed to get structure name and seed info
            parts = filename.split('_seed-')
            if len(parts) == 2:
                struct_name = parts[0]
                seed_part = parts[1]  # e.g., "1_sample-0"
                seed_name = f'seed-{seed_part}'
            else:
                print(f"  SKIP: Cannot parse filename: {npz_path.name}")
                skipped += 1
                continue
            
            # Map structure to domains
            matched_domains = map_structure_to_domains(struct_name, domains)
            
            if len(matched_domains) < 2:
                # Need at least 2 domains for domain pairs
                skipped += 1
                continue
            
            # Calculate domain distances
            domain_distances = calculate_domain_distances_from_matrix(
                dist_matrix, residues, domains, matched_domains
            )
            
            # Store results
            for pair, dist in domain_distances.items():
                all_domain_distances.append({
                    'structure': struct_name,
                    'seed': seed_name,
                    'domain_pair': pair,
                    'distance': dist
                })
            
            processed += 1
            
            if processed % 500 == 0:
                print(f"  Processed: {processed}/{len(npz_files)} ({100*processed/len(npz_files):.1f}%)")
                # Save intermediate results
                if all_domain_distances:
                    df_temp = pd.DataFrame(all_domain_distances)
                    df_temp.to_csv(output_file, index=False)
                    print(f"    Saved {len(all_domain_distances)} measurements so far...")
        
        except Exception as e:
            print(f"  ERROR processing {npz_path.name}: {e}")
            skipped += 1
            continue
    
    # Save results
    if all_domain_distances:
        df_distances = pd.DataFrame(all_domain_distances)
        df_distances.to_csv(output_file, index=False)
        
        print(f"\n{'='*70}")
        print(f"✓ Successfully generated domain_distances.csv")
        print(f"{'='*70}")
        print(f"Processed: {processed} distance matrices")
        print(f"Skipped: {skipped} files")
        print(f"Total measurements: {len(all_domain_distances)}")
        print(f"Unique structures: {df_distances['structure'].nunique()}")
        print(f"Unique domain pairs: {df_distances['domain_pair'].nunique()}")
        print(f"\nOutput: {output_file}")
        print(f"{'='*70}")
    else:
        print("\n⚠️  No domain distances calculated!")
        print("Check that the distance matrices contain multi-domain structures.")


if __name__ == '__main__':
    main()
