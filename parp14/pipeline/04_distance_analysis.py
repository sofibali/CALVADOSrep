#!/usr/bin/env python3
"""
Distance analysis of AlphaFold3 structures.

Outputs:
- Contact maps (PDF visualizations)
- Distance matrices (NPZ format for re-analysis)
- Domain distances (CSV format)

Usage:
    python 04_distance_analysis.py --input alphafold_outputs/ --domains Domain_boundries.csv --output analysis_output/

Requires: biopython, scipy
    pip install biopython scipy
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser
from scipy.spatial.distance import pdist, squareform
import csv
from collections import defaultdict
import warnings
import time
from datetime import datetime, timedelta
import sys
warnings.filterwarnings('ignore')


def read_domains(domain_file):
    """Read domain boundaries from CSV file.
    
    Expected format (comma-separated):
    Domain,Start,End
    RRM1,1,145
    MD1,738,1006
    ...
    """
    domains = {}
    domain_order = []
    
    try:
        with open(domain_file, 'r') as f:
            reader = csv.DictReader(f)
            
            # Validate headers
            if reader.fieldnames is None or not all(h in reader.fieldnames for h in ['Domain', 'Start', 'End']):
                print(f"\nERROR: Invalid CSV format in {domain_file}")
                print("Expected format (comma-separated):")
                print("  Domain,Start,End")
                print("  RRM1,1,145")
                print("  MD1,738,1006")
                print("  ...")
                print(f"\nFound headers: {reader.fieldnames}")
                sys.exit(1)
            
            for line_num, row in enumerate(reader, start=2):  # start=2 because header is line 1
                try:
                    domain_name = row['Domain'].strip()
                    start = int(row['Start'].strip())
                    end = int(row['End'].strip())
                    
                    if not domain_name:
                        print(f"WARNING: Empty domain name on line {line_num}, skipping")
                        continue
                    
                    if start >= end:
                        print(f"ERROR: Invalid range for {domain_name} on line {line_num}: Start ({start}) >= End ({end})")
                        sys.exit(1)
                    
                    domains[domain_name] = (start, end)
                    domain_order.append(domain_name)
                    
                except ValueError as e:
                    print(f"\nERROR: Invalid data on line {line_num} of {domain_file}")
                    print(f"  Start and End must be integers")
                    print(f"  Row: {row}")
                    print(f"  Error: {e}")
                    sys.exit(1)
                except KeyError as e:
                    print(f"\nERROR: Missing required column on line {line_num}: {e}")
                    print(f"  Required columns: Domain, Start, End")
                    sys.exit(1)
    
    except FileNotFoundError:
        print(f"\nERROR: Domain boundaries file not found: {domain_file}")
        sys.exit(1)
    
    if not domains:
        print(f"\nERROR: No valid domains found in {domain_file}")
        sys.exit(1)
    
    return domains, domain_order


def extract_ca_coords(structure):
    """Extract CA atom coordinates from BioPython structure."""
    ca_coords = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_num = residue.get_id()[1]
                if residue.has_id('CA'):
                    ca_coords[res_num] = residue['CA'].coord
    
    return ca_coords


def calculate_distance_matrix(ca_coords):
    """Calculate pairwise CA distance matrix using vectorized operations."""
    residues = sorted(ca_coords.keys())
    n = len(residues)
    
    # Stack coordinates into array
    coords_array = np.array([ca_coords[res] for res in residues])
    
    # Use scipy's pdist for fast pairwise distance calculation
    # pdist returns condensed distance matrix, squareform converts to square matrix
    dist_matrix = squareform(pdist(coords_array))
    
    return dist_matrix, residues


def map_structure_to_domains(structure_name, domains):
    """Map structure filename to actual domain boundaries.
    
    Structure filenames contain lowercase domain names (e.g., kh1_kh2_rrm1_md1)
    but Domain_boundries.csv may have compound names (e.g., KH1-KH6).
    
    This function identifies which domains from the CSV are present in the structure.
    """
    # Parse structure name into individual tokens
    tokens = structure_name.lower().replace('-', '_').split('_')
    tokens_set = set(tokens)
    
    matched_domains = []
    
    for domain_name, (start, end) in domains.items():
        domain_lower = domain_name.lower()
        
        # Check if domain name appears in tokens
        # Handle compound domains like "KH1-KH6" or "KHb-KH8"
        if '-' in domain_name:
            # Compound domain - check if any component is in tokens
            components = domain_lower.replace('-', '_').split('_')
            if any(comp in tokens_set for comp in components):
                matched_domains.append(domain_name)
        else:
            # Simple domain - direct match
            if domain_lower in tokens_set:
                matched_domains.append(domain_name)
    
    return matched_domains


def calculate_domain_distances(ca_coords, domains, matched_domains):
    """Calculate inter-domain center-of-mass distances."""
    domain_coms = {}
    
    # Calculate COM for each domain
    for domain_name in matched_domains:
        if domain_name in domains:
            start, end = domains[domain_name]
            domain_cas = [ca_coords[rn] for rn in range(start, end + 1) if rn in ca_coords]
            if domain_cas:
                domain_coms[domain_name] = np.mean(domain_cas, axis=0)
    
    # Calculate pairwise distances
    distances = {}
    domain_list = sorted(domain_coms.keys())
    
    for i in range(len(domain_list)):
        for j in range(i + 1, len(domain_list)):
            dom1, dom2 = domain_list[i], domain_list[j]
            dist = np.linalg.norm(domain_coms[dom1] - domain_coms[dom2])
            distances[f'{dom1}-{dom2}'] = dist
    
    return distances


def plot_contact_map(dist_matrix, residues, output_path, title="Contact Map"):
    """Generate and save contact map visualization."""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot distance matrix with color scale
    im = ax.imshow(dist_matrix, cmap='viridis_r', vmin=0, vmax=50, origin='lower')
    
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('Residue Index')
    ax.set_title(title)
    
    plt.colorbar(im, ax=ax, label='Distance (Å)')
    plt.tight_layout()
    plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Analyze distances in AlphaFold3 structures')
    parser.add_argument('--input', required=True, help='AlphaFold3 output directory')
    parser.add_argument('--domains', required=True, help='Domain boundaries CSV file')
    parser.add_argument('--output', default='distance_analysis', help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directories
    os.makedirs(args.output, exist_ok=True)
    contact_map_dir = os.path.join(args.output, 'contact_maps')
    distance_data_dir = os.path.join(args.output, 'distance_matrices')
    os.makedirs(contact_map_dir, exist_ok=True)
    os.makedirs(distance_data_dir, exist_ok=True)
    
    # Read domain boundaries
    print(f"Reading domain boundaries: {args.domains}")
    domains, domain_order = read_domains(args.domains)
    print(f"  Domains: {', '.join(domain_order)}\n")
    
    # Initialize BioPython CIF parser
    parser = MMCIFParser(QUIET=True)
    
    # Find all structure directories
    print(f"Scanning AlphaFold3 outputs: {args.input}")
    structure_dirs = sorted([d for d in os.listdir(args.input) 
                            if os.path.isdir(os.path.join(args.input, d))])
    
    # Count total models (structure × seeds)
    total_models = 0
    for struct_name in structure_dirs:
        struct_dir = os.path.join(args.input, struct_name)
        seed_dirs = [d for d in os.listdir(struct_dir) 
                    if os.path.isdir(os.path.join(struct_dir, d)) and 'seed-' in d]
        total_models += len(seed_dirs)
    
    print(f"  Found {len(structure_dirs)} structures with {total_models} total models\n")
    
    # Store domain distances for summary
    all_domain_distances = []
    
    processed = 0
    skipped = 0
    start_time = time.time()
    last_report_time = start_time
    
    print(f"Starting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Progress updates every 50 models...\n")
    sys.stdout.flush()
    
    # Track skip reasons
    skip_reasons = defaultdict(int)
    
    model_idx = 0
    for struct_name in structure_dirs:
        struct_dir = os.path.join(args.input, struct_name)
        
        # Find all seed directories
        seed_dirs = sorted([d for d in os.listdir(struct_dir) 
                           if os.path.isdir(os.path.join(struct_dir, d)) and 'seed-' in d])
        
        if not seed_dirs:
            skipped += 1
            skip_reasons['no_seed_dirs'] += 1
            if skip_reasons['no_seed_dirs'] <= 5:
                print(f"  SKIP: {struct_name} - No seed directories found")
                sys.stdout.flush()
            continue
        
        # Process each seed
        for seed_name in seed_dirs:
            model_idx += 1
            seed_dir = os.path.join(struct_dir, seed_name)
            cif_path = os.path.join(seed_dir, 'model.cif')
            model_name = f"{struct_name}/{seed_name}"
            
            if not os.path.isfile(cif_path):
                skipped += 1
                skip_reasons['no_cif_file'] += 1
                
                # Debug: show first 5 skipped structures
                if skip_reasons['no_cif_file'] <= 5:
                    print(f"  SKIP: {model_name} - CIF file not found at: {cif_path}")
                    sys.stdout.flush()
                continue
            
            try:
                # Load structure with BioPython
                structure = parser.get_structure('model', cif_path)
                
                # Extract CA coordinates
                ca_coords = extract_ca_coords(structure)
                
                if not ca_coords:
                    skipped += 1
                    skip_reasons['no_ca_atoms'] += 1
                    if skip_reasons['no_ca_atoms'] <= 5:
                        print(f"  SKIP: {model_name} - No CA atoms found in structure")
                        sys.stdout.flush()
                    continue
                
                # Calculate distance matrix
                dist_matrix, residues = calculate_distance_matrix(ca_coords)
                
                # Save distance matrix as NPZ
                npz_name = f'{struct_name}_{seed_name}_distances.npz'
                npz_path = os.path.join(distance_data_dir, npz_name)
                np.savez_compressed(npz_path, 
                                   distance_matrix=dist_matrix,
                                   residue_numbers=residues,
                                   structure_name=struct_name,
                                   seed_name=seed_name)
                
                # Generate contact map
                contact_map_name = f'{struct_name}_{seed_name}_contact_map.pdf'
                contact_map_path = os.path.join(contact_map_dir, contact_map_name)
                plot_contact_map(dist_matrix, residues, contact_map_path, 
                               title=f'Contact Map: {struct_name} ({seed_name})')
                
                # Calculate domain distances
                matched_domains = map_structure_to_domains(struct_name, domains)
                domain_distances = calculate_domain_distances(ca_coords, domains, matched_domains)
                
                for pair, dist in domain_distances.items():
                    all_domain_distances.append({
                        'structure': struct_name,
                        'seed': seed_name,
                        'domain_pair': pair,
                        'distance': dist
                    })
                
                processed += 1
                
                # Progress reporting every 50 models or every 5 minutes
                current_time = time.time()
                if processed % 50 == 0 or (current_time - last_report_time) >= 300:
                    elapsed_time = current_time - start_time
                    avg_time_per_model = elapsed_time / processed
                    remaining = total_models - model_idx
                    eta_seconds = avg_time_per_model * remaining
                    eta = timedelta(seconds=int(eta_seconds))
                    
                    print(f"[{datetime.now().strftime('%H:%M:%S')}] Progress: {model_idx}/{total_models} models")
                    print(f"  Processed: {processed} | Skipped: {skipped}")
                    print(f"  Remaining: {remaining} models")
                    print(f"  Rate: {processed/elapsed_time:.2f} models/sec ({avg_time_per_model:.2f} sec/model)")
                    print(f"  ETA: {eta} (estimated completion: {(datetime.now() + timedelta(seconds=eta_seconds)).strftime('%H:%M:%S')})")
                    print()
                    sys.stdout.flush()
                    last_report_time = current_time
            
            except Exception as e:
                skipped += 1
                error_type = type(e).__name__
                skip_reasons[f'error_{error_type}'] += 1
                if skip_reasons[f'error_{error_type}'] <= 5:
                    print(f"  ERROR: {model_name} - {error_type}: {str(e)[:100]}")
                    sys.stdout.flush()
                continue
    
    # Save domain distances summary
    if all_domain_distances:
        df_distances = pd.DataFrame(all_domain_distances)
        csv_path = os.path.join(args.output, 'domain_distances.csv')
        df_distances.to_csv(csv_path, index=False)
        print(f"\n✓ Domain distances saved to: {csv_path}")
    
    total_time = time.time() - start_time
    
    print(f"\n{'='*70}")
    print(f"Distance Analysis Complete")
    print(f"{'='*70}")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total time: {timedelta(seconds=int(total_time))}")
    print(f"Structures: {len(structure_dirs)}")
    print(f"Models processed: {processed}")
    print(f"Models skipped: {skipped}")
    if processed > 0:
        print(f"Average rate: {processed/total_time:.2f} models/second")
    
    if skip_reasons:
        print(f"\nSkip reasons:")
        for reason, count in sorted(skip_reasons.items(), key=lambda x: x[1], reverse=True):
            print(f"  {reason}: {count}")
    
    if skipped > 0 and processed == 0:
        print(f"\n⚠️  WARNING: All models were skipped!")
        print(f"  Check the directory structure and file names in: {args.input}")
        print(f"  Expected: <structure_name>/seed-N_sample-M/model.cif")
    
    print(f"\nOutput files:")
    print(f"  Contact maps: {contact_map_dir}/")
    print(f"  Distance matrices: {distance_data_dir}/")
    print(f"  Domain distances CSV: {args.output}/domain_distances.csv")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
