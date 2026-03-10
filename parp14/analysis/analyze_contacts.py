#!/usr/bin/env python3
"""
Generate contact maps from AlphaFold3 predicted structures.
Analyzes inter-residue contacts and domain-domain interfaces.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import MMCIFParser, NeighborSearch
from pathlib import Path


def parse_structure(cif_file):
    """Parse CIF structure file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)
    return structure


def calculate_contact_map(structure, distance_cutoff=8.0):
    """
    Calculate contact map based on CA-CA distances.
    
    Args:
        structure: Bio.PDB structure object
        distance_cutoff: Maximum distance (Angstroms) to consider as contact
    
    Returns:
        contact_matrix: NxN numpy array of contacts
        residue_info: List of residue information
    """
    # Get all CA atoms
    ca_atoms = []
    residue_info = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    for atom in residue:
                        if atom.name == 'CA':
                            ca_atoms.append(atom)
                            residue_info.append({
                                'chain': chain.id,
                                'resname': residue.resname,
                                'resid': residue.id[1]
                            })
                            break
    
    n_residues = len(ca_atoms)
    contact_matrix = np.zeros((n_residues, n_residues))
    distance_matrix = np.zeros((n_residues, n_residues))
    
    # Calculate pairwise distances
    for i in range(n_residues):
        for j in range(i, n_residues):
            dist = ca_atoms[i] - ca_atoms[j]  # Distance between atoms
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
            
            if dist <= distance_cutoff:
                contact_matrix[i, j] = 1
                contact_matrix[j, i] = 1
    
    return contact_matrix, distance_matrix, residue_info


def plot_contact_map(contact_matrix, residue_info, output_file, title=None):
    """Plot and save contact map."""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Plot contact map
    im = ax.imshow(contact_matrix, cmap='Blues', aspect='auto', origin='lower')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Contact')
    
    # Labels
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('Residue Index')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Contact Map (8Å cutoff)')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved contact map to {output_file}")


def plot_distance_map(distance_matrix, output_file, title=None, vmax=30):
    """Plot distance matrix."""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Mask diagonal and distances > vmax
    masked_dist = np.ma.masked_where(distance_matrix > vmax, distance_matrix)
    
    # Plot distance map
    im = ax.imshow(masked_dist, cmap='viridis_r', aspect='auto', 
                   origin='lower', vmin=0, vmax=vmax)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Distance (Å)')
    
    # Labels
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('Residue Index')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title(f'Distance Map (max {vmax}Å)')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved distance map to {output_file}")


def analyze_domain_interfaces(contact_matrix, domain_boundaries):
    """
    Analyze contacts between domains.
    
    Args:
        contact_matrix: NxN contact matrix
        domain_boundaries: Dict mapping domain names to (start, end) indices
    
    Returns:
        DataFrame with inter-domain contact counts
    """
    interface_data = []
    
    domain_list = list(domain_boundaries.keys())
    
    for i, domain1 in enumerate(domain_list):
        start1, end1 = domain_boundaries[domain1]
        
        for domain2 in domain_list[i+1:]:
            start2, end2 = domain_boundaries[domain2]
            
            # Extract sub-matrix for domain pair
            sub_matrix = contact_matrix[start1:end1+1, start2:end2+1]
            num_contacts = np.sum(sub_matrix)
            
            interface_data.append({
                'domain1': domain1,
                'domain2': domain2,
                'num_contacts': num_contacts,
                'domain1_size': end1 - start1 + 1,
                'domain2_size': end2 - start2 + 1,
                'normalized_contacts': num_contacts / ((end1-start1+1) * (end2-start2+1))
            })
    
    return pd.DataFrame(interface_data)


def process_all_structures(output_dir, results_dir, distance_cutoff=8.0):
    """Process all structures in the output directory."""
    
    # Create subdirectories
    contact_dir = os.path.join(results_dir, 'contact_maps')
    distance_dir = os.path.join(results_dir, 'distance_maps')
    os.makedirs(contact_dir, exist_ok=True)
    os.makedirs(distance_dir, exist_ok=True)
    
    # Find all prediction directories
    pred_dirs = glob.glob(os.path.join(output_dir, '*'))
    pred_dirs = [d for d in pred_dirs if os.path.isdir(d)]
    
    print(f"Processing {len(pred_dirs)} structures...")
    
    summary_data = []
    
    for pred_dir in pred_dirs:
        domain_name = os.path.basename(pred_dir)
        print(f"\nProcessing: {domain_name}")
        
        # Find CIF file
        cif_files = glob.glob(os.path.join(pred_dir, '*.cif'))
        if not cif_files:
            print(f"  Warning: No CIF file found")
            continue
        
        cif_file = cif_files[0]
        
        try:
            # Parse structure
            structure = parse_structure(cif_file)
            
            # Calculate contact map
            contact_matrix, distance_matrix, residue_info = calculate_contact_map(
                structure, distance_cutoff
            )
            
            # Count contacts
            total_contacts = np.sum(np.triu(contact_matrix, k=1))
            n_residues = len(residue_info)
            contact_density = total_contacts / (n_residues * (n_residues - 1) / 2)
            
            # Plot contact map
            contact_file = os.path.join(contact_dir, f'{domain_name}_contact.png')
            plot_contact_map(
                contact_matrix, 
                residue_info, 
                contact_file,
                title=f'{domain_name} - Contact Map'
            )
            
            # Plot distance map
            dist_file = os.path.join(distance_dir, f'{domain_name}_distance.png')
            plot_distance_map(
                distance_matrix,
                dist_file,
                title=f'{domain_name} - Distance Map'
            )
            
            summary_data.append({
                'domain_combination': domain_name,
                'num_residues': n_residues,
                'total_contacts': total_contacts,
                'contact_density': contact_density,
                'avg_contact_per_residue': total_contacts / n_residues
            })
            
            print(f"  Residues: {n_residues}, Contacts: {total_contacts}, Density: {contact_density:.4f}")
            
        except Exception as e:
            print(f"  Error processing {domain_name}: {e}")
            continue
    
    # Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(results_dir, 'contact_analysis_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSaved contact analysis summary to: {summary_file}")
    
    return summary_df


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze contacts in AlphaFold3 structures')
    parser.add_argument('--output_dir', default='alphafold_outputs',
                        help='AlphaFold3 output directory')
    parser.add_argument('--results_dir', default='analysis_results',
                        help='Directory to save analysis results')
    parser.add_argument('--distance_cutoff', type=float, default=8.0,
                        help='Distance cutoff for contacts (Angstroms)')
    
    args = parser.parse_args()
    
    print("="*70)
    print("PARP14 Contact Map Analysis")
    print("="*70)
    print(f"\nInput directory: {args.output_dir}")
    print(f"Results directory: {args.results_dir}")
    print(f"Distance cutoff: {args.distance_cutoff}Å\n")
    
    # Process all structures
    summary_df = process_all_structures(
        args.output_dir,
        args.results_dir,
        args.distance_cutoff
    )
    
    print("\n" + "="*70)
    print("Contact analysis complete!")
    print("="*70)
    print(f"\nProcessed {len(summary_df)} structures")
    print(f"Contact maps saved to: {os.path.join(args.results_dir, 'contact_maps')}")
    print(f"Distance maps saved to: {os.path.join(args.results_dir, 'distance_maps')}")


if __name__ == '__main__':
    main()
