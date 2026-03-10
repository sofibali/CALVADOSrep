#!/usr/bin/env python3
"""
Analyze domain-domain interfaces in multi-domain PARP14 structures.
Identifies interface residues and calculates interface properties.
"""

import os
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import MMCIFParser
from collections import defaultdict


def parse_structure(cif_file):
    """Parse CIF structure file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)
    return structure


def get_domain_boundaries_from_name(domain_name_str):
    """
    Parse domain combination name and map to approximate residue ranges.
    Returns dict of domain_name: (start_idx, end_idx) in the structure.
    
    Note: This is approximate - actual boundaries depend on the specific
    sequence extracted for each combination.
    """
    # Standard domain sizes (approximate)
    domain_sizes = {
        'rrm1': 145, 'rrm2': 79, 'rrm3': 90,
        'kh1': 82, 'kh2': 69, 'kh3': 65, 'kh4': 66, 'kh5': 73, 'kh6': 68,
        'kh7a': 52,
        'md1': 192, 'md2': 190, 'md3': 182,
        'khb': 65, 'kh8': 80,
        'wwe': 69, 'art': 199
    }
    
    # Parse domain names from the combination string
    domains = domain_name_str.lower().split('_')
    
    boundaries = {}
    current_pos = 0
    
    for domain in domains:
        if domain in domain_sizes:
            size = domain_sizes[domain]
            boundaries[domain] = (current_pos, current_pos + size - 1)
            current_pos += size
    
    return boundaries


def calculate_interface_residues(structure, domain_boundaries, distance_cutoff=5.0):
    """
    Identify interface residues between domains.
    
    Args:
        structure: Bio.PDB structure
        domain_boundaries: Dict of domain_name: (start, end) residue indices
        distance_cutoff: Distance cutoff for interface (Angstroms)
    
    Returns:
        DataFrame with interface residue pairs
    """
    # Get all CA atoms with domain labels
    ca_atoms = []
    
    for model in structure:
        for chain in model:
            res_idx = 0
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    # Find which domain this residue belongs to
                    domain_label = None
                    for domain, (start, end) in domain_boundaries.items():
                        if start <= res_idx <= end:
                            domain_label = domain
                            break
                    
                    for atom in residue:
                        if atom.name == 'CA':
                            ca_atoms.append({
                                'atom': atom,
                                'domain': domain_label,
                                'resid': residue.id[1],
                                'resname': residue.resname,
                                'res_idx': res_idx
                            })
                            break
                    res_idx += 1
    
    # Find interface residues
    interface_pairs = []
    
    for i, res1 in enumerate(ca_atoms):
        for res2 in ca_atoms[i+1:]:
            # Only consider inter-domain contacts
            if res1['domain'] != res2['domain'] and res1['domain'] and res2['domain']:
                dist = res1['atom'] - res2['atom']
                
                if dist <= distance_cutoff:
                    interface_pairs.append({
                        'domain1': res1['domain'],
                        'domain2': res2['domain'],
                        'resid1': res1['resid'],
                        'resname1': res1['resname'],
                        'resid2': res2['resid'],
                        'resname2': res2['resname'],
                        'distance': dist
                    })
    
    return pd.DataFrame(interface_pairs)


def summarize_interfaces(interface_df):
    """Summarize interface statistics."""
    if len(interface_df) == 0:
        return pd.DataFrame()
    
    # Group by domain pairs
    summary = []
    
    for (dom1, dom2), group in interface_df.groupby(['domain1', 'domain2']):
        summary.append({
            'domain1': dom1,
            'domain2': dom2,
            'num_interface_contacts': len(group),
            'avg_distance': group['distance'].mean(),
            'min_distance': group['distance'].min(),
            'unique_residues_dom1': group['resid1'].nunique(),
            'unique_residues_dom2': group['resid2'].nunique()
        })
    
    return pd.DataFrame(summary)


def plot_interface_network(interface_summary, output_file):
    """Plot domain interface network."""
    if len(interface_summary) == 0:
        print("No interfaces to plot")
        return
    
    import networkx as nx
    
    # Create network graph
    G = nx.Graph()
    
    # Add edges weighted by number of contacts
    for _, row in interface_summary.iterrows():
        G.add_edge(
            row['domain1'], 
            row['domain2'],
            weight=row['num_interface_contacts']
        )
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=3000, node_color='lightblue',
                           edgecolors='black', linewidths=2, ax=ax)
    
    # Draw edges with width proportional to contacts
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    max_weight = max(weights) if weights else 1
    
    nx.draw_networkx_edges(G, pos, width=[w/max_weight*10 for w in weights],
                           alpha=0.6, ax=ax)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)
    
    # Add edge labels (number of contacts)
    edge_labels = {(u, v): f"{G[u][v]['weight']}" for u, v in G.edges()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=8, ax=ax)
    
    ax.set_title('Domain Interface Network\n(Edge width ~ contact strength)', 
                 fontsize=14, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved interface network to {output_file}")


def process_all_structures(output_dir, results_dir, distance_cutoff=5.0):
    """Process all multi-domain structures."""
    
    interface_dir = os.path.join(results_dir, 'interface_analysis')
    os.makedirs(interface_dir, exist_ok=True)
    
    # Find all prediction directories
    pred_dirs = glob.glob(os.path.join(output_dir, '*'))
    pred_dirs = [d for d in pred_dirs if os.path.isdir(d)]
    
    print(f"Processing {len(pred_dirs)} structures...")
    
    all_interfaces = []
    
    for pred_dir in pred_dirs:
        domain_name = os.path.basename(pred_dir)
        
        # Skip single domain structures
        if '_' not in domain_name:
            continue
        
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
            
            # Get domain boundaries
            boundaries = get_domain_boundaries_from_name(domain_name)
            
            if len(boundaries) < 2:
                print(f"  Skipping: only one domain detected")
                continue
            
            # Calculate interfaces
            interface_df = calculate_interface_residues(
                structure, boundaries, distance_cutoff
            )
            
            if len(interface_df) > 0:
                interface_df['structure'] = domain_name
                all_interfaces.append(interface_df)
                
                # Summarize interfaces for this structure
                summary = summarize_interfaces(interface_df)
                
                print(f"  Found {len(summary)} domain interfaces:")
                for _, row in summary.iterrows():
                    print(f"    {row['domain1']} - {row['domain2']}: "
                          f"{row['num_interface_contacts']} contacts, "
                          f"avg dist {row['avg_distance']:.2f}Å")
                
                # Plot interface network for this structure
                network_file = os.path.join(interface_dir, f'{domain_name}_interface_network.png')
                plot_interface_network(summary, network_file)
            else:
                print(f"  No interfaces found")
            
        except Exception as e:
            print(f"  Error processing {domain_name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Combine all interfaces
    if all_interfaces:
        combined_df = pd.concat(all_interfaces, ignore_index=True)
        
        # Save detailed interface data
        detail_file = os.path.join(results_dir, 'interface_residues_detailed.csv')
        combined_df.to_csv(detail_file, index=False)
        print(f"\nSaved detailed interface data to: {detail_file}")
        
        # Create global summary
        global_summary = []
        for structure, group in combined_df.groupby('structure'):
            summary = summarize_interfaces(group)
            summary['structure'] = structure
            global_summary.append(summary)
        
        if global_summary:
            global_summary_df = pd.concat(global_summary, ignore_index=True)
            summary_file = os.path.join(results_dir, 'interface_analysis_summary.csv')
            global_summary_df.to_csv(summary_file, index=False)
            print(f"Saved interface summary to: {summary_file}")
            
            return global_summary_df
    
    return pd.DataFrame()


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze domain interfaces in structures')
    parser.add_argument('--output_dir', default='alphafold_outputs',
                        help='AlphaFold3 output directory')
    parser.add_argument('--results_dir', default='analysis_results',
                        help='Directory to save analysis results')
    parser.add_argument('--distance_cutoff', type=float, default=5.0,
                        help='Distance cutoff for interface (Angstroms)')
    
    args = parser.parse_args()
    
    print("="*70)
    print("PARP14 Domain Interface Analysis")
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
    print("Interface analysis complete!")
    print("="*70)
    if len(summary_df) > 0:
        print(f"\nAnalyzed interfaces in {summary_df['structure'].nunique()} structures")
        print(f"Total interface pairs: {len(summary_df)}")


if __name__ == '__main__':
    main()
