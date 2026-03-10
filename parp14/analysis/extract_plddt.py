#!/usr/bin/env python3
"""
Extract pLDDT scores from AlphaFold3 predictions and generate visualizations.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import MMCIFParser


def extract_plddt_from_cif(cif_file):
    """Extract pLDDT scores from CIF file (stored in B-factor column)."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)
    
    plddt_data = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    for atom in residue:
                        if atom.name == 'CA':  # Use CA atom
                            plddt_data.append({
                                'chain': chain.id,
                                'resid': residue.id[1],
                                'resname': residue.resname,
                                'plddt': atom.bfactor
                            })
                            break
    
    return pd.DataFrame(plddt_data)


def plot_plddt_per_residue(plddt_df, output_file, title=None):
    """Plot pLDDT scores along the sequence."""
    fig, ax = plt.subplots(figsize=(14, 5))
    
    # Color regions by confidence
    colors = []
    for plddt in plddt_df['plddt']:
        if plddt >= 90:
            colors.append('#0053D6')  # Very high (dark blue)
        elif plddt >= 70:
            colors.append('#65CBF3')  # Confident (light blue)
        elif plddt >= 50:
            colors.append('#FFDB13')  # Low (yellow)
        else:
            colors.append('#FF7D45')  # Very low (orange)
    
    # Plot bars
    ax.bar(plddt_df['resid'], plddt_df['plddt'], color=colors, width=1.0, edgecolor='none')
    
    # Add confidence threshold lines
    ax.axhline(90, color='darkblue', linestyle='--', linewidth=1, alpha=0.5, label='Very high (≥90)')
    ax.axhline(70, color='royalblue', linestyle='--', linewidth=1, alpha=0.5, label='Confident (≥70)')
    ax.axhline(50, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Low (≥50)')
    
    # Labels and formatting
    ax.set_xlabel('Residue Index', fontsize=12)
    ax.set_ylabel('pLDDT Score', fontsize=12)
    ax.set_ylim(0, 100)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved pLDDT plot to {output_file}")


def categorize_confidence(plddt):
    """Categorize pLDDT score into confidence levels."""
    if plddt >= 90:
        return 'very_high'
    elif plddt >= 70:
        return 'confident'
    elif plddt >= 50:
        return 'low'
    else:
        return 'very_low'


def process_all_structures(output_dir, results_dir):
    """Extract pLDDT scores from all structures."""
    
    plddt_dir = os.path.join(results_dir, 'plddt_plots')
    os.makedirs(plddt_dir, exist_ok=True)
    
    # Find all prediction directories
    pred_dirs = glob.glob(os.path.join(output_dir, '*'))
    pred_dirs = [d for d in pred_dirs if os.path.isdir(d)]
    
    print(f"Processing {len(pred_dirs)} structures...")
    
    all_plddt_data = []
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
            # Extract pLDDT scores
            plddt_df = extract_plddt_from_cif(cif_file)
            plddt_df['structure'] = domain_name
            
            # Add confidence categories
            plddt_df['confidence'] = plddt_df['plddt'].apply(categorize_confidence)
            
            all_plddt_data.append(plddt_df)
            
            # Calculate statistics
            mean_plddt = plddt_df['plddt'].mean()
            median_plddt = plddt_df['plddt'].median()
            
            # Count residues by confidence level
            conf_counts = plddt_df['confidence'].value_counts()
            
            summary_data.append({
                'structure': domain_name,
                'num_residues': len(plddt_df),
                'mean_plddt': mean_plddt,
                'median_plddt': median_plddt,
                'min_plddt': plddt_df['plddt'].min(),
                'max_plddt': plddt_df['plddt'].max(),
                'std_plddt': plddt_df['plddt'].std(),
                'very_high': conf_counts.get('very_high', 0),
                'confident': conf_counts.get('confident', 0),
                'low': conf_counts.get('low', 0),
                'very_low': conf_counts.get('very_low', 0)
            })
            
            print(f"  Mean pLDDT: {mean_plddt:.2f}")
            print(f"  Very high (≥90): {conf_counts.get('very_high', 0)} residues")
            print(f"  Confident (70-90): {conf_counts.get('confident', 0)} residues")
            print(f"  Low (50-70): {conf_counts.get('low', 0)} residues")
            print(f"  Very low (<50): {conf_counts.get('very_low', 0)} residues")
            
            # Plot pLDDT per residue
            plot_file = os.path.join(plddt_dir, f'{domain_name}_plddt.png')
            plot_plddt_per_residue(
                plddt_df,
                plot_file,
                title=f'{domain_name} - pLDDT per Residue'
            )
            
        except Exception as e:
            print(f"  Error processing {domain_name}: {e}")
            continue
    
    # Combine all data
    if all_plddt_data:
        combined_df = pd.concat(all_plddt_data, ignore_index=True)
        
        # Save detailed pLDDT data
        detail_file = os.path.join(results_dir, 'plddt_per_residue.csv')
        combined_df.to_csv(detail_file, index=False)
        print(f"\nSaved detailed pLDDT data to: {detail_file}")
    
    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(results_dir, 'plddt_summary.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"Saved pLDDT summary to: {summary_file}")
        
        # Plot summary statistics
        plot_summary_statistics(summary_df, results_dir)
        
        return summary_df
    
    return pd.DataFrame()


def plot_summary_statistics(summary_df, results_dir):
    """Plot summary statistics across all structures."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Sort by mean pLDDT for better visualization
    summary_df = summary_df.sort_values('mean_plddt', ascending=False)
    
    # 1. Mean pLDDT distribution
    axes[0, 0].hist(summary_df['mean_plddt'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].axvline(70, color='r', linestyle='--', linewidth=2, label='Confidence threshold')
    axes[0, 0].set_xlabel('Mean pLDDT Score', fontsize=11)
    axes[0, 0].set_ylabel('Number of Structures', fontsize=11)
    axes[0, 0].set_title('Distribution of Mean pLDDT Scores', fontsize=12, fontweight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # 2. pLDDT vs structure size
    axes[0, 1].scatter(summary_df['num_residues'], summary_df['mean_plddt'], 
                      alpha=0.6, s=50, color='steelblue')
    axes[0, 1].axhline(70, color='r', linestyle='--', linewidth=2, alpha=0.5)
    axes[0, 1].set_xlabel('Number of Residues', fontsize=11)
    axes[0, 1].set_ylabel('Mean pLDDT Score', fontsize=11)
    axes[0, 1].set_title('pLDDT vs Structure Size', fontsize=12, fontweight='bold')
    axes[0, 1].grid(alpha=0.3)
    
    # 3. Confidence distribution (stacked bar for top structures)
    top_n = min(20, len(summary_df))
    top_structures = summary_df.head(top_n)
    
    conf_data = top_structures[['very_high', 'confident', 'low', 'very_low']].values
    x_pos = np.arange(len(top_structures))
    
    axes[1, 0].bar(x_pos, conf_data[:, 0], label='Very high (≥90)', color='#0053D6')
    axes[1, 0].bar(x_pos, conf_data[:, 1], bottom=conf_data[:, 0], 
                   label='Confident (70-90)', color='#65CBF3')
    axes[1, 0].bar(x_pos, conf_data[:, 2], 
                   bottom=conf_data[:, 0] + conf_data[:, 1],
                   label='Low (50-70)', color='#FFDB13')
    axes[1, 0].bar(x_pos, conf_data[:, 3], 
                   bottom=conf_data[:, 0] + conf_data[:, 1] + conf_data[:, 2],
                   label='Very low (<50)', color='#FF7D45')
    
    axes[1, 0].set_xlabel(f'Top {top_n} Structures (by mean pLDDT)', fontsize=11)
    axes[1, 0].set_ylabel('Number of Residues', fontsize=11)
    axes[1, 0].set_title('Confidence Distribution', fontsize=12, fontweight='bold')
    axes[1, 0].legend(loc='upper right', fontsize=8)
    axes[1, 0].set_xticks([])
    
    # 4. Percentage of high-confidence residues
    summary_df['pct_confident'] = (summary_df['very_high'] + summary_df['confident']) / summary_df['num_residues'] * 100
    
    axes[1, 1].hist(summary_df['pct_confident'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    axes[1, 1].set_xlabel('% Confident Residues (pLDDT ≥ 70)', fontsize=11)
    axes[1, 1].set_ylabel('Number of Structures', fontsize=11)
    axes[1, 1].set_title('Distribution of Confident Residues', fontsize=12, fontweight='bold')
    axes[1, 1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plot_file = os.path.join(results_dir, 'plddt_summary_plots.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved summary plots to: {plot_file}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract pLDDT scores from structures')
    parser.add_argument('--output_dir', default='alphafold_outputs',
                        help='AlphaFold3 output directory')
    parser.add_argument('--results_dir', default='analysis_results',
                        help='Directory to save analysis results')
    
    args = parser.parse_args()
    
    print("="*70)
    print("PARP14 pLDDT Score Analysis")
    print("="*70)
    print(f"\nInput directory: {args.output_dir}")
    print(f"Results directory: {args.results_dir}\n")
    
    # Process all structures
    summary_df = process_all_structures(args.output_dir, args.results_dir)
    
    print("\n" + "="*70)
    print("pLDDT analysis complete!")
    print("="*70)
    if len(summary_df) > 0:
        print(f"\nProcessed {len(summary_df)} structures")
        print(f"Average mean pLDDT: {summary_df['mean_plddt'].mean():.2f}")
        print(f"High confidence structures (mean pLDDT ≥70): {len(summary_df[summary_df['mean_plddt'] >= 70])}")


if __name__ == '__main__':
    main()
