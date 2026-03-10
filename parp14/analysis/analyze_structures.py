#!/usr/bin/env python3
"""
Analyze AlphaFold3 predicted structures for PARP14 domain combinations.
Extracts and analyzes structural quality metrics including pLDDT scores,
PAE (Predicted Aligned Error), and ranking scores.
"""

import os
import json
import glob
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from Bio.PDB import MMCIFParser, PDBIO
import matplotlib.pyplot as plt
import seaborn as sns


def parse_cif_file(cif_file):
    """Parse CIF file and extract structure information."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)
    
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    residues.append({
                        'chain': chain.id,
                        'resname': residue.resname,
                        'resid': residue.id[1],
                        'residue': residue
                    })
    
    return structure, residues


def extract_plddt_scores(cif_file):
    """Extract pLDDT scores from CIF file."""
    structure, residues = parse_cif_file(cif_file)
    
    plddt_scores = []
    for res_info in residues:
        residue = res_info['residue']
        # pLDDT is stored in the B-factor column
        for atom in residue:
            if atom.name == 'CA':  # Use CA atom
                plddt_scores.append({
                    'resid': res_info['resid'],
                    'resname': res_info['resname'],
                    'plddt': atom.bfactor
                })
                break
    
    return pd.DataFrame(plddt_scores)


def analyze_confidence_json(confidence_file):
    """Analyze summary_confidences JSON file."""
    with open(confidence_file, 'r') as f:
        data = json.load(f)
    
    metrics = {
        'ptm': data.get('ptm', None),
        'iptm': data.get('iptm', None),
        'ranking_score': data.get('ranking_score', None),
        'fraction_disordered': data.get('fraction_disordered', None),
    }
    
    # Extract PAE if available
    if 'pae' in data:
        pae_matrix = np.array(data['pae'])
        metrics['mean_pae'] = np.mean(pae_matrix)
        metrics['median_pae'] = np.median(pae_matrix)
    
    return metrics


def analyze_structure_directory(output_dir):
    """Analyze all structures in the AlphaFold output directory."""
    results = []
    
    # Find all prediction directories
    pred_dirs = glob.glob(os.path.join(output_dir, '*'))
    pred_dirs = [d for d in pred_dirs if os.path.isdir(d)]
    
    print(f"Found {len(pred_dirs)} prediction directories")
    
    for pred_dir in pred_dirs:
        domain_name = os.path.basename(pred_dir)
        
        # Find CIF file (typically fold_model_0.cif or similar)
        cif_files = glob.glob(os.path.join(pred_dir, '*.cif'))
        if not cif_files:
            print(f"Warning: No CIF file found in {domain_name}")
            continue
        
        cif_file = cif_files[0]
        
        # Extract pLDDT scores
        try:
            plddt_df = extract_plddt_scores(cif_file)
            mean_plddt = plddt_df['plddt'].mean()
            median_plddt = plddt_df['plddt'].median()
            min_plddt = plddt_df['plddt'].min()
            max_plddt = plddt_df['plddt'].max()
        except Exception as e:
            print(f"Error extracting pLDDT from {domain_name}: {e}")
            continue
        
        # Find confidence JSON
        conf_files = glob.glob(os.path.join(pred_dir, 'summary_confidences_*.json'))
        
        result = {
            'domain_combination': domain_name,
            'num_residues': len(plddt_df),
            'mean_plddt': mean_plddt,
            'median_plddt': median_plddt,
            'min_plddt': min_plddt,
            'max_plddt': max_plddt,
        }
        
        # Add confidence metrics if available
        if conf_files:
            try:
                conf_metrics = analyze_confidence_json(conf_files[0])
                result.update(conf_metrics)
            except Exception as e:
                print(f"Error reading confidence JSON for {domain_name}: {e}")
        
        results.append(result)
    
    return pd.DataFrame(results)


def plot_plddt_distribution(results_df, output_file='plddt_distribution.png'):
    """Plot distribution of pLDDT scores across all predictions."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Mean pLDDT histogram
    axes[0, 0].hist(results_df['mean_plddt'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Mean pLDDT Score')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Mean pLDDT Scores')
    axes[0, 0].axvline(70, color='r', linestyle='--', label='Confidence threshold (70)')
    axes[0, 0].legend()
    
    # pLDDT vs number of residues
    axes[0, 1].scatter(results_df['num_residues'], results_df['mean_plddt'], alpha=0.6)
    axes[0, 1].set_xlabel('Number of Residues')
    axes[0, 1].set_ylabel('Mean pLDDT Score')
    axes[0, 1].set_title('pLDDT vs Domain Size')
    axes[0, 1].axhline(70, color='r', linestyle='--', alpha=0.5)
    
    # PTM scores (if available)
    if 'ptm' in results_df.columns and results_df['ptm'].notna().any():
        axes[1, 0].hist(results_df['ptm'].dropna(), bins=30, edgecolor='black', alpha=0.7)
        axes[1, 0].set_xlabel('PTM Score')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Distribution of PTM Scores')
    
    # Ranking scores (if available)
    if 'ranking_score' in results_df.columns and results_df['ranking_score'].notna().any():
        axes[1, 1].hist(results_df['ranking_score'].dropna(), bins=30, edgecolor='black', alpha=0.7)
        axes[1, 1].set_xlabel('Ranking Score')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].set_title('Distribution of Ranking Scores')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {output_file}")
    plt.close()


def generate_summary_report(results_df, output_file='structure_analysis_summary.txt'):
    """Generate a text summary report."""
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("PARP14 ALPHAFOLD3 STRUCTURE ANALYSIS SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total structures analyzed: {len(results_df)}\n\n")
        
        f.write("pLDDT Score Statistics:\n")
        f.write(f"  Mean: {results_df['mean_plddt'].mean():.2f}\n")
        f.write(f"  Median: {results_df['median_plddt'].median():.2f}\n")
        f.write(f"  Std Dev: {results_df['mean_plddt'].std():.2f}\n")
        f.write(f"  Min: {results_df['mean_plddt'].min():.2f}\n")
        f.write(f"  Max: {results_df['mean_plddt'].max():.2f}\n\n")
        
        # High confidence structures (pLDDT > 70)
        high_conf = results_df[results_df['mean_plddt'] > 70]
        f.write(f"High confidence structures (pLDDT > 70): {len(high_conf)} ({len(high_conf)/len(results_df)*100:.1f}%)\n\n")
        
        # Top 10 structures by pLDDT
        f.write("Top 10 structures by mean pLDDT:\n")
        top_10 = results_df.nlargest(10, 'mean_plddt')
        for idx, row in top_10.iterrows():
            f.write(f"  {row['domain_combination']}: {row['mean_plddt']:.2f}\n")
        
        f.write("\n" + "="*70 + "\n")
    
    print(f"Saved summary report to {output_file}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze AlphaFold3 structure predictions')
    parser.add_argument('--output_dir', default='alphafold_outputs',
                        help='AlphaFold3 output directory')
    parser.add_argument('--results_dir', default='analysis_results',
                        help='Directory to save analysis results')
    
    args = parser.parse_args()
    
    # Create results directory
    os.makedirs(args.results_dir, exist_ok=True)
    
    print("="*70)
    print("PARP14 AlphaFold3 Structure Analysis")
    print("="*70)
    print(f"\nAnalyzing structures in: {args.output_dir}")
    print(f"Saving results to: {args.results_dir}\n")
    
    # Analyze all structures
    results_df = analyze_structure_directory(args.output_dir)
    
    if len(results_df) == 0:
        print("No structures found to analyze!")
        return
    
    # Save results to CSV
    csv_file = os.path.join(args.results_dir, 'structure_metrics.csv')
    results_df.to_csv(csv_file, index=False)
    print(f"\nSaved structure metrics to: {csv_file}")
    
    # Generate plots
    plot_file = os.path.join(args.results_dir, 'plddt_distribution.png')
    plot_plddt_distribution(results_df, plot_file)
    
    # Generate summary report
    report_file = os.path.join(args.results_dir, 'structure_analysis_summary.txt')
    generate_summary_report(results_df, report_file)
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)


if __name__ == '__main__':
    main()
