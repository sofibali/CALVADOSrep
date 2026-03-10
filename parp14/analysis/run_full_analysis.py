#!/usr/bin/env python3
"""
Master script to run comprehensive structural analysis on all AlphaFold3 predictions.
Runs all analysis modules and generates a comprehensive report.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import pandas as pd


def run_command(command, description):
    """Run a command and report status."""
    print("\n" + "="*70)
    print(f"Running: {description}")
    print("="*70)
    
    result = subprocess.run(command, shell=True)
    
    if result.returncode == 0:
        print(f"✓ {description} completed successfully")
        return True
    else:
        print(f"✗ {description} failed with return code {result.returncode}")
        return False


def generate_master_report(results_dir):
    """Generate a comprehensive master report combining all analyses."""
    
    report_file = os.path.join(results_dir, 'MASTER_ANALYSIS_REPORT.md')
    
    with open(report_file, 'w') as f:
        f.write("# PARP14 AlphaFold3 Structure Analysis Report\n\n")
        f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Overview\n\n")
        f.write("This report summarizes the comprehensive structural analysis of ")
        f.write("PARP14 domain combinations predicted by AlphaFold3.\n\n")
        
        # Structure Quality Metrics
        f.write("## 1. Structure Quality Metrics\n\n")
        
        metrics_file = os.path.join(results_dir, 'structure_metrics.csv')
        if os.path.exists(metrics_file):
            df = pd.read_csv(metrics_file)
            f.write(f"- **Total structures analyzed**: {len(df)}\n")
            f.write(f"- **Average mean pLDDT**: {df['mean_plddt'].mean():.2f}\n")
            f.write(f"- **Structures with pLDDT > 70**: {len(df[df['mean_plddt'] > 70])} ")
            f.write(f"({len(df[df['mean_plddt'] > 70])/len(df)*100:.1f}%)\n\n")
            
            f.write("### Top 10 Highest Quality Structures\n\n")
            f.write("| Domain Combination | Mean pLDDT | Residues |\n")
            f.write("|-------------------|------------|----------|\n")
            for _, row in df.nlargest(10, 'mean_plddt').iterrows():
                f.write(f"| {row['domain_combination']} | {row['mean_plddt']:.2f} | {row['num_residues']} |\n")
            f.write("\n")
        else:
            f.write("*Structure metrics file not found*\n\n")
        
        # pLDDT Analysis
        f.write("## 2. pLDDT Score Distribution\n\n")
        
        plddt_file = os.path.join(results_dir, 'plddt_summary.csv')
        if os.path.exists(plddt_file):
            df = pd.read_csv(plddt_file)
            
            f.write("### Confidence Level Distribution (across all structures)\n\n")
            total_very_high = df['very_high'].sum()
            total_confident = df['confident'].sum()
            total_low = df['low'].sum()
            total_very_low = df['very_low'].sum()
            total_residues = total_very_high + total_confident + total_low + total_very_low
            
            f.write(f"- **Very high confidence (≥90)**: {total_very_high} residues ")
            f.write(f"({total_very_high/total_residues*100:.1f}%)\n")
            f.write(f"- **Confident (70-90)**: {total_confident} residues ")
            f.write(f"({total_confident/total_residues*100:.1f}%)\n")
            f.write(f"- **Low (50-70)**: {total_low} residues ")
            f.write(f"({total_low/total_residues*100:.1f}%)\n")
            f.write(f"- **Very low (<50)**: {total_very_low} residues ")
            f.write(f"({total_very_low/total_residues*100:.1f}%)\n\n")
        else:
            f.write("*pLDDT summary file not found*\n\n")
        
        # Contact Analysis
        f.write("## 3. Contact Map Analysis\n\n")
        
        contact_file = os.path.join(results_dir, 'contact_analysis_summary.csv')
        if os.path.exists(contact_file):
            df = pd.read_csv(contact_file)
            f.write(f"- **Structures analyzed**: {len(df)}\n")
            f.write(f"- **Average contact density**: {df['contact_density'].mean():.4f}\n")
            f.write(f"- **Average contacts per residue**: {df['avg_contact_per_residue'].mean():.2f}\n\n")
            
            f.write("### Structures with Highest Contact Density\n\n")
            f.write("| Domain Combination | Contact Density | Total Contacts |\n")
            f.write("|-------------------|-----------------|----------------|\n")
            for _, row in df.nlargest(10, 'contact_density').iterrows():
                f.write(f"| {row['domain_combination']} | {row['contact_density']:.4f} | {row['total_contacts']:.0f} |\n")
            f.write("\n")
        else:
            f.write("*Contact analysis file not found*\n\n")
        
        # Interface Analysis
        f.write("## 4. Domain Interface Analysis\n\n")
        
        interface_file = os.path.join(results_dir, 'interface_analysis_summary.csv')
        if os.path.exists(interface_file):
            df = pd.read_csv(interface_file)
            f.write(f"- **Total interface pairs analyzed**: {len(df)}\n")
            f.write(f"- **Average interface contacts**: {df['num_interface_contacts'].mean():.2f}\n")
            f.write(f"- **Average interface distance**: {df['avg_distance'].mean():.2f} Å\n\n")
            
            f.write("### Most Common Domain Interfaces\n\n")
            interface_counts = df.groupby(['domain1', 'domain2']).size().sort_values(ascending=False)
            f.write("| Domain 1 | Domain 2 | Frequency |\n")
            f.write("|----------|----------|----------|\n")
            for (dom1, dom2), count in interface_counts.head(10).items():
                f.write(f"| {dom1} | {dom2} | {count} |\n")
            f.write("\n")
        else:
            f.write("*Interface analysis file not found*\n\n")
        
        # Output Files Summary
        f.write("## 5. Generated Output Files\n\n")
        f.write("### Data Files\n")
        f.write("- `structure_metrics.csv` - Overall structure quality metrics\n")
        f.write("- `plddt_summary.csv` - pLDDT score summaries\n")
        f.write("- `plddt_per_residue.csv` - Detailed per-residue pLDDT scores\n")
        f.write("- `contact_analysis_summary.csv` - Contact map statistics\n")
        f.write("- `interface_analysis_summary.csv` - Domain interface statistics\n")
        f.write("- `interface_residues_detailed.csv` - Detailed interface residue pairs\n\n")
        
        f.write("### Visualization Directories\n")
        f.write("- `contact_maps/` - Contact map images\n")
        f.write("- `distance_maps/` - Distance matrix images\n")
        f.write("- `plddt_plots/` - Per-structure pLDDT plots\n")
        f.write("- `interface_analysis/` - Domain interface network plots\n\n")
        
        f.write("### Summary Plots\n")
        f.write("- `plddt_distribution.png` - pLDDT score distributions\n")
        f.write("- `plddt_summary_plots.png` - Comprehensive pLDDT analysis\n\n")
        
        f.write("## Notes\n\n")
        f.write("- **pLDDT**: Per-residue confidence score (0-100). Values >70 indicate high confidence.\n")
        f.write("- **Contact**: Defined as CA-CA distance ≤ 8Å\n")
        f.write("- **Interface**: Defined as CA-CA distance ≤ 5Å between different domains\n")
        f.write("- All distance measurements are in Ångströms (Å)\n")
    
    print(f"\n✓ Generated master report: {report_file}")
    return report_file


def main():
    parser = argparse.ArgumentParser(
        description='Run comprehensive structural analysis on AlphaFold3 predictions'
    )
    parser.add_argument('--output_dir', default='alphafold_outputs',
                        help='AlphaFold3 output directory')
    parser.add_argument('--results_dir', default='analysis_results',
                        help='Directory to save analysis results')
    parser.add_argument('--skip_structures', action='store_true',
                        help='Skip structure quality analysis')
    parser.add_argument('--skip_plddt', action='store_true',
                        help='Skip pLDDT extraction')
    parser.add_argument('--skip_contacts', action='store_true',
                        help='Skip contact map analysis')
    parser.add_argument('--skip_interfaces', action='store_true',
                        help='Skip interface analysis')
    
    args = parser.parse_args()
    
    # Get the analysis directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    print("="*70)
    print("PARP14 COMPREHENSIVE STRUCTURAL ANALYSIS")
    print("="*70)
    print(f"\nInput directory: {args.output_dir}")
    print(f"Results directory: {args.results_dir}")
    print(f"Analysis scripts: {script_dir}\n")
    
    # Create results directory
    os.makedirs(args.results_dir, exist_ok=True)
    
    results = []
    
    # Run structure quality analysis
    if not args.skip_structures:
        cmd = f"python {os.path.join(script_dir, 'analyze_structures.py')} "
        cmd += f"--output_dir {args.output_dir} --results_dir {args.results_dir}"
        results.append(run_command(cmd, "Structure Quality Analysis"))
    
    # Run pLDDT extraction
    if not args.skip_plddt:
        cmd = f"python {os.path.join(script_dir, 'extract_plddt.py')} "
        cmd += f"--output_dir {args.output_dir} --results_dir {args.results_dir}"
        results.append(run_command(cmd, "pLDDT Score Extraction"))
    
    # Run contact analysis
    if not args.skip_contacts:
        cmd = f"python {os.path.join(script_dir, 'analyze_contacts.py')} "
        cmd += f"--output_dir {args.output_dir} --results_dir {args.results_dir}"
        results.append(run_command(cmd, "Contact Map Analysis"))
    
    # Run interface analysis
    if not args.skip_interfaces:
        cmd = f"python {os.path.join(script_dir, 'analyze_interfaces.py')} "
        cmd += f"--output_dir {args.output_dir} --results_dir {args.results_dir}"
        results.append(run_command(cmd, "Domain Interface Analysis"))
    
    # Generate master report
    print("\n" + "="*70)
    print("Generating master report...")
    print("="*70)
    report_file = generate_master_report(args.results_dir)
    
    # Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    
    successful = sum(results)
    total = len(results)
    
    print(f"\nCompleted {successful}/{total} analysis modules")
    print(f"\nResults saved to: {args.results_dir}/")
    print(f"Master report: {report_file}")
    
    if successful < total:
        print("\n⚠ Some analysis modules failed. Check logs above for details.")
        sys.exit(1)
    else:
        print("\n✓ All analyses completed successfully!")
        sys.exit(0)


if __name__ == '__main__':
    main()
