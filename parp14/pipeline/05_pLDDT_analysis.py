#!/usr/bin/env python3
"""
Comprehensive confidence analysis of AlphaFold3 structures.

Extracts both per-residue (pLDDT from B-factor) and per-structure 
(PTM, iPTM, ranking from JSON) confidence metrics in a single efficient pass.

pLDDT confidence levels:
  - Very High: ≥90
  - Confident: 70-90
  - Low: 50-70
  - Very Low: <50

Output: 
  - Per-residue CSV with pLDDT and confidence level
  - Per-structure CSV with PTM, iPTM, ranking, and pLDDT statistics

Usage:
    python 05_pLDDT_analysis.py --input alphafold_outputs/ --output-residue pLDDT_per_residue.csv --output-structure confidence_scores.csv
"""

import argparse
import os
import json
import pandas as pd
from Bio.PDB import MMCIFParser
import warnings
import time
from datetime import datetime, timedelta
from collections import defaultdict
import sys
warnings.filterwarnings('ignore')


def classify_confidence(plddt):
    """Classify confidence level based on pLDDT score."""
    if plddt >= 90:
        return 'Very High'
    elif plddt >= 70:
        return 'Confident'
    elif plddt >= 50:
        return 'Low'
    else:
        return 'Very Low'


def extract_plddt(structure):
    """
    Extract pLDDT scores from B-factor column.
    
    AlphaFold stores pLDDT (confidence) in the B-factor column.
    Returns dictionary of residue -> pLDDT and confidence level.
    """
    try:
        model = structure[0]
        plddt_dict = {}
        
        for chain in model:
            for residue in chain:
                # Skip heteroatoms
                if residue.id[0] != ' ':
                    continue
                
                res_id = residue.id[1]
                res_name = residue.resname
                
                # Get pLDDT from B-factor (average across atoms in residue)
                plddt_values = []
                for atom in residue:
                    plddt_values.append(atom.bfactor)
                
                if plddt_values:
                    avg_plddt = sum(plddt_values) / len(plddt_values)
                    
                    plddt_dict[res_id] = {
                        'residue_name': res_name,
                        'plddt': avg_plddt,
                        'confidence_level': classify_confidence(avg_plddt)
                    }
        
        return plddt_dict
    
    except Exception as e:
        print(f"    pLDDT extraction failed: {e}")
        return None


def extract_confidence_scores(seed_dir):
    """
    Extract PTM, iPTM, and ranking scores from JSON files.
    
    Returns dict with ptm, iptm, ranking_score, has_clash.
    """
    summary_json = os.path.join(seed_dir, 'summary_confidences.json')
    confidences_json = os.path.join(seed_dir, 'confidences.json')
    
    scores = {
        'ptm': None,
        'iptm': None,
        'ranking_score': None,
        'has_clash': None
    }
    
    # Read summary_confidences.json for PTM/iPTM/ranking
    if os.path.exists(summary_json):
        try:
            with open(summary_json, 'r') as f:
                summary_data = json.load(f)
                scores['ptm'] = summary_data.get('ptm', None)
                scores['iptm'] = summary_data.get('iptm', None)
                scores['ranking_score'] = summary_data.get('ranking_score', None)
        except (json.JSONDecodeError, Exception) as e:
            pass
    
    # Read confidences.json for additional info
    if os.path.exists(confidences_json):
        try:
            with open(confidences_json, 'r') as f:
                conf_data = json.load(f)
                scores['has_clash'] = conf_data.get('has_clash', None)
        except (json.JSONDecodeError, Exception) as e:
            pass
    
    return scores


def main():
    parser = argparse.ArgumentParser(description='Comprehensive confidence analysis for AlphaFold3 structures')
    parser.add_argument('--input', required=True, help='AlphaFold3 output directory')
    parser.add_argument('--output-residue', default='pLDDT_per_residue.csv', help='Output CSV file for per-residue data')
    parser.add_argument('--output-structure', default='confidence_scores.csv', help='Output CSV file for per-structure data')
    
    args = parser.parse_args()
    
    # Initialize parser
    parser_pdb = MMCIFParser(QUIET=True)
    
    # Find all structure directories
    print(f"Scanning AlphaFold3 outputs: {args.input}")
    structure_dirs = sorted([d for d in os.listdir(args.input) 
                            if os.path.isdir(os.path.join(args.input, d))])
    
    # Count total models
    total_models = 0
    for struct_name in structure_dirs:
        struct_dir = os.path.join(args.input, struct_name)
        seed_dirs = [d for d in os.listdir(struct_dir) 
                    if os.path.isdir(os.path.join(struct_dir, d)) and 'seed-' in d]
        total_models += len(seed_dirs)
    
    print(f"  Found {len(structure_dirs)} structures with {total_models} total models\n")
    
    # Store all data
    all_residue_data = []
    all_structure_data = []
    
    processed = 0
    skipped = 0
    start_time = time.time()
    last_report_time = start_time
    skip_reasons = defaultdict(int)
    
    print(f"Starting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Progress updates every 250 models...\n")
    sys.stdout.flush()
    
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
                if skip_reasons['no_cif_file'] <= 5:
                    print(f"  SKIP: {model_name} - CIF file not found")
                    sys.stdout.flush()
                continue
            
            try:
                # Parse structure
                structure = parser_pdb.get_structure('model', cif_path)
                
                # Extract pLDDT from structure
                plddt_data = extract_plddt(structure)
                
                # Extract confidence scores from JSON
                confidence_scores = extract_confidence_scores(seed_dir)
                
                if plddt_data:
                    # Calculate structure-level statistics
                    plddt_values = [data['plddt'] for data in plddt_data.values()]
                    mean_plddt = sum(plddt_values) / len(plddt_values) if plddt_values else None
                    
                    # Count confidence levels
                    confidence_counts = defaultdict(int)
                    for data in plddt_data.values():
                        confidence_counts[data['confidence_level']] += 1
                    
                    total_residues = len(plddt_values)
                    
                    # Store per-residue data
                    for res_num, data in plddt_data.items():
                        all_residue_data.append({
                            'structure': struct_name,
                            'seed': seed_name,
                            'residue_number': res_num,
                            'residue_name': data['residue_name'],
                            'plddt': data['plddt'],
                            'confidence_level': data['confidence_level']
                        })
                    
                    # Store per-structure data
                    all_structure_data.append({
                        'structure': struct_name,
                        'seed': seed_name,
                        'ptm': confidence_scores['ptm'],
                        'iptm': confidence_scores['iptm'],
                        'ranking_score': confidence_scores['ranking_score'],
                        'has_clash': confidence_scores['has_clash'],
                        'mean_plddt': mean_plddt,
                        'num_residues': total_residues,
                        'very_high_conf': confidence_counts.get('Very High', 0),
                        'confident': confidence_counts.get('Confident', 0),
                        'low_conf': confidence_counts.get('Low', 0),
                        'very_low_conf': confidence_counts.get('Very Low', 0),
                        'fraction_very_high': confidence_counts.get('Very High', 0) / total_residues if total_residues > 0 else 0,
                        'fraction_confident': confidence_counts.get('Confident', 0) / total_residues if total_residues > 0 else 0,
                        'fraction_low': confidence_counts.get('Low', 0) / total_residues if total_residues > 0 else 0,
                        'fraction_very_low': confidence_counts.get('Very Low', 0) / total_residues if total_residues > 0 else 0
                    })
                    
                    processed += 1
                else:
                    skipped += 1
                    skip_reasons['analysis_failed'] += 1
                    if skip_reasons['analysis_failed'] <= 5:
                        print(f"  SKIP: {model_name} - pLDDT extraction failed")
                        sys.stdout.flush()
                
                # Progress reporting every 250 models or every 5 minutes
                current_time = time.time()
                if processed % 250 == 0 or (current_time - last_report_time) >= 300:
                    elapsed_time = current_time - start_time
                    avg_time_per_model = elapsed_time / processed if processed > 0 else 0
                    remaining = total_models - model_idx
                    eta_seconds = avg_time_per_model * remaining if avg_time_per_model > 0 else 0
                    eta = timedelta(seconds=int(eta_seconds))
                    
                    print(f"[{datetime.now().strftime('%H:%M:%S')}] Progress: {model_idx}/{total_models} models")
                    print(f"  Processed: {processed} | Skipped: {skipped}")
                    print(f"  Remaining: {remaining} models")
                    if processed > 0:
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
    
    # Save to CSV
    total_time = time.time() - start_time
    
    if all_residue_data and all_structure_data:
        df_residue = pd.DataFrame(all_residue_data)
        df_structure = pd.DataFrame(all_structure_data)
        
        df_residue.to_csv(args.output_residue, index=False)
        df_structure.to_csv(args.output_structure, index=False)
        
        print(f"\n{'='*70}")
        print(f"Confidence Analysis Complete")
        print(f"{'='*70}")
        print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total time: {timedelta(seconds=int(total_time))}")
        print(f"Structures: {len(structure_dirs)}")
        print(f"Models processed: {processed}")
        print(f"Models skipped: {skipped}")
        print(f"Total residues analyzed: {len(all_residue_data)}")
        if processed > 0:
            print(f"Average rate: {processed/total_time:.2f} models/second")
        
        if skip_reasons:
            print(f"\nSkip reasons:")
            for reason, count in sorted(skip_reasons.items(), key=lambda x: x[1], reverse=True):
                print(f"  {reason}: {count}")
        
        print(f"\nOutput files:")
        print(f"  Per-residue: {args.output_residue}")
        print(f"  Per-structure: {args.output_structure}")
        print(f"{'='*70}")
        
        # Summary statistics
        print(f"\nConfidence level distribution (per-residue):")
        conf_counts = df_residue['confidence_level'].value_counts()
        total = len(df_residue)
        for conf_level in ['Very High', 'Confident', 'Low', 'Very Low']:
            if conf_level in conf_counts.index:
                count = conf_counts[conf_level]
                pct = 100 * count / total
                print(f"  {conf_level}: {count:,} residues ({pct:.1f}%)")
        
        print(f"\npLDDT statistics (per-residue):")
        print(f"  Mean: {df_residue['plddt'].mean():.2f}")
        print(f"  Median: {df_residue['plddt'].median():.2f}")
        print(f"  Min: {df_residue['plddt'].min():.2f}")
        print(f"  Max: {df_residue['plddt'].max():.2f}")
        
        # Structure-level statistics
        print(f"\nStructure-level confidence scores:")
        if df_structure['ptm'].notna().any():
            print(f"  PTM: {df_structure['ptm'].mean():.3f} ± {df_structure['ptm'].std():.3f}")
        if df_structure['iptm'].notna().any():
            print(f"  iPTM: {df_structure['iptm'].mean():.3f} ± {df_structure['iptm'].std():.3f}")
        if df_structure['ranking_score'].notna().any():
            print(f"  Ranking: {df_structure['ranking_score'].mean():.3f} ± {df_structure['ranking_score'].std():.3f}")
        
        print(f"\nMean structure-level pLDDT:")
        print(f"  Overall: {df_structure['mean_plddt'].mean():.2f} ± {df_structure['mean_plddt'].std():.2f}")
        
        print(f"\nAverage confidence distribution across structures:")
        print(f"  Very High (≥90): {df_structure['fraction_very_high'].mean()*100:.1f}%")
        print(f"  Confident (70-90): {df_structure['fraction_confident'].mean()*100:.1f}%")
        print(f"  Low (50-70): {df_structure['fraction_low'].mean()*100:.1f}%")
        print(f"  Very Low (<50): {df_structure['fraction_very_low'].mean()*100:.1f}%")
        
    else:
        print(f"\n{'='*70}")
        print(f"Confidence Analysis Complete")
        print(f"{'='*70}")
        print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total time: {timedelta(seconds=int(total_time))}")
        print(f"Models processed: {processed}")
        print(f"Models skipped: {skipped}")
        
        if skip_reasons:
            print(f"\nSkip reasons:")
            for reason, count in sorted(skip_reasons.items(), key=lambda x: x[1], reverse=True):
                print(f"  {reason}: {count}")
        
        print(f"\n⚠️  WARNING: No confidence data collected")
        print(f"{'='*70}")


if __name__ == '__main__':
    main()
