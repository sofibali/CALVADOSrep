#!/usr/bin/env python3
"""
Solvent Accessible Surface Area (SASA) analysis per residue.

Output: Per-residue SASA values in CSV format.

Usage:
    python 06_sasa_analysis.py --input alphafold_outputs/ --output sasa_per_residue.csv
"""

import argparse
import os
import pandas as pd
from Bio.PDB import MMCIFParser
import freesasa
import time
from datetime import datetime, timedelta
from collections import defaultdict
import sys


def calculate_sasa(structure):
    """
    Calculate SASA for all residues using freesasa.
    
    Returns dictionary of residue_number -> SASA value.
    """
    try:
        # Convert BioPython structure to freesasa structure
        struct_sasa = freesasa.structureFromBioPDB(structure)
        
        # Calculate SASA
        result = freesasa.calc(struct_sasa)
        
        # Extract per-residue SASA
        sasa_dict = {}
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_num = residue.get_id()[1]
                    res_name = residue.get_resname()
                    
                    try:
                        # Get SASA for this residue
                        selections = freesasa.selectArea(
                            [f'res, resi {res_num}'],
                            struct_sasa, result
                        )
                        sasa_dict[res_num] = {
                            'residue_name': res_name,
                            'sasa': selections['res']
                        }
                    except:
                        pass
        
        return sasa_dict
    
    except Exception as e:
        print(f"    SASA calculation failed: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description='Calculate SASA per residue for AlphaFold3 structures')
    parser.add_argument('--input', required=True, help='AlphaFold3 output directory')
    parser.add_argument('--output', default='sasa_per_residue.csv', help='Output CSV file')
    parser.add_argument('--domains', help='Optional: Domain boundaries CSV for domain annotation')
    
    args = parser.parse_args()
    
    # Initialize parser
    parser_pdb = MMCIFParser(QUIET=True)
    
    # Find all structure directories and seed subdirectories
    print(f"Scanning AlphaFold3 outputs: {args.input}")
    structure_dirs = sorted([d for d in os.listdir(args.input) 
                            if os.path.isdir(os.path.join(args.input, d))])
    
    # Count total models (structures × seeds)
    total_models = 0
    model_list = []
    for struct_name in structure_dirs:
        struct_dir = os.path.join(args.input, struct_name)
        seed_dirs = sorted([d for d in os.listdir(struct_dir)
                           if os.path.isdir(os.path.join(struct_dir, d)) and d.startswith('seed-')])
        if seed_dirs:
            for seed_name in seed_dirs:
                model_list.append((struct_name, seed_name))
                total_models += 1
    
    print(f"  Found {len(structure_dirs)} structures with {total_models} total models\n")
    
    print(f"Starting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Progress updates every 50 models...\n")
    sys.stdout.flush()
    
    # Optional: Read domain boundaries for annotation
    domain_annotation = None
    if args.domains:
        import csv
        domain_annotation = {}
        with open(args.domains, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                domain_name = row['Domain'].strip()
                start = int(row['Start'])
                end = int(row['End'])
                for res_num in range(start, end + 1):
                    domain_annotation[res_num] = domain_name
    
    # Store all SASA data
    all_sasa_data = []
    
    processed = 0
    skipped = 0
    skip_reasons = defaultdict(int)
    skip_examples = defaultdict(list)
    
    # Progress tracking
    start_time = time.time()
    last_report_time = start_time
    
    # Process each model
    for model_idx, (struct_name, seed_name) in enumerate(model_list, 1):
        struct_dir = os.path.join(args.input, struct_name)
        seed_dir = os.path.join(struct_dir, seed_name)
        cif_path = os.path.join(seed_dir, 'model.cif')
        model_name = f"{struct_name}/{seed_name}"
        
        if not os.path.isfile(cif_path):
            skipped += 1
            skip_reasons['no_cif_file'] += 1
            if len(skip_examples['no_cif_file']) < 5:
                skip_examples['no_cif_file'].append(model_name)
                print(f"  SKIP: {model_name} - CIF file not found")
                sys.stdout.flush()
            continue
        
        try:
            # Parse structure
            structure = parser_pdb.get_structure('model', cif_path)
            
            # Calculate SASA
            sasa_data = calculate_sasa(structure)
            
            if sasa_data:
                # Store data for each residue
                for res_num, data in sasa_data.items():
                    row = {
                        'structure': struct_name,
                        'seed': seed_name,
                        'residue_number': res_num,
                        'residue_name': data['residue_name'],
                        'sasa': data['sasa']
                    }
                    
                    # Add domain annotation if available
                    if domain_annotation and res_num in domain_annotation:
                        row['domain'] = domain_annotation[res_num]
                    
                    all_sasa_data.append(row)
                
                processed += 1
            else:
                skipped += 1
                skip_reasons['sasa_failed'] += 1
                if len(skip_examples['sasa_failed']) < 5:
                    skip_examples['sasa_failed'].append(model_name)
        
        except Exception as e:
            skipped += 1
            error_type = type(e).__name__
            skip_reasons[f'error_{error_type}'] += 1
            if len(skip_examples[f'error_{error_type}']) < 5:
                skip_examples[f'error_{error_type}'].append(f"{model_name}: {str(e)[:50]}")
                print(f"  ERROR: {model_name} - {error_type}: {str(e)[:100]}")
                sys.stdout.flush()
            continue
        
        # Progress reporting every 50 models or every 5 minutes
        current_time = time.time()
        if model_idx % 50 == 0 or (current_time - last_report_time) >= 300:
            elapsed_time = current_time - start_time
            avg_time_per_model = elapsed_time / model_idx
            remaining = total_models - model_idx
            eta_seconds = avg_time_per_model * remaining
            eta = timedelta(seconds=int(eta_seconds))
            
            print(f"[{datetime.now().strftime('%H:%M:%S')}] Progress: {model_idx}/{total_models} models")
            print(f"  Processed: {processed} | Skipped: {skipped}")
            print(f"  Remaining: {remaining} models")
            print(f"  Rate: {model_idx/elapsed_time:.2f} models/sec ({avg_time_per_model:.2f} sec/model)")
            print(f"  ETA: {eta} (estimated completion: {(datetime.now() + timedelta(seconds=eta_seconds)).strftime('%H:%M:%S')})")
            print()
            sys.stdout.flush()
            last_report_time = current_time
    
    # Save to CSV
    total_time = time.time() - start_time
    
    if all_sasa_data:
        df_sasa = pd.DataFrame(all_sasa_data)
        df_sasa.to_csv(args.output, index=False)
        
        print(f"\n{'='*70}")
        print(f"SASA Analysis Complete")
        print(f"{'='*70}")
        print(f"Total models: {total_models}")
        print(f"Models processed: {processed}")
        print(f"Models skipped: {skipped}")
        print(f"Total residues analyzed: {len(all_sasa_data)}")
        print(f"Total time: {timedelta(seconds=int(total_time))}")
        print(f"\nOutput file: {args.output}")
        print(f"{'='*70}")
        
        # Skip reasons summary
        if skip_reasons:
            print(f"\nSkip reasons:")
            for reason, count in sorted(skip_reasons.items(), key=lambda x: x[1], reverse=True):
                print(f"  {reason}: {count}")
                if skip_examples[reason]:
                    print(f"    Examples: {', '.join(skip_examples[reason][:3])}")
        
        # Summary statistics
        print(f"\nSASA statistics:")
        print(f"  Mean SASA: {df_sasa['sasa'].mean():.2f} Ų")
        print(f"  Median SASA: {df_sasa['sasa'].median():.2f} Ų")
        print(f"  Min SASA: {df_sasa['sasa'].min():.2f} Ų")
        print(f"  Max SASA: {df_sasa['sasa'].max():.2f} Ų")
        
        if 'domain' in df_sasa.columns:
            print(f"\nMean SASA by domain:")
            domain_means = df_sasa.groupby('domain')['sasa'].mean().sort_values(ascending=False)
            for domain, mean_sasa in domain_means.items():
                print(f"  {domain}: {mean_sasa:.2f} Ų")
    else:
        print(f"\n✗ No SASA data collected")
        print(f"Total time: {timedelta(seconds=int(total_time))}")


if __name__ == '__main__':
    main()
