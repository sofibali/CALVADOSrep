#!/usr/bin/env python3
"""
Generate AlphaFold3 JSON input files for all domain combinations.
Creates JSON files optimized for batch processing with AlphaFold3.
"""

import os
import json
from collections import defaultdict

def read_fasta(fasta_file):
    """Read FASTA file and return header and sequence."""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip().lstrip('>')
    sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence

def create_alphafold_json(name, sequence, output_file):
    """Create AlphaFold3 JSON input file."""
    data = {
        "name": name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": sequence
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    return output_file

def main():
    # Configuration
    fasta_dir = "domain_combinations"
    json_dir = "alphafold_inputs"
    
    # Create JSON directory
    os.makedirs(json_dir, exist_ok=True)
    
    # Get all FASTA files
    fasta_files = sorted([f for f in os.listdir(fasta_dir) if f.endswith('.fasta')])
    total_files = len(fasta_files)
    
    print("="*70)
    print("Generate AlphaFold3 JSON Input Files")
    print("="*70)
    print(f"\nTotal FASTA files: {total_files}")
    print(f"Input directory: {fasta_dir}/")
    print(f"Output directory: {json_dir}/")
    print()
    
    # Process each FASTA file
    for idx, fasta_file in enumerate(fasta_files, 1):
        fasta_path = os.path.join(fasta_dir, fasta_file)
        name = os.path.splitext(fasta_file)[0]
        
        # Read FASTA
        header, sequence = read_fasta(fasta_path)
        
        # Create JSON input file
        json_file = os.path.join(json_dir, f"{name}.json")
        create_alphafold_json(name, sequence, json_file)
        
        print(f"[{idx}/{total_files}] {name}: {len(sequence)} aa -> {json_file}")
    
    # Group by sequence length for bucket optimization
    length_groups = defaultdict(list)
    for idx, fasta_file in enumerate(fasta_files):
        fasta_path = os.path.join(fasta_dir, fasta_file)
        _, sequence = read_fasta(fasta_path)
        length_groups[len(sequence)].append(os.path.splitext(fasta_file)[0])
    
    print(f"\n{'='*70}")
    print(f"Created {total_files} JSON files in '{json_dir}/'")
    print(f"\nSequence length distribution:")
    for length in sorted(length_groups.keys()):
        print(f"  {length} aa: {len(length_groups[length])} files")
    print(f"{'='*70}")
    print(f"\nTo run ALL predictions efficiently using AlphaFold3's batch mode:")
    print(f"  CUDA_VISIBLE_DEVICES=0 run_alphafold.py \\")
    print(f"    --db_dir /mnt/alphafold3 \\")
    print(f"    --model_dir /mnt/alphafold3 \\")
    print(f"    --output_dir alphafold_outputs \\")
    print(f"    --input_dir {json_dir} \\")
    print(f"    --jax_compilation_cache_dir ./jax_cache")
    print(f"\nOr use the automated batch script:")
    print(f"  python3 run_alphafold_batch.py")

if __name__ == "__main__":
    main()
