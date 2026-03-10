#!/usr/bin/env python3
"""
Generate AlphaFold3 JSON input files from FASTA files.

Usage:
    python 02_generate_af3_inputs.py --fasta-dir fasta_files/ --output alphafold_inputs/ --seed 5
"""

import argparse
import os
import json


def read_fasta(fasta_file):
    """Read FASTA file and return header and sequence."""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip()[1:]  # Remove '>'
    sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence


def generate_af3_json(name, sequence, seeds):
    """Generate AlphaFold3 JSON input structure with multiple seeds."""
    return {
        "name": name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": sequence
                }
            }
        ],
        "modelSeeds": seeds,
        "dialect": "alphafold3",
        "version": 1
    }


def main():
    parser = argparse.ArgumentParser(description='Generate AlphaFold3 JSON input files')
    parser.add_argument('--fasta-dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output', default='alphafold_inputs', help='Output directory for JSON files')
    parser.add_argument('--seed', type=int, default=5, help='Number of model seeds (default: 5)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Generate seed list
    seeds = list(range(1, args.seed + 1))
    
    # Find all FASTA files
    fasta_files = sorted([f for f in os.listdir(args.fasta_dir) if f.endswith('.fasta')])
    
    print(f"Found {len(fasta_files)} FASTA files in {args.fasta_dir}")
    print(f"Generating {args.seed} seeds per structure: {seeds}")
    print(f"Output directory: {args.output}\n")
    
    for i, fasta_file in enumerate(fasta_files, 1):
        fasta_path = os.path.join(args.fasta_dir, fasta_file)
        base_name = os.path.splitext(fasta_file)[0]
        
        # Read sequence
        header, sequence = read_fasta(fasta_path)
        
        # Generate single JSON with all seeds
        json_data = generate_af3_json(base_name, sequence, seeds)
        
        # Write JSON file
        json_filename = f"{base_name}.json"
        json_path = os.path.join(args.output, json_filename)
        
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        if i % 100 == 0 or i == len(fasta_files):
            print(f"  Processed {i}/{len(fasta_files)} FASTA files...")
    
    print(f"\n✓ Generated {len(fasta_files)} JSON input files")
    print(f"  (Each with {args.seed} seeds: {seeds})")


if __name__ == '__main__':
    main()
