#!/usr/bin/env python3
"""
Generate FASTA files for domain combinations.

Usage:
    python 01_generate_fasta.py --fasta PARP14.fasta --domains Domain_boundries.csv --output fasta_files/
"""

import argparse
import os
from itertools import combinations
import csv


def read_fasta(fasta_file):
    """Read FASTA file and return sequence."""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence


def read_domains(domain_file):
    """Read domain boundaries from CSV file."""
    domains = {}
    domain_order = []
    
    with open(domain_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            domain_name = row['Domain'].strip()
            start = int(row['Start'])
            end = int(row['End'])
            domains[domain_name] = (start, end)
            domain_order.append(domain_name)
    
    return domains, domain_order


def generate_domain_combinations(domains, domain_order, max_size=4):
    """Generate all sequential domain combinations up to max_size."""
    all_combinations = []
    
    for size in range(1, max_size + 1):
        # Generate all possible sequential combinations
        for i in range(len(domain_order) - size + 1):
            combo = tuple(domain_order[i:i + size])
            all_combinations.append(combo)
    
    return all_combinations


def extract_domain_sequence(sequence, domains, domain_combo):
    """Extract concatenated sequence for a domain combination."""
    combined_seq = ''
    
    for domain_name in domain_combo:
        if domain_name in domains:
            start, end = domains[domain_name]
            # Convert to 0-based indexing for Python slicing
            domain_seq = sequence[start-1:end]
            combined_seq += domain_seq
    
    return combined_seq


def main():
    parser = argparse.ArgumentParser(description='Generate FASTA files for domain combinations')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--domains', required=True, help='Domain boundaries CSV file')
    parser.add_argument('--output', default='fasta_files', help='Output directory')
    parser.add_argument('--max-size', type=int, default=4, help='Maximum domain combination size')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Read input files
    print(f"Reading FASTA file: {args.fasta}")
    header, sequence = read_fasta(args.fasta)
    print(f"  Sequence length: {len(sequence)} aa")
    
    print(f"\nReading domain boundaries: {args.domains}")
    domains, domain_order = read_domains(args.domains)
    print(f"  Domains found: {len(domains)}")
    for domain in domain_order:
        start, end = domains[domain]
        print(f"    {domain}: {start}-{end}")
    
    # Generate combinations
    print(f"\nGenerating domain combinations (max size: {args.max_size})...")
    combinations_list = generate_domain_combinations(domains, domain_order, args.max_size)
    print(f"  Total combinations: {len(combinations_list)}")
    
    # Generate FASTA files
    print(f"\nGenerating FASTA files in: {args.output}")
    
    for combo in combinations_list:
        # Create filename from domain names
        combo_name = '_'.join(combo)
        fasta_path = os.path.join(args.output, f"{combo_name}.fasta")
        
        # Extract sequence
        combo_seq = extract_domain_sequence(sequence, domains, combo)
        
        # Write FASTA file
        with open(fasta_path, 'w') as f:
            f.write(f">PARP14_{combo_name}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(combo_seq), 80):
                f.write(combo_seq[i:i+80] + '\n')
    
    print(f"✓ Generated {len(combinations_list)} FASTA files")
    
    # Summary by size
    from collections import Counter
    size_counts = Counter(len(combo) for combo in combinations_list)
    print(f"\nCombination size breakdown:")
    for size in sorted(size_counts.keys()):
        print(f"  {size} domains: {size_counts[size]} combinations")


if __name__ == '__main__':
    main()
