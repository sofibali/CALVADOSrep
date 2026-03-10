#!/usr/bin/env python3
"""
Generate FASTA files for all domain combinations of PARP14.
Reads domain boundaries from Domain_boundries.csv and the full sequence from PARP14.fasta.
Treats KH1-6 as a single group and KH7a-KHb-KH8 as another group.
"""

import csv
import os
from itertools import combinations

def read_fasta(fasta_file):
    """Read FASTA file and return the sequence."""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence

def read_domains(domain_file):
    """Read domain boundaries from CSV file, maintaining order."""
    domains = {}
    domain_order = []
    with open(domain_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header line and parse domain lines
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            domain_name = parts[0].strip()
            start = int(parts[1].strip())
            end = int(parts[2].strip())
            domains[domain_name] = (start, end)
            domain_order.append(domain_name)
    return domains, domain_order

def create_domain_groups(domain_order):
    """
    Create domain groups treating KH1-6 and KH7a-KHb-KH8 as units.
    MD1, MD2, MD3 remain as independent groups.
    Returns list of groups where each group is a list of domain names.
    """
    groups = []
    i = 0
    while i < len(domain_order):
        domain = domain_order[i]
        
        # Check if this is KH1 (start of KH1-6 group)
        if domain == 'KH1':
            kh_group = []
            # Collect KH1 through KH6
            while i < len(domain_order) and domain_order[i] in ['KH1', 'KH2', 'KH3', 'KH4', 'KH5', 'KH6']:
                kh_group.append(domain_order[i])
                i += 1
            groups.append(kh_group)
        # Check if this is KH7a (start of KH7a-KHb-KH8 group, but not MD domains)
        elif domain == 'KH7a':
            # Only KH7a goes in this group, MD domains will be separate
            groups.append([domain])
            i += 1
        # Check if this is KHb (continue the KH7a group with KHb and KH8)
        elif domain == 'KHb':
            kh7_group = []
            # Collect KHb and KH8
            while i < len(domain_order) and domain_order[i] in ['KHb', 'KH8']:
                kh7_group.append(domain_order[i])
                i += 1
            groups.append(kh7_group)
        else:
            # Single domain as its own group (including MD1, MD2, MD3)
            groups.append([domain])
            i += 1
    
    return groups

def generate_sequential_combinations(groups, combo_size):
    """
    Generate combinations of domain groups respecting their sequential order.
    Enforces constraint that KH7a and KHb+KH8 must always appear together.
    """
    result = []
    n = len(groups)
    
    def backtrack(start, current_combo):
        if len(current_combo) == combo_size:
            # Validate KH7a and KHb constraint
            if is_valid_kh7_constraint(groups, current_combo):
                result.append(tuple(current_combo))
            return
        
        for i in range(start, n):
            current_combo.append(i)  # Store group index
            backtrack(i + 1, current_combo)
            current_combo.pop()
    
    backtrack(0, [])
    return result

def is_valid_kh7_constraint(groups, group_indices):
    """
    Check if KH7a and KHb+KH8 appear together.
    If one is present, the other must also be present.
    """
    has_kh7a = False
    has_khb_kh8 = False
    
    for idx in group_indices:
        group = groups[idx]
        if 'KH7a' in group:
            has_kh7a = True
        if 'KHb' in group or 'KH8' in group:
            has_khb_kh8 = True
    
    # Both present or both absent = valid
    return has_kh7a == has_khb_kh8

def extract_domain_sequence(sequence, start, end):
    """Extract domain sequence (1-indexed, inclusive)."""
    # Convert from 1-indexed to 0-indexed for Python slicing
    return sequence[start-1:end]

def write_fasta(filename, header_prefix, sequence):
    """Write FASTA file."""
    with open(filename, 'w') as f:
        f.write(f">{header_prefix}\n")
        # Write sequence in lines of 70 characters
        for i in range(0, len(sequence), 70):
            f.write(sequence[i:i+70] + '\n')

def main():
    # File paths
    fasta_file = 'PARP14.fasta'
    domain_file = 'Domain_boundries.csv'
    output_dir = 'domain_combinations'
    
    # Maximum number of domain groups in a combination
    max_combo_size = None
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read input files
    header, sequence = read_fasta(fasta_file)
    domains, domain_order = read_domains(domain_file)
    
    # Create domain groups
    groups = create_domain_groups(domain_order)
    
    print(f"Total sequence length: {len(sequence)}")
    print(f"Number of individual domains: {len(domains)}")
    print(f"\nDomain groups created:")
    for i, group in enumerate(groups):
        group_name = '+'.join(group)
        print(f"  Group {i+1}: {group_name}")
    
    print(f"\nTotal number of groups: {len(groups)}")
    print(f"\nConstraints:")
    print(f"  - KH1-6 are treated as a single unit (all or nothing)")
    print(f"  - KH7a and KHb+KH8 must always appear together")
    print(f"  - MD1, MD2, MD3 are independent and can appear individually")
    print(f"  - Domain order is preserved (sequential combinations only)")
    print(f"\nGenerating domain group combinations...")
    
    # Set maximum combination size
    if max_combo_size is None:
        max_size = len(groups)
    else:
        max_size = min(max_combo_size, len(groups))
    
    total_files = 0
    
    # Generate files for all combination sizes
    for combo_size in range(1, max_size + 1):
        combinations_for_size = generate_sequential_combinations(groups, combo_size)
        print(f"\nGenerating combinations of {combo_size} group(s): {len(combinations_for_size)} files")
        
        for group_indices in combinations_for_size:
            # Build combined sequence from selected groups
            combined_sequence = ""
            position_info = []
            domain_names = []
            
            for group_idx in group_indices:
                group = groups[group_idx]
                for domain in group:
                    start, end = domains[domain]
                    seq = extract_domain_sequence(sequence, start, end)
                    combined_sequence += seq
                    position_info.append(f"{domain}({start}-{end})")
                    domain_names.append(domain)
            
            # Create filename - use underscore to separate domains
            filename = os.path.join(output_dir, "_".join(domain_names) + ".fasta")
            
            # Create header with position information
            header_info = "_".join(position_info)
            
            # Write file
            write_fasta(filename, header_info, combined_sequence)
            total_files += 1
        
        print(f"  Created {len(combinations_for_size)} files")
    
    print(f"\n{'='*60}")
    print(f"Total domain combination files generated: {total_files}")
    print(f"Location: '{output_dir}/' directory")
    print(f"Max combination size: {max_size} groups")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
