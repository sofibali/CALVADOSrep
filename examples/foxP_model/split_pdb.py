#!/usr/bin/env python3
"""
Split multi-chain PDB file into separate single-chain PDB files for CALVADOS.
"""

import os
from argparse import ArgumentParser

def split_pdb_by_chain(input_pdb, output_dir, chain_names=None):
    """
    Split a multi-chain PDB file into separate files by chain.
    
    Parameters:
    -----------
    input_pdb : str
        Path to input PDB file with multiple chains
    output_dir : str
        Directory to write output PDB files
    chain_names : dict, optional
        Dictionary mapping chain IDs to output names
        e.g., {'A': 'FOXP4', 'B': 'FOX'}
        If None, uses chain IDs as names
    """
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Read PDB file and separate by chain
    chains = {}
    current_chain = None
    
    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract chain ID (column 22, 0-indexed column 21)
                chain_id = line[21]
                
                if chain_id not in chains:
                    chains[chain_id] = []
                
                chains[chain_id].append(line)
                current_chain = chain_id
            
            elif line.startswith('TER') and current_chain is not None:
                chains[current_chain].append(line)
            
            elif line.startswith('END'):
                # Add END record to all chains
                for chain_id in chains:
                    chains[chain_id].append(line)
    
    # Write separate PDB files for each chain
    output_files = {}
    for chain_id, lines in chains.items():
        # Determine output filename
        if chain_names and chain_id in chain_names:
            output_name = chain_names[chain_id]
        else:
            output_name = f'chain_{chain_id}'
        
        output_path = os.path.join(output_dir, f'{output_name}.pdb')
        
        with open(output_path, 'w') as f:
            f.writelines(lines)
        
        output_files[chain_id] = output_path
        
        # Count atoms
        n_atoms = sum(1 for line in lines if line.startswith('ATOM'))
        print(f"Chain {chain_id} -> {output_path} ({n_atoms} atoms)")
    
    return output_files

def main():
    parser = ArgumentParser(description='Split multi-chain PDB file for CALVADOS')
    parser.add_argument('--input', required=True, help='Input PDB file')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--chain-a-name', default='FOXP4', help='Name for chain A')
    parser.add_argument('--chain-b-name', default='FOX', help='Name for chain B')
    parser.add_argument('--chain-c-name', default='chain_C', help='Name for chain C')
    parser.add_argument('--chain-d-name', default='chain_D', help='Name for chain D')  
    parser.add_argument('--chain-e-name', default='chain_E', help='Name for chain E') 
    args = parser.parse_args()
    
    chain_names = {
        'A': args.chain_a_name,
        'B': args.chain_b_name,
        'C': args.chain_c_name,
        'D': args.chain_d_name,
        'E': args.chain_e_name
    }
    
    print(f"\n{'='*60}")
    print("Splitting PDB file by chain")
    print(f"{'='*60}")
    print(f"Input: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print(f"{'='*60}\n")
    
    output_files = split_pdb_by_chain(args.input, args.output_dir, chain_names)
    
    print(f"\n{'='*60}")
    print("Splitting complete!")
    print(f"{'='*60}\n")

if __name__ == '__main__':
    main()
