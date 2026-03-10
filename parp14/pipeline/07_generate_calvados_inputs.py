#!/usr/bin/env python3
"""
Generate CALVADOS simulation inputs from AlphaFold3 structures.

Creates directory structure using CALVADOS Config and Components classes.

Usage:
    python 07_generate_calvados_inputs.py \
        --input alphafold_outputs/ \
        --domains Domain_boundries.csv \
        --output calvados_simulations/ \
        --residues-csv /path/to/residues_CALVADOS3.csv
"""

import argparse
import os
import subprocess
from Bio.PDB import MMCIFParser, PDBIO, Select
import csv
from collections import defaultdict
from calvados.cfg import Config, Components
import yaml


class BackboneSelect(Select):
    """Select only backbone atoms (N, CA, C, O) for coarse-grained simulations."""
    def accept_atom(self, atom):
        return atom.get_name() in ['N', 'CA', 'C', 'O']


def read_domains(domains_csv):
    """
    Read domain boundaries from CSV file.
    
    Returns dict: {domain_name: (start, end)}
    """
    domains = {}
    with open(domains_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            domain_name = row['Domain'].strip().lower()
            start = int(row['Start'])
            end = int(row['End'])
            domains[domain_name] = (start, end)
    return domains


def get_structure_domains(structure_name, all_domains):
    """
    Extract domain boundaries for a given structure name.
    
    Returns list of [start, end] pairs in order of appearance.
    """
    # Parse domain names from structure (e.g., "kh1_kh2_kh3_art" -> [kh1, kh2, kh3, art])
    domain_names = structure_name.split('_')
    
    structure_domains = []
    for domain_name in domain_names:
        if domain_name in all_domains:
            start, end = all_domains[domain_name]
            structure_domains.append([start, end])
    
    return structure_domains


def convert_cif_to_pdb(cif_path, pdb_path):
    """
    Convert mmCIF file to PDB format (backbone only for CALVADOS).
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('model', cif_path)
    
    # Save as PDB with only backbone atoms
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path, BackboneSelect())


def create_domains_yaml(structure_name, domains_list, output_file):
    """
    Create domains.yaml file with domain boundaries.
    
    Format:
      structure_name:
        - [start1, end1]
        - [start2, end2]
        ...
    """
    domains_dict = {structure_name: domains_list}
    
    with open(output_file, 'w') as f:
        yaml.dump(domains_dict, f, default_flow_style=False, sort_keys=False)


def main():
    parser = argparse.ArgumentParser(
        description='Generate CALVADOS simulation inputs from AlphaFold3 structures'
    )
    parser.add_argument('--input', required=True, help='AlphaFold3 output directory')
    parser.add_argument('--domains', required=True, help='Domain boundaries CSV file')
    parser.add_argument('--output', required=True, help='Output directory for CALVADOS inputs')
    parser.add_argument('--residues-csv', required=True, 
                       help='Path to residues_CALVADOS3.csv file')
    parser.add_argument('--seed', default='seed-1_sample-0',
                       help='Which seed/sample to use (default: seed-1_sample-0)')
    
    args = parser.parse_args()
    
    # Read domain boundaries
    print(f"Reading domain boundaries from {args.domains}")
    all_domains = read_domains(args.domains)
    print(f"  Loaded {len(all_domains)} domains\n")
    
    # Find all structure directories
    print(f"Scanning AlphaFold3 outputs: {args.input}")
    structure_dirs = sorted([d for d in os.listdir(args.input) 
                            if os.path.isdir(os.path.join(args.input, d))])
    print(f"  Found {len(structure_dirs)} structures\n")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    processed = 0
    skipped = 0
    skip_reasons = defaultdict(int)
    
    for struct_name in structure_dirs:
        struct_dir = os.path.join(args.input, struct_name)
        seed_dir = os.path.join(struct_dir, args.seed)
        cif_path = os.path.join(seed_dir, 'model.cif')
        
        # Check if CIF file exists
        if not os.path.isfile(cif_path):
            skipped += 1
            skip_reasons['no_cif_file'] += 1
            if skip_reasons['no_cif_file'] <= 3:
                print(f"  SKIP: {struct_name} - CIF file not found at {cif_path}")
            continue
        
        try:
            # Get domain boundaries for this structure
            structure_domains = get_structure_domains(struct_name, all_domains)
            
            if not structure_domains:
                skipped += 1
                skip_reasons['no_domains'] += 1
                if skip_reasons['no_domains'] <= 3:
                    print(f"  SKIP: {struct_name} - No domain boundaries found")
                continue
            
            # Create directory structure
            output_struct_dir = os.path.join(args.output, struct_name)
            input_dir = os.path.join(output_struct_dir, 'input')
            subprocess.run(f'mkdir -p {output_struct_dir}', shell=True)
            subprocess.run(f'mkdir -p {input_dir}', shell=True)
            
            # Convert CIF to PDB
            pdb_path = os.path.join(input_dir, f'{struct_name}.pdb')
            convert_cif_to_pdb(cif_path, pdb_path)
            
            # Create domains.yaml first
            domains_yaml = os.path.join(output_struct_dir, 'domains.yaml')
            create_domains_yaml(struct_name, structure_domains, domains_yaml)
            
            # Prepare analysis script
            analyses = f"""
from calvados.analysis import save_conf_prop

save_conf_prop(path="{output_struct_dir}",name="{struct_name}",residues_file="{args.residues_csv}",output_path="{args.output}/data",start=100,is_idr=False,select='all')
"""
            
            # Create config.yaml using Config class
            config = Config(
                sysname = struct_name,
                box = [123, 123, 123],  # nm, based on single_MDP example
                temp = 293,  # K
                ionic = 0.19,  # molar
                pH = 7.0,
                topol = 'center',
                wfreq = 8000,  # dcd writing interval
                steps = 32000000,  # 4000 frames * 8000 steps
                runtime = 0,
                platform = 'CPU',
                threads = 4,
                restart = 'checkpoint',
                frestart = 'restart.chk',
                verbose = True,
            )
            config.write(output_struct_dir, name='config.yaml', analyses=analyses)
            
            # Create components.yaml using Components class
            components = Components(
                molecule_type = 'protein',
                nmol = 1,
                restraint = True,
                charge_termini = 'both',
                fresidues = args.residues_csv,
                fdomains = os.path.join(output_struct_dir, 'domains.yaml'),
                pdb_folder = input_dir,
                restraint_type = 'harmonic',
                use_com = True,
                colabfold = 1,
                k_harmonic = 700.0,
            )
            components.add(name=struct_name)
            components.write(output_struct_dir, name='components.yaml')
            
            processed += 1
            
            if processed % 100 == 0:
                print(f"  Processed {processed} structures...")
        
        except Exception as e:
            skipped += 1
            error_type = type(e).__name__
            skip_reasons[f'error_{error_type}'] += 1
            if skip_reasons[f'error_{error_type}'] <= 3:
                print(f"  ERROR: {struct_name} - {error_type}: {str(e)[:100]}")
            continue
    
    # Summary
    print(f"\n{'='*70}")
    print(f"CALVADOS Input Generation Complete")
    print(f"{'='*70}")
    print(f"Structures processed: {processed}")
    print(f"Structures skipped: {skipped}")
    print(f"\nOutput directory: {args.output}")
    print(f"{'='*70}")
    
    if skip_reasons:
        print(f"\nSkip reasons:")
        for reason, count in sorted(skip_reasons.items(), key=lambda x: x[1], reverse=True):
            print(f"  {reason}: {count}")
    
    print(f"\nEach structure directory contains:")
    print(f"  - config.yaml: Simulation parameters")
    print(f"  - components.yaml: Molecule definitions")
    print(f"  - domains.yaml: Domain boundaries")
    print(f"  - input/{struct_name}.pdb: Structure file")


if __name__ == '__main__':
    main()
