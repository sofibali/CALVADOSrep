#!/usr/bin/env python
"""
Prepare CALVADOS simulation for FOXP complex as a SINGLE merged protein.

This approach treats all 3 chains as one continuous protein with:
- All domain restraints preserved
- Chain connectivity maintained through the merged sequence
- No need for inter-chain restraints (handled by bonded interactions)

This script:
1. Merges sequences from all chains into one
2. Creates merged PDB with correct coordinates
3. Creates merged domains file
4. Generates config.yaml and components.yaml for simulation
"""

import os
import sys
import numpy as np
import yaml
import subprocess

# Add CALVADOS to path
calvados_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, calvados_root)

from calvados.cfg import Config, Components

# ============================================================================
# CONFIGURATION
# ============================================================================

cwd = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(cwd)
input_dir = f'{parent_dir}/input'

# System name
sysname = 'foxp_single_complex'

# Chain order for merging (this defines the order in the merged sequence)
# Note: We'll add short linkers between chains to allow some flexibility
chain_order = ['FOXP4', 'FOX', 'chain_C']

# Chain mapping: CALVADOS name -> PDB chain ID
chain_mapping = {
    'FOXP4': 'A',
    'FOX': 'B',
    'chain_C': 'C'
}

# Original domain definitions (from domains.yaml)
original_domains = {
    'FOXP4': [(118, 199), (299, 375), (451, 551)],
    'FOX': [(118, 299), (319, 375), (450, 539)],
    'chain_C': [(232, 321)]
}

# Linker between chains (flexible glycine-serine linker)
LINKER_SEQUENCE = "GSGSGSGSGS"  # 10 residue flexible linker
LINKER_LENGTH = len(LINKER_SEQUENCE)

# Restraint parameters
INTRA_DOMAIN_K = 2000.0   # kJ/mol/nm^2 - force constant for domain restraints

# Simulation parameters
BOX_SIZE = 80  # nm - larger box for single long protein
N_SAVE = 5000  # frames between saves
N_FRAMES = 1000  # total frames to save

# ============================================================================
# FUNCTIONS
# ============================================================================

def read_fasta_sequences(fasta_file):
    """Read sequences from FASTA file (without BioPython)."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Get ID (first word after >)
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def parse_pdb_ca_atoms(pdb_file, chain_id=None):
    """Parse CA atom coordinates from PDB file."""
    ca_atoms = {}  # {resnum: (x, y, z, resname)}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                chain = line[21]
                if chain_id and chain != chain_id:
                    continue
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_atoms[resnum] = (x, y, z, resname)

    return ca_atoms


def calculate_domain_offset(chain_order, sequences, linker_length):
    """
    Calculate the residue offset for each chain in the merged sequence.
    Returns dict: {chain_name: offset}
    """
    offsets = {}
    current_offset = 0

    for chain in chain_order:
        offsets[chain] = current_offset
        chain_length = len(sequences[chain])
        current_offset += chain_length + linker_length

    return offsets


def create_merged_sequence(sequences, chain_order, linker_seq):
    """Create merged sequence from all chains with linkers."""
    merged = ""
    for i, chain in enumerate(chain_order):
        merged += sequences[chain]
        if i < len(chain_order) - 1:  # Don't add linker after last chain
            merged += linker_seq
    return merged


def create_merged_domains(original_domains, chain_order, offsets):
    """
    Create merged domain definitions with correct offsets.
    Returns list of (start, end) tuples for the merged protein.
    """
    merged_domains = []

    for chain in chain_order:
        offset = offsets[chain]
        for start, end in original_domains[chain]:
            # Adjust for 1-based indexing in original domains
            merged_start = start + offset
            merged_end = end + offset
            merged_domains.append((merged_start, merged_end))

    return merged_domains


def create_merged_pdb(complex_pdb, chain_mapping, chain_order, sequences,
                      linker_seq, output_pdb):
    """
    Create a merged PDB file with coordinates from the complex.
    Linker residues are placed at interpolated positions.
    """
    # Parse all CA atoms from complex
    all_ca_atoms = {}
    for name, pdb_chain in chain_mapping.items():
        all_ca_atoms[name] = parse_pdb_ca_atoms(complex_pdb, pdb_chain)

    # Calculate offsets
    offsets = calculate_domain_offset(chain_order, sequences, len(linker_seq))

    with open(output_pdb, 'w') as f:
        atom_num = 1
        merged_resnum = 1

        for chain_idx, chain in enumerate(chain_order):
            ca_atoms = all_ca_atoms[chain]
            seq = sequences[chain]

            # Handle chain_C special case (starts at residue 232)
            if chain == 'chain_C':
                start_resnum = 232
            else:
                start_resnum = 1

            # Write CA atoms for this chain
            for i, aa in enumerate(seq):
                orig_resnum = start_resnum + i
                if orig_resnum in ca_atoms:
                    x, y, z, resname = ca_atoms[orig_resnum]
                else:
                    # Residue not in PDB (disordered region) - interpolate
                    # Find nearest residues with coordinates
                    x, y, z = interpolate_position(ca_atoms, orig_resnum, start_resnum)
                    resname = aa_one_to_three(aa)

                f.write(f"ATOM  {atom_num:5d}  CA  {resname:3s} A{merged_resnum:4d}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                atom_num += 1
                merged_resnum += 1

            # Add linker between chains (except after last chain)
            if chain_idx < len(chain_order) - 1:
                # Get end position of current chain and start of next
                next_chain = chain_order[chain_idx + 1]
                next_ca_atoms = all_ca_atoms[next_chain]

                # Get last coord of current chain
                last_resnum = max(ca_atoms.keys())
                if last_resnum in ca_atoms:
                    end_coord = np.array(ca_atoms[last_resnum][:3])
                else:
                    end_coord = np.array([0, 0, 0])

                # Get first coord of next chain
                if next_chain == 'chain_C':
                    first_resnum = 232
                else:
                    first_resnum = 1

                # Find first available residue in next chain
                next_keys = sorted(next_ca_atoms.keys())
                if next_keys:
                    start_coord = np.array(next_ca_atoms[next_keys[0]][:3])
                else:
                    start_coord = end_coord

                # Interpolate linker positions
                for j, aa in enumerate(linker_seq):
                    t = (j + 1) / (len(linker_seq) + 1)
                    coord = end_coord + t * (start_coord - end_coord)
                    resname = aa_one_to_three(aa)
                    f.write(f"ATOM  {atom_num:5d}  CA  {resname:3s} A{merged_resnum:4d}    "
                           f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00           C\n")
                    atom_num += 1
                    merged_resnum += 1

        f.write("TER\n")
        f.write("END\n")

    print(f"Created merged PDB with {atom_num-1} CA atoms: {output_pdb}")
    return merged_resnum - 1  # Total residues


def interpolate_position(ca_atoms, target_resnum, start_resnum):
    """Interpolate position for missing residue."""
    sorted_resnums = sorted(ca_atoms.keys())

    # Find flanking residues
    prev_res = None
    next_res = None

    for r in sorted_resnums:
        if r < target_resnum:
            prev_res = r
        elif r > target_resnum and next_res is None:
            next_res = r
            break

    if prev_res and next_res:
        # Interpolate between flanking residues
        x1, y1, z1, _ = ca_atoms[prev_res]
        x2, y2, z2, _ = ca_atoms[next_res]
        t = (target_resnum - prev_res) / (next_res - prev_res)
        return (x1 + t*(x2-x1), y1 + t*(y2-y1), z1 + t*(z2-z1))
    elif prev_res:
        x, y, z, _ = ca_atoms[prev_res]
        return (x, y, z + 0.38)  # Extend by ~bond length
    elif next_res:
        x, y, z, _ = ca_atoms[next_res]
        return (x, y, z - 0.38)
    else:
        return (0, 0, 0)


def aa_one_to_three(one_letter):
    """Convert one-letter amino acid code to three-letter."""
    mapping = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    return mapping.get(one_letter.upper(), 'UNK')


def write_merged_domains_yaml(merged_domains, merged_name, output_file):
    """Write merged domains to YAML file."""
    domains_dict = {merged_name: [[s, e] for s, e in merged_domains]}

    with open(output_file, 'w') as f:
        yaml.dump(domains_dict, f, default_flow_style=False)

    print(f"Wrote {len(merged_domains)} domains to {output_file}")


def write_merged_fasta(merged_sequence, merged_name, output_file):
    """Write merged sequence to FASTA file (without BioPython)."""
    with open(output_file, 'w') as f:
        f.write(f">{merged_name}\n")
        # Write sequence in lines of 60 characters
        for i in range(0, len(merged_sequence), 60):
            f.write(merged_sequence[i:i+60] + "\n")
    print(f"Wrote merged sequence ({len(merged_sequence)} residues) to {output_file}")


# ============================================================================
# MAIN SCRIPT
# ============================================================================

if __name__ == '__main__':

    # Paths
    complex_pdb = f'{input_dir}/COMPLEX_WITHDNA.pdb'
    fasta_file = f'{input_dir}/proteins.fasta'
    output_dir = cwd
    sim_path = f'{output_dir}/{sysname}'

    # Create directories
    subprocess.run(f'mkdir -p {sim_path}', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/input', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/data', shell=True)

    # Copy residues file
    subprocess.run(f'cp {input_dir}/residues_CALVADOS3.csv {output_dir}/input/', shell=True)

    print("="*60)
    print("Creating merged FOXP complex...")
    print("="*60)

    # Read original sequences
    sequences = read_fasta_sequences(fasta_file)

    # Calculate offsets for domain mapping
    offsets = calculate_domain_offset(chain_order, sequences, LINKER_LENGTH)

    print(f"\nChain lengths and offsets:")
    for chain in chain_order:
        print(f"  {chain}: {len(sequences[chain])} residues, offset = {offsets[chain]}")

    # Create merged sequence
    merged_name = "FOXP_complex"
    merged_sequence = create_merged_sequence(sequences, chain_order, LINKER_SEQUENCE)
    print(f"\nMerged sequence length: {len(merged_sequence)} residues")
    print(f"  (includes {len(chain_order)-1} linkers of {LINKER_LENGTH} residues each)")

    # Create merged domains
    merged_domains = create_merged_domains(original_domains, chain_order, offsets)

    print(f"\nMerged domains (with offsets applied):")
    domain_idx = 0
    for chain in chain_order:
        print(f"  From {chain} (offset {offsets[chain]}):")
        for start, end in original_domains[chain]:
            new_start = start + offsets[chain]
            new_end = end + offsets[chain]
            print(f"    Domain {domain_idx+1}: [{start}, {end}] -> [{new_start}, {new_end}]")
            domain_idx += 1

    # Write merged files
    merged_fasta = f'{output_dir}/input/merged_complex.fasta'
    merged_domains_file = f'{output_dir}/input/merged_domains.yaml'
    merged_pdb = f'{output_dir}/input/FOXP_complex.pdb'

    write_merged_fasta(merged_sequence, merged_name, merged_fasta)
    write_merged_domains_yaml(merged_domains, merged_name, merged_domains_file)

    # Create merged PDB with coordinates
    total_residues = create_merged_pdb(
        complex_pdb,
        chain_mapping,
        chain_order,
        sequences,
        LINKER_SEQUENCE,
        merged_pdb
    )

    print("\n" + "="*60)
    print("Creating CALVADOS configuration files...")
    print("="*60)

    # Create Config
    residues_file = f'{output_dir}/input/residues_CALVADOS3.csv'

    config = Config(
        # GENERAL
        sysname=sysname,
        box=[BOX_SIZE, BOX_SIZE, BOX_SIZE],
        temp=293,
        ionic=0.15,
        pH=7.0,
        topol='center',  # Single protein at center

        # RUNTIME SETTINGS
        wfreq=N_SAVE,
        steps=N_FRAMES * N_SAVE,
        platform='CPU',
        restart='checkpoint',
        frestart='restart.chk',
        verbose=True,
    )

    # Analysis script
    analyses = f'''
from calvados.analysis import save_conf_prop

save_conf_prop(path="{sim_path}", name="{sysname}",
               residues_file="{residues_file}",
               output_path="{output_dir}/data",
               start=100, is_idr=False, select='all')
'''

    config.write(sim_path, name='config.yaml', analyses=analyses)

    # Create Components
    components = Components(
        # Defaults
        molecule_type='protein',
        nmol=1,
        restraint=True,
        charge_termini='both',

        # INPUT
        fresidues=residues_file,
        ffasta=merged_fasta,
        fdomains=merged_domains_file,
        pdb_folder=f'{output_dir}/input',

        # RESTRAINTS - High value for rigid domains
        restraint_type='harmonic',
        k_harmonic=INTRA_DOMAIN_K,
        cutoff_restr=0.9,
        use_com=True,
        colabfold=1,
        ext_restraint=True,
    )

    # Add the merged complex as single protein
    components.add(name=merged_name, restraint=True)

    components.write(sim_path, name='components.yaml')

    print(f"\nConfiguration files written to: {sim_path}/")
    print(f"  - config.yaml")
    print(f"  - components.yaml")

    print("\n" + "="*60)
    print("SETUP COMPLETE")
    print("="*60)
    print(f"\nTo run the simulation:")
    print(f"  cd {sim_path}")
    print(f"  python run.py")
    print(f"\nMerged complex: {merged_name}")
    print(f"Total residues: {len(merged_sequence)}")
    print(f"Number of domains: {len(merged_domains)}")
    print(f"Intra-domain restraint strength: k = {INTRA_DOMAIN_K} kJ/mol/nm^2")
    print(f"\nNote: Chain connectivity is maintained through bonded interactions.")
    print(f"Linker regions ({LINKER_SEQUENCE}) allow flexibility between original chains.")
