#!/usr/bin/env python
"""
Prepare CALVADOS simulation for FOXP4 and FOXP1 (2 chains).
Uses domain restraints only - NO inter-chain restraints.

This script:
1. Creates chain PDB files for CALVADOS to read domain restraints
2. Creates an initial positions PDB with chains positioned in the box
3. Uses restart='pdb' to load the correct initial coordinates
4. Generates config.yaml and components.yaml for simulation

KEY: Uses topol='random' + restart='pdb' with correctly positioned chains.
     Intra-domain restraints are applied, but NO inter-chain restraints.
"""

import os
import sys
import numpy as np
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
sysname = 'two_chains_coils_interchain'

# Chain mapping: CALVADOS name -> PDB chain ID
chain_mapping = {
    'FOXP4': 'A',
    'FOXP1': 'B',  # Renamed from FOX
}

# Domain definitions (from domains.yaml)
# Domain 3 is the third entry for each protein
domains = {
    'FOXP4': [(118, 199), (299, 375), (451, 551)],   # domain3: 451-551
    'FOXP1': [(118, 199), (299, 375), (451, 551)],   # domain3: 451-551
}

# Restraint parameters
INTRA_DOMAIN_K = 2000.0   # kJ/mol/nm^2 - force constant for intra-domain restraints
INTER_CHAIN_K = 700.0     # kJ/mol/nm^2 - force constant for inter-chain restraints changed to 0.0 for testing

# Simulation parameters
BOX_SIZE = 50  # nm
N_SAVE = 500  # frames between saves
N_FRAMES = 1000  # total frames to save

# ============================================================================
# FUNCTIONS
# ============================================================================

def parse_pdb_ca_atoms(pdb_file):
    """Parse CA atom coordinates from PDB file."""
    ca_atoms = {}  # {(chain, resnum): (x, y, z)}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                chain = line[21]
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_atoms[(chain, resnum)] = np.array([x, y, z])

    return ca_atoms


def calculate_complex_com(pdb_file, chain_ids):
    """Calculate center of mass of the complex (CA atoms only)."""
    ca_atoms = parse_pdb_ca_atoms(pdb_file)

    coords = []
    for (chain, resnum), coord in ca_atoms.items():
        if chain in chain_ids:
            coords.append(coord)

    coords = np.array(coords)
    com = np.mean(coords, axis=0)
    return com


def extract_chain_pdb(complex_pdb, chain_id, output_pdb):
    """
    Extract a single chain from complex PDB (original coordinates).
    This is used by CALVADOS to read the internal structure for restraints.
    """
    with open(complex_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM') and line[21] == chain_id:
                f_out.write(line)
        f_out.write("TER\n")
        f_out.write("END\n")

    print(f"Extracted chain {chain_id} -> {output_pdb}")


def create_initial_positions_pdb(complex_pdb, chain_mapping, box_size, output_pdb):
    """
    Create a coarse-grained PDB with chains positioned in the box.

    This file will be used with restart='pdb' to set initial coordinates.
    Coordinates are:
    1. Extracted from complex PDB (CA atoms only)
    2. Centered at complex COM
    3. Shifted to box center

    Format matches CALVADOS top.pdb: CA atoms only, sequential residue numbering.
    """
    # Parse CA atoms from complex
    ca_atoms = parse_pdb_ca_atoms(complex_pdb)

    # Calculate complex COM
    protein_chains = list(chain_mapping.values())
    coords = [coord for (chain, _), coord in ca_atoms.items() if chain in protein_chains]
    complex_com = np.mean(coords, axis=0)

    # Box center in Angstrom
    box_center = np.array([box_size * 10 / 2] * 3)  # nm to Angstrom

    with open(output_pdb, 'w') as f:
        atom_num = 1
        res_num = 1

        for name in chain_mapping.keys():  # Maintain order: FOXP4, FOXP1
            chain_id = chain_mapping[name]

            # Get all CA atoms for this chain, sorted by residue number
            chain_cas = [(rnum, coord) for (c, rnum), coord in ca_atoms.items() if c == chain_id]
            chain_cas.sort(key=lambda x: x[0])

            for orig_resnum, coord in chain_cas:
                # Center at complex COM and shift to box center
                new_coord = coord - complex_com + box_center

                # Get residue name from original PDB
                resname = "ALA"  # Default, will be overwritten
                with open(complex_pdb, 'r') as pdb_in:
                    for line in pdb_in:
                        if (line.startswith('ATOM') and ' CA ' in line and
                            line[21] == chain_id and int(line[22:26].strip()) == orig_resnum):
                            resname = line[17:20].strip()
                            break

                # Write CA atom with sequential numbering
                f.write(f"ATOM  {atom_num:5d}  CA  {resname:3s} A{res_num:4d}    "
                       f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}"
                       f"  1.00  0.00           C\n")
                atom_num += 1
                res_num += 1

        f.write("END\n")

    print(f"Created initial positions PDB with {atom_num-1} atoms -> {output_pdb}")
    print(f"  Complex COM: [{complex_com[0]:.2f}, {complex_com[1]:.2f}, {complex_com[2]:.2f}] A")
    print(f"  Box center: [{box_center[0]:.2f}, {box_center[1]:.2f}, {box_center[2]:.2f}] A")


# ============================================================================
# MAIN SCRIPT
# ============================================================================

if __name__ == '__main__':

    # Paths
    complex_pdb = f'{input_dir}/COMPLEX_WITHDNA.pdb'
    output_dir = cwd
    restraints_file = f'{output_dir}/input/inter_chain_restraints.txt'
    sim_path = f'{output_dir}/{sysname}'
    initial_pdb = f'{sim_path}/restart.pdb'  # Initial positions file

    # Create directories
    subprocess.run(f'mkdir -p {sim_path}', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/input', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/data', shell=True)

    # Copy required input files
    subprocess.run(f'cp {input_dir}/residues_CALVADOS3.csv {output_dir}/input/', shell=True)
    subprocess.run(f'cp {input_dir}/proteins.fasta {output_dir}/input/', shell=True)
    subprocess.run(f'cp {input_dir}/domains.yaml {output_dir}/input/', shell=True)

    print("="*60)
    print("Step 1: Extract individual chain PDBs (for restraint calculation)")
    print("="*60)

    # Extract each chain with original coordinates (for CALVADOS to read structure)
    for name, pdb_chain in chain_mapping.items():
        output_pdb = f'{output_dir}/input/{name}.pdb'
        extract_chain_pdb(complex_pdb, pdb_chain, output_pdb)

    print("\n" + "="*60)
    print("Step 2: Create initial positions PDB (for restart)")
    print("="*60)

    # Create the initial positions PDB with both chains correctly positioned
    create_initial_positions_pdb(complex_pdb, chain_mapping, BOX_SIZE, initial_pdb)

    print("\n" + "="*60)
    print("Step 3: Creating CALVADOS configuration files")
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

        # Use 'random' topology (required for multiple molecules)
        # The actual positions will be loaded from restart.pdb
        topol='random',

        # RUNTIME SETTINGS
        wfreq=N_SAVE,
        steps=N_FRAMES * N_SAVE,
        platform='CPU',

        # KEY: Use restart='pdb' to load positions from restart.pdb
        # This overrides the random topology placement!
        restart='pdb',
        frestart='restart.pdb',

        verbose=True,

        # INTER-CHAIN RESTRAINTS
        custom_restraints=True,
        custom_restraint_type='harmonic',
        fcustom_restraints=restraints_file,
    )

    # Analysis script
    analyses = f'''
from calvados.analysis import calc_com_traj, calc_contact_map

chainid_dict = dict(FOXP4=0, FOXP1=1)
calc_com_traj(path="{sim_path}",sysname="{sysname}",output_path="{output_dir}/data",
              residues_file="{residues_file}",chainid_dict=chainid_dict,start=100)
calc_contact_map(path="{sim_path}",sysname="{sysname}",output_path="{output_dir}/data",
                 chainid_dict=chainid_dict)
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
        ffasta=f'{output_dir}/input/proteins.fasta',
        fdomains=f'{output_dir}/input/domains.yaml',
        pdb_folder=f'{output_dir}/input',

        # RESTRAINTS - High value for rigid domains
        restraint_type='harmonic',
        k_harmonic=INTRA_DOMAIN_K,
        cutoff_restr=0.9,
        use_com=True,
        colabfold=1,
        ext_restraint=True,
    )

    # Add each protein
    components.add(name='FOXP4', restraint=True)
    components.add(name='FOXP1', restraint=True)

    components.write(sim_path, name='components.yaml')

    print(f"\nConfiguration files written to: {sim_path}/")
    print(f"  - config.yaml")
    print(f"  - components.yaml")
    print(f"  - restart.pdb (initial positions from complex)")

    print("\n" + "="*60)
    print("SETUP COMPLETE")
    print("="*60)
    print(f"\nTo run the simulation:")
    print(f"  cd {sim_path}")
    print(f"  python run.py")
    print(f"\nKey settings:")
    print(f"  - topol: 'random' (required for multiple molecules)")
    print(f"  - restart: 'pdb' (loads positions from restart.pdb)")
    print(f"  - frestart: 'restart.pdb' (contains chain positions)")
    print(f"  - Intra-domain restraints: k = {INTRA_DOMAIN_K} kJ/mol/nm^2")
    print(f"  - NO inter-chain restraints applied")
    print(f"\nNote: Chains will start at their positions from the complex!")
    print(f"      The restart.pdb file contains both chains centered in the {BOX_SIZE} nm box.")
