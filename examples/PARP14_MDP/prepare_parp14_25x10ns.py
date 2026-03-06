#!/usr/bin/env python3
"""
Prepare 25 x 10 ns CALVADOS production simulations for PARP14 construct.

Construct: kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art (1474 residues)
Initial structures: 25 AlphaFold3 seed/sample models (5 seeds x 5 samples)

Usage:
    python prepare_parp14_25x10ns.py
"""

import os
import subprocess
import shutil
import json
import yaml
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from calvados.cfg import Config, Components

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.abspath(__file__))

# AlphaFold3 structure directory
AF3_DIR = '/home/sbali/parp14/alphafold_outputs/kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art'

# Component name used in CALVADOS
COMP_NAME = 'parp14_construct'

# Residues file
RESIDUES_FILE = os.path.join(CWD, 'input/residues_CALVADOS3.csv')

# Shared input directory for domains.yaml (same for all 25 sims)
SHARED_INPUT = os.path.join(CWD, 'input_construct')

# Simulation parameters
N_STEPS = 1_000_000     # 10 ns at 0.01 ps timestep
WFREQ = 1_000           # save every 1000 steps -> 1000 frames
BOX_L = 150             # nm cubic box
TEMP = 293              # K
IONIC = 0.19            # M
PH = 7.0
THREADS = 4
PLATFORM = 'CPU'

# Restraint parameters
K_HARMONIC = 700.0      # kJ/mol/nm^2
COLABFOLD = 2           # AF3 PAE format: 'pae' key in JSON

# Construct domain boundaries (construct numbering, 1-indexed)
# Full-length 315-1801, skipping linker 1194-1206
# Offset: FL->construct = -314 (before gap), -327 (after gap)
CONSTRUCT_DOMAINS = [
    [1, 423],       # KH1-KH6 (FL 315-737)
    [424, 475],     # KH7a (FL 738-789)
    [476, 690],     # MD1L1 (FL 790-1004)
    [690, 879],     # MD2 (FL 1004-1193)
    [880, 1061],    # MD3 (FL 1207-1388)
    [1062, 1206],   # KHb-KH8 (FL 1389-1533)
    [1207, 1275],   # WWE (FL 1534-1602)
    [1276, 1474],   # ART (FL 1603-1801)
]

DOMAIN_LABELS = ['KH1-KH6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

# Seeds and samples
SEEDS = range(1, 6)     # seed-1 to seed-5
SAMPLES = range(0, 5)   # sample-0 to sample-4


# ============================================================
# Helper classes and functions
# ============================================================

class BackboneSelect(Select):
    """Select only backbone atoms for PDB output."""
    def accept_atom(self, atom):
        return atom.get_name() in ['N', 'CA', 'C', 'O']


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB format (backbone only)."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('model', cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path, BackboneSelect())


def prepare_pae_json(confidences_json, output_json):
    """
    Convert AF3 confidences.json to CALVADOS-compatible PAE JSON.
    AF3 uses 'pae' key; CALVADOS with colabfold=2 reads 'pae' key.
    """
    shutil.copy2(confidences_json, output_json)


def create_run_script(sim_dir):
    """Write the standard CALVADOS run.py script."""
    run_py = os.path.join(sim_dir, 'run.py')
    with open(run_py, 'w') as f:
        f.write("""from calvados import sim
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--path', nargs='?', default='.', const='.', type=str)
    parser.add_argument('--config', nargs='?', default='config.yaml', const='config.yaml', type=str)
    parser.add_argument('--components', nargs='?', default='components.yaml', const='components.yaml', type=str)

    args = parser.parse_args()
    sim.run(path=args.path, fconfig=args.config, fcomponents=args.components)
""")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("PARP14 Construct: Preparing 25 x 10 ns Production Simulations")
    print("=" * 70)

    # Create shared input directory
    os.makedirs(SHARED_INPUT, exist_ok=True)

    # Write construct domains.yaml
    domains_yaml_path = os.path.join(SHARED_INPUT, 'domains.yaml')
    domains_dict = {COMP_NAME: CONSTRUCT_DOMAINS}
    with open(domains_yaml_path, 'w') as f:
        yaml.dump(domains_dict, f, default_flow_style=False, sort_keys=False)
    print(f"\nDomains YAML: {domains_yaml_path}")
    for label, bounds in zip(DOMAIN_LABELS, CONSTRUCT_DOMAINS):
        print(f"  {label}: {bounds}")

    # Copy residues file to shared input
    residues_dest = os.path.join(SHARED_INPUT, 'residues_CALVADOS3.csv')
    if not os.path.exists(residues_dest):
        shutil.copy2(RESIDUES_FILE, residues_dest)

    # Iterate over all 25 seed/sample combinations
    n_prepared = 0
    n_skipped = 0

    for seed in SEEDS:
        for sample in SAMPLES:
            seed_sample = f'seed-{seed}_sample-{sample}'
            sim_name = f'parp14_{seed_sample}'
            sim_dir = os.path.join(CWD, sim_name)
            input_dir = os.path.join(sim_dir, 'input')

            # Source files
            cif_path = os.path.join(AF3_DIR, seed_sample, 'model.cif')
            conf_json = os.path.join(AF3_DIR, seed_sample, 'confidences.json')

            if not os.path.isfile(cif_path):
                print(f"  SKIP: {seed_sample} - model.cif not found")
                n_skipped += 1
                continue

            print(f"\n  Preparing: {sim_name}")

            # Create directories
            os.makedirs(sim_dir, exist_ok=True)
            os.makedirs(input_dir, exist_ok=True)

            # Convert CIF -> PDB (backbone only)
            pdb_path = os.path.join(input_dir, f'{COMP_NAME}.pdb')
            convert_cif_to_pdb(cif_path, pdb_path)
            print(f"    PDB: {pdb_path}")

            # Copy PAE JSON
            pae_path = os.path.join(input_dir, f'{COMP_NAME}.json')
            prepare_pae_json(conf_json, pae_path)
            print(f"    PAE: {pae_path}")

            # Create config.yaml
            config = Config(
                sysname=COMP_NAME,
                box=[BOX_L, BOX_L, BOX_L],
                temp=TEMP,
                ionic=IONIC,
                pH=PH,
                topol='center',
                wfreq=WFREQ,
                steps=N_STEPS,
                runtime=0,
                platform=PLATFORM,
                threads=THREADS,
                restart='checkpoint',
                frestart='restart.chk',
                verbose=True,
            )
            config.write(sim_dir, name='config.yaml')
            print(f"    Config: {N_STEPS} steps, wfreq={WFREQ} -> {N_STEPS//WFREQ} frames")

            # Create components.yaml
            components = Components(
                molecule_type='protein',
                nmol=1,
                restraint=True,
                charge_termini='both',
                fresidues=os.path.join(SHARED_INPUT, 'residues_CALVADOS3.csv'),
                fdomains=domains_yaml_path,
                pdb_folder=input_dir,
                restraint_type='harmonic',
                use_com=True,
                colabfold=COLABFOLD,
                k_harmonic=K_HARMONIC,
            )
            components.add(name=COMP_NAME)
            components.write(sim_dir, name='components.yaml')

            # Write run.py
            create_run_script(sim_dir)

            n_prepared += 1

    # Summary
    print("\n" + "=" * 70)
    print("PREPARATION COMPLETE")
    print("=" * 70)
    print(f"  Simulations prepared: {n_prepared}")
    print(f"  Simulations skipped:  {n_skipped}")
    print(f"\n  Simulation parameters:")
    print(f"    Steps:       {N_STEPS:,} ({N_STEPS * 0.01 / 1000:.1f} ns)")
    print(f"    Save freq:   {WFREQ} ({N_STEPS // WFREQ} frames)")
    print(f"    Box:         {BOX_L} nm cubic")
    print(f"    Temperature: {TEMP} K")
    print(f"    Ionic:       {IONIC} M")
    print(f"    Restraints:  harmonic, k={K_HARMONIC} kJ/mol/nm^2")
    print(f"    Platform:    {PLATFORM}")
    print(f"\n  To run a simulation:")
    print(f"    cd parp14_seed-1_sample-0 && python run.py")
    print(f"\n  To run all 25 in parallel (e.g. with SLURM):")
    print(f"    for d in parp14_seed-*_sample-*; do")
    print(f"      cd $d && python run.py &")
    print(f"      cd ..")
    print(f"    done")
    print("=" * 70)


if __name__ == '__main__':
    main()
