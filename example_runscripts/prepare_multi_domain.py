#!/usr/bin/env python
"""
prepare_multi_domain.py - Prepare simulation for a multi-domain protein with restraints

Usage:
    python prepare_multi_domain.py --name ProteinName

Requirements:
    - input/ProteinName.pdb       (3D structure)
    - input/domains.yaml          (structured domain definitions)
    - input/residues_CALVADOS3.csv (force field parameters)

This script creates:
    - config.yaml
    - components.yaml
    - run.py

Then run with:
    python ProteinName/run.py --path ProteinName
"""

import os
from argparse import ArgumentParser
from calvados.cfg import Config, Components

def prepare_multi_domain(name, pdb_folder, domains_file, residues_file,
                         box_size=100, steps=10000000, wfreq=10000,
                         platform='CPU', threads=4,
                         restraint_type='harmonic', k_harmonic=700.0):
    """
    Prepare simulation files for a multi-domain protein.

    Parameters
    ----------
    name : str
        Protein name (must match PDB filename and domains.yaml entry)
    pdb_folder : str
        Directory containing PDB file
    domains_file : str
        Path to domains.yaml
    residues_file : str
        Path to residues parameter file
    box_size : float
        Cubic box side length in nm
    steps : int
        Total simulation steps
    wfreq : int
        Trajectory write frequency
    platform : str
        'CPU' or 'CUDA'
    threads : int
        CPU threads
    restraint_type : str
        'harmonic' or 'go'
    k_harmonic : float
        Harmonic restraint force constant (kJ/mol/nm²)
    """

    cwd = os.getcwd()

    # Verify input files exist
    pdb_file = f'{pdb_folder}/{name}.pdb'
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    if not os.path.exists(domains_file):
        raise FileNotFoundError(f"Domains file not found: {domains_file}")

    # Simulation configuration
    config = Config(
        sysname=name,
        box=[box_size, box_size, box_size],
        temp=293,
        ionic=0.15,
        pH=7.0,
        topol='center',
        steps=steps,
        wfreq=wfreq,
        platform=platform,
        threads=threads,
        restart='checkpoint',
        verbose=True,
    )

    # Component definition with restraints
    components = Components(
        molecule_type='protein',
        nmol=1,
        charge_termini='both',
        fresidues=residues_file,
        fdomains=domains_file,
        pdb_folder=pdb_folder,
        restraint_type=restraint_type,
        k_harmonic=k_harmonic,
        use_com=True,           # Restraints on center of mass
        cutoff_restr=0.9,       # Distance cutoff for restraints (nm)
    )

    # Add the protein with restraints
    components.add(name=name, restraint=True)

    # Create output directory
    path = f'{cwd}/{name}'
    os.makedirs(path, exist_ok=True)

    # Analysis code
    analyses = f'''
from calvados.analysis import save_conf_prop
import os

output_path = "{cwd}/analysis"
os.makedirs(output_path, exist_ok=True)

save_conf_prop(
    path="{path}",
    name="{name}",
    residues_file="{residues_file}",
    output_path=output_path,
    start=100,
    is_idr=False,  # Multi-domain protein
    select='all'
)
print(f"Analysis saved to {{output_path}}")
'''

    config.write(path, analyses=analyses)
    components.write(path)

    print(f"Prepared multi-domain simulation for {name}")
    print(f"  PDB: {pdb_file}")
    print(f"  Domains: {domains_file}")
    print(f"  Restraint type: {restraint_type}")
    print(f"  k_harmonic: {k_harmonic} kJ/mol/nm²")
    print(f"  Box size: {box_size} nm")
    print(f"  Steps: {steps} ({steps * 0.01 / 1000:.1f} ns)")
    print(f"\nTo run:")
    print(f"  python {name}/run.py --path {name}")

    return path


if __name__ == "__main__":
    parser = ArgumentParser(description="Prepare multi-domain protein simulation")
    parser.add_argument('--name', required=True, type=str,
                        help='Protein name (must match PDB filename)')
    parser.add_argument('--pdb_folder', default='input', type=str,
                        help='Directory containing PDB file')
    parser.add_argument('--domains', default='input/domains.yaml', type=str,
                        help='Path to domains.yaml')
    parser.add_argument('--residues', default='input/residues_CALVADOS3.csv', type=str,
                        help='Path to residues parameter file')
    parser.add_argument('--box', default=100, type=float,
                        help='Box size in nm (default: 100)')
    parser.add_argument('--steps', default=10000000, type=int,
                        help='Simulation steps (default: 10M = 100 ns)')
    parser.add_argument('--wfreq', default=10000, type=int,
                        help='Write frequency (default: 10000 = 100 ps)')
    parser.add_argument('--platform', default='CPU', choices=['CPU', 'CUDA'],
                        help='Compute platform')
    parser.add_argument('--threads', default=4, type=int,
                        help='CPU threads')
    parser.add_argument('--restraint_type', default='harmonic',
                        choices=['harmonic', 'go'],
                        help='Restraint type')
    parser.add_argument('--k_harmonic', default=700.0, type=float,
                        help='Harmonic force constant (default: 700)')

    args = parser.parse_args()

    prepare_multi_domain(
        name=args.name,
        pdb_folder=args.pdb_folder,
        domains_file=args.domains,
        residues_file=args.residues,
        box_size=args.box,
        steps=args.steps,
        wfreq=args.wfreq,
        platform=args.platform,
        threads=args.threads,
        restraint_type=args.restraint_type,
        k_harmonic=args.k_harmonic
    )
