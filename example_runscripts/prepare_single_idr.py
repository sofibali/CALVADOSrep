#!/usr/bin/env python
"""
prepare_single_idr.py - Prepare simulation for a single intrinsically disordered protein

Usage:
    python prepare_single_idr.py --name ProteinName --fasta input/sequences.fasta

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

def prepare_single_idr(name, fasta_file, residues_file, box_size=50, steps=10000000,
                       wfreq=10000, platform='CPU', threads=4):
    """
    Prepare simulation files for a single IDR.

    Parameters
    ----------
    name : str
        Protein name (must match entry in FASTA file)
    fasta_file : str
        Path to FASTA file containing sequence
    residues_file : str
        Path to residues parameter file (e.g., residues_CALVADOS3.csv)
    box_size : float
        Cubic box side length in nm (default: 50)
    steps : int
        Total simulation steps (default: 10M = 100 ns)
    wfreq : int
        Trajectory write frequency (default: 10000 = 100 ps)
    platform : str
        'CPU' or 'CUDA'
    threads : int
        CPU threads (only used if platform='CPU')
    """

    cwd = os.getcwd()

    # Simulation configuration
    config = Config(
        sysname=name,
        box=[box_size, box_size, box_size],
        temp=293,           # K
        ionic=0.15,         # M
        pH=7.0,
        topol='center',     # Single molecule at box center
        steps=steps,
        wfreq=wfreq,
        platform=platform,
        threads=threads,
        restart='checkpoint',
        verbose=True,
    )

    # Component definition
    components = Components(
        molecule_type='protein',
        nmol=1,
        charge_termini='both',
        fresidues=residues_file,
        ffasta=fasta_file,
    )

    # Add the protein (no restraints for IDR)
    components.add(name=name, restraint=False)

    # Create output directory
    path = f'{cwd}/{name}'
    os.makedirs(path, exist_ok=True)

    # Optional: Add analysis code to run after simulation
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
    is_idr=True,
    select='all'
)
print(f"Analysis saved to {{output_path}}")
'''

    # Write configuration files
    config.write(path, analyses=analyses)
    components.write(path)

    print(f"Prepared simulation for {name}")
    print(f"  Output directory: {path}")
    print(f"  Box size: {box_size} nm")
    print(f"  Steps: {steps} ({steps * 0.01 / 1000:.1f} ns)")
    print(f"  Platform: {platform}")
    print(f"\nTo run:")
    print(f"  python {name}/run.py --path {name}")

    return path


if __name__ == "__main__":
    parser = ArgumentParser(description="Prepare single IDR simulation")
    parser.add_argument('--name', required=True, type=str,
                        help='Protein name (must match FASTA entry)')
    parser.add_argument('--fasta', default='input/sequences.fasta', type=str,
                        help='Path to FASTA file')
    parser.add_argument('--residues', default='input/residues_CALVADOS3.csv', type=str,
                        help='Path to residues parameter file')
    parser.add_argument('--box', default=50, type=float,
                        help='Box size in nm (default: 50)')
    parser.add_argument('--steps', default=10000000, type=int,
                        help='Simulation steps (default: 10M = 100 ns)')
    parser.add_argument('--wfreq', default=10000, type=int,
                        help='Write frequency (default: 10000 = 100 ps)')
    parser.add_argument('--platform', default='CPU', choices=['CPU', 'CUDA'],
                        help='Compute platform')
    parser.add_argument('--threads', default=4, type=int,
                        help='CPU threads')

    args = parser.parse_args()

    prepare_single_idr(
        name=args.name,
        fasta_file=args.fasta,
        residues_file=args.residues,
        box_size=args.box,
        steps=args.steps,
        wfreq=args.wfreq,
        platform=args.platform,
        threads=args.threads
    )
