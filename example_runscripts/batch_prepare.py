#!/usr/bin/env python
"""
batch_prepare.py - Prepare simulations for multiple proteins or variants

Usage:
    # Prepare from list of names
    python batch_prepare.py --names protein1 protein2 protein3 --fasta sequences.fasta

    # Prepare from CSV file
    python batch_prepare.py --csv compositions.csv --pdb_dir structures/

    # Prepare specific range (for parallel job preparation)
    python batch_prepare.py --csv compositions.csv --start 1 --end 100

CSV format:
    name,domains,structure_file
    comp_0001,"[1,3,5]",structures/comp_0001.pdb
    comp_0002,"[2,4,6]",structures/comp_0002.pdb
"""

import os
import sys
import pandas as pd
from argparse import ArgumentParser

# Add CALVADOS to path
try:
    from calvados.cfg import Config, Components
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from calvados.cfg import Config, Components


def prepare_idr(name, fasta_file, residues_file, output_base, **config_kwargs):
    """Prepare simulation for a single IDR."""

    path = f'{output_base}/{name}'
    os.makedirs(path, exist_ok=True)

    config = Config(
        sysname=name,
        box=config_kwargs.get('box', [50, 50, 50]),
        temp=config_kwargs.get('temp', 293),
        ionic=config_kwargs.get('ionic', 0.15),
        pH=config_kwargs.get('pH', 7.0),
        steps=config_kwargs.get('steps', 5000000),
        wfreq=config_kwargs.get('wfreq', 10000),
        platform=config_kwargs.get('platform', 'CUDA'),
        restart=None,
        verbose=False,
    )

    components = Components(
        fresidues=residues_file,
        ffasta=fasta_file,
    )

    components.add(name=name, restraint=False)

    config.write(path)
    components.write(path)

    return path


def prepare_mdp(name, pdb_folder, domains_file, residues_file, output_base, **config_kwargs):
    """Prepare simulation for a multi-domain protein."""

    path = f'{output_base}/{name}'
    os.makedirs(path, exist_ok=True)

    config = Config(
        sysname=name,
        box=config_kwargs.get('box', [100, 100, 100]),
        temp=config_kwargs.get('temp', 293),
        ionic=config_kwargs.get('ionic', 0.19),
        pH=config_kwargs.get('pH', 7.0),
        steps=config_kwargs.get('steps', 5000000),
        wfreq=config_kwargs.get('wfreq', 10000),
        platform=config_kwargs.get('platform', 'CUDA'),
        restart=None,
        verbose=False,
    )

    components = Components(
        fresidues=residues_file,
        fdomains=domains_file,
        pdb_folder=pdb_folder,
        restraint_type='harmonic',
        k_harmonic=config_kwargs.get('k_harmonic', 700.),
        use_com=True,
    )

    components.add(name=name, restraint=True)

    config.write(path)
    components.write(path)

    return path


def prepare_from_names(names, fasta_file, residues_file, output_base, **config_kwargs):
    """Prepare simulations for a list of protein names."""

    prepared = []

    for name in names:
        try:
            path = prepare_idr(name, fasta_file, residues_file, output_base, **config_kwargs)
            prepared.append({'name': name, 'path': path, 'status': 'success'})
            print(f"  Prepared: {name}")
        except Exception as e:
            prepared.append({'name': name, 'error': str(e), 'status': 'failed'})
            print(f"  Failed: {name} - {e}")

    return prepared


def prepare_from_csv(csv_file, residues_file, output_base, start=None, end=None,
                     mode='mdp', **config_kwargs):
    """
    Prepare simulations from CSV file.

    CSV columns:
        name : str - Protein/composition name
        domains : str - Domain list (optional, for MDP mode)
        pdb_folder : str - Path to PDB folder (optional)
        domains_file : str - Path to domains.yaml (optional)
    """

    df = pd.read_csv(csv_file)

    # Apply range filter
    if start is not None:
        df = df.iloc[start-1:]
    if end is not None:
        df = df.iloc[:end-start+1] if start else df.iloc[:end]

    print(f"Preparing {len(df)} simulations from {csv_file}")

    prepared = []

    for idx, row in df.iterrows():
        name = row['name']

        try:
            if mode == 'idr':
                path = prepare_idr(
                    name=name,
                    fasta_file=row.get('fasta_file', 'input/sequences.fasta'),
                    residues_file=residues_file,
                    output_base=output_base,
                    **config_kwargs
                )
            else:  # mdp
                path = prepare_mdp(
                    name=name,
                    pdb_folder=row.get('pdb_folder', 'input'),
                    domains_file=row.get('domains_file', f'input/domains/{name}.yaml'),
                    residues_file=residues_file,
                    output_base=output_base,
                    **config_kwargs
                )

            prepared.append({'name': name, 'path': path, 'status': 'success'})
            print(f"  [{idx+1}/{len(df)}] Prepared: {name}")

        except Exception as e:
            prepared.append({'name': name, 'error': str(e), 'status': 'failed'})
            print(f"  [{idx+1}/{len(df)}] Failed: {name} - {e}")

    return prepared


def prepare_parp14_library(compositions_file, structures_dir, domains_dir,
                           residues_file, output_base, n_structures=25,
                           start_comp=1, end_comp=None, **config_kwargs):
    """
    Prepare PARP14 domain deletion library.

    For each composition:
        - Creates simulations for n_structures initial structures
        - Uses composition-specific domains.yaml

    Directory structure expected:
        structures_dir/comp_XXXX/s01.pdb, s02.pdb, ...
        domains_dir/domains_XXXX.yaml
    """

    df = pd.read_csv(compositions_file)

    if end_comp is None:
        end_comp = len(df)

    df_subset = df.iloc[start_comp-1:end_comp]

    print(f"Preparing compositions {start_comp} to {end_comp}")
    print(f"Structures per composition: {n_structures}")

    prepared = []

    for _, row in df_subset.iterrows():
        comp_id = row['composition_id']
        comp_str = f'comp_{comp_id:04d}'

        for struct_id in range(1, n_structures + 1):
            name = f'{comp_str}_s{struct_id:02d}'
            path = f'{output_base}/{comp_str}/rep_{struct_id:02d}'

            try:
                os.makedirs(path, exist_ok=True)

                config = Config(
                    sysname=name,
                    box=config_kwargs.get('box', [100, 100, 100]),
                    temp=config_kwargs.get('temp', 293),
                    ionic=config_kwargs.get('ionic', 0.19),
                    pH=7.0,
                    steps=config_kwargs.get('steps', 5000000),
                    wfreq=config_kwargs.get('wfreq', 10000),
                    platform='CUDA',
                    restart=None,
                    verbose=False,
                )

                components = Components(
                    fresidues=residues_file,
                    fdomains=f'{domains_dir}/domains_{comp_id:04d}.yaml',
                    pdb_folder=f'{structures_dir}/{comp_str}',
                    restraint_type='harmonic',
                    k_harmonic=700.,
                    use_com=True,
                )

                components.add(name=f's{struct_id:02d}', restraint=True)

                config.write(path)
                components.write(path)

                prepared.append({
                    'composition': comp_str,
                    'structure': struct_id,
                    'path': path,
                    'status': 'success'
                })

            except Exception as e:
                prepared.append({
                    'composition': comp_str,
                    'structure': struct_id,
                    'error': str(e),
                    'status': 'failed'
                })

        print(f"  Prepared: {comp_str} ({n_structures} structures)")

    return prepared


def main():
    parser = ArgumentParser(description="Batch prepare CALVADOS simulations")

    # Input sources (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--names', nargs='+', type=str,
                             help='List of protein names')
    input_group.add_argument('--csv', type=str,
                             help='CSV file with protein definitions')
    input_group.add_argument('--parp14', type=str,
                             help='PARP14 library mode: path to compositions CSV')

    # Common arguments
    parser.add_argument('--residues', default='input/residues_CALVADOS3.csv',
                        help='Residues parameter file')
    parser.add_argument('--output', default='simulations',
                        help='Output base directory')
    parser.add_argument('--fasta', default='input/sequences.fasta',
                        help='FASTA file (for IDR mode)')

    # Range selection
    parser.add_argument('--start', type=int, default=None,
                        help='Start index (1-based)')
    parser.add_argument('--end', type=int, default=None,
                        help='End index (inclusive)')

    # Simulation parameters
    parser.add_argument('--mode', choices=['idr', 'mdp'], default='mdp',
                        help='Simulation mode')
    parser.add_argument('--box', type=float, default=100,
                        help='Box size (nm)')
    parser.add_argument('--steps', type=int, default=5000000,
                        help='Simulation steps')
    parser.add_argument('--wfreq', type=int, default=10000,
                        help='Write frequency')
    parser.add_argument('--platform', choices=['CPU', 'CUDA'], default='CUDA',
                        help='Compute platform')

    # PARP14-specific
    parser.add_argument('--structures_dir', default='input/structures',
                        help='Directory with structure PDBs (PARP14 mode)')
    parser.add_argument('--domains_dir', default='input/domains',
                        help='Directory with domains.yaml files (PARP14 mode)')
    parser.add_argument('--n_structures', type=int, default=25,
                        help='Structures per composition (PARP14 mode)')

    args = parser.parse_args()

    config_kwargs = {
        'box': [args.box, args.box, args.box],
        'steps': args.steps,
        'wfreq': args.wfreq,
        'platform': args.platform,
    }

    if args.names:
        results = prepare_from_names(
            names=args.names,
            fasta_file=args.fasta,
            residues_file=args.residues,
            output_base=args.output,
            **config_kwargs
        )

    elif args.csv:
        results = prepare_from_csv(
            csv_file=args.csv,
            residues_file=args.residues,
            output_base=args.output,
            start=args.start,
            end=args.end,
            mode=args.mode,
            **config_kwargs
        )

    elif args.parp14:
        results = prepare_parp14_library(
            compositions_file=args.parp14,
            structures_dir=args.structures_dir,
            domains_dir=args.domains_dir,
            residues_file=args.residues,
            output_base=args.output,
            n_structures=args.n_structures,
            start_comp=args.start or 1,
            end_comp=args.end,
            **config_kwargs
        )

    # Summary
    df_results = pd.DataFrame(results)
    n_success = len(df_results[df_results['status'] == 'success'])
    n_total = len(df_results)

    print(f"\n{'='*50}")
    print(f"Preparation complete: {n_success}/{n_total} successful")
    print(f"{'='*50}")

    # Save summary
    summary_file = f'{args.output}/preparation_summary.csv'
    df_results.to_csv(summary_file, index=False)
    print(f"Summary saved to: {summary_file}")


if __name__ == "__main__":
    main()
