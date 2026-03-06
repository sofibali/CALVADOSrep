#!/usr/bin/env python
"""
compile_trajectories.py - Merge multiple replicate trajectories into one

Usage:
    python compile_trajectories.py --input_dir simulations/comp_0001 --output compiled.dcd

Features:
    - Merges multiple DCD files from replicate simulations
    - Skips equilibration frames from each replicate
    - Handles different trajectory lengths
    - Saves compiled trajectory and merged topology

Useful for:
    - Combining multiple short simulations into ensemble trajectory
    - Creating longer effective sampling from parallel runs
    - PARP14 library: compile 25 × 40ns = 1 μs per composition
"""

import os
import sys
import glob
import numpy as np
from argparse import ArgumentParser

import MDAnalysis as mda
from MDAnalysis import Writer


def find_trajectories(input_dir, pattern='*.dcd'):
    """Find all trajectory files in directory structure."""

    # Try different directory structures
    patterns = [
        f'{input_dir}/{pattern}',              # Direct
        f'{input_dir}/*/{pattern}',            # One level down
        f'{input_dir}/rep_*/{pattern}',        # rep_XX subdirs
        f'{input_dir}/*/*{pattern}',           # Two levels down
    ]

    traj_files = []
    for p in patterns:
        found = sorted(glob.glob(p))
        if found:
            traj_files = found
            break

    return traj_files


def compile_trajectories(input_dir, output_file, top_file=None,
                         skip_ns=10, wfreq=10000, pattern='*.dcd'):
    """
    Compile multiple replicate trajectories.

    Parameters
    ----------
    input_dir : str
        Directory containing replicate simulations
    output_file : str
        Output DCD file path
    top_file : str
        Topology file (PDB). If None, auto-detect from first replicate
    skip_ns : float
        Nanoseconds to skip at start of each trajectory (equilibration)
    wfreq : int
        Write frequency used in simulations (to calculate skip frames)
    pattern : str
        Glob pattern for trajectory files
    """

    # Find trajectory files
    traj_files = find_trajectories(input_dir, pattern)

    if not traj_files:
        raise FileNotFoundError(f"No trajectory files found in {input_dir}")

    print(f"Found {len(traj_files)} trajectory files:")
    for f in traj_files:
        print(f"  {f}")

    # Find or verify topology
    if top_file is None:
        # Try to find top.pdb
        possible_tops = [
            f'{input_dir}/top.pdb',
            f'{input_dir}/rep_01/top.pdb',
            os.path.join(os.path.dirname(traj_files[0]), 'top.pdb'),
        ]
        for t in possible_tops:
            if os.path.exists(t):
                top_file = t
                break

        if top_file is None:
            raise FileNotFoundError("Could not find topology file. Specify with --top")

    print(f"\nUsing topology: {top_file}")

    # Calculate frames to skip
    # 1 step = 0.01 ps, so wfreq steps = wfreq * 0.01 ps per frame
    ps_per_frame = wfreq * 0.01
    skip_frames = int(skip_ns * 1000 / ps_per_frame)

    print(f"Skipping {skip_ns} ns = {skip_frames} frames per trajectory")

    # Get total atom count
    u_test = mda.Universe(top_file, traj_files[0])
    n_atoms = u_test.atoms.n_atoms

    # Compile trajectories
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    total_frames = 0
    n_replicates = 0

    with Writer(output_file, n_atoms) as W:
        for traj_file in traj_files:
            try:
                u = mda.Universe(top_file, traj_file)
                n_frames_total = len(u.trajectory)

                if n_frames_total <= skip_frames:
                    print(f"  Skipping {traj_file}: only {n_frames_total} frames")
                    continue

                frames_written = 0
                for ts in u.trajectory[skip_frames:]:
                    W.write(u.atoms)
                    frames_written += 1
                    total_frames += 1

                n_replicates += 1
                print(f"  Added {frames_written} frames from {os.path.basename(traj_file)}")

            except Exception as e:
                print(f"  Error processing {traj_file}: {e}")
                continue

    # Copy topology
    output_top = output_file.replace('.dcd', '_top.pdb')
    u_top = mda.Universe(top_file)
    u_top.atoms.write(output_top)

    # Summary
    total_time_ns = total_frames * ps_per_frame / 1000

    print(f"\n{'='*50}")
    print(f"Compilation complete!")
    print(f"{'='*50}")
    print(f"  Replicates compiled: {n_replicates}")
    print(f"  Total frames: {total_frames}")
    print(f"  Total time: {total_time_ns:.1f} ns ({total_time_ns/1000:.3f} μs)")
    print(f"  Output trajectory: {output_file}")
    print(f"  Output topology: {output_top}")

    return output_file, total_frames


def compile_library(library_dir, output_dir, n_compositions,
                    skip_ns=10, wfreq=10000):
    """
    Compile trajectories for an entire library of compositions.

    Parameters
    ----------
    library_dir : str
        Base directory containing simulations/comp_XXXX/ subdirs
    output_dir : str
        Directory for compiled trajectories
    n_compositions : int
        Number of compositions to process
    skip_ns : float
        Equilibration time to skip
    wfreq : int
        Write frequency
    """

    os.makedirs(output_dir, exist_ok=True)

    results = []

    for comp_id in range(1, n_compositions + 1):
        comp_str = f'comp_{comp_id:04d}'
        input_dir = f'{library_dir}/simulations/{comp_str}'
        output_file = f'{output_dir}/{comp_str}.dcd'

        if not os.path.exists(input_dir):
            print(f"Skipping {comp_str}: directory not found")
            continue

        try:
            out_file, n_frames = compile_trajectories(
                input_dir=input_dir,
                output_file=output_file,
                skip_ns=skip_ns,
                wfreq=wfreq
            )
            results.append({
                'composition': comp_str,
                'output_file': out_file,
                'n_frames': n_frames,
                'status': 'success'
            })
        except Exception as e:
            print(f"Error compiling {comp_str}: {e}")
            results.append({
                'composition': comp_str,
                'error': str(e),
                'status': 'failed'
            })

    # Summary
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(f'{output_dir}/compilation_summary.csv', index=False)

    n_success = len(df[df['status'] == 'success'])
    print(f"\nLibrary compilation: {n_success}/{n_compositions} successful")
    print(f"Summary saved to: {output_dir}/compilation_summary.csv")


if __name__ == "__main__":
    parser = ArgumentParser(description="Compile replicate trajectories")
    parser.add_argument('--input_dir', required=True, type=str,
                        help='Directory containing replicate simulations')
    parser.add_argument('--output', default='compiled.dcd', type=str,
                        help='Output DCD file')
    parser.add_argument('--top', default=None, type=str,
                        help='Topology file (auto-detect if not specified)')
    parser.add_argument('--skip_ns', default=10, type=float,
                        help='Nanoseconds to skip (equilibration)')
    parser.add_argument('--wfreq', default=10000, type=int,
                        help='Write frequency used in simulations')
    parser.add_argument('--pattern', default='*.dcd', type=str,
                        help='Glob pattern for trajectory files')

    # Library mode
    parser.add_argument('--library', action='store_true',
                        help='Process entire library')
    parser.add_argument('--n_compositions', default=0, type=int,
                        help='Number of compositions (library mode)')

    args = parser.parse_args()

    if args.library:
        if args.n_compositions == 0:
            print("Error: --n_compositions required for library mode")
            sys.exit(1)
        compile_library(
            library_dir=args.input_dir,
            output_dir=os.path.dirname(args.output) or 'trajectories',
            n_compositions=args.n_compositions,
            skip_ns=args.skip_ns,
            wfreq=args.wfreq
        )
    else:
        compile_trajectories(
            input_dir=args.input_dir,
            output_file=args.output,
            top_file=args.top,
            skip_ns=args.skip_ns,
            wfreq=args.wfreq,
            pattern=args.pattern
        )
