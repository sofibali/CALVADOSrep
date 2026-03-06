#!/usr/bin/env python
"""
Example: Analyze a single CALVADOS trajectory

Usage:
    python analyze_single.py path/to/simulation/

This script demonstrates how to use parp14_tools functions
to analyze a simulation trajectory.
"""

import os
import sys

# Import functions from our package
from parp14_tools import (
    load_trajectory,
    get_trajectory_info,
    compute_rg,
    compute_ree,
    compute_all_metrics,
    plot_rg_distribution,
    plot_ree_distribution,
    plot_analysis_summary,
    export_frames,
    generate_vmd_script,
    block_error,
)


def main(simulation_dir):
    """Analyze a single simulation."""

    # File paths
    top_file = f'{simulation_dir}/top.pdb'
    traj_file = None

    # Find trajectory file
    for f in os.listdir(simulation_dir):
        if f.endswith('.dcd') and not f.startswith('equilibration'):
            traj_file = f'{simulation_dir}/{f}'
            break

    if traj_file is None:
        print(f"No .dcd file found in {simulation_dir}")
        return

    # Create output directory
    output_dir = f'{simulation_dir}/analysis'
    os.makedirs(output_dir, exist_ok=True)

    #---------------------------------------------------------------------------
    # Step 1: Get trajectory info (without loading full trajectory)
    #---------------------------------------------------------------------------
    print("=" * 60)
    print("Trajectory Information")
    print("=" * 60)

    info = get_trajectory_info(top_file, traj_file)
    print(f"  Frames: {info['n_frames']}")
    print(f"  Atoms: {info['n_atoms']}")
    print(f"  Time: {info['time_ns']:.1f} ns")
    print(f"  Time step: {info['dt_ps']:.1f} ps")

    #---------------------------------------------------------------------------
    # Step 2: Load trajectory (skip equilibration)
    #---------------------------------------------------------------------------
    print("\nLoading trajectory...")

    # Skip first 10% of frames as equilibration
    skip_frames = info['n_frames'] // 10

    traj = load_trajectory(top_file, traj_file, skip=skip_frames)
    print(f"  Loaded {traj.n_frames} frames (skipped {skip_frames})")

    #---------------------------------------------------------------------------
    # Step 3: Compute properties
    #---------------------------------------------------------------------------
    print("\nComputing properties...")

    # Radius of gyration
    rg = compute_rg(traj)
    rg_error = block_error(rg)
    print(f"  Rg = {rg_error['mean']:.3f} ± {rg_error['block_sem']:.3f} nm")

    # End-to-end distance
    ree = compute_ree(traj)
    ree_error = block_error(ree)
    print(f"  Ree = {ree_error['mean']:.3f} ± {ree_error['block_sem']:.3f} nm")

    #---------------------------------------------------------------------------
    # Step 4: Create plots
    #---------------------------------------------------------------------------
    print("\nCreating plots...")

    # Individual plots
    plot_rg_distribution(rg, f'{output_dir}/rg_distribution.png')
    plot_ree_distribution(ree, f'{output_dir}/ree_distribution.png')

    # All metrics + summary plot
    metrics = compute_all_metrics(traj)
    plot_analysis_summary(metrics, f'{output_dir}/summary.png')

    #---------------------------------------------------------------------------
    # Step 5: Export frames for visualization
    #---------------------------------------------------------------------------
    print("\nExporting frames...")

    export_frames(traj, f'{output_dir}/structures')
    generate_vmd_script(top_file, traj_file, f'{output_dir}/visualize.vmd')

    #---------------------------------------------------------------------------
    # Summary
    #---------------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)
    print(f"Output saved to: {output_dir}/")
    print("\nTo visualize in VMD:")
    print(f"  vmd -e {output_dir}/visualize.vmd")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_single.py path/to/simulation/")
        print("\nExample:")
        print("  python analyze_single.py ../examples/single_IDR/FUSRGG3/")
        sys.exit(1)

    simulation_dir = sys.argv[1]

    if not os.path.isdir(simulation_dir):
        print(f"Error: Directory not found: {simulation_dir}")
        sys.exit(1)

    main(simulation_dir)
