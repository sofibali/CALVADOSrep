#!/usr/bin/env python
"""
analyze_trajectory.py - Comprehensive trajectory analysis

Usage:
    python analyze_trajectory.py --top simulation/top.pdb --traj simulation/protein.dcd

Computes:
    - Radius of gyration (Rg) with error estimation
    - End-to-end distance (Ree)
    - Scaling exponent (nu) for IDRs
    - Contact map
    - RMSF per residue

Outputs:
    - analysis/rg_timeseries.npy
    - analysis/ree_timeseries.npy
    - analysis/contact_map.npy
    - analysis/rmsf.npy
    - analysis/metrics.csv
    - analysis/rg_distribution.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser

import MDAnalysis as mda
from MDAnalysis.analysis import rms

# Add CALVADOS to path if needed
try:
    from calvados.analysis import calc_rg, calc_ete, calc_cmap, fit_scaling_exp
    from calvados.BLOCKING.main import BlockAnalysis
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from calvados.analysis import calc_rg, calc_ete, calc_cmap, fit_scaling_exp
    from calvados.BLOCKING.main import BlockAnalysis


def analyze_trajectory(top_file, traj_file, residues_file, output_dir='analysis',
                       skip_frames=100, is_idr=True):
    """
    Comprehensive trajectory analysis.

    Parameters
    ----------
    top_file : str
        Path to topology (PDB) file
    traj_file : str
        Path to trajectory (DCD) file
    residues_file : str
        Path to residues CSV file
    output_dir : str
        Output directory for results
    skip_frames : int
        Number of initial frames to skip (equilibration)
    is_idr : bool
        If True, compute scaling exponent
    """

    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading trajectory: {traj_file}")
    u = mda.Universe(top_file, traj_file)
    ag = u.select_atoms('all')

    n_frames = len(u.trajectory)
    n_atoms = ag.n_atoms
    print(f"  Frames: {n_frames}")
    print(f"  Atoms: {n_atoms}")
    print(f"  Skipping first {skip_frames} frames")

    # Load residue parameters
    residues = pd.read_csv(residues_file).set_index('three')

    results = {}

    # 1. Radius of gyration
    print("\nComputing radius of gyration...")
    rgs = calc_rg(u, ag, ag.resnames.tolist(), residues, start=skip_frames)
    np.save(f'{output_dir}/rg_timeseries.npy', rgs)

    block_rg = BlockAnalysis(rgs)
    block_rg.SEM()

    results['rg_mean'] = np.mean(rgs)
    results['rg_std'] = np.std(rgs)
    results['rg_sem'] = block_rg.sem
    print(f"  Rg = {results['rg_mean']:.3f} ± {results['rg_sem']:.3f} nm")

    # 2. End-to-end distance
    print("\nComputing end-to-end distance...")
    rees, ree_mean, ree_sem = calc_ete(u, ag, start=skip_frames)
    np.save(f'{output_dir}/ree_timeseries.npy', rees)

    results['ree_mean'] = ree_mean
    results['ree_sem'] = ree_sem
    print(f"  Ree = {results['ree_mean']:.3f} ± {results['ree_sem']:.3f} nm")

    # 3. Scaling exponent (IDRs only)
    if is_idr:
        print("\nFitting scaling exponent...")
        try:
            ij, dij, r0, nu, nu_err = fit_scaling_exp(u, ag, start=skip_frames)
            np.save(f'{output_dir}/internal_distances.npy', np.array([ij, dij]))
            results['nu'] = nu
            results['nu_err'] = nu_err
            results['r0'] = r0
            print(f"  ν = {nu:.3f} ± {nu_err:.3f}")
        except Exception as e:
            print(f"  Could not fit scaling exponent: {e}")

    # 4. Contact map
    print("\nComputing contact map...")
    cmap_sum = np.zeros((n_atoms, n_atoms))
    n_analyzed = 0

    for ts in u.trajectory[skip_frames:]:
        cmap_sum += calc_cmap(ag, ag, cutoff=1.0)
        n_analyzed += 1

    cmap_avg = cmap_sum / n_analyzed
    np.save(f'{output_dir}/contact_map.npy', cmap_avg)
    print(f"  Averaged over {n_analyzed} frames")

    # 5. RMSF
    print("\nComputing RMSF...")
    rmsf_analysis = rms.RMSF(ag).run(start=skip_frames)
    rmsf = rmsf_analysis.results.rmsf / 10.0  # Convert to nm
    np.save(f'{output_dir}/rmsf.npy', rmsf)
    print(f"  Mean RMSF = {np.mean(rmsf):.3f} nm")

    # Save summary metrics
    df = pd.DataFrame([results])
    df.to_csv(f'{output_dir}/metrics.csv', index=False)
    print(f"\nMetrics saved to {output_dir}/metrics.csv")

    # Generate plots
    print("\nGenerating plots...")
    generate_plots(rgs, rees, cmap_avg, rmsf, output_dir)

    return results


def generate_plots(rgs, rees, cmap, rmsf, output_dir):
    """Generate analysis plots."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Rg distribution
    ax = axes[0, 0]
    ax.hist(rgs, bins=50, density=True, alpha=0.7, color='steelblue')
    ax.axvline(np.mean(rgs), color='red', linestyle='--',
               label=f'Mean: {np.mean(rgs):.2f} nm')
    ax.set_xlabel('Rg (nm)')
    ax.set_ylabel('Probability density')
    ax.set_title('Radius of Gyration Distribution')
    ax.legend()

    # Ree distribution
    ax = axes[0, 1]
    ax.hist(rees, bins=50, density=True, alpha=0.7, color='forestgreen')
    ax.axvline(np.mean(rees), color='red', linestyle='--',
               label=f'Mean: {np.mean(rees):.2f} nm')
    ax.set_xlabel('Ree (nm)')
    ax.set_ylabel('Probability density')
    ax.set_title('End-to-End Distance Distribution')
    ax.legend()

    # Contact map
    ax = axes[1, 0]
    im = ax.imshow(cmap, cmap='hot', origin='lower', aspect='equal')
    ax.set_xlabel('Residue')
    ax.set_ylabel('Residue')
    ax.set_title('Average Contact Map')
    plt.colorbar(im, ax=ax, label='Contact probability')

    # RMSF
    ax = axes[1, 1]
    ax.plot(np.arange(1, len(rmsf) + 1), rmsf, color='purple', lw=1)
    ax.fill_between(np.arange(1, len(rmsf) + 1), 0, rmsf, alpha=0.3, color='purple')
    ax.set_xlabel('Residue')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title('Root Mean Square Fluctuation')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/analysis_summary.png', dpi=150)
    plt.close()

    print(f"  Saved: {output_dir}/analysis_summary.png")


if __name__ == "__main__":
    parser = ArgumentParser(description="Analyze CALVADOS trajectory")
    parser.add_argument('--top', required=True, type=str,
                        help='Topology file (PDB)')
    parser.add_argument('--traj', required=True, type=str,
                        help='Trajectory file (DCD)')
    parser.add_argument('--residues', default='input/residues_CALVADOS3.csv', type=str,
                        help='Residues parameter file')
    parser.add_argument('--output', default='analysis', type=str,
                        help='Output directory')
    parser.add_argument('--skip', default=100, type=int,
                        help='Frames to skip (equilibration)')
    parser.add_argument('--idr', action='store_true',
                        help='Compute scaling exponent (for IDRs)')

    args = parser.parse_args()

    analyze_trajectory(
        top_file=args.top,
        traj_file=args.traj,
        residues_file=args.residues,
        output_dir=args.output,
        skip_frames=args.skip,
        is_idr=args.idr
    )
