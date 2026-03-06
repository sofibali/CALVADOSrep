#!/usr/bin/env python
"""
visualize_trajectory.py - Trajectory visualization and export utilities

Usage:
    python visualize_trajectory.py --top simulation/top.pdb --traj simulation/protein.dcd

Features:
    - Export first/last frames as PDB
    - Create subsampled trajectory
    - Generate Rg/Ree time series plot
    - Create contact map heatmap
    - Export trajectory to other formats (XTC, XYZ)
    - Optional: Launch NGLView in Jupyter

Outputs:
    - viz/first_frame.pdb
    - viz/last_frame.pdb
    - viz/subsampled.dcd (every 10th frame)
    - viz/timeseries.png
    - viz/contact_map.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

import mdtraj as md
import MDAnalysis as mda


def export_frames(traj, output_dir):
    """Export first, last, and representative frames."""

    # First and last
    traj[0].save_pdb(f'{output_dir}/first_frame.pdb')
    traj[-1].save_pdb(f'{output_dir}/last_frame.pdb')
    print(f"  Exported: first_frame.pdb, last_frame.pdb")

    # Middle frame
    mid = len(traj) // 2
    traj[mid].save_pdb(f'{output_dir}/middle_frame.pdb')
    print(f"  Exported: middle_frame.pdb (frame {mid})")

    # Representative frame (closest to mean Rg)
    rg = md.compute_rg(traj)
    mean_rg = np.mean(rg)
    closest_idx = np.argmin(np.abs(rg - mean_rg))
    traj[closest_idx].save_pdb(f'{output_dir}/representative_frame.pdb')
    print(f"  Exported: representative_frame.pdb (frame {closest_idx}, Rg={rg[closest_idx]:.3f} nm)")


def subsample_trajectory(traj, output_dir, stride=10):
    """Create subsampled trajectory."""

    subsampled = traj[::stride]
    subsampled.save_dcd(f'{output_dir}/subsampled.dcd')
    print(f"  Exported: subsampled.dcd ({len(subsampled)} frames, stride={stride})")

    return subsampled


def convert_formats(traj, output_dir):
    """Export trajectory to other formats."""

    # XTC (GROMACS)
    traj.save_xtc(f'{output_dir}/trajectory.xtc')
    print(f"  Exported: trajectory.xtc")

    # XYZ (simple text format)
    traj[::100].save_xyz(f'{output_dir}/trajectory.xyz')  # Subsample for XYZ
    print(f"  Exported: trajectory.xyz (every 100th frame)")


def plot_timeseries(traj, output_dir):
    """Plot Rg and Ree time series."""

    # Compute properties
    rg = md.compute_rg(traj)

    # End-to-end distance (first to last CA)
    pairs = [[0, traj.n_atoms - 1]]
    ree = md.compute_distances(traj, pairs).flatten()

    time_ns = traj.time / 1000  # Convert ps to ns

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Rg time series
    ax = axes[0, 0]
    ax.plot(time_ns, rg, lw=0.5, color='steelblue')
    ax.axhline(np.mean(rg), color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg (nm)')
    ax.set_title(f'Radius of Gyration (mean: {np.mean(rg):.3f} nm)')

    # Rg histogram
    ax = axes[0, 1]
    ax.hist(rg, bins=50, density=True, alpha=0.7, color='steelblue')
    ax.axvline(np.mean(rg), color='red', linestyle='--')
    ax.set_xlabel('Rg (nm)')
    ax.set_ylabel('Probability')
    ax.set_title('Rg Distribution')

    # Ree time series
    ax = axes[1, 0]
    ax.plot(time_ns, ree, lw=0.5, color='forestgreen')
    ax.axhline(np.mean(ree), color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Ree (nm)')
    ax.set_title(f'End-to-End Distance (mean: {np.mean(ree):.3f} nm)')

    # Ree histogram
    ax = axes[1, 1]
    ax.hist(ree, bins=50, density=True, alpha=0.7, color='forestgreen')
    ax.axvline(np.mean(ree), color='red', linestyle='--')
    ax.set_xlabel('Ree (nm)')
    ax.set_ylabel('Probability')
    ax.set_title('Ree Distribution')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/timeseries.png', dpi=150)
    plt.close()

    print(f"  Saved: timeseries.png")

    return rg, ree


def plot_contact_map(traj, output_dir, skip_fraction=0.2, cutoff=0.8):
    """Compute and plot average contact map."""

    n_skip = int(skip_fraction * traj.n_frames)
    traj_prod = traj[n_skip:]

    # Compute contact map using MDTraj
    # Using CA-CA distances
    n_residues = traj.n_atoms
    pairs = []
    for i in range(n_residues):
        for j in range(i + 2, n_residues):  # Skip neighbors
            pairs.append([i, j])

    if len(pairs) == 0:
        print("  Not enough residues for contact map")
        return

    pairs = np.array(pairs)
    distances = md.compute_distances(traj_prod, pairs)

    # Average contact probability
    contacts = (distances < cutoff).astype(float)
    avg_contacts = np.mean(contacts, axis=0)

    # Reconstruct matrix
    cmap = np.zeros((n_residues, n_residues))
    for idx, (i, j) in enumerate(pairs):
        cmap[i, j] = avg_contacts[idx]
        cmap[j, i] = avg_contacts[idx]

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(cmap, cmap='hot', origin='lower', aspect='equal',
                   vmin=0, vmax=1)
    ax.set_xlabel('Residue')
    ax.set_ylabel('Residue')
    ax.set_title(f'Average Contact Map (cutoff={cutoff} nm)')
    plt.colorbar(im, ax=ax, label='Contact probability')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/contact_map.png', dpi=150)
    plt.close()

    # Save data
    np.save(f'{output_dir}/contact_map.npy', cmap)

    print(f"  Saved: contact_map.png, contact_map.npy")


def generate_vmd_script(top_file, traj_file, output_dir):
    """Generate VMD visualization script."""

    vmd_script = f'''# VMD visualization script for CALVADOS trajectory
# Usage: vmd -e {output_dir}/visualize.vmd

# Load structure and trajectory
mol new {top_file} type pdb
mol addfile {traj_file} type dcd waitfor all

# Visualization style for coarse-grained
mol modstyle 0 0 VDW 0.8 12
mol modcolor 0 0 ResType
mol modmaterial 0 0 AOChalky

# Alternative: connected beads
# mol modstyle 0 0 Licorice 0.3 12 12

# Color by residue index (rainbow)
# mol modcolor 0 0 Index

# Center view
display resetview

# Smooth trajectory playback
animate speed 0.9
animate forward

puts "Trajectory loaded. Use animation controls to play."
puts "Frames: [molinfo 0 get numframes]"
'''

    with open(f'{output_dir}/visualize.vmd', 'w') as f:
        f.write(vmd_script)

    print(f"  Saved: visualize.vmd")
    print(f"    Run with: vmd -e {output_dir}/visualize.vmd")


def main(top_file, traj_file, output_dir='viz', stride=10, convert=False):
    """Main visualization routine."""

    os.makedirs(output_dir, exist_ok=True)

    print(f"\nLoading trajectory...")
    traj = md.load(traj_file, top=top_file)
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    print(f"  Time range: {traj.time[0]:.1f} - {traj.time[-1]:.1f} ps")

    print(f"\nExporting frames...")
    export_frames(traj, output_dir)

    #print(f"\nSubsampling trajectory...")
    #subsample_trajectory(traj, output_dir, stride=stride)

    if convert:
        print(f"\nConverting formats...")
        convert_formats(traj, output_dir)

    print(f"\nPlotting time series...")
    plot_timeseries(traj, output_dir)

    #print(f"\nComputing contact map...")
    #plot_contact_map(traj, output_dir)

    print(f"\nGenerating VMD script...")
    generate_vmd_script(top_file, traj_file, output_dir)

    print(f"\n{'='*50}")
    print(f"Visualization outputs saved to: {output_dir}/")
    print(f"{'='*50}")


if __name__ == "__main__":
    parser = ArgumentParser(description="Visualize CALVADOS trajectory")
    parser.add_argument('--top', required=True, type=str,
                        help='Topology file (PDB)')
    parser.add_argument('--traj', required=True, type=str,
                        help='Trajectory file (DCD)')
    parser.add_argument('--output', default='viz', type=str,
                        help='Output directory')
    parser.add_argument('--stride', default=10, type=int,
                        help='Stride for subsampling')
    parser.add_argument('--convert', action='store_true',
                        help='Convert to XTC/XYZ formats')

    args = parser.parse_args()

    main(
        top_file=args.top,
        traj_file=args.traj,
        output_dir=args.output,
        stride=args.stride,
        convert=args.convert
    )
