"""
Visualization functions for CALVADOS trajectories.

Functions for creating plots and exporting structures.
All plotting functions save to files (no interactive display).
"""

import numpy as np
import os


def plot_rg_distribution(rg_values, output_file, title=None):
    """
    Plot radius of gyration distribution.

    Parameters
    ----------
    rg_values : numpy.ndarray
        Rg values from compute_rg()
    output_file : str
        Path to save the plot (e.g., 'rg_distribution.png')
    title : str, optional
        Custom title for the plot

    Example
    -------
    >>> rg = compute_rg(traj)
    >>> plot_rg_distribution(rg, 'analysis/rg_distribution.png')
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.hist(rg_values, bins=50, density=True, alpha=0.7, color='steelblue',
            edgecolor='white', linewidth=0.5)

    mean_rg = np.mean(rg_values)
    std_rg = np.std(rg_values)

    ax.axvline(mean_rg, color='red', linestyle='--', linewidth=2,
               label=f'Mean: {mean_rg:.2f} nm')

    ax.set_xlabel('Radius of Gyration (nm)', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f'Rg = {mean_rg:.2f} ± {std_rg:.2f} nm', fontsize=14)

    ax.legend()
    plt.tight_layout()

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def plot_ree_distribution(ree_values, output_file, title=None):
    """
    Plot end-to-end distance distribution.

    Parameters
    ----------
    ree_values : numpy.ndarray
        Ree values from compute_ree()
    output_file : str
        Path to save the plot

    Example
    -------
    >>> ree = compute_ree(traj)
    >>> plot_ree_distribution(ree, 'analysis/ree_distribution.png')
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.hist(ree_values, bins=50, density=True, alpha=0.7, color='forestgreen',
            edgecolor='white', linewidth=0.5)

    mean_ree = np.mean(ree_values)
    std_ree = np.std(ree_values)

    ax.axvline(mean_ree, color='red', linestyle='--', linewidth=2,
               label=f'Mean: {mean_ree:.2f} nm')

    ax.set_xlabel('End-to-End Distance (nm)', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f'Ree = {mean_ree:.2f} ± {std_ree:.2f} nm', fontsize=14)

    ax.legend()
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def plot_timeseries(values, output_file, ylabel='Value', dt_ns=0.1, title=None):
    """
    Plot a time series of values.

    Parameters
    ----------
    values : numpy.ndarray
        Values to plot (one per frame)
    output_file : str
        Path to save the plot
    ylabel : str
        Label for y-axis
    dt_ns : float
        Time step in nanoseconds (default: 0.1 ns = 100 ps)
    title : str, optional
        Plot title

    Example
    -------
    >>> rg = compute_rg(traj)
    >>> plot_timeseries(rg, 'analysis/rg_timeseries.png', ylabel='Rg (nm)')
    """
    import matplotlib.pyplot as plt

    time_ns = np.arange(len(values)) * dt_ns

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(time_ns, values, lw=0.5, color='steelblue', alpha=0.8)
    ax.axhline(np.mean(values), color='red', linestyle='--', lw=1.5,
               label=f'Mean: {np.mean(values):.2f}')

    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)

    ax.legend()
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def plot_contact_map(contact_map, output_file, title='Contact Map', cmap='hot'):
    """
    Plot a contact map heatmap.

    Parameters
    ----------
    contact_map : numpy.ndarray
        Contact map from compute_contact_map()
    output_file : str
        Path to save the plot
    title : str
        Plot title
    cmap : str
        Matplotlib colormap name (default: 'hot')

    Example
    -------
    >>> cmap = compute_contact_map(traj)
    >>> plot_contact_map(cmap, 'analysis/contact_map.png')
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 7))

    im = ax.imshow(contact_map, cmap=cmap, origin='lower', aspect='equal',
                   vmin=0, vmax=1)

    ax.set_xlabel('Residue', fontsize=12)
    ax.set_ylabel('Residue', fontsize=12)
    ax.set_title(title, fontsize=14)

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Contact Probability', fontsize=11)

    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def plot_rmsf(rmsf_values, output_file, title='RMSF per Residue'):
    """
    Plot RMSF values.

    Parameters
    ----------
    rmsf_values : numpy.ndarray
        RMSF values from compute_rmsf()
    output_file : str
        Path to save the plot

    Example
    -------
    >>> rmsf = compute_rmsf(traj)
    >>> plot_rmsf(rmsf, 'analysis/rmsf.png')
    """
    import matplotlib.pyplot as plt

    residues = np.arange(1, len(rmsf_values) + 1)

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.fill_between(residues, 0, rmsf_values, alpha=0.4, color='purple')
    ax.plot(residues, rmsf_values, lw=1, color='purple')

    ax.axhline(np.mean(rmsf_values), color='red', linestyle='--', lw=1,
               label=f'Mean: {np.mean(rmsf_values):.3f} nm')

    ax.set_xlabel('Residue', fontsize=12)
    ax.set_ylabel('RMSF (nm)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(1, len(rmsf_values))
    ax.set_ylim(0, None)

    ax.legend()
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def plot_analysis_summary(metrics, output_file):
    """
    Create a 4-panel summary figure.

    Parameters
    ----------
    metrics : dict
        Output from compute_all_metrics()
    output_file : str
        Path to save the plot

    Example
    -------
    >>> metrics = compute_all_metrics(traj, skip=100)
    >>> plot_analysis_summary(metrics, 'analysis/summary.png')
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Rg distribution
    ax = axes[0, 0]
    rg = metrics['rg_values']
    ax.hist(rg, bins=50, density=True, alpha=0.7, color='steelblue')
    ax.axvline(np.mean(rg), color='red', linestyle='--',
               label=f'Mean: {np.mean(rg):.2f} nm')
    ax.set_xlabel('Rg (nm)')
    ax.set_ylabel('Probability')
    ax.set_title(f"Rg = {metrics['rg_mean']:.2f} ± {metrics['rg_std']:.2f} nm")
    ax.legend()

    # Ree distribution
    ax = axes[0, 1]
    ree = metrics['ree_values']
    ax.hist(ree, bins=50, density=True, alpha=0.7, color='forestgreen')
    ax.axvline(np.mean(ree), color='red', linestyle='--',
               label=f'Mean: {np.mean(ree):.2f} nm')
    ax.set_xlabel('Ree (nm)')
    ax.set_ylabel('Probability')
    ax.set_title(f"Ree = {metrics['ree_mean']:.2f} ± {metrics['ree_std']:.2f} nm")
    ax.legend()

    # Contact map
    ax = axes[1, 0]
    im = ax.imshow(metrics['contact_map'], cmap='hot', origin='lower', aspect='equal')
    ax.set_xlabel('Residue')
    ax.set_ylabel('Residue')
    ax.set_title('Contact Map')
    plt.colorbar(im, ax=ax, fraction=0.046)

    # RMSF
    ax = axes[1, 1]
    rmsf = metrics['rmsf']
    ax.fill_between(range(1, len(rmsf)+1), 0, rmsf, alpha=0.4, color='purple')
    ax.plot(range(1, len(rmsf)+1), rmsf, lw=1, color='purple')
    ax.set_xlabel('Residue')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title('Root Mean Square Fluctuation')

    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def export_frames(traj, output_dir, frames='key'):
    """
    Export specific frames as PDB files.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory
    output_dir : str
        Directory to save PDB files
    frames : str or list
        Which frames to export:
        - 'key': first, middle, last, and representative (default)
        - 'all': every frame (warning: many files!)
        - list of ints: specific frame indices

    Returns
    -------
    saved_files : list
        List of saved file paths

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd')
    >>> files = export_frames(traj, 'structures/')
    """
    import mdtraj as md

    os.makedirs(output_dir, exist_ok=True)
    saved_files = []

    if frames == 'key':
        # First frame
        traj[0].save_pdb(f'{output_dir}/first_frame.pdb')
        saved_files.append(f'{output_dir}/first_frame.pdb')

        # Middle frame
        mid = len(traj) // 2
        traj[mid].save_pdb(f'{output_dir}/middle_frame.pdb')
        saved_files.append(f'{output_dir}/middle_frame.pdb')

        # Last frame
        traj[-1].save_pdb(f'{output_dir}/last_frame.pdb')
        saved_files.append(f'{output_dir}/last_frame.pdb')

        # Representative frame (closest to mean Rg)
        rg = md.compute_rg(traj)
        rep_idx = np.argmin(np.abs(rg - np.mean(rg)))
        traj[rep_idx].save_pdb(f'{output_dir}/representative_frame.pdb')
        saved_files.append(f'{output_dir}/representative_frame.pdb')

    elif frames == 'all':
        for i in range(len(traj)):
            fname = f'{output_dir}/frame_{i:06d}.pdb'
            traj[i].save_pdb(fname)
            saved_files.append(fname)

    elif isinstance(frames, (list, tuple)):
        for i in frames:
            fname = f'{output_dir}/frame_{i:06d}.pdb'
            traj[i].save_pdb(fname)
            saved_files.append(fname)

    print(f"Exported {len(saved_files)} frames to {output_dir}/")
    return saved_files


def generate_vmd_script(topology_file, trajectory_file, output_file='visualize.vmd'):
    """
    Generate a VMD visualization script.

    Parameters
    ----------
    topology_file : str
        Path to PDB topology
    trajectory_file : str
        Path to DCD trajectory
    output_file : str
        Path to save VMD script

    Example
    -------
    >>> generate_vmd_script('top.pdb', 'simulation.dcd', 'viz/view.vmd')
    >>> # Then run: vmd -e viz/view.vmd
    """
    script = f'''# VMD visualization script
# Run with: vmd -e {output_file}

# Load structure and trajectory
mol new {topology_file} type pdb
mol addfile {trajectory_file} type dcd waitfor all

# Coarse-grained visualization style
mol modstyle 0 0 VDW 0.8 12
mol modcolor 0 0 ResType
mol modmaterial 0 0 AOChalky

# Center view
display resetview

# Playback settings
animate speed 0.9

puts "Loaded trajectory with [molinfo 0 get numframes] frames"
puts "Use animation controls or type 'animate forward' to play"
'''

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(script)

    print(f"Saved VMD script: {output_file}")
    print(f"  Run with: vmd -e {output_file}")
