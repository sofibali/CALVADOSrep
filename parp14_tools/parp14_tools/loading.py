"""
Loading functions for trajectories and structures.

Functions for loading DCD trajectories with PDB topology files.
Supports both MDTraj and MDAnalysis backends.
"""

import os


def load_trajectory(topology_file, trajectory_file, skip=0, stride=1):
    """
    Load a trajectory using MDTraj.

    Parameters
    ----------
    topology_file : str
        Path to PDB topology file (e.g., 'top.pdb')
    trajectory_file : str
        Path to DCD trajectory file (e.g., 'simulation.dcd')
    skip : int
        Number of initial frames to skip (default: 0)
    stride : int
        Load every Nth frame (default: 1 = all frames)

    Returns
    -------
    traj : mdtraj.Trajectory
        Loaded trajectory object

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd', skip=100)
    >>> print(f"Loaded {traj.n_frames} frames")
    """
    import mdtraj as md

    # Load trajectory
    traj = md.load(trajectory_file, top=topology_file, stride=stride)

    # Skip equilibration if requested
    if skip > 0:
        traj = traj[skip:]

    return traj


def load_trajectory_mda(topology_file, trajectory_file):
    """
    Load a trajectory using MDAnalysis.

    MDAnalysis is better for large trajectories (doesn't load all into memory).

    Parameters
    ----------
    topology_file : str
        Path to PDB topology file
    trajectory_file : str
        Path to DCD trajectory file

    Returns
    -------
    universe : MDAnalysis.Universe
        MDAnalysis Universe object

    Example
    -------
    >>> u = load_trajectory_mda('top.pdb', 'simulation.dcd')
    >>> ag = u.select_atoms('all')
    >>> for ts in u.trajectory[100:]:  # Skip first 100 frames
    ...     # Process each frame
    ...     pass
    """
    import MDAnalysis as mda

    universe = mda.Universe(topology_file, trajectory_file)
    return universe


def get_trajectory_info(topology_file, trajectory_file):
    """
    Get basic information about a trajectory.

    Parameters
    ----------
    topology_file : str
        Path to PDB topology file
    trajectory_file : str
        Path to DCD trajectory file

    Returns
    -------
    info : dict
        Dictionary with trajectory information:
        - n_frames: number of frames
        - n_atoms: number of atoms/residues
        - time_ns: total time in nanoseconds
        - dt_ps: time between frames in picoseconds

    Example
    -------
    >>> info = get_trajectory_info('top.pdb', 'simulation.dcd')
    >>> print(f"Trajectory has {info['n_frames']} frames over {info['time_ns']:.1f} ns")
    """
    import mdtraj as md

    # Load just metadata (not full trajectory)
    traj = md.load(trajectory_file, top=topology_file)

    info = {
        'n_frames': traj.n_frames,
        'n_atoms': traj.n_atoms,
        'n_residues': traj.n_residues,
        'time_ps': traj.time[-1] if len(traj.time) > 0 else 0,
        'time_ns': traj.time[-1] / 1000 if len(traj.time) > 0 else 0,
        'dt_ps': traj.time[1] - traj.time[0] if len(traj.time) > 1 else 0,
        'topology_file': os.path.abspath(topology_file),
        'trajectory_file': os.path.abspath(trajectory_file),
    }

    return info
