"""
Analysis functions for CALVADOS trajectories.

Functions for computing structural properties:
- Radius of gyration (Rg)
- End-to-end distance (Ree)
- Contact maps
- RMSF
- Scaling exponent
"""

import numpy as np


def compute_rg(traj):
    """
    Compute radius of gyration for each frame.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory (from load_trajectory)

    Returns
    -------
    rg : numpy.ndarray
        Radius of gyration values in nanometers, shape (n_frames,)

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd')
    >>> rg = compute_rg(traj)
    >>> print(f"Mean Rg: {np.mean(rg):.2f} nm")
    """
    import mdtraj as md

    rg = md.compute_rg(traj)
    return rg


def compute_ree(traj, first_atom=0, last_atom=None):
    """
    Compute end-to-end distance for each frame.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory
    first_atom : int
        Index of first atom (default: 0)
    last_atom : int
        Index of last atom (default: last atom in trajectory)

    Returns
    -------
    ree : numpy.ndarray
        End-to-end distance in nanometers, shape (n_frames,)

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd')
    >>> ree = compute_ree(traj)
    >>> print(f"Mean Ree: {np.mean(ree):.2f} nm")
    """
    import mdtraj as md

    if last_atom is None:
        last_atom = traj.n_atoms - 1

    pairs = [[first_atom, last_atom]]
    distances = md.compute_distances(traj, pairs)

    return distances.flatten()


def compute_contact_map(traj, cutoff=0.8, skip_neighbors=2):
    """
    Compute average contact map over trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory
    cutoff : float
        Distance cutoff in nm (default: 0.8 nm)
    skip_neighbors : int
        Skip contacts between residues this close in sequence (default: 2)

    Returns
    -------
    contact_map : numpy.ndarray
        Contact probability matrix, shape (n_atoms, n_atoms)

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd', skip=100)
    >>> cmap = compute_contact_map(traj)
    >>> plot_contact_map(cmap, 'contacts.png')
    """
    import mdtraj as md

    n_atoms = traj.n_atoms

    # Create atom pairs (excluding neighbors)
    pairs = []
    for i in range(n_atoms):
        for j in range(i + skip_neighbors + 1, n_atoms):
            pairs.append([i, j])

    if len(pairs) == 0:
        return np.zeros((n_atoms, n_atoms))

    pairs = np.array(pairs)

    # Compute distances for all frames
    distances = md.compute_distances(traj, pairs)

    # Compute contact probability (soft cutoff)
    contacts = 0.5 - 0.5 * np.tanh((distances - cutoff) / 0.3)
    avg_contacts = np.mean(contacts, axis=0)

    # Reconstruct symmetric matrix
    contact_map = np.zeros((n_atoms, n_atoms))
    for idx, (i, j) in enumerate(pairs):
        contact_map[i, j] = avg_contacts[idx]
        contact_map[j, i] = avg_contacts[idx]

    return contact_map


def compute_rmsf(traj, reference='mean'):
    """
    Compute root mean square fluctuation per atom.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory
    reference : str or int
        Reference structure: 'mean' for average, or frame index

    Returns
    -------
    rmsf : numpy.ndarray
        RMSF values in nanometers, shape (n_atoms,)

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd', skip=100)
    >>> rmsf = compute_rmsf(traj)
    >>> plot_rmsf(rmsf, 'rmsf.png')
    """
    import mdtraj as md

    # Superpose to reference
    if reference == 'mean':
        # Use first frame as initial reference, then iterate
        traj = traj.superpose(traj, 0)
        mean_xyz = np.mean(traj.xyz, axis=0)
    else:
        traj = traj.superpose(traj, reference)
        mean_xyz = np.mean(traj.xyz, axis=0)

    # Compute RMSF
    deviations = traj.xyz - mean_xyz
    rmsf = np.sqrt(np.mean(np.sum(deviations**2, axis=2), axis=0))

    return rmsf


def compute_scaling_exponent(traj, min_separation=5):
    """
    Compute polymer scaling exponent (nu) from internal distances.

    For an ideal chain, nu = 0.5
    For a self-avoiding walk, nu ~ 0.588
    For a collapsed globule, nu ~ 0.33

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory
    min_separation : int
        Minimum sequence separation for fitting (default: 5)

    Returns
    -------
    result : dict
        Dictionary with:
        - nu: scaling exponent
        - nu_err: fitting error
        - r0: prefactor
        - separations: sequence separations used
        - distances: mean distances at each separation

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd', skip=100)
    >>> result = compute_scaling_exponent(traj)
    >>> print(f"Scaling exponent: {result['nu']:.3f} ± {result['nu_err']:.3f}")
    """
    import mdtraj as md
    from scipy.optimize import curve_fit

    n_atoms = traj.n_atoms

    # Compute all pairwise distances
    all_pairs = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            all_pairs.append([i, j])
    all_pairs = np.array(all_pairs)

    distances = md.compute_distances(traj, all_pairs)

    # Average over frames, then group by sequence separation
    mean_distances = np.mean(distances, axis=0)

    # Group by separation |i-j|
    separations = []
    mean_dist_by_sep = []

    for sep in range(1, n_atoms):
        indices = [idx for idx, (i, j) in enumerate(all_pairs) if j - i == sep]
        if len(indices) > 0:
            separations.append(sep)
            mean_dist_by_sep.append(np.mean(mean_distances[indices]))

    separations = np.array(separations)
    mean_dist_by_sep = np.array(mean_dist_by_sep)

    # Fit power law: R = r0 * N^nu
    def power_law(n, r0, nu):
        return r0 * n**nu

    # Fit only for separations >= min_separation
    mask = separations >= min_separation
    if np.sum(mask) < 3:
        return {'nu': np.nan, 'nu_err': np.nan, 'r0': np.nan,
                'separations': separations, 'distances': mean_dist_by_sep}

    try:
        popt, pcov = curve_fit(power_law, separations[mask], mean_dist_by_sep[mask],
                               p0=[0.5, 0.5])
        perr = np.sqrt(np.diag(pcov))

        result = {
            'nu': popt[1],
            'nu_err': perr[1],
            'r0': popt[0],
            'separations': separations,
            'distances': mean_dist_by_sep,
        }
    except Exception:
        result = {
            'nu': np.nan,
            'nu_err': np.nan,
            'r0': np.nan,
            'separations': separations,
            'distances': mean_dist_by_sep,
        }

    return result


def compute_all_metrics(traj, skip=0):
    """
    Compute all standard metrics at once.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory
    skip : int
        Frames to skip (equilibration)

    Returns
    -------
    metrics : dict
        Dictionary with all computed metrics:
        - rg_mean, rg_std, rg_values
        - ree_mean, ree_std, ree_values
        - nu, nu_err (scaling exponent)
        - contact_map
        - rmsf

    Example
    -------
    >>> traj = load_trajectory('top.pdb', 'simulation.dcd')
    >>> metrics = compute_all_metrics(traj, skip=100)
    >>> print(f"Rg = {metrics['rg_mean']:.2f} ± {metrics['rg_std']:.2f} nm")
    """
    if skip > 0:
        traj = traj[skip:]

    # Compute all properties
    rg = compute_rg(traj)
    ree = compute_ree(traj)
    cmap = compute_contact_map(traj)
    rmsf = compute_rmsf(traj)
    scaling = compute_scaling_exponent(traj)

    metrics = {
        # Rg
        'rg_mean': np.mean(rg),
        'rg_std': np.std(rg),
        'rg_values': rg,
        # Ree
        'ree_mean': np.mean(ree),
        'ree_std': np.std(ree),
        'ree_values': ree,
        # Scaling
        'nu': scaling['nu'],
        'nu_err': scaling['nu_err'],
        # Maps
        'contact_map': cmap,
        'rmsf': rmsf,
    }

    return metrics
