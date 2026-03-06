"""
Utility functions for CALVADOS workflows.

Helper functions for:
- Sequence extraction from PDB
- FASTA file handling
- Error estimation (block averaging)
"""

import numpy as np


def extract_sequence_from_pdb(pdb_file):
    """
    Extract amino acid sequence from a PDB file.

    Uses CALVADOS built-in function if available, otherwise uses MDAnalysis.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file

    Returns
    -------
    sequence : str
        One-letter amino acid sequence

    Example
    -------
    >>> seq = extract_sequence_from_pdb('structure.pdb')
    >>> print(f"Sequence: {seq[:50]}...")
    >>> print(f"Length: {len(seq)} residues")
    """
    try:
        # Try CALVADOS function first
        from calvados.sequence import seq_from_pdb
        sequence, _, _ = seq_from_pdb(pdb_file)
        return sequence
    except ImportError:
        pass

    # Fallback to MDAnalysis
    import MDAnalysis as mda
    from MDAnalysis.lib.util import convert_aa_code

    u = mda.Universe(pdb_file)

    # Get unique residues (by residue number)
    sequence = ''
    for residue in u.residues:
        resname = residue.resname
        try:
            one_letter = convert_aa_code(resname)
            sequence += one_letter
        except Exception:
            sequence += 'X'  # Unknown residue

    return sequence


def write_fasta(sequences, output_file, names=None):
    """
    Write sequences to a FASTA file.

    Parameters
    ----------
    sequences : str or list
        Single sequence or list of sequences
    output_file : str
        Path to output FASTA file
    names : str or list, optional
        Sequence name(s). If not provided, uses 'seq_1', 'seq_2', etc.

    Example
    -------
    >>> seq = extract_sequence_from_pdb('protein.pdb')
    >>> write_fasta(seq, 'protein.fasta', names='MyProtein')

    >>> seqs = [extract_sequence_from_pdb(f) for f in pdb_files]
    >>> write_fasta(seqs, 'all_proteins.fasta', names=protein_names)
    """
    # Handle single sequence
    if isinstance(sequences, str):
        sequences = [sequences]
    if names is None:
        names = [f'seq_{i+1}' for i in range(len(sequences))]
    elif isinstance(names, str):
        names = [names]

    with open(output_file, 'w') as f:
        for name, seq in zip(names, sequences):
            f.write(f'>{name}\n')
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(f'{seq[i:i+60]}\n')

    print(f"Wrote {len(sequences)} sequence(s) to {output_file}")


def block_error(data, n_blocks=10):
    """
    Estimate standard error using block averaging.

    Block averaging accounts for correlation in time series data,
    giving more accurate error estimates than simple std/sqrt(N).

    Parameters
    ----------
    data : numpy.ndarray
        Time series data (e.g., Rg values per frame)
    n_blocks : int
        Number of blocks to divide data into (default: 10)

    Returns
    -------
    result : dict
        Dictionary with:
        - mean: mean value
        - std: standard deviation
        - sem: simple standard error (std / sqrt(N))
        - block_sem: block-averaged standard error
        - n_blocks: number of blocks used

    Example
    -------
    >>> rg = compute_rg(traj)
    >>> error = block_error(rg)
    >>> print(f"Rg = {error['mean']:.3f} ± {error['block_sem']:.3f} nm")
    """
    data = np.asarray(data)
    n_points = len(data)

    # Simple statistics
    mean = np.mean(data)
    std = np.std(data)
    sem = std / np.sqrt(n_points)

    # Block averaging
    block_size = n_points // n_blocks
    if block_size < 2:
        # Not enough data for blocking, return simple SEM
        return {
            'mean': mean,
            'std': std,
            'sem': sem,
            'block_sem': sem,
            'n_blocks': 1,
        }

    # Compute block means
    block_means = []
    for i in range(n_blocks):
        start = i * block_size
        end = start + block_size
        block_means.append(np.mean(data[start:end]))

    block_means = np.array(block_means)
    block_sem = np.std(block_means) / np.sqrt(n_blocks)

    return {
        'mean': mean,
        'std': std,
        'sem': sem,
        'block_sem': block_sem,
        'n_blocks': n_blocks,
    }


def save_metrics_csv(metrics, output_file):
    """
    Save analysis metrics to a CSV file.

    Parameters
    ----------
    metrics : dict
        Dictionary of metrics (from compute_all_metrics or custom)
    output_file : str
        Path to output CSV file

    Example
    -------
    >>> metrics = compute_all_metrics(traj, skip=100)
    >>> save_metrics_csv(metrics, 'analysis/metrics.csv')
    """
    import pandas as pd
    import os

    # Filter out array values (can't save to CSV easily)
    scalar_metrics = {}
    for key, value in metrics.items():
        if isinstance(value, (int, float, np.integer, np.floating)):
            scalar_metrics[key] = value
        elif isinstance(value, np.ndarray) and value.ndim == 0:
            scalar_metrics[key] = float(value)

    df = pd.DataFrame([scalar_metrics])

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    df.to_csv(output_file, index=False)

    print(f"Saved metrics to {output_file}")
