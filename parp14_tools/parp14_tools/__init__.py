"""
PARP14 Tools - Reusable functions for CALVADOS simulation analysis

This package provides easy-to-use functions for:
- Loading and analyzing trajectories
- Visualizing simulation results
- Compiling multiple trajectories
- Extracting sequences from PDB files

Installation:
    cd parp14_tools
    pip install -e .

Usage:
    from parp14_tools import load_trajectory, compute_rg, plot_rg_distribution

    # Load trajectory
    traj = load_trajectory('top.pdb', 'trajectory.dcd')

    # Compute properties
    rg_values = compute_rg(traj)

    # Plot
    plot_rg_distribution(rg_values, 'output.png')
"""

# Import main functions so users can do:
# from parp14_tools import load_trajectory, compute_rg, etc.

from .loading import (
    load_trajectory,
    load_trajectory_mda,
    get_trajectory_info,
)

from .analysis import (
    compute_rg,
    compute_ree,
    compute_contact_map,
    compute_rmsf,
    compute_scaling_exponent,
    compute_all_metrics,
)

from .visualization import (
    plot_rg_distribution,
    plot_ree_distribution,
    plot_timeseries,
    plot_contact_map,
    plot_rmsf,
    plot_analysis_summary,
    export_frames,
    generate_vmd_script,
)

from .utils import (
    extract_sequence_from_pdb,
    write_fasta,
    block_error,
)

# Version
__version__ = '0.1.0'

# What shows up with "from parp14_tools import *"
__all__ = [
    # Loading
    'load_trajectory',
    'load_trajectory_mda',
    'get_trajectory_info',
    # Analysis
    'compute_rg',
    'compute_ree',
    'compute_contact_map',
    'compute_rmsf',
    'compute_scaling_exponent',
    'compute_all_metrics',
    # Visualization
    'plot_rg_distribution',
    'plot_ree_distribution',
    'plot_timeseries',
    'plot_contact_map',
    'plot_rmsf',
    'plot_analysis_summary',
    'export_frames',
    'generate_vmd_script',
    # Utilities
    'extract_sequence_from_pdb',
    'write_fasta',
    'block_error',
]
