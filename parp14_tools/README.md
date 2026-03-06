# PARP14 Tools

A simple Python package for analyzing and visualizing CALVADOS simulation trajectories.

## Installation

```bash
cd parp14_tools
pip install -e .
```

The `-e` flag installs in "editable" mode - any changes you make to the code take effect immediately.

## Quick Start

```python
from parp14_tools import (
    load_trajectory,
    compute_rg,
    compute_all_metrics,
    plot_rg_distribution,
    plot_analysis_summary,
)

# Load trajectory (skip first 100 frames for equilibration)
traj = load_trajectory('simulation/top.pdb', 'simulation/protein.dcd', skip=100)

# Compute radius of gyration
rg = compute_rg(traj)
print(f"Mean Rg: {rg.mean():.2f} nm")

# Create a plot
plot_rg_distribution(rg, 'analysis/rg.png')

# Or compute everything at once
metrics = compute_all_metrics(traj)
plot_analysis_summary(metrics, 'analysis/summary.png')
```

## Available Functions

### Loading
- `load_trajectory(top, traj)` - Load trajectory with MDTraj
- `load_trajectory_mda(top, traj)` - Load with MDAnalysis (for large files)
- `get_trajectory_info(top, traj)` - Get frame count, time, etc.

### Analysis
- `compute_rg(traj)` - Radius of gyration
- `compute_ree(traj)` - End-to-end distance
- `compute_contact_map(traj)` - Average contact map
- `compute_rmsf(traj)` - Root mean square fluctuation
- `compute_scaling_exponent(traj)` - Polymer scaling (nu)
- `compute_all_metrics(traj)` - All of the above at once

### Visualization
- `plot_rg_distribution(rg, output)` - Rg histogram
- `plot_ree_distribution(ree, output)` - Ree histogram
- `plot_timeseries(values, output)` - Time series plot
- `plot_contact_map(cmap, output)` - Contact map heatmap
- `plot_rmsf(rmsf, output)` - RMSF bar plot
- `plot_analysis_summary(metrics, output)` - 4-panel summary
- `export_frames(traj, output_dir)` - Save PDB frames
- `generate_vmd_script(top, traj, output)` - VMD script

### Utilities
- `extract_sequence_from_pdb(pdb)` - Get amino acid sequence
- `write_fasta(sequences, output)` - Write FASTA file
- `block_error(data)` - Block averaging for error estimation

## Examples

See `examples/` folder for complete workflows:
- `analyze_single.py` - Analyze one trajectory
- `analyze_library.py` - Batch analyze many trajectories
- `compare_variants.py` - Compare multiple variants

## Package Structure

```
parp14_tools/
├── setup.py              # Installation script
├── README.md             # This file
├── parp14_tools/         # The actual package
│   ├── __init__.py       # Makes functions importable
│   ├── loading.py        # Trajectory loading
│   ├── analysis.py       # Compute properties
│   ├── visualization.py  # Create plots
│   └── utils.py          # Helper functions
└── examples/             # Example scripts
```
