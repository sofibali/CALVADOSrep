# CALVADOS Quickstart Guide

A step-by-step guide to set up, run, and analyze coarse-grained protein simulations.

---

## 1. Installation

### Create conda environment

```bash
conda create -n calvados python=3.10
conda activate calvados
```

### (Optional) For GPU support

```bash
conda install -c conda-forge openmm=8.2.0 cudatoolkit=11.8
```

### Install CALVADOS

```bash
git clone https://github.com/KULL-Centre/CALVADOS.git
cd CALVADOS
pip install .
```

### Verify installation

```bash
python -m pytest
```

---

## 2. Prepare Input Files

### Required files

| Simulation Type | Required Files |
|-----------------|----------------|
| IDR (disordered) | `sequences.fasta`, `residues.csv` |
| MDP (structured) | `protein.pdb`, `domains.yaml`, `residues.csv` |
| Mixed IDR + MDP | All of the above |

### Create FASTA from PDB (if needed)

```python
from calvados.sequence import seq_from_pdb, write_fasta, record_from_seq

seq, _, _ = seq_from_pdb('input/protein.pdb')
record = record_from_seq(seq, 'ProteinName')
write_fasta([record], 'input/sequences.fasta')
```

### domains.yaml format

```yaml
ProteinName:
  - [6, 82]      # Domain 1: residues 6-82
  - [95, 172]    # Domain 2: residues 95-172
```

---

## 3. Create prepare.py Script

```python
import os
from calvados.cfg import Config, Components

cwd = os.getcwd()

# Configure simulation
config = Config(
    sysname = 'my_simulation',
    box = [30, 30, 30],       # Box size in nm
    temp = 293,                # Temperature in K
    ionic = 0.15,              # Ionic strength in M
    pH = 7.0,
    steps = 1000000,           # Total steps (1 step = 0.01 ps)
    wfreq = 1000,              # Save every N steps
    platform = 'CPU',          # 'CPU' or 'CUDA'
    topol = 'random',          # Molecule placement
)

# Define components
components = Components(
    fresidues = f'{cwd}/input/residues.csv',
    ffasta = f'{cwd}/input/sequences.fasta',
    fdomains = f'{cwd}/input/domains.yaml',
    pdb_folder = f'{cwd}/input',
)

# Add molecules
components.add(name='MyIDR', restraint=False)           # Disordered protein
components.add(name='MyMDP', restraint=True)            # Structured protein

# Create output directory and write files
path = f'{cwd}/simulation'
os.makedirs(path, exist_ok=True)

config.write(path)
components.write(path)
```

---

## 4. Run Simulation

### Step 1: Prepare system

```bash
python prepare.py
```

This creates:
```
simulation/
├── config.yaml
├── components.yaml
└── run.py
```

### Step 2: Execute simulation

```bash
python simulation/run.py --path simulation
```

### Output files

| File | Description |
|------|-------------|
| `top.pdb` | Coarse-grained topology |
| `my_simulation.dcd` | Trajectory |
| `my_simulation.log` | Energy log |
| `restart.chk` | Checkpoint (for restart) |
| `checkpoint.pdb` | Final structure |

---

## 5. Restart a Simulation

Simulations automatically restart from checkpoint:

```bash
# Just run again - will continue from restart.chk
python simulation/run.py --path simulation
```

---

## 6. Analyze Results

### Basic analysis

```python
import numpy as np
import pandas as pd
import MDAnalysis as mda
from calvados.analysis import calc_rg, calc_ete, calc_cmap, fit_scaling_exp

# Load trajectory
u = mda.Universe('simulation/top.pdb', 'simulation/my_simulation.dcd')
ag = u.select_atoms('all')

# Load residue parameters
residues = pd.read_csv('input/residues_CALVADOS3.csv').set_index('three')

# Radius of gyration
rgs = calc_rg(u, ag, ag.resnames.tolist(), residues, start=100)
print(f"Rg = {np.mean(rgs):.2f} ± {np.std(rgs):.2f} nm")

# End-to-end distance
rees, ree_mean, ree_sem = calc_ete(u, ag, start=100)
print(f"Ree = {ree_mean:.2f} ± {ree_sem:.2f} nm")

# Scaling exponent (for IDRs)
ij, dij, r0, nu, nu_err = fit_scaling_exp(u, ag, start=100)
print(f"ν = {nu:.3f} ± {nu_err:.3f}")

# Contact map
cmap = calc_cmap(ag, ag, cutoff=1.0)
np.save('contact_map.npy', cmap)
```

### Center-of-mass trajectory

```python
from calvados.analysis import calc_com_traj

calc_com_traj(
    path='simulation',
    sysname='my_simulation',
    output_path='analysis',
    residues_file='input/residues.csv',
    chainid_dict={'MyIDR': 0, 'MyMDP': 1},
    start=100
)
```

### Block analysis for error estimation

```python
from calvados.BLOCKING.main import BlockAnalysis

block = BlockAnalysis(rgs)
block.SEM()
print(f"Rg = {np.mean(rgs):.2f} ± {block.sem:.2f} nm")
```

---

## 7. Common Simulation Types

### Single IDR

```python
config = Config(
    sysname='single_IDR',
    box=[50, 50, 50],
    steps=10000000,
)

components = Components(
    fresidues=f'{cwd}/input/residues.csv',
    ffasta=f'{cwd}/input/sequences.fasta',
)

components.add(name='MyProtein', restraint=False)
```

### Multi-domain protein with restraints

```python
components = Components(
    fresidues=f'{cwd}/input/residues.csv',
    fdomains=f'{cwd}/input/domains.yaml',
    pdb_folder=f'{cwd}/input',
    restraint_type='harmonic',
    k_harmonic=700.,
)

components.add(name='MyMDP', restraint=True)
```

### Two interacting proteins

```python
config = Config(
    sysname='two_proteins',
    box=[30, 30, 30],
    topol='random',
)

components.add(name='Protein1', restraint=False)
components.add(name='Protein2', restraint=True)
```

### Phase separation (slab geometry)

```python
config = Config(
    sysname='slab_simulation',
    box=[20, 20, 100],
    topol='slab',
    slab_eq=True,
    steps_eq=100000,
)

components = Components(nmol=100)  # 100 copies
components.add(name='MyIDR', restraint=False)
```

---

## 8. GPU Acceleration

```python
config = Config(
    platform='CUDA',
    gpu_id=0,  # GPU device ID
)
```

---

## 9. HPC Job Submission

```python
from calvados.cfg import Job

job = Job(
    batch_sys='SLURM',
    envname='calvados',
)

job.write(path, config, components, name='job.sh')
job.submit(path)
```

---

## 10. Quick Reference

### Simulation time conversion

- 1 step = 0.01 ps = 10 fs
- 1,000,000 steps = 10 ns
- `wfreq=1000` saves every 10 ps

### Typical box sizes

| System | Box Size (nm) |
|--------|---------------|
| Single protein | 50 × 50 × 50 |
| Two proteins | 30 × 30 × 30 |
| Phase separation | 20 × 20 × 100 |

### Key parameters to tune

| Parameter | Effect |
|-----------|--------|
| `ionic` | Electrostatic screening (higher = weaker) |
| `k_harmonic` | Restraint stiffness (higher = more rigid) |
| `cutoff_restr` | Max distance for restraints |
| `temp` | Temperature (affects dynamics and phase behavior) |

---

## 11. Visualizing Trajectories

CALVADOS outputs DCD trajectory files with PDB topology. Several tools can visualize these coarse-grained simulations.

### Visualization Tools Comparison

| Tool | Best For | Install |
|------|----------|---------|
| **VMD** | Interactive 3D, movies, publication figures | [Download](https://www.ks.uiuc.edu/Research/vmd/) |
| **PyMOL** | High-quality static images | `conda install -c conda-forge pymol-open-source` |
| **NGLView** | Jupyter notebook integration | `pip install nglview` |
| **MDTraj** | Python scripting, format conversion | Included with CALVADOS |
| **OVITO** | Coarse-grained/particle systems | [Download](https://www.ovito.org/) |

### VMD (Recommended for trajectories)

VMD natively reads DCD files and is ideal for coarse-grained simulations.

```bash
# Open trajectory in VMD
vmd simulation/top.pdb simulation/my_simulation.dcd
```

**VMD visualization tips for coarse-grained:**
```tcl
# In VMD TkConsole:
mol modstyle 0 0 VDW 1.0 12       # Show as spheres
mol modcolor 0 0 ResType          # Color by residue type
mol modmaterial 0 0 AOChalky      # Nice material for CG

# Or show as connected beads:
mol modstyle 0 0 Licorice 0.3 12 12
```

### NGLView (Jupyter notebooks)

```python
import nglview as nv
import mdtraj as md

# Load trajectory
traj = md.load('simulation/my_simulation.dcd', top='simulation/top.pdb')

# Create interactive view
view = nv.show_mdtraj(traj)
view.add_representation('spacefill', selection='all', color='residueindex')
view.center()
view
```

### MDTraj (Python scripting)

```python
import mdtraj as md

# Load trajectory
traj = md.load('simulation/my_simulation.dcd', top='simulation/top.pdb')

# Basic info
print(f"Frames: {traj.n_frames}")
print(f"Atoms: {traj.n_atoms}")
print(f"Time: {traj.time[-1]} ps")

# Save specific frames as PDB
traj[0].save_pdb('first_frame.pdb')
traj[-1].save_pdb('last_frame.pdb')

# Save subsampled trajectory
traj[::10].save_dcd('subsampled.dcd')  # Every 10th frame

# Convert to other formats
traj.save_xtc('trajectory.xtc')        # GROMACS format
traj.save_xyz('trajectory.xyz')        # XYZ format
```

### MDAnalysis (Alternative Python approach)

```python
import MDAnalysis as mda

# Load trajectory
u = mda.Universe('simulation/top.pdb', 'simulation/my_simulation.dcd')

# Write specific frames
u.trajectory[0]  # Go to first frame
u.atoms.write('first_frame.pdb')

u.trajectory[-1]  # Go to last frame
u.atoms.write('last_frame.pdb')

# Write subsampled trajectory
with mda.Writer('subsampled.dcd', u.atoms.n_atoms) as W:
    for ts in u.trajectory[::10]:  # Every 10th frame
        W.write(u.atoms)

# Center trajectory in box
from MDAnalysis import transformations
ag = u.atoms
transform = transformations.center_in_box(ag, center='geometry')
u.trajectory.add_transformations(transform)
```

### PyMOL (Static images)

```python
# In PyMOL or via pymol script
import pymol
from pymol import cmd

cmd.load('simulation/top.pdb', 'protein')
cmd.load_traj('simulation/my_simulation.dcd', 'protein')

# Visualization settings for coarse-grained
cmd.show('spheres', 'all')
cmd.set('sphere_scale', 0.4)
cmd.spectrum('count', 'rainbow', 'all')

# Navigate frames
cmd.frame(1)    # First frame
cmd.frame(100)  # Frame 100

# Save image
cmd.ray(1920, 1080)
cmd.png('snapshot.png', dpi=300)
```

### Create Movie with VMD

```tcl
# VMD script: make_movie.tcl
# Usage: vmd -e make_movie.tcl

mol new simulation/top.pdb
mol addfile simulation/my_simulation.dcd waitfor all

mol modstyle 0 0 VDW 0.8 12
mol modcolor 0 0 ResType

# Rotate and render frames
set nframes [molinfo 0 get numframes]
for {set i 0} {$i < $nframes} {incr i 10} {
    animate goto $i
    render TachyonInternal frame_[format %04d $i].tga
}

# Then use ffmpeg to combine:
# ffmpeg -framerate 30 -i frame_%04d.tga -c:v libx264 -pix_fmt yuv420p movie.mp4
```

### Quick Visualization Script

```python
#!/usr/bin/env python
"""visualize_trajectory.py - Quick trajectory visualization"""

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

def quick_viz(top_file, traj_file, output_prefix='viz'):
    """Generate quick visualization of trajectory."""

    traj = md.load(traj_file, top=top_file)

    # 1. Save first/last frames
    traj[0].save_pdb(f'{output_prefix}_first.pdb')
    traj[-1].save_pdb(f'{output_prefix}_last.pdb')

    # 2. Compute and plot Rg over time
    rg = md.compute_rg(traj)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(traj.time / 1000, rg, lw=0.5)  # Time in ns
    axes[0].set_xlabel('Time (ns)')
    axes[0].set_ylabel('Rg (nm)')
    axes[0].set_title('Radius of Gyration')

    axes[1].hist(rg, bins=50, density=True, alpha=0.7)
    axes[1].set_xlabel('Rg (nm)')
    axes[1].set_ylabel('Probability')
    axes[1].set_title(f'Rg Distribution: {np.mean(rg):.2f} ± {np.std(rg):.2f} nm')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_rg.png', dpi=150)
    plt.close()

    # 3. Contact map (last 10% of trajectory)
    n_skip = int(0.9 * traj.n_frames)
    contacts = md.compute_contacts(traj[n_skip:], scheme='closest-heavy')
    avg_distances = np.mean(contacts[0], axis=0)

    print(f"Saved: {output_prefix}_first.pdb, {output_prefix}_last.pdb, {output_prefix}_rg.png")

    return traj

if __name__ == "__main__":
    import sys
    top = sys.argv[1] if len(sys.argv) > 1 else 'top.pdb'
    traj = sys.argv[2] if len(sys.argv) > 2 else 'trajectory.dcd'
    quick_viz(top, traj)
```

---

## 13. Troubleshooting

### Simulation crashes with overlap

- Increase box size
- Use `topol='random'` for better initial placement
- Run energy minimization (automatic)

### Slow performance

- Use `platform='CUDA'` for GPU
- Increase `threads` for CPU
- Reduce `wfreq` (save less frequently)

### Restart not working

- Check `restart.chk` exists
- Verify `restart='checkpoint'` in config
- Ensure `.dcd` file exists to append to

---

## 14. Examples

Run example simulations:

```bash
cd examples/single_IDR
python prepare.py --name FUSRGG3
python FUSRGG3/run.py --path FUSRGG3
```

```bash
cd examples/two_IDR_MDP
python prepare.py --name_1 Tau35 --name_2 TIA1
python Tau35_TIA1/run.py --path Tau35_TIA1
```

---

---

## 15. Example Run Scripts

Ready-to-use scripts are available in `example_runscripts/`:

| Script | Purpose |
|--------|---------|
| `prepare_single_idr.py` | Prepare single IDR simulation |
| `prepare_multi_domain.py` | Prepare MDP with restraints |
| `batch_prepare.py` | Batch prepare multiple simulations |
| `analyze_trajectory.py` | Comprehensive analysis (Rg, Ree, contacts) |
| `visualize_trajectory.py` | Export frames, plots, VMD scripts |
| `compile_trajectories.py` | Merge replicate trajectories |
| `extract_fasta_from_pdb.py` | Extract FASTA from PDB |
| `slurm_array_template.sh` | SLURM array job template |

```bash
# Example: Quick analysis workflow
python example_runscripts/analyze_trajectory.py \
    --top simulation/top.pdb \
    --traj simulation/protein.dcd \
    --output analysis/ \
    --idr

python example_runscripts/visualize_trajectory.py \
    --top simulation/top.pdb \
    --traj simulation/protein.dcd \
    --output viz/
```

See `example_runscripts/README.md` for detailed usage.

---

## References

- Repository: https://github.com/KULL-Centre/CALVADOS
- Documentation: https://doi.org/10.48550/arXiv.2504.10408
- Issues: https://github.com/KULL-Centre/CALVADOS/issues
