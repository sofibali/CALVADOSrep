# CALVADOS Example Run Scripts

Ready-to-use Python scripts for common CALVADOS workflows.

## Scripts Overview

| Script | Purpose |
|--------|---------|
| `prepare_single_idr.py` | Prepare simulation for a single IDR |
| `prepare_multi_domain.py` | Prepare simulation for multi-domain protein with restraints |
| `batch_prepare.py` | Batch prepare multiple simulations from list or CSV |
| `analyze_trajectory.py` | Comprehensive trajectory analysis (Rg, Ree, contacts) |
| `visualize_trajectory.py` | Export frames, create plots, generate VMD scripts |
| `compile_trajectories.py` | Merge replicate trajectories into single file |
| `extract_fasta_from_pdb.py` | Extract FASTA sequences from PDB files |
| `slurm_array_template.sh` | SLURM array job template for HPC |

---

## Quick Usage Examples

### 1. Single IDR Simulation

```bash
# Prepare
python prepare_single_idr.py --name MyProtein --fasta input/sequences.fasta

# Run
python MyProtein/run.py --path MyProtein

# Analyze
python analyze_trajectory.py --top MyProtein/top.pdb --traj MyProtein/MyProtein.dcd --idr
```

### 2. Multi-Domain Protein

```bash
# Prepare (requires input/MyMDP.pdb and input/domains.yaml)
python prepare_multi_domain.py --name MyMDP --pdb_folder input --domains input/domains.yaml

# Run
python MyMDP/run.py --path MyMDP

# Analyze
python analyze_trajectory.py --top MyMDP/top.pdb --traj MyMDP/MyMDP.dcd
```

### 3. Batch Processing

```bash
# From list of names
python batch_prepare.py --names prot1 prot2 prot3 --fasta input/sequences.fasta --mode idr

# From CSV file
python batch_prepare.py --csv proteins.csv --mode mdp --output simulations/

# PARP14 library
python batch_prepare.py --parp14 compositions.csv --structures_dir input/structures --n_structures 25
```

### 4. Visualization

```bash
# Generate visualizations
python visualize_trajectory.py --top simulation/top.pdb --traj simulation/protein.dcd

# Open in VMD
vmd -e viz/visualize.vmd
```

### 5. Compile Replicates

```bash
# Compile 25 replicates, skip 10 ns equilibration
python compile_trajectories.py --input_dir simulations/comp_0001 --output compiled/comp_0001.dcd --skip_ns 10
```

### 6. Extract Sequences

```bash
# From single PDB
python extract_fasta_from_pdb.py --pdb input/protein.pdb --output sequences.fasta --stats

# From directory
python extract_fasta_from_pdb.py --pdb_dir input/structures/ --output all_sequences.fasta
```

### 7. HPC Job Submission

```bash
# Edit slurm_array_template.sh to set BASE_DIR and other options
# Then submit:
sbatch --array=1-100 slurm_array_template.sh
```

---

## Output Files

### From `analyze_trajectory.py`
- `analysis/rg_timeseries.npy` - Rg values per frame
- `analysis/ree_timeseries.npy` - Ree values per frame
- `analysis/contact_map.npy` - Average contact map
- `analysis/rmsf.npy` - RMSF per residue
- `analysis/metrics.csv` - Summary statistics
- `analysis/analysis_summary.png` - Combined plots

### From `visualize_trajectory.py`
- `viz/first_frame.pdb` - First frame structure
- `viz/last_frame.pdb` - Last frame structure
- `viz/representative_frame.pdb` - Frame closest to mean Rg
- `viz/subsampled.dcd` - Every 10th frame
- `viz/timeseries.png` - Rg/Ree time series
- `viz/contact_map.png` - Contact map heatmap
- `viz/visualize.vmd` - VMD script

---

## Dependencies

These scripts use:
- CALVADOS (calvados package)
- MDAnalysis
- MDTraj
- NumPy, Pandas
- Matplotlib

All are installed with CALVADOS:
```bash
pip install .
```
