# PARP14 Structure-to-Dynamics Pipeline

## Overview

This documentation covers the complete computational pipeline for studying PARP14 domain variants, from initial structure prediction through coarse-grained molecular dynamics simulations.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        PARP14 COMPUTATIONAL PIPELINE                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  PHASE 1: Structure Generation & Analysis    (/home/sbali/parp14)   │   │
│  │                                                                      │   │
│  │  • Generate domain combinations from full PARP14 sequence           │   │
│  │  • Predict structures with AlphaFold3 (1,285 combinations × 5 seeds)│   │
│  │  • Analyze confidence (pLDDT), distances, and solvent accessibility │   │
│  │  • Select high-confidence structures for MD simulations             │   │
│  │  • Generate CALVADOS-compatible input files                          │   │
│  └──────────────────────────────┬──────────────────────────────────────┘   │
│                                 │                                           │
│                                 ▼                                           │
│                     ┌───────────────────────┐                               │
│                     │   PDB structures +    │                               │
│                     │   domain definitions  │                               │
│                     │   (YAML configs)      │                               │
│                     └───────────┬───────────┘                               │
│                                 │                                           │
│                                 ▼                                           │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  PHASE 2: Coarse-Grained MD Simulations    (/home/sbali/CALVADOS)   │   │
│  │                                                                      │   │
│  │  • Run CALVADOS simulations on selected domain variants             │   │
│  │  • 25 replicates per composition (conformational sampling)          │   │
│  │  • Compile 1 μs effective trajectories per composition              │   │
│  │  • Analyze structural dynamics: Rg, Ree, contacts, RMSF             │   │
│  │  • Compare domain deletion effects on structural ensemble           │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Project Locations

| Phase | Directory | Purpose |
|-------|-----------|---------|
| **Phase 1** | `/home/sbali/parp14` | AlphaFold3 structure prediction and preliminary analysis |
| **Phase 2** | `/home/sbali/CALVADOS` | Coarse-grained MD simulations and dynamics analysis |

---

## PARP14 Protein Architecture

Full-length PARP14 is an 1801-residue multi-domain protein involved in ADP-ribosylation.

| Domain | Residues | Function |
|--------|----------|----------|
| RRM1-3 | 1-314 | RNA Recognition Motifs |
| KH1-KH6 | 315-737 | K-Homology domains (RNA binding) |
| KH7a | 738-789 | K-Homology domain |
| MD1 | 790-981 | Macrodomain 1 (ADP-ribose binding) |
| MD2 | 1005-1193 | Macrodomain 2 |
| MD3 | 1207-1388 | Macrodomain 3 |
| KHb-KH8 | 1389-1533 | K-Homology domains |
| WWE | 1534-1602 | WWE domain |
| ART | 1603-1801 | ADP-ribosyltransferase (catalytic) |

**IDR Linkers**: Disordered regions connecting structured domains at positions 89-149, 302-790, 979-1002, 1191-1215, 1388-1522.

---

## Phase 1: Structure Generation (parp14)

### Purpose

Generate and analyze AlphaFold3 structure predictions for all biologically relevant domain combinations of PARP14 to:

1. Understand structural confidence across domain compositions
2. Identify stable domain combinations for MD simulations
3. Prepare input structures for CALVADOS simulations

### Pipeline Scripts

| Stage | Script | Output |
|-------|--------|--------|
| 1 | `01_generate_fasta.py` | Domain combination FASTA files |
| 2 | `02_generate_af3_inputs.py` | AlphaFold3 JSON inputs |
| 3 | `03_run_alphafold.sh` | Predicted structures (mmCIF) |
| 4 | `04_distance_analysis.py` | Inter-domain distance metrics |
| 5 | `05_pLDDT_analysis.py` | Confidence scores |
| 6 | `06_sasa_analysis.py` | Solvent accessibility |
| 7 | `07_generate_calvados_inputs.py` | **CALVADOS-ready configs** |

### Key Outputs for Phase 2

The script `07_generate_calvados_inputs.py` generates CALVADOS-compatible files:

```
calvados_simulations/
├── rrm1_rrm2_rrm3/
│   ├── config.yaml           # Simulation parameters
│   ├── components.yaml       # Molecule definitions
│   ├── domains.yaml          # Domain boundaries for restraints
│   └── input/
│       └── structure.pdb     # Backbone-only PDB from AlphaFold3
└── ... (484 domain combinations)
```

### Structure Selection Criteria

Before passing to CALVADOS, structures are filtered by:

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| Mean pLDDT | >70 | Reasonable structural confidence |
| No clashes | True | Valid backbone geometry |
| PTM score | >0.5 | Good predicted fold quality |

---

## Phase 2: MD Simulations (CALVADOS)

### Purpose

Run coarse-grained molecular dynamics simulations to characterize the structural dynamics of PARP14 domain variants:

1. Sample conformational ensembles for each domain composition
2. Compare flexibility, compaction, and domain contacts across variants
3. Understand how domain deletion affects structural dynamics

### Using parp14 Structures in CALVADOS

**Option A: Use pre-generated CALVADOS inputs**

```bash
# Structures are already prepared in parp14/calvados_simulations/
cd /home/sbali/parp14/calvados_simulations/rrm1_rrm2_rrm3

# Run simulation directly
python /home/sbali/CALVADOS/calvados/sim.py --path .
```

**Option B: Manual CALVADOS setup from AlphaFold3 structures**

```python
from calvados.cfg import Config, Components
import os

# Point to parp14 AlphaFold3 outputs
parp14_dir = '/home/sbali/parp14'
structure_name = 'rrm1_rrm2_rrm3'

config = Config(
    sysname = f'parp14_{structure_name}',
    box = [100, 100, 100],
    temp = 293,
    ionic = 0.19,
    pH = 7.0,
    steps = 5000000,
    wfreq = 10000,
    platform = 'CUDA',
)

components = Components(
    fresidues = '/home/sbali/CALVADOS/calvados/data/residues_CALVADOS3.csv',
    fdomains = f'{parp14_dir}/calvados_simulations/{structure_name}/domains.yaml',
    pdb_folder = f'{parp14_dir}/calvados_simulations/{structure_name}/input',
    restraint_type = 'harmonic',
    k_harmonic = 700.,
    use_com = True,
)

components.add(name='structure', restraint=True)

sim_path = f'./simulations/{structure_name}'
os.makedirs(sim_path, exist_ok=True)
config.write(sim_path)
components.write(sim_path)
```

### Simulation Protocol

| Parameter | Value | Notes |
|-----------|-------|-------|
| Simulation length | 50 ns | 5,000,000 steps |
| Replicates per composition | 25 | Different starting seeds |
| Equilibration (discarded) | 10 ns | First 100 frames |
| Production (analyzed) | 40 ns | Last 400 frames |
| Compiled trajectory | 1 μs | 25 × 40 ns |

### Analysis Metrics

For each domain composition:

| Metric | Description | CALVADOS Function |
|--------|-------------|-------------------|
| Rg | Radius of gyration | `calc_rg()` |
| Ree | End-to-end distance | `calc_ete()` |
| Contact map | Inter-residue contacts | `calc_cmap()` |
| Scaling exponent | Polymer ν parameter | `fit_scaling_exp()` |
| Domain distances | Pairwise domain COM | Custom analysis |

---

## Complete Workflow Example

### Step 1: Generate structures (parp14)

```bash
cd /home/sbali/parp14/pipeline

# Generate domain combinations
python 01_generate_fasta.py --fasta ../PARP14.fasta --domains ../Domain_boundries.csv --output ../fasta_files/

# Create AlphaFold3 inputs
python 02_generate_af3_inputs.py --fasta-dir ../fasta_files/ --output ../alphafold_inputs/

# Run AlphaFold3 (takes hours/days)
bash 03_run_alphafold.sh

# Analyze structures
python 05_pLDDT_analysis.py --input ../alphafold_outputs/ --output-residue ../pLDDT_per_residue.csv

# Generate CALVADOS inputs for high-confidence structures
python 07_generate_calvados_inputs.py --input ../alphafold_outputs/ --output ../calvados_simulations/
```

### Step 2: Run CALVADOS simulations

```bash
cd /home/sbali/CALVADOS

# Activate environment
conda activate calvados

# Run simulation for a specific domain composition
cd /home/sbali/parp14/calvados_simulations/rrm1_rrm2_rrm3
python run.py --path .

# Or use the library preparation script
python /home/sbali/CALVADOS/scripts/prepare_variant.py 1 1 /home/sbali/parp14/calvados_simulations
```

### Step 3: Analyze trajectories

```python
import MDAnalysis as mda
import pandas as pd
from calvados.analysis import calc_rg, calc_ete

# Load compiled trajectory
u = mda.Universe('trajectories/rrm1_rrm2_rrm3_top.pdb',
                 'trajectories/rrm1_rrm2_rrm3.dcd')
ag = u.select_atoms('all')

# Load residue parameters
residues = pd.read_csv('input/residues_CALVADOS3.csv').set_index('three')

# Calculate metrics
rgs = calc_rg(u, ag, ag.resnames.tolist(), residues)
rees, ree_mean, ree_sem = calc_ete(u, ag)

print(f"Rg = {np.mean(rgs):.2f} ± {np.std(rgs):.2f} nm")
print(f"Ree = {ree_mean:.2f} ± {ree_sem:.2f} nm")
```

---

## Data Flow Summary

```
Phase 1 (parp14)                         Phase 2 (CALVADOS)
================                         ==================

PARP14.fasta
     │
     ▼
Domain_boundries.csv ──► 01_generate_fasta.py
     │
     ▼
domain_combinations/*.fasta ──► 02_generate_af3_inputs.py
     │
     ▼
alphafold_inputs/*.json ──► AlphaFold3
     │
     ▼
alphafold_outputs/
├── model.cif
├── confidences.json
└── summary_confidences.json
     │
     ├──► 04_distance_analysis.py ──► distance_analysis/
     ├──► 05_pLDDT_analysis.py ──► pLDDT_per_residue.csv
     ├──► 06_sasa_analysis.py ──► sasa_per_residue.csv
     │
     └──► 07_generate_calvados_inputs.py
               │
               ▼
     calvados_simulations/
     ├── config.yaml        ─────────►  CALVADOS Sim
     ├── components.yaml                      │
     ├── domains.yaml                         ▼
     └── input/*.pdb                    Trajectories (DCD)
                                              │
                                              ▼
                                        Analysis
                                        ├── Rg, Ree
                                        ├── Contact maps
                                        ├── Scaling exponent
                                        └── Domain dynamics
```

---

## Documentation Index

### Phase 1: Structure Generation (parp14)

| Document | Description |
|----------|-------------|
| [PARP14_ARCHITECTURE.md](PARP14_ARCHITECTURE.md) | Technical reference for pipeline stages and data formats |
| [PARP14_QUICKSTART.md](PARP14_QUICKSTART.md) | Step-by-step guide for running the structure pipeline |
| [PARP14_ANALYSIS_GUIDE.md](PARP14_ANALYSIS_GUIDE.md) | Detailed analysis and visualization examples |

### Phase 2: MD Simulations (CALVADOS)

| Document | Description |
|----------|-------------|
| [CALVADOS_ARCHITECTURE.md](../../CALVADOS/docs/CALVADOS_ARCHITECTURE.md) | CALVADOS classes, force fields, and analysis tools |
| [CALVADOS_QUICKSTART.md](../../CALVADOS/docs/CALVADOS_QUICKSTART.md) | Running and analyzing CALVADOS simulations |
| [PARP14_PROJECT.md](../../CALVADOS/docs/PARP14_PROJECT.md) | PARP14-specific simulation protocol and HPC scripts |

---

## Resource Requirements

### Phase 1 (Structure Generation)

| Resource | Requirement |
|----------|-------------|
| GPU | NVIDIA with ≥16 GB VRAM (AlphaFold3) |
| Storage | ~615 GB total |
| Time | ~200-500 GPU-hours |

### Phase 2 (MD Simulations)

| Resource | Requirement |
|----------|-------------|
| GPU | NVIDIA CUDA-capable |
| Storage | ~1.5-2 TB (all trajectories) |
| Time | ~12,500-37,500 GPU-hours (full library) |

---

## Key Files Reference

| File | Location | Purpose |
|------|----------|---------|
| `PARP14.fasta` | `/home/sbali/parp14/` | Full protein sequence |
| `Domain_boundries.csv` | `/home/sbali/parp14/` | Domain definitions |
| `pLDDT_per_residue.csv` | `/home/sbali/parp14/` | AlphaFold3 confidence scores |
| `calvados_simulations/` | `/home/sbali/parp14/` | CALVADOS-ready inputs |
| `residues_CALVADOS3.csv` | `/home/sbali/CALVADOS/calvados/data/` | Force field parameters |
| `examples/PARP14_MDP/` | `/home/sbali/CALVADOS/` | Reference PARP14 simulation |
