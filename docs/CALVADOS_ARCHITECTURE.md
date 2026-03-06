# CALVADOS Architecture

Coarse-grained implicit-solvent simulations of biomolecules using the OpenMM framework.

## Repository Structure

```
CALVADOS/
├── calvados/                 # Main Python package
│   ├── __init__.py
│   ├── cfg.py               # Config, Components, Job classes
│   ├── sim.py               # Sim class - simulation engine
│   ├── build.py             # System building utilities
│   ├── components.py        # Molecule type definitions
│   ├── interactions.py      # Force field implementations
│   ├── sequence.py          # FASTA/PDB sequence utilities
│   ├── analysis.py          # Post-simulation analysis tools
│   ├── utilities.py         # Helper functions
│   ├── BLOCKING/            # Block analysis for error estimation
│   └── data/                # Default configurations & templates
│       ├── default_config.yaml
│       ├── default_component.yaml
│       ├── default_job.yaml
│       ├── residues.csv
│       └── templates/       # HPC job templates (SLURM, PBS)
├── examples/                 # 18+ example simulations
├── tests/                    # Unit tests
└── docs/                     # Documentation
```

---

## Core Classes

### 1. Config (`calvados/cfg.py`)

Manages simulation parameters. All settings have sensible defaults.

**Key Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sysname` | `'default_simulation'` | Simulation name |
| `box` | `[L, L, L]` | Box dimensions in nm |
| `temp` | `293` | Temperature in K |
| `ionic` | `0.15` | Ionic strength in M |
| `pH` | `7.0` | pH value |
| `steps` | `100000000` | Total simulation steps |
| `wfreq` | `100000` | Trajectory write frequency |
| `platform` | `'CPU'` | `'CPU'` or `'CUDA'` |
| `threads` | `1` | CPU threads |
| `topol` | `'center'` | Topology: `'random'`, `'center'`, `'grid'`, `'slab'` |
| `restart` | `'checkpoint'` | Restart mode: `'checkpoint'`, `'pdb'`, or `None` |

**Force Field Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps_lj` | `0.2` | LJ energy scale (kcal/mol, converted to kJ/mol) |
| `cutoff_lj` | `2.0` | LJ cutoff (nm) |
| `cutoff_yu` | `4.0` | Yukawa electrostatic cutoff (nm) |
| `fixed_lambda` | `0` | Fixed hydrophobicity for crowder interactions |

**Equilibration Options:**

| Parameter | Description |
|-----------|-------------|
| `slab_eq` | Equilibrate slab geometry with restraints |
| `bilayer_eq` | Equilibrate lipid bilayer under zero lateral tension |
| `box_eq` | Equilibrate with Monte Carlo barostat |
| `k_eq` | Equilibration restraint force constant |
| `steps_eq` | Equilibration steps |

---

### 2. Components (`calvados/cfg.py`)

Defines molecular components in the system.

**Default Component Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `molecule_type` | `'protein'` | Type: `protein`, `rna`, `lipid`, `crowder`, `cyclic`, `seastar`, `ptm_protein` |
| `nmol` | `1` | Number of molecules |
| `charge_termini` | `'both'` | Charge termini: `'N'`, `'C'`, `'both'`, or `None` |
| `restraint` | `False` | Apply structural restraints |
| `restraint_type` | `'harmonic'` | Restraint type: `'harmonic'` or `'go'` |
| `k_harmonic` | `700.0` | Harmonic restraint force constant |
| `k_go` | `15.0` | Go-model restraint force constant |
| `cutoff_restr` | `0.9` | Restraint distance cutoff (nm) |
| `use_com` | `True` | Use center of mass for restraints |
| `kb` | `8033.0` | Bond force constant (kJ/mol/nm²) |

**AlphaFold Integration Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `colabfold` | `0` | PAE format: 0=EBI, 1/2=ColabFold |
| `bfac_shift` | `0.8` | B-factor sigmoid shift |
| `bfac_width` | `50.0` | B-factor sigmoid width |
| `pae_shift` | `0.3` | PAE sigmoid shift |
| `pae_width` | `15.0` | PAE sigmoid width |

---

### 3. Sim (`calvados/sim.py`)

Main simulation engine. Handles system construction and MD execution.

**Key Methods:**

```python
sim = Sim(path, config, components)
sim.build_system()   # Initialize OpenMM system
sim.simulate()       # Run molecular dynamics
```

**Workflow:**
1. Parse configuration and component files
2. Initialize Debye-Hückel parameters
3. Create component objects (Protein, RNA, Lipid, etc.)
4. Build particle coordinates based on topology
5. Add forces (bonds, Ashbaugh-Hatch, Yukawa, restraints)
6. Run Langevin dynamics with energy minimization

---

### 4. Job (`calvados/cfg.py`)

HPC job submission management.

```python
job = Job(
    batch_sys='SLURM',      # or 'PBS'
    envname='calvados',     # conda environment
    template='robust.jinja' # job template
)
job.write(path, config, components)
job.submit(path, njobs=1)
```

---

## Molecule Types (`calvados/components.py`)

### Protein
- Standard amino acid chains
- IDRs (intrinsically disordered regions) or MDPs (multi-domain proteins)
- Supports harmonic or Go-model restraints for structured domains
- AlphaFold PAE integration for confidence-weighted restraints

### RNA
- Two-bead per nucleotide model (phosphate + base)
- Supports structured and unstructured regions
- Phosphate-phosphate and phosphate-base bonds
- Base-base stacking interactions

### Lipid / Cooke_Lipid
- Coarse-grained lipid models for membrane simulations
- Cosine attractive potential for tail-tail interactions
- Bilayer self-assembly capability

### Crowder
- Spherical molecular crowders
- Excluded volume effects
- PEG and other polymeric crowders

### Cyclic
- Cyclic peptides with N-C terminal bond

### Seastar
- Branched/star-shaped peptides
- Multi-arm architectures

### PTMProtein
- Proteins with post-translational modifications
- Attached PTM chains at specified locations

---

## Force Field (`calvados/interactions.py`)

### Ashbaugh-Hatch (AH) Potential
Hydrophobic/hydrophilic interactions between residues.

```
U_AH = ε × [LJ(r) - λ × LJ(rc)]  for r > 2^(1/6)σ
     = ε × [LJ(r) + (1-λ)]       for r ≤ 2^(1/6)σ
```

- `λ = 0`: Fully hydrophilic (repulsive only)
- `λ = 1`: Fully hydrophobic (attractive)
- Parameters from CALVADOS2/CALVADOS3 force fields

### Yukawa (Debye-Hückel) Potential
Screened electrostatic interactions.

```
U_YU = ε_yu × q1 × q2 × exp(-κr)/r
```

- `ε_yu`: Bjerrum length × kT
- `κ`: Inverse Debye length (depends on ionic strength)

### Bonded Interactions
- Harmonic bonds between consecutive residues
- RNA: additional angle terms for backbone stiffness

### Restraints
- **Harmonic**: Simple harmonic springs between residue pairs
- **Go-model**: 12-10 potential for native contacts

---

## Sequence Utilities (`calvados/sequence.py`)

### Input/Output Functions

```python
from calvados.sequence import seq_from_pdb, read_fasta, write_fasta, record_from_seq

# Extract sequence from PDB
seq, n_termini, c_termini = seq_from_pdb('protein.pdb')

# Read FASTA file
records = read_fasta('sequences.fasta')

# Create and write FASTA
record = record_from_seq(sequence, 'ProteinName')
write_fasta([record], 'output.fasta')
```

### Sequence Analysis

| Function | Description |
|----------|-------------|
| `get_qs()` | Calculate residue charges |
| `calc_SCD()` | Sequence charge decoration |
| `calc_SHD()` | Sequence hydropathy decoration |
| `calc_kappa()` | Charge patterning parameter (κ) |
| `calc_aromatics()` | Fraction of aromatic residues |
| `calc_mw()` | Molecular weight |

---

## Analysis Tools (`calvados/analysis.py`)

### Structural Properties

| Function | Description |
|----------|-------------|
| `calc_rg()` | Radius of gyration |
| `calc_ete()` | End-to-end distance |
| `calc_dmap()` | Distance map |
| `calc_cmap()` | Contact map |
| `calc_rmsd()` | RMSD to reference |
| `calc_fnc()` | Fraction of native contacts |
| `calc_ocf()` | Orientational correlation function |
| `fit_scaling_exp()` | Polymer scaling exponent (ν) |

### Trajectory Processing

| Function | Description |
|----------|-------------|
| `center_traj()` | Center trajectory |
| `subsample_traj()` | Subsample trajectory |
| `calc_com_traj()` | Calculate center-of-mass trajectory |
| `calc_contact_map()` | Average contact map over trajectory |

### Slab Analysis

`SlabAnalysis` class for phase separation simulations:
- Density profiles along z-axis
- Dilute/dense phase concentrations
- Partition coefficients

### Block Analysis

Error estimation using block averaging (from `BLOCKING` module):

```python
from calvados.BLOCKING.main import BlockAnalysis

block = BlockAnalysis(data_array)
block.SEM()  # Standard error of mean
```

---

## Residue Parameters

Located in `calvados/data/residues.csv` or custom files.

| Column | Description |
|--------|-------------|
| `one` | One-letter amino acid code |
| `three` | Three-letter code |
| `MW` | Molecular weight (Da) |
| `lambdas` | Hydrophobicity (0-1) |
| `sigmas` | Effective size (nm) |
| `q` | Formal charge |
| `bondlength` | Bond length to next residue (nm) |

CALVADOS versions:
- **CALVADOS2**: Original parameterization
- **CALVADOS3**: Updated hydrophobicity parameters

---

## Input File Requirements

### For IDRs (Intrinsically Disordered Regions)

| File | Required | Description |
|------|----------|-------------|
| `sequences.fasta` | Yes | Protein sequences |
| `residues.csv` | Yes | Amino acid parameters |

### For MDPs (Multi-Domain Proteins)

| File | Required | Description |
|------|----------|-------------|
| `protein.pdb` | Yes | 3D structure (AlphaFold or experimental) |
| `domains.yaml` | Yes | Structured domain residue ranges |
| `protein.json` | Optional | PAE matrix for Go restraints |
| `residues.csv` | Yes | Amino acid parameters |

### domains.yaml Format

```yaml
ProteinName:
  - [start1, end1]   # Domain 1 residue range
  - [start2, end2]   # Domain 2 residue range
```

---

## Examples Overview

| Example | Description |
|---------|-------------|
| `single_IDR` | Single intrinsically disordered protein |
| `single_MDP` | Single multi-domain protein with restraints |
| `single_RNA` | Single-stranded RNA |
| `single_dsRNA` | Double-stranded RNA |
| `single_pIDR` | Phosphorylated IDR |
| `two_IDR` | Two IDRs interacting |
| `two_IDR_MDP` | IDR + multi-domain protein |
| `slab_IDR` | IDR phase separation (slab geometry) |
| `slab_MDP` | MDP phase separation |
| `slab_IDR_MDP` | Mixed IDR/MDP phase separation |
| `slab_IDR_PEG` | IDR with PEG crowders |
| `slab_mixed` | Multiple component phase separation |
| `ten_IDR_cyl` | 10 IDRs in cylindrical geometry |
| `custom_restraints` | User-defined inter-molecular restraints |
| `single_AF_CALVADOS` | AlphaFold structure integration |
| `PARP14_MDP` | PARP14 multi-domain protein |
| `foxP_model` | FOXP transcription factors |
| `single_IDR_box_eq` | Box equilibration example |

---

## Topology Options

| Topology | Description |
|----------|-------------|
| `'center'` | Place single molecule at box center |
| `'random'` | Random placement avoiding overlaps |
| `'grid'` | Regular grid placement |
| `'slab'` | Slab geometry for phase separation |
| `'shift_ref_bead'` | Center on reference bead |

---

## Platform Options

| Platform | Description |
|----------|-------------|
| `'CPU'` | CPU execution (multi-threaded) |
| `'CUDA'` | NVIDIA GPU acceleration |

---

## Output Files

| File | Description |
|------|-------------|
| `top.pdb` | Coarse-grained topology |
| `{sysname}.dcd` | MD trajectory |
| `{sysname}.log` | Energy/speed log |
| `{sysname}.xml` | Serialized OpenMM system |
| `restart.chk` | Checkpoint for restart |
| `checkpoint.pdb` | Final structure |
| `bonds_{name}.txt` | Bond list |
| `restr_{name}.txt` | Restraint list |

---

## References

1. Tesei et al. PNAS (2021) - CALVADOS model
2. Tesei & Lindorff-Larsen. Open Res Europe (2022) - CALVADOS2
3. Cao et al. Protein Science (2024) - Multi-domain proteins
4. von Bülow et al. arXiv (2025) - Software package description
