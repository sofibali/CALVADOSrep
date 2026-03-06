# FOXP Transcription Factor Simulations

Coarse-grained CALVADOS simulations of FOXP1/FOXP4 transcription factor complexes. Folded domains are restrained to maintain structure while disordered linkers explore conformational space.

## Domain Architecture

### FOXP4 (551 residues, chain A)

| Region | Residues | Type |
|--------|----------|------|
| N-terminal IDR | 1-117 | Disordered |
| Domain 1 | 118-199 | Folded |
| Linker | 200-298 | Disordered |
| Domain 2 | 299-375 | Folded |
| Linker | 376-450 | Disordered |
| Domain 3 | 451-551 | Folded |

### FOX/FOXP1 (539 residues, chain B)

| Region | Residues | Type |
|--------|----------|------|
| N-terminal IDR | 1-117 | Disordered |
| Domain 1 | 118-299 | Folded |
| Linker | 300-318 | Disordered |
| Domain 2 | 319-375 | Folded |
| Linker | 376-449 | Disordered |
| Domain 3 | 450-539 | Folded |

### Other chains

| Chain | Residues | Description |
|-------|----------|-------------|
| chain_C | 90 res, domain [232,321] | Small protein |
| chain_D | 2 res | Spacer |
| chain_E | 21 nt (modeled as RNA) | DNA strand |
| chain_F | 21 nt (modeled as RNA) | DNA strand |

## Directory Structure

```
foxP_model/
├── input/                    # Shared PDBs, FASTA, domains.yaml
├── proteins_only/            # 3 proteins (FOXP4, FOX, chain_C) - start here
├── proteins_with_DNA/        # 3 proteins + 2 DNA chains (as RNA)
├── individual_chains/        # Multi-chain variants with different setups
├── single_complex/           # Single complex setup
├── FOX_FOXP4/                # Two-protein setup
├── FOXP4_FOX_chain_C/        # Original 3-protein setup
├── simulations/              # Simulation outputs
├── data/                     # Analysis data
├── analysis/                 # Analysis scripts/outputs
├── viz/                      # Visualization outputs
├── rendered_movies/          # PyMOL rendered movies (.pse, .mp4)
├── prepare.py                # Main preparation script
├── prepare_complex.py        # Complex preparation
├── split_pdb.py              # Split multi-chain PDB into individual chains
├── recunstruct_all_atom.py   # Reconstruct all-atom from CG trajectory
├── render_frames_pymol.py    # Render trajectory frames with PyMOL
├── render_sample_movie.py    # Sample movie rendering
├── save_pse_movies.py        # Save PyMOL session movies
├── export_movie.py           # Export movie files
├── prepare_movie_trajectories.py  # Prepare trajectories for movie rendering
├── pymol_trajectory_movie.pml     # PyMOL movie script
├── pymol_movie_coils_interchain.pml  # Inter-chain coils movie
├── pymol_movie_option1_pull.pml      # Pull visualization movie
└── CHANGELOG.md              # Change log
```

## Quick Start

### Proteins Only (recommended first run)

```bash
cd proteins_only
conda activate calvados
python run.py
```

### Proteins + DNA

```bash
cd proteins_with_DNA
conda activate calvados
python run.py
```

### From Scratch

```bash
python prepare.py          # generates config.yaml, components.yaml, run.py
cd <output_dir>
python run.py
```

### Analyzing a Running Simulation

```bash
python /home/sbali/CALVADOS/example_runscripts/analyze_trajectory.py \
    --top proteins_only/top.pdb \
    --traj proteins_only/proteins_only.dcd \
    --skip 0
```

## Simulation Setups

| Setup | Chains | Description |
|-------|--------|-------------|
| `proteins_only/` | FOXP4, FOX, chain_C | Three proteins, no DNA. Best for testing. |
| `proteins_with_DNA/` | FOXP4, FOX, chain_C, chain_E, chain_F | Full complex with DNA modeled as RNA. |
| `individual_chains/` | Various | Multi-chain variants exploring different combinations |
| `single_complex/` | Full complex | All chains in single complex |
| `FOX_FOXP4/` | FOX, FOXP4 | Two-protein interaction study |

## Key Parameters

- **Restraint strength:** k_harmonic = 2000 kJ/mol/nm^2 (strongly restrained domains)
- **Box size:** 30 nm cubic
- **Steps:** 5M for initial testing
- **Platform:** CPU (GPU optional for longer runs)

## Visualization Pipeline

The rendering pipeline converts CG trajectories into publication-quality movies:

1. **`recunstruct_all_atom.py`** - Reconstruct all-atom coordinates from CG beads
2. **`prepare_movie_trajectories.py`** - Prepare trajectory frames for rendering
3. **`render_frames_pymol.py`** - Render individual frames with PyMOL
4. **`save_pse_movies.py`** - Save PyMOL session files with trajectory
5. **`export_movie.py`** - Export final movie files

PyMOL scripts (`.pml`) define visualization styles for different movie types.

## Input Files

- **PDBs:** `input/FOXP4.pdb`, `input/FOX.pdb`, `input/chain_C.pdb`, etc.
- **FASTA:** `input/proteins.fasta`, `input/proteins2.fasta`
- **Domains:** `input/domains.yaml` - defines folded domain residue ranges per chain
- **Complex:** `input/COMPLEX_WITHDNA.pdb` - full complex structure
