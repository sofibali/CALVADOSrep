# CALVADOS Computational Projects

Personal fork of [KULL-Centre/CALVADOS](https://github.com/KULL-Centre/CALVADOS) used as a computational project record for coarse-grained biomolecular simulations.

## Projects

### 1. FOXP1/FOXP4 Transcription Factor Complex

Coarse-grained simulations of FOXP transcription factor complexes: multi-domain proteins with folded domains restrained from AlphaFold/experimental structures and disordered linkers free to explore conformational space.

- **Location:** [`examples/foxP_model/`](examples/foxP_model/)
- **Components:** FOXP4 (551 res), FOX/FOXP1 (539 res), chain_C (90 res), DNA chains
- **Simulation types:** proteins-only, proteins+DNA, individual chain variants
- **Visualization:** PyMOL movie rendering pipeline, trajectory analysis

### 2. PARP14 Domain Deletion Library

Two-phase pipeline studying how domain composition affects the structural dynamics of PARP14 (1801 residues, 17 domains).

| Phase | Description | Location |
|-------|-------------|----------|
| Phase 1 | AlphaFold3 structure predictions for 1,285 domain combinations | [`parp14/`](parp14/) |
| Phase 2 | CALVADOS coarse-grained MD (25 replicates x 50 ns per structure) | [`examples/PARP14_MDP/`](examples/PARP14_MDP/) |

## Installation

```bash
conda create -n calvados python=3.10
conda activate calvados

# Optional: GPU support
conda install -c conda-forge openmm=8.2.0 cudatoolkit=11.8

# Install CALVADOS
pip install -e .

# Verify
python -m pytest
```

## Running a Simulation

```bash
cd examples/foxP_model/proteins_only
python run.py
```

Or from scratch:
```bash
python prepare.py          # generates config.yaml, components.yaml, run.py
cd <output_dir>
python run.py
```

See [`docs/CALVADOS_QUICKSTART.md`](docs/CALVADOS_QUICKSTART.md) for detailed usage.

## Documentation

| Document | Description |
|----------|-------------|
| [`docs/CALVADOS_QUICKSTART.md`](docs/CALVADOS_QUICKSTART.md) | How to set up, run, and analyze simulations |
| [`docs/CALVADOS_ARCHITECTURE.md`](docs/CALVADOS_ARCHITECTURE.md) | Code architecture, classes, force fields |
| [`docs/PARP14_PROJECT.md`](docs/PARP14_PROJECT.md) | PARP14 simulation protocol and analysis |
| [`examples/foxP_model/README.md`](examples/foxP_model/README.md) | FOXP project overview and setups |
| [`examples/PARP14_MDP/README.md`](examples/PARP14_MDP/README.md) | PARP14 Phase 2 simulations |
| [`parp14/README.md`](parp14/README.md) | PARP14 Phase 1 AlphaFold3 pipeline |

## Upstream CALVADOS

This fork is based on [KULL-Centre/CALVADOS](https://github.com/KULL-Centre/CALVADOS). Please cite the following when using the CALVADOS software:

- Tesei et al. PNAS (2021), 118(44):e2111696118. [DOI: 10.1073/pnas.2111696118](https://doi.org/10.1073/pnas.2111696118)
- Tesei & Lindorff-Larsen. Open Research Europe (2022), 2(94). [DOI: 10.12688/openreseurope.14967.2](https://doi.org/10.12688/openreseurope.14967.2)
- Cao et al. Protein Science (2024), 33(11):e5172. [DOI: 10.1002/pro.5172](https://doi.org/10.1002/pro.5172)
- von Bulow et al. arXiv (2025). [DOI: 10.48550/arXiv.2504.10408](https://doi.org/10.48550/arXiv.2504.10408)
