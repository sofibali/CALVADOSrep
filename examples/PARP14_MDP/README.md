# PARP14 Domain Deletion - CALVADOS Simulations (Phase 2)

Coarse-grained MD simulations of PARP14 domain deletion variants using CALVADOS. This is **Phase 2** of the pipeline; Phase 1 (AlphaFold3 structure generation) is in [`../../parp14/`](../../parp14/).

## Overview

- **Construct:** KH1-KH6, KH7a, MD1L1, MD2, MD3, KHb-KH8, WWE, ART (1474 residues)
- **Initial structures:** 25 AlphaFold3 models (5 seeds x 5 samples)
- **Protocol:** 10 ns per simulation, 1000 frames, harmonic domain restraints

## Quick Start

### Prepare all 25 simulations

```bash
conda activate calvados
python prepare_parp14_25x10ns.py
```

This converts AF3 mmCIF structures to PDB, sets up domains.yaml and CALVADOS configs for each seed/sample combination.

### Run a single simulation

```bash
cd parp14_seed-1_sample-0
python run.py
```

### Run all 25 in parallel

```bash
for d in parp14_seed-*_sample-*; do
    (cd "$d" && python run.py) &
done
wait
```

### Analyze results

```bash
python analyze_parp14.py
```

## Directory Structure

```
PARP14_MDP/
├── prepare_parp14_25x10ns.py   # Prepare 25 simulations from AF3 structures
├── analyze_parp14.py           # Per-domain analysis (Rg, RMSD, FNC, contacts)
├── input/                      # Shared input (domains.yaml, residues CSV)
├── input_construct/            # Construct-specific domains.yaml
├── parp14_seed-{1-5}_sample-{0-4}/  # 25 individual simulation directories
│   ├── config.yaml
│   ├── components.yaml
│   ├── run.py
│   └── input/
│       ├── parp14_construct.pdb    # Backbone PDB from AF3
│       └── parp14_construct.json   # PAE confidence data
├── data/                       # Analysis outputs
│   ├── domain_rgs.npz
│   ├── domain_fnc.npz
│   ├── interdomain_com_dist.npy
│   ├── interdomain_contact_freq.npy
│   └── *.png                  # Time series and heatmap plots
└── parp14/                    # Compiled trajectory directory
```

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Steps | 1,000,000 (10 ns) |
| Write frequency | 1,000 (1000 frames) |
| Box | 150 nm cubic |
| Temperature | 293 K |
| Ionic strength | 0.19 M |
| Restraints | harmonic, k=700 kJ/mol/nm^2 |
| Platform | CPU (set PLATFORM='CUDA' for GPU) |

## Domain Definitions (Construct Numbering)

| Domain | Construct residues | Full-length residues |
|--------|-------------------|---------------------|
| KH1-KH6 | 1-423 | 315-737 |
| KH7a | 424-475 | 738-789 |
| MD1L1 | 476-690 | 790-1004 |
| MD2 | 690-879 | 1004-1193 |
| MD3 | 880-1061 | 1207-1388 |
| KHb-KH8 | 1062-1206 | 1389-1533 |
| WWE | 1207-1275 | 1534-1602 |
| ART | 1276-1474 | 1603-1801 |

## Analysis Outputs

`analyze_parp14.py` computes:

- **Per-domain Rg** - radius of gyration time series
- **Per-domain RMSD** - vs reference structure and vs mean
- **Per-domain RMSF** - per-residue fluctuations
- **Per-domain FNC** - fraction of native contacts
- **Interdomain COM distances** - pairwise distance matrix heatmap
- **Interdomain contact frequency** - contact map heatmap
- **Global Rg and Ree** - whole-protein properties

All plots saved to `data/` as PNG files, numerical data as NPZ/NPY.

## Phase 1 Connection

Initial structures come from AlphaFold3 predictions in:
```
/home/sbali/parp14/alphafold_outputs/kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art/
```

See [`../../parp14/README.md`](../../parp14/README.md) for the structure generation pipeline, and [`../../docs/PARP14_PROJECT.md`](../../docs/PARP14_PROJECT.md) for the full protocol.
