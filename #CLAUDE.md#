# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is CALVADOS?

Coarse-grained implicit-solvent simulation framework for biomolecules (proteins, RNA, lipids) built on OpenMM. One bead per amino acid (CA-level), with Ashbaugh-Hatch hydrophobic interactions and Debye-Huckel electrostatics.

## Build & Test

```bash
conda create -n calvados python=3.10
conda activate calvados
pip install -e .          # editable install
python -m pytest          # run all tests (test_potentials, RNA bond order, custom restraints)
```

GPU support requires: `conda install -c conda-forge openmm=8.2.0 cudatoolkit=11.8` before pip install.

## Running a Simulation

Simulations are configured via two YAML files and a prepare script:

1. **`prepare.py`**: Uses `calvados.cfg.Config` and `calvados.cfg.Components` to generate `config.yaml`, `components.yaml`, and `run.py`
2. **`run.py`**: Calls `sim.run(path, fconfig, fcomponents)` which instantiates `Sim`, calls `build_system()`, then `simulate()`

```bash
python prepare.py          # generates config.yaml, components.yaml, run.py in output dir
cd <output_dir> && python run.py
```

## Architecture

### Core module (`calvados/`)

- **`sim.py`** — `Sim` class: orchestrates system setup. Key methods:
  - `build_system()`: creates OpenMM System, places molecules, adds all forces
  - `simulate()`: creates Langevin integrator, handles restart (checkpoint/pdb), runs MD
  - `make_components()`: instantiates component objects from YAML, computes properties
  - `add_forces_to_system()`: adds AH, YU, bonded, restraint, and external forces

- **`components.py`** — Component hierarchy: `Component` (base) → `Protein`, `RNA`, `Lipid`, `Crowder`, `Cyclic`, `Seastar`, `PTMProtein`. Each component manages its own sequence, bond force, and restraint force.

- **`interactions.py`** — Force initialization:
  - `init_ah_interactions()`: Ashbaugh-Hatch (CustomNonbondedForce)
  - `init_yu_interactions()`: Yukawa/Debye-Huckel electrostatics
  - `init_bonded_interactions()`: harmonic bonds
  - `init_restraints()`: harmonic or Go-model restraints between domain residue pairs

- **`build.py`** — Coordinate generation (linear, spiral, compact, random placement), PDB parsing, box construction, clash checking

- **`cfg.py`** — `Config` (simulation parameters), `Components` (system composition), `Job` (SLURM/PBS submission). These read defaults from `calvados/data/default_config.yaml` and `default_component.yaml`.

### Key configuration parameters

**config.yaml** (simulation):
- `topol`: molecule placement strategy (`center`, `random`, `slab`, `grid`)
- `restart`/`frestart`: restart mode (`checkpoint`+`.chk`, `pdb`+`.pdb`, or `null`)
- `box`: simulation box dimensions in nm
- `steps`, `wfreq`: total steps and write frequency

**components.yaml** (system composition):
- `system`: dict of named components, each with `molecule_type`, `nmol`, `ffasta`/`restraint`/`pdb_folder`, `fdomains`, etc.
- `defaults`: fallback values from `default_component.yaml`
- Restraints require: `restraint: true`, a PDB in `pdb_folder/<name>.pdb`, and `fdomains` pointing to a domains YAML

### Multi-chain with domain restraints

```yaml
# components.yaml system section
system:
  FOXP4:
    molecule_type: protein
    restraint: true
    pdb_folder: input
    fdomains: input/domains.yaml
  FOXP1:
    molecule_type: protein
    restraint: true
    pdb_folder: input
    fdomains: input/domains.yaml
```

Chain order in the system follows the order in the YAML `system` section. For restart from PDB (`topol: random`, `restart: pdb`, `frestart: restart.pdb`), the PDB must have coordinates in Angstroms matching the chain/residue count.

### Adding custom forces

After `build_system()` and before `simulate()`:
```python
mysim = sim.Sim(path, config, components)
mysim.build_system()
custom_force = openmm.CustomExternalForce(...)
mysim.system.addForce(custom_force)
mysim.simulate()
```

## Examples

Examples are in `examples/`. Each has a `prepare.py` that generates simulation input. The `foxP_model/` example demonstrates multi-domain protein modeling with domain restraints, inter-chain interactions, and visualization pipelines.

## PARP14 Project (`parp14/`)

Two-phase pipeline for studying PARP14 (1801 residues) domain deletion variants:

**Phase 1 — Structure Generation** (`parp14/`):
- AlphaFold3 predictions for all 2^11 - 1 = 2,047 domain combinations (5 seeds × 5 samples = 25 models each)
- 774 completed, 1,288 remaining (inputs in `alphafold_inputs_missing/`, run via `run_af3_missing.sh`)
- Pipeline scripts in `parp14/pipeline/`: FASTA generation → AF3 inputs → predictions → analysis
- Outputs: structures in `alphafold_outputs/` (gitignored), CALVADOS-ready inputs in `calvados_simulations/`
- Quality metrics: `pLDDT_per_residue.csv`, `sasa_per_residue.csv`, `confidence_scores.csv`
- AF3 run: 2-phase (MSA on CPU parallel, then inference with `--input_dir` for single model load per GPU)

**Phase 2 — MD Simulations** (`examples/PARP14_MDP/`):
- CALVADOS CG-MD on AF3 structures, 6 construct sets x 25 replicates x 20 ns = 500 ns per set
- Protocol: 20 ns per replicate, discard first 0.5 ns equilibration, 25 x 19.5 ns = 487.5 ns effective per set
- **Simulation sets** (prepared by `prepare_and_run_all.py`):
  - `fl` — Full-length (1801 res, 300 nm box, EBI AF2, `colabfold=1`)
  - `md` — Macrodomains only: MD1L1+MD2+MD3 (586 res, 80 nm box)
  - `core` — KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE+ART (1051 res, 120 nm box)
  - `mka` — MD1L1+MD2+MD3+KHb-KH8+WWE+ART (999 res, 100 nm box)
  - `norrm` — No-RRM: KH1-6+KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE+ART (1474 res, 250 nm box)
  - `noart` — No-ART: KH1-6+KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE (1275 res, 200 nm box)
- Domain restraints use FL structured-core boundaries (from `input/domains.yaml`), mapped to construct numbering; inter-domain linkers left flexible (unrestrained)
- Directory layout (nested): `{set}/seed-{1-5}_sample-{0-4}/` (e.g. `md/seed-1_sample-0/`)
- Shared inputs per construct: `{set}/input/domains.yaml` + `residues_CALVADOS3.csv`
- Launch: `bash run_all.sh parallel` (or `run_all.sh fl`, `run_all.sh norrm`, etc.)
- New constructs: `bash run_new.sh parallel` (norrm + noart only)
- **Analysis** (`analyze_all.py`): parallelized with ProcessPoolExecutor, 8 modules (conf-prop, dmap, cmap, fnc, energy, wcn, active-sites, accessibility). Figures in `figures/` (PNG + SVG). Energy analysis uses trimmed domain boundaries (3 residues from each zero/small-gap border) to avoid steric-clash artifacts; all plots on canonical 11-domain grid with gray for missing domains, ±1 kJ/mol fixed scale, plus difference maps vs FL and all-vs-all.
- **Accessibility figures** (`figure_accessibility.py`): 7 publication figures (heatmap, per-site SAA, metrics by construct, domain context, docking feasibility, schematic, ranking). Violin + jitter points + SD style.
- **Analysis glossary** (`analysis_glossary.txt`): detailed description of every metric, what it measures physically, units, typical value ranges, and which figures it produces.
- See `docs/PARP14_PROJECT.md` for full protocol

**Phase 2.5 — Restraint Optimization** (`examples/PARP14_MDP/test_restraints.py`):
- Systematic boundary trim sweep (0 / 5 / 10 / 15 / 20 / 25 / 30 residues from each side of full domain extents), each tested with/without KH7a-KHb custom inter-domain restraints
- 16 test configurations (2 baselines + 7 trims × 2 with/without custom)
- Scoring uses experimental Rg references from PDBs 3VFQ (MD1, MD2), 3Q6Z (MD1), 3GOY (WWE, ART), 1X4R (WWE NMR) plus KH homolog scaling from FUBP1 6Y2D / Nova1 2ANR (Saleh 2024 JMB). Crystal PDBs in `input/xtal_refs/`.
- Custom KH7a-KHb restraints in `input/custom_restraints.txt` (141 pairs at k=350 kJ/mol/nm²) via CALVADOS `custom_restraints: true` feature
- Boundary schematic with serine positions (`PARP14_domains_atS.fasta` boundary serines marked) + simulation summary CSV/MD in `restraint_tests/summary_plots/`
- Comparison + PyMOL session: `compare_restraints.py` produces RMSD-vs-crystal-or-AF2 boxplots, Rg/CMap heatmap, KH-domain RMSD distributions, best-config-per-domain, and PyMOL scenes (`load_kh_comparison.pml`)

**Phase 2.6 — Optimized FL + Fragment Library**:
- `prepare_fl_optimized.py` — 25-replicate FL with per-domain custom trim values (RRM1 trim 10 from structured core 6-88; RRM/RRM trim 10; KH1/2/3/5 trim 5; KH4/6 trim 10; KH7a, KHb, KH8 trim 0; MDs trim 10; WWE trim 15; ART trim 10) + KH7a-KHb custom restraints
- `prepare_all_fragments.py` — All 66 contiguous PARP14 fragments (48 prepared, 17 missing AF3 structures); per-fragment in `fragments/{name}/seed-{1-5}_sample-{0-4}/`. Box size scales with construct length. Launch via `fragments/run_all_fragments.sh N_PARALLEL` (240-core machine → ~60 parallel sims).

**Phase 3 — Per-residue analysis**:
- `figure_sasa_faces.py` — RSA (relative SASA, Chothia Gly-X-Gly normalized per Wu 2017) heatmap + per-domain active-site-face vs back-face annotation. Outputs: `figures/rsa_heatmap_FL.png`, `figures/rsa_per_domain_matrix_FL.png`, PyMOL session, per-residue CSV.
- `figure_md_distances.py` — Inter-domain COM-COM (or min CA-CA) distance violin plots across all 25 replicates pooled. Default pairs: MD1L1-MD2/MD3/ART, MD2-MD3, MD3-ART. 1.0 nm contact cutoff line.
- `cluster_states.py` — K-means clustering of frames using inter-domain distances (optional Rg), with silhouette/elbow sweep, PCA scatter, per-state distance profile heatmap, replicate-by-state distribution heatmap, and extraction of centroid-representative PDB structures per state. Requires `scikit-learn`.

**PARP14 Domain Architecture (11 grouped units for combinatorial library):**

| Domain | Residues | Function |
|--------|----------|----------|
| RRM1 | 1-145 | RNA Recognition Motif 1 |
| RRM2 | 146-224 | RNA Recognition Motif 2 |
| RRM3 | 225-314 | RNA Recognition Motif 3 |
| KH1-KH6 | 315-737 | K-Homology domains (grouped) |
| KH7a | 738-789 | K-Homology domain |
| MD1L1 | 790-1004 | Macrodomain 1 + linker to MD2 |
| MD2 | 1004-1193 | Macrodomain 2 |
| MD3 | 1207-1388 | Macrodomain 3 |
| KHb-KH8 | 1389-1533 | K-Homology domains (grouped) |
| WWE | 1534-1602 | WWE domain |
| ART | 1603-1801 | ADP-ribosyltransferase (catalytic) |

**Active Sites** (defined in `parp14/input/active_sites.yaml`):
- MD1: D831, N923, D962 (ADP-ribose hydrolase)
- MD2: 1035, 1046, 1134, 1171 (ADP-ribose reader)
- MD3: 1248, 1259, 1330, 1371 (ADP-ribose binding)
- ART: H1684, Y1705, E1706, 1722 (H-Y-E catalytic triad)

**Visualization:** `parp14/make_pse.py` generates PyMOL PSE with domain coloring + active site labels (uses `pymol-render` conda env).

**Key tools:** `parp14_tools/` contains helper scripts and examples for the pipeline.
