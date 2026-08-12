# TICA → state-selection → full-atom pipeline

Three small, single-purpose scripts (plus a wrapper) that turn CALVADOS CG
trajectories into a handful of **full-atom, clash-free representative structures**.
Each stage is runnable on its own with clear flags; `run_pipeline.sh` chains them.

```
01_featurize_tica.py   features -> TICA projection        (-> data/tica_<tag>.npz)
02_select_states.py    TICA -> k representative CG states  (-> states/<tag>/state_*.pdb)
03_backmap_minimize.py CG states -> full-atom, clash-free  (-> states/<tag>_allatom/)
run_pipeline.sh        runs 1 -> 2 -> 3 from one config
```

## Setup

```bash
conda env create -f ../environment.yml      # creates calvados-tica
conda activate calvados-tica
pip install -e /path/to/CALVADOS            # only needed to run new simulations
```
The two dependencies added beyond CALVADOS's `setup.py` are **scikit-learn**
(k-means microstates) and **pdbfixer** (terminal/missing-atom fixing for backmap).

## Quick start

```bash
# md_full, rigid-body pose, 4 evenly-spread states, full-atom:
bash tica_pipeline/run_pipeline.sh

# full-length set, CA25 feature, 4 states:
SET=fl_optimized FEATURES=ca K=4 \
  REFERENCE=input/parp14.pdb FDOMAINS=input/domains.yaml \
  bash tica_pipeline/run_pipeline.sh
```

Run stages individually (all paths relative to `examples/PARP14_MDP/`):

```bash
python tica_pipeline/01_featurize_tica.py --set md_full --features pose \
    --tica-lag 200 --out data/tica_md_full_pose.npz
python tica_pipeline/02_select_states.py --tica-npz data/tica_md_full_pose.npz \
    --k 4 --rep-mode spread --out-dir states/md_full_pose
python tica_pipeline/03_backmap_minimize.py --states-dir states/md_full_pose \
    --reference md_full/input/ref_allatom.pdb --fdomains md_full/input/domains.yaml
```

## Key choices

**Featurization** (`--features`, stage 1) — `pose` (rigid-body relative position +
orientation per domain pair) best matches the arrangements you see by eye, because
plain distances (`ca`, `com`) are rotation-invariant. `pose_iface` adds the interface
contact register. The feature *definition* is identical across construct sizes, so the
same components are computed for md_full, fl_optimized, fragments, etc.

**State selection** (`--rep-mode`, stage 2):
- `spread` — farthest-point sampling + Voronoi tiling → evenly-distributed,
  non-overlapping states. The honest choice when implied timescales don't plateau
  (diffusive system): "distinct sampled conformations", not Markovian basins.
- `pcca` — MSM + PCCA+ metastable macrostates. Only when the timescales converge.

**Lag** (`--tica-lag`, frames; 1 frame = 0.01 ns) — use the shortest lag where the
ITS plateau. For a diffusive system they don't, so a short lag (200 = 2 ns) keeps
statistics/resolution; 8 ns just bleeds counts. (See `../WHAT_DOES_THIS_ALL_MEAN.md`.)

**Reference** (stage 3) — an all-atom structure with side chains, in the construct's
numbering: `md_full/input/ref_allatom.pdb` (FL 790–1388 excised) or `input/parp14.pdb`
(full length). Back-mapping grafts these domains onto each state's CG frame; clashes
are then removed by restrained minimization, with a rigid-body declash fallback for
deep inter-domain overlap.

## Reuse / provenance
The engine modules live in this folder (`cluster_states.py`,
`fullatom_minimize_states.py`, `backmap_states.py`), so the pipeline is
self-contained; the shared `sim_registry.py` / `_fig_layout.py` stay in the project
root and are imported via a root-bootstrap at the top of each script.
- Stages 1–2 import the heavy functions from `cluster_states.py`
  (`collect_features`, `compute_tica`, `cluster_kmeans`, `pcca_plus`,
  `find_spread_frames`, `extract_pdb_for_frame`).
- Stage 3 wraps `fullatom_minimize_states.py` (which uses `backmap_states.py`).
- `cluster_states.py` remains usable standalone (it does all three at once); this
  pipeline is the decomposed, shareable version.
- Paths (`data/`, `figures/`, `representative_frames/`, sim folders) always resolve
  to the project root regardless of where you invoke from.
