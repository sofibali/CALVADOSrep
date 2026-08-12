# PARP14 CALVADOS analysis suite

All analysis scripts here share **`sim_registry.py`** and run on **any** simulation
set or folder via `--set` / `--sim-folder` — no code edits to add a new construct.
Same environment as the rest of the project (`environment.yml`, `calvados-tica`):
matplotlib / numpy / MDAnalysis / pandas / CALVADOS — all already included.

## The shared pattern (`sim_registry.py`)

Every script accepts:
- `--set NAME` — a named set (`fl`, `fl_optimized`, `md`, `md_full`, `core`, `mka`,
  `norrm`, `noart`, `md3art`) or a fragment name.
- `--sim-folder PATH [PATH ...]` — any folder of replicate trajectories. Normally
  `PATH/seed-{1-5}_sample-{0-4}/{top.pdb, <sysname>.dcd}`, but every script also
  falls back to a **flat scan** (any immediate subdirectory containing a `.dcd`,
  in sorted order) for folders that don't follow that grid — e.g. `fl_go`'s
  `state-N_tica_seed-*_sample-*_frF/` dirs (a handful of independent long runs,
  not a 5x5 seed/sample grid). Topology is read from `top.pdb`, falling back to
  `restart.pdb` then `checkpoint.pdb` (fresh/still-running sims never write
  `top.pdb`).
- `--units u1 u2 ...` — the FL domain units in that construct (e.g. `md1l1 md2 md3`),
  **only needed if** the folder has no `metadata.json`.

`sim_registry` resolves the construct's domain units (from `--units`, the folder's
`metadata.json`, or a known set name), then remaps FL-numbered domains / active sites
into the construct's local numbering automatically. Drop a `metadata.json` into new
folders with `write_metadata()` (or `stamp_metadata.py`) so `--units` isn't needed.

## One-shot wrapper

Run from inside `sim_analysis/` (or prefix with the folder from the project root,
e.g. `bash sim_analysis/run_analysis.sh ...`). Scripts always read/write the project
root's `data/` and `figures/` regardless of where they live or are invoked from.

```bash
cd sim_analysis
bash run_analysis.sh --set fl_optimized          # all modules on a named set
bash run_analysis.sh --sim-folder ../md_full --units md1l1 md2 md3
WORKERS=24 bash run_analysis.sh --sim-folder ../fl_go     # any new sim folder
WORKERS=24 bash run_analysis.sh --sim-folder ../fl_go --force   # force recompute
```
Runs `analyze_all.py` (8 modules) + `analyze_lys_contacts.py`, writing arrays to
`data/`. Figures are separate (below) since they read `data/`.

## Scripts

**Trajectory analysis (sim_registry, `--set`/`--sim-folder`):**
| Script | What it computes | Example (fl_go) |
|--------|------------------|------------------|
| `analyze_all.py` | master; 8 modules — `--conf-prop --dmap --cmap --fnc --energy --wcn --active-sites --accessibility` (run a subset by flag, all by default). `--workers N`, `--include-fragments`, `--force` to recompute. | `python analyze_all.py --sim-folder ../fl_go --workers 24` |
| `analyze_lys_contacts.py` | lysine <-> acidic / lysine <-> lysine contact maps + domain-sum heatmaps (plain + `_v2` inter-domain-only color scale). | `python analyze_lys_contacts.py --sim-folder ../fl_go` |
| `analyze_active_sites.py` | standalone, more detailed than `analyze_all.py --active-sites`: contact number, inter-site distances, radial position, RMSF, domain contacts, one combined ensemble figure. `--seed/--sample/--ensemble/--skip-frames`. | `python analyze_active_sites.py --sim-folder ../fl_go` |
| `analyze_accessibility.py` | standalone, same underlying method as `analyze_all.py --accessibility` (writes the same `data/accessibility_stats.npz`) but with its own SAA/cone/summary figures and `--nrays`. | `python analyze_accessibility.py --sim-folder ../fl_go --nrays 200` |

**Figures that talk to trajectories directly (sim_registry, `--set`/`--sim-folder`):**
| Script | What it computes | Example (fl_go) |
|--------|------------------|------------------|
| `figure_md_distances.py` | inter-domain COM-COM (or `--metric min` CA-CA) distance violins, all replicates pooled. `--pairs` to override the default MD1L1-MD2/MD3/ART, MD2-MD3 set. | `python figure_md_distances.py --set fl_go` (or `--sim-folder ../fl_go`) |
| `figure_lys_exposed_persistence.py` | per-residue exposed Lys<->acidic contact persistence bars (needs `data/sasa_face_per_residue.csv` from `figure_sasa_faces.py` first). | `python figure_lys_exposed_persistence.py --sets fl_go` (or `--sim-folder ../fl_go`) |

**Figures that read cached `data/` only (run `analyze_all.py` first; not directly
`--sim-folder`-aware — they operate on the FL reference structure or compare fixed
named sets):** `figure_accessibility.py`, `figure_sasa_faces.py`,
`figure_domain_comparisons.py`, `figure_kh_domains.py`, `figure_schematic.py`,
`make_analysis_pse.py`. They share `_fig_layout.py` (style).

**Conformational states / TICA clustering:** `../tica_pipeline/cluster_states.py`
(its own README lives in that folder). For an arbitrary folder like `fl_go`, call
it directly with `--sim-folder` — the staged `01_featurize_tica.py` /
`02_select_states.py` / `03_backmap_minimize.py` scripts only accept `--set` (named
sets), not `--sim-folder`.
```bash
cd ../tica_pipeline
python cluster_states.py --sim-folder ../fl_go --features com \
    --reduce tica --tica-lag 200 --k 5 --rep-mode spread --workers 8
```
`--features` options: `com` (COM-COM domain-distance pairs, default — treats each
restrained domain as a rigid body via its centroid), `pose` (adds a full relative
rotation matrix per domain pair — high-dimensional, 660 features for 11 domains),
`site_orient` (COM-COM distances + a rotationally-aware but low-dimensional
complement: for each active-site domain, the cosine angle between its
COM->pocket-centroid "site vector" and the direction to every other domain, plus
site-vector-vs-site-vector angles between active sites — invariant to the whole
molecule's tumbling, sensitive to real reorientation). Restrict to just the
active-site domains with `--domains MD1L1 MD2 MD3 ART`. A `--domains` subset gets
a short hash suffix on its cache file / output dir so it never collides with a
different `--domains` run at the same `--features` on the same set.

**Legacy / special-purpose (not CG-MD trajectory analysis):**
`analyze_kh_domains.py` (AF3 PAE/pLDDT domain-boundary analysis on structure
predictions, not MD trajectories).

## Outputs
- `data/<set>_<analysis>.npy|.npz` — per-set arrays (modules skip if present; `--force` to redo).
- `data/*_per_residue.csv` — per-residue tables (RSA, SASA, contacts).
- `figures/` — PNG + SVG, dated subdirectories per analysis category.
- `../figures/05_clustering/<date>/<set>_<features>[-d<hash>]_<reduce>/` — TICA landscape,
  ITS, silhouette, PCCA plots.
- `../representative_frames/<date>/<set>_<features>_<reduce>/state_*.pdb` — TICA
  representative-frame structures.

## Adding a new construct
1. Run the sim (`prepare_*` -> `run.py`), producing trajectories anywhere under
   `PATH/` -- either the `seed-*_sample-*/<sysname>.dcd` grid or any flat layout
   of replicate subdirectories, each with a `.dcd` and a `top.pdb`/`restart.pdb`/
   `checkpoint.pdb`.
2. (Recommended) drop a `metadata.json` with its `units`.
3. `bash run_analysis.sh --sim-folder PATH` — no code edits.
