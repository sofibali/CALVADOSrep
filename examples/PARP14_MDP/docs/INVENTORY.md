# PARP14 project inventory — where everything is and what it contains

Last surveyed 2026-09-21. Everything below is under
`/home/sbali/CALVADOS/examples/PARP14_MDP/` unless noted, with Phase-1
structure generation in the sibling `/home/sbali/CALVADOS/parp14/`.

**~60 GB total.** Most of it — every trajectory, every cache in `data/`, every
figure — is **gitignored and exists only on this filesystem**. Only inputs,
scripts and docs are in version control. Treat this machine as the sole copy.

---

## 1. Start here

| path | what it is |
|---|---|
| `docs/QUESTIONS.md` | the science entry point: each question the simulations ask, whether it is answered, and the answer with caveats |
| `docs/questions/Q1–Q10_*.md` | one file per question. Q3 accessibility, Q4 writer-vs-eraser, Q9 self-association/puncta are the ones with live analyses |
| `docs/NUMBERING_AUDIT.md` | the residue-numbering audit (found 2026-09-18, fixed 2026-09-21) and what it did/didn't invalidate |
| `analysis_glossary.txt` | every metric: what it measures physically, units, typical ranges, which figures it produces |
| `slab/README.md` | multi-chain campaign: measured cost, GPU-sharing policy, gotchas |
| `slab/SESSION_2026-09-18_lyra.md` | what was measured on lyra, tree state, what to check when moving machines |

---

## 2. Single-chain trajectories (~30 GB)

**290 replicate directories** across 22 sets. Layout is always
`{set}/seed-{1-5}_sample-{0-4}/`, one CALVADOS run per directory.

| tier | sets | reps each |
|---|---|---|
| primary | `fl`, `fl_optimized`, `norrm`, `noart`, `core`, `mka`, `md`, `md3art`, `md_full` | 25 |
| extended | `mka_full`, `core_full_go`, `kh1_art_full`, `kh1_wwe_full`, `md2_art_full`, `md2_wwe_full`, `md3_art_full`, `md3_wwe_full`, `core_wwe_full_go`, `mka_wwe_full`, `fl_wwe_full_go`, `md1_md2`, `md2_md3` | 5 |
| fragments | `fragments/{name}/` — 48 contiguous PARP14 fragments (7.2 GB) | 25 |

Largest on disk: `mka_full` 5.7 G, `fl_optimized` 5.1 G, `md_full` 4.8 G,
`fl` 1.9 G, `norrm` 1.9 G, `noart` 1.6 G.

### What a replicate directory holds

```
parp14.dcd            the trajectory -- the actual result. Frame every wfreq steps
parp14.log            StateDataReporter: Step / Speed(ns/day) / Elapsed(s).
                      Production steps only -- equilibration is not in it.
restart.chk           OpenMM checkpoint; resuming runs `steps` ADDITIONAL steps
checkpoint.pdb        last frame as PDB
parp14_<timestamp>.pdb final frame, one per completed leg
parp14.xml            serialised OpenMM System (~40 MB, regenerated every run)
bonds_parp14.txt      bonded terms written at setup
restr_parp14.txt      the harmonic/Go restraint list actually applied
config.yaml           simulation parameters as run
components.yaml       system composition (absolute paths)
input/                starting PDB, domains.yaml, residues_CALVADOS3.csv
run.py                entry point
```

Sets, their domain units and their FL residue spans are declared in
`sim_registry.py` — that is the authority, not directory names.

---

## 3. Analysis caches — `data/` (5.9 GB, gitignored)

Key grammar is `{set}_{site}_{metric}`, and each array is **one value per
replicate**, so `shape == (n_reps,)` — 25 for primary sets, 5 for extended.

| file | keys | contents |
|---|---|---|
| `active_site_stats.npz` | **275** | per-site `radial` (distance from protein COM / Rg) and `rmsf`; plus inter-site COM distances keyed by pair (`MD1-MD2`, `MD1-ART`, `MD2-MD3`, `WWE-ART`, …). **Refreshed 2026-09-21 with the corrected residue numbers.** |
| `accessibility_stats.npz` | 222 | solid-angle accessibility per site: `saa` (fraction of 200 rays escaping to 5 nm), `cone`, `density`. Feeds Q3 and Q4. |
| `cluster_assignments_*_ca25_tica.npz` | ~30 files | per-frame state labels from TICA/k-means. The bulk of the 5.9 GB (`fl` alone is 895 MB). |
| `cluster_features_*_ca25.npz` | ~30 files | the feature matrices those assignments were computed from |
| `puncta_features.csv` | 2047 × 39 | every domain combination: `combo`, `combo_key`, `barcode`, `n_residues`, `rna_valence`, `adpr_reader_valence`, `has_writer_ART`, `has_eraser_MD1`, plus sequence features (FCR, NCPR, κ, SCD, aromatics). Keyed for joining to the sort-seq screen. |
| `art_steric_effect_saa.csv` | — | ±ART matched-pair steric-coupling result |

A note on the older `cluster_*` filenames: they still spell the MD1 unit
`md1l1`. That name was retired on 2026-09-21 in favour of `md1`, but
`sim_registry.UNIT_ALIASES` resolves the old spelling, so these files stay
readable and were deliberately not renamed.

---

## 4. Figures — `figures/` (1.4 GB, ~2,900 files)

Every figure is written as both `.png` and `.svg` unless noted.

| directory | files | contents |
|---|---|---|
| `by_sim/{set}/` | 946 | per-set output, 23 sets. Each has `02_main_analysis/`, `03_accessibility/`, `04_md_distances/`, `05_clustering/`, `07_lysine_contacts/`, plus **`dashboard.html`** and `dashboard_montage.png` |
| `02_main_analysis/` | 468 | conformational properties, distance/contact maps, energy maps, Rg convergence, active sites |
| `05_clustering/` | 743 | TICA/k-means states, silhouette/elbow sweeps, PCA scatters, per-state profiles |
| `09_surface_gallery/` | 237 | surface renders per state/face |
| `comparisons/` | 124 | cross-set comparisons; subdirectory name is the concatenated set list |
| `03_accessibility/` | 67 | the 7 publication accessibility figures |
| `07_lysine_contacts/` | 42 | Q6 |
| `01_static_FL/` | 18 | static full-length reference panels |
| `04_md_distances/` | 13 | inter-domain COM–COM violins |
| `_archive/`, `figures_5ns/` | 249 | superseded 5 ns-era output, kept for provenance |

**23 `dashboard.html` files** — one per set, the quickest way to eyeball a set.

---

## 5. Multi-chain slab campaign — `slab/` (294 MB)

The only part of the project that can answer *does this construct
self-associate*. Everything single-chain above cannot, by construction.

```
slab/
  homotypic/{construct}/     N copies of one construct
  rna/{construct}/           the same + polyU40 RNA
  benchmark/md2_md3/         throughput calibration (its job is done)
  logs/                      launcher logs, one per arm/construct
  README.md                  cost, GPU policy, gotchas
  SESSION_2026-09-18_lyra.md what was measured and what state the tree is in
  run_slab_queue.sh          launcher for a dedicated GPU
  run_slab_opportunistic.sh  launcher for a shared GPU (yields to other users)
  monitor_slab.py            progress / ns-per-day / ETA / GPU load
  slab_steps.py              remaining-step arithmetic
  submit_slab.slurm          for a SLURM cluster; NOT usable on lyra
```

10 constructs × 2 arms = 20 runs, all validated (build + run on CUDA with
correct bead counts). Each construct directory holds `config.yaml`,
`components.yaml`, `run.py`, `slab_meta.yaml` and `input/` (starting PDB,
`domains.yaml`, `residues_CALVADOS3.csv`, plus `rna.fasta` in the rna arm).

### Current output

Only `core_full_go` has data; everything else is at 0 steps.

| path | state |
|---|---|
| `slab/homotypic/core_full_go/parp14_core_full_go.dcd` | 186 MB, ~26M/200M steps |
| `slab/homotypic/core_full_go/restart.chk` | resumable, equilibration done |
| `slab/rna/core_full_go/` | equilibration only; that arm is waiting on a contended GPU |

**The result this campaign exists to produce is c_sat**, which lands in
`{construct}/{sysname}_ps_results.csv` as the `c_dilute` column (mM), via
`calvados.analysis.SlabAnalysis`. No run has reached that point yet.

`components.yaml` carries **absolute paths**, so the tree is not portable —
re-run `prepare_slab.py` at a new location rather than copying.

---

## 6. Tools

| path | role |
|---|---|
| `sim_registry.py` | **single source of truth**: sets, domain units, FL↔construct residue mapping, active sites. Everything imports from here. |
| `parp14/input/active_sites.yaml` | the corrected catalytic/pocket residue definitions `sim_registry` loads |
| `sim_analysis/` | 32 scripts. Key: `analyze_all.py` (8 parallel modules), `analyze_accessibility.py`, `analyze_active_sites.py`, `predict_enzyme_dominance.py`, `predict_puncta_propensity.py`, `figure_sasa_faces.py`, `figure_md_distances.py` |
| `sim_analysis/verify_numbering.py` | checks every residue definition against the **live UniProt Q460N5 record**; exits non-zero on failure, so it works in CI |
| `tica_pipeline/` | featurisation, clustering, state extraction, surface galleries |
| `prepare_*.py` | input generation per campaign (`prepare_slab.py`, `prepare_fl_optimized.py`, `prepare_all_fragments.py`, …) |

---

## 7. Adjacent work

| path | size | what |
|---|---|---|
| `complex/` | 18 GB | PARP9/DTX3L complexes — clustering, docked states, binding analysis |
| `bindcraft_md/` | 4.1 GB | binder design campaign vs the macrodomains; `sweep/run_queue.sh` drives 9 targets sequentially |
| `representative_frames/` | 18 MB | extracted centroid structures per TICA state |
| `states/` | 2 MB | backmapped state structures |
| `archive_analysis/` | — | superseded scripts, kept so old figures reproduce. **Contains the retired active-site numbers; flagged in-file, do not copy them.** |
| `../parp14/` | — | Phase 1: AF3 structure generation, `PARP14.fasta` (UniProt Q460N5), `input/active_sites.yaml`, `alphafold_outputs/` (gitignored) |

---

## 8. Things worth knowing before you trust an artifact

- **Nothing large is in git.** `data/`, all trajectory directories and
  `figures/` are gitignored. There is no second copy.
- **Anything produced before 2026-09-21 that keyed on *catalytic* residues used
  the wrong numbers.** `data/active_site_stats.npz` has been refreshed; the
  per-set `figures/by_sim/*/02_main_analysis/*_active_sites.png` were being
  regenerated at the time of writing. Accessibility results (Q3, Q4) were never
  affected — they key on `pocket_residues`, which were already correct.
- **`figures/_archive/` and `figures_5ns/` are from the 5 ns protocol**, long
  superseded by the 20 ns / 1 µs runs. Don't mix them into current comparisons.
- **Older filenames spell the MD1 unit `md1l1`.** Resolved by an alias; not a
  sign of staleness.
- **Trajectories written by the opportunistic slab runner can contain
  overlapping frames** where a run was preempted between checkpoints. This does
  not bias c_sat but does break frame-index-to-time. See `slab/README.md`.
