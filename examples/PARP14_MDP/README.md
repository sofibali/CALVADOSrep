# PARP14 CALVADOS Simulations (Phase 2)

Coarse-grained MD simulations of PARP14 (1801 residues) and domain-deletion/
domain-truncation variants, using CALVADOS. This is **Phase 2** of the
pipeline; Phase 1 (AlphaFold3 structure generation for the full domain
combinatorial library) is in [`../../parp14/`](../../parp14/).

## Project shape

Two families of constructs exist side by side:

1. **Original "short" sets** (`fl`, `fl_optimized`, `norrm`, `noart`, `core`,
   `mka`, `md`, `md3art`) — deletion-style constructs built by
   `prepare_and_run_all.py` from the AF3 combinatorial library, with
   inter-domain linkers excised. 25 replicates each (5 seeds x 5 samples),
   originally 10-100 ns, later checkpoint-extended once (see per-set length
   below). `fl`/`fl_optimized`/`fl_go` are full-length (nothing excised, no
   sequence is missing) but differ in restraint scheme.
2. **Newer "full/contiguous" sets** (`md_full`, `mka_full`, `core_full_go`,
   `kh1_art_full`, `kh1_wwe_full`, `md2_art_full`, `md2_wwe_full`,
   `md3_art_full`, `md3_wwe_full`, `core_wwe_full_go`, `mka_wwe_full`,
   `fl_wwe_full_go`) — sliced directly from the shared FL AF2 model
   (`input/parp14.pdb`), keeping **every** residue in their span including
   inter-domain linkers (nothing excised), 5 replicates x 1 us each. Built by
   `prepare_md_full.py` / `prepare_mka_full.py` / `prepare_core_full_go.py` /
   `prepare_extra_full_constructs.py`. Constructs spanning **both** KH7a
   (738-789) and KHb-KH8 (1389-1533) together also carry a Go-model
   KH7a-KHb custom restraint (141 CA-CA pairs, k=15 kJ/mol/nm^2, reused from
   `fl_go`) since those two KH pieces are sequence-split but spatially
   adjacent.

`fl_go` is a third variant: full-length, Go-model KH7a-KHb restraints, but
run as **11 independent 2 us extensions** seeded from representative TICA
cluster-state frames rather than a uniform 25x-replicate grid (directory
names are `state-N_tica_seed-X_sample-Y_frZZZZ/`, not `seed-N_sample-0/`).

**Redundancy note:** 11 of the 12 new full/contiguous constructs have the
same domain-unit composition as an existing entry in the fragment library
(`fragments/`, 20 ns each) — the new runs extend those to 1 us rather than
covering new sequence territory (`fl_wwe_full_go` is the one exception, no
fragment counterpart). The old short `md`/`mka`/`md3art` sets are similarly
superseded on sampling depth by `md_full`/`mka_full`/`md3_art_full` (same
restraint style, same span, far deeper sampling) — their trajectory files
have been removed to reclaim disk space (config.yaml/components.yaml/
metadata.json/input/ retained for reproducibility). `core`/`norrm`/`noart`
differ from their `_full`/`_go` counterparts by the added Go-model restraint,
so they remain a narrower but genuine comparison point.

## Current simulation inventory

| Set | Span (FL) | Residues | Reps x length | Restraints |
|---|---|---|---|---|
| `fl` | 1-1801 | 1801 | 25 x ~35 ns | harmonic |
| `fl_optimized` | 1-1801 | 1801 | 25 x 100 ns | harmonic + Go-model KH7a-KHb (141 pairs, k=350) |
| `fl_go` | 1-1801 | 1801 | 11 x 2000 ns (TICA-seeded) | harmonic + Go-model KH7a-KHb (141 pairs, k=15) |
| `norrm` | 315-1801 | 1474 | 25 x ~42 ns | harmonic |
| `noart` | 315-1602 | 1275 | 25 x ~42 ns | harmonic |
| `core` | 738-1801 | 1051 | 25 x ~30 ns | harmonic |
| `mka` | 790-1801 | 999 | 25 x (traj. removed, superseded) | harmonic |
| `md` | 790-1388 | 586 | 25 x (traj. removed, superseded) | harmonic |
| `md3art` | 1207-1801 | 595 | 25 x (traj. removed, superseded) | harmonic |
| `md_full` | 790-1388 | 599 | 5 x 1000 ns | harmonic |
| `mka_full` | 790-1801 | 1012 | 5 x 1000 ns | harmonic |
| `core_full_go` | 738-1801 | 1064 | 5 x 1000 ns | harmonic + Go-model KH7a-KHb |
| `kh1_art_full` | 315-1801 | 1487 | 5 x 1000 ns | harmonic + Go-model KH7a-KHb |
| `kh1_wwe_full` | 315-1602 | 1288 | 5 x 1000 ns | harmonic + Go-model KH7a-KHb |
| `md2_art_full` | 1004-1801 | 798 | 5 x 1000 ns | harmonic |
| `md2_wwe_full` | 1004-1602 | 599 | 5 x 1000 ns | harmonic |
| `md3_art_full` | 1207-1801 | 595 | 5 x 1000 ns | harmonic |
| `md3_wwe_full` | 1207-1602 | 396 | 5 x 1000 ns | harmonic |
| `core_wwe_full_go` | 738-1602 | 865 | 5 x 1000 ns | harmonic + Go-model KH7a-KHb |
| `mka_wwe_full` | 790-1602 | 813 | 5 x 1000 ns | harmonic |
| `fl_wwe_full_go` | 1-1602 | 1602 | 5 x 1000 ns | harmonic + Go-model KH7a-KHb |

Plus a **fragment library** (`fragments/`) of 49 (of 66 possible) contiguous
PARP14 sub-sequences at 20-25 replicates x 20 ns each — see
`prepare_all_fragments.py`.

Source of truth for every set's span/units/sysname/color is
[`sim_registry.py`](sim_registry.py)'s `SETS` dict — read from there (or via
`sim_registry.get_units(set_key)` etc.) rather than hardcoding, since several
analysis scripts have had bugs from keeping stale local copies.

## Domain units (FL numbering)

| Unit | FL range |
|---|---|
| RRM1 | 1-145 |
| RRM2 | 146-224 |
| RRM3 | 225-314 |
| KH1-KH6 | 315-737 |
| KH7a | 738-789 |
| MD1L1 | 790-1004 |
| MD2 | 1004-1193 |
| MD3 | 1207-1388 |
| KHb-KH8 | 1389-1533 |
| WWE | 1534-1602 |
| ART | 1603-1801 |

Note the real, un-excised gaps between some units (e.g. 1194-1206 between
MD2 and MD3) — deletion-style sets compress these out; contiguous `_full`
sets deliberately keep them (see `sim_registry.CONTIGUOUS_FL_RANGE`).

## Active sites

Defined in `parp14/input/active_sites.yaml` and mirrored in
`sim_registry.ACTIVE_SITES_FL`:

| Site | Catalytic residues | Basis |
|---|---|---|
| MD1 | 831, 923, 962 | ADP-ribose hydrolase |
| MD2 | 1035, 1046, 1134, 1171 | ADP-ribose reader |
| MD3 | 1248, 1259, 1330, 1371 | ADP-ribose binding |
| WWE | 1539, 1548, 1570, 1576 | iso-ADP-ribose binding (2 Tyr + 1 Phe aromatics, 1 Lys charged) |
| ART | 1684, 1705, 1706, 1722 | H-Y-E catalytic triad |

WWE was added by structural alignment (PyMOL `cealign`, RMSD 1.83 A over 64
residues) of PARP14's WWE domain onto RNF146's WWE domain bound to
iso-ADP-ribose (PDB 3V3L; Wang et al. 2012 *Genes Dev*, PMID 22267412) — no
PARP14-specific ligand-bound structure exists, so treat this as a
structure-derived working hypothesis, not a direct experimental measurement.

## Quick start

### Prepare + launch a new contiguous construct

```bash
conda activate calvados
python prepare_extra_full_constructs.py     # edit CONSTRUCTS list for a new span
bash <construct_key>/run_all.sh parallel
```

Or for the original short-sim style: `python prepare_and_run_all.py` (see
that script for the full construct list and CLI).

### Run the analysis pipeline on a finished set

```bash
cd sim_analysis
bash run_analysis.sh --set <construct_key>        # analyze_all.py (8 modules) + lysine contacts
python figure_md_distances.py --set <construct_key>
python analyze_active_sites.py --set <construct_key>
python analyze_accessibility.py --set <construct_key>
python figure_lys_exposed_persistence.py --sets <construct_key>
python build_sim_dashboard.py --sim <construct_key>   # figures/by_sim/<construct_key>/dashboard.html
```

Works on **any** construct with no code changes via `--set NAME` (registered
in `sim_registry.py`) or `--sim-folder PATH` (arbitrary folder with a
`metadata.json` — see `sim_registry.write_metadata`), including irregular
layouts like `fl_go`'s TICA-seeded state directories.

### Cross-set comparison figures

```bash
# Any subset of sim_registry sets, one shared distance/Rg cache:
python figure_md_distances.py --set md_full mka_full core_full_go \
    --sim-folder-as fl_go_5rep1us:/path/to/fl_go \
    --max-reps fl_go_5rep1us:5 --max-ns fl_go_5rep1us:1000   # length/replicate-matching for fl_go

# Domain x domain distance-distribution grid (reads that cache):
python compare_full_sims_grid.py --sets <sets...>

# % frames in contact bar chart, one bracket per sim-pair vs a reference
# (use --reference past a handful of sets -- C(n,2) brackets is unreadable):
python compare_full_sims_contact_bars_pairwise_fdr.py --sets <sets...> --reference fl_go_5rep1us

# Single metric across sets (rg, ree, or a DOMAIN_DOMAIN pair):
python compare_sims.py --metric rg --sims <sets...>
```

## Repository layout

Scripts are grouped by role. Build/simulation scripts and shared helpers stay
in this root; analysis and the TICA state pipeline live in labeled folders.
**Every script resolves `data/`, `figures/`, and the simulation folders to
this root**, no matter where it lives or is invoked from (a root-bootstrap
handles the shared imports).

| Location | Contents |
|---|---|
| **root** (build/simulation) | `prepare_and_run_all.py` (original 8 short sets), `prepare_md_full.py` / `prepare_mka_full.py` / `prepare_core_full_go.py` / `prepare_extra_full_constructs.py` (contiguous `_full`/`_full_go` sets), `prepare_all_fragments.py` (fragment library), `prepare_fl_optimized.py`, `stamp_metadata.py` |
| **root** (shared) | `sim_registry.py` (sets / domain units / FL<->construct remap / `--sim-folder`), `_fig_layout.py` (figure styling/output routing), `environment.yml` (env `calvados-tica`) |
| **`sim_analysis/`** | `run_analysis.sh` + `analyze_all.py` (8 modules) + `analyze_lys_contacts.py` + `analyze_active_sites.py` + `analyze_accessibility.py` + `figure_md_distances.py` + `figure_lys_exposed_persistence.py` + `figure_sasa_faces.py` + `build_sim_dashboard.py` + `compare_sims.py` + `compare_full_sims_grid.py` + `compare_full_sims_contact_bars.py` / `compare_full_sims_contact_bars_pairwise_fdr.py` + `make_analysis_pse.py` + `ANALYSIS_README.md` |
| **`tica_pipeline/`** | `run_pipeline.sh` + `01/02/03_*.py` + engine (`cluster_states.py`, `fullatom_minimize_states.py`, `backmap_states.py`) + `sweep_its.py` + `make_cluster_pse.py` + `make_surface_gallery.py` (electrostatic/hydrophobic/pocket surface galleries via APBS) + `README.md` |
| **`archive_analysis/`** | superseded / one-off scripts |
| **shared outputs** | `data/` (npz caches, incl. per-replicate arrays for error-bar stats), `figures/by_sim/<set>/` (per-set dashboards + figures) and `figures/comparisons/<category>/<sets-joined>/` (cross-set figures), `input/`, `states/`, `fragments/`; simulation set directories named per the inventory table above |

Setup: `conda env create -f environment.yml && conda activate calvados-tica`
(then `pip install -e /path/to/CALVADOS` only to run new simulations).

### Typical flow

```bash
python prepare_md_full.py; bash md_full/run_all.sh parallel      # 1. build+run   (root)
cd sim_analysis && bash run_analysis.sh --set md_full; cd ..     # 2. analyze     (sim_analysis/)
SET=md_full FEATURES=pose K=4 bash tica_pipeline/run_pipeline.sh # 3. states      (tica_pipeline/)
```

## Statistics conventions (contact-fraction comparisons)

`compare_full_sims_contact_bars*.py` compute "% frames in contact" against
two thresholds (Rg_A+Rg_B "surfaces touching", and +2 nm "AH potential could
still reach") **per replicate**, so the across-replicate SEM is a real error
estimate rather than a pooled-frame statistic. With only N=5 replicates/set,
a Mann-Whitney U test's smallest possible two-sided p-value is fixed at
2/C(10,5) = 0.0079 (perfect separation) — Holm-Bonferroni correcting that
across many repeated panels for one sim-pair (the same 5+5 replicates feed
every domain-pair/threshold combination, so those tests are correlated, not
independent) can push everything to "ns" regardless of true effect size.
Two scripts exist for two different correction scopes:

- `compare_full_sims_contact_bars.py` — corrects the sim-pairs compared
  *within one panel* (Holm-Bonferroni); good for "which sims differ at this
  specific domain pair."
- `compare_full_sims_contact_bars_pairwise_fdr.py` — corrects one sim-pair
  *across all the panels it's repeated in* (Benjamini-Hochberg FDR, since
  Holm at that breadth is unusable at N=5); good for "is sim A different from
  sim B overall." Use `--reference SET` once comparing more than ~4 sets, or
  the bracket count (C(n,2) per panel) becomes unreadable.

---

Phase 1 (AlphaFold3 structure generation) lives in
[`../../parp14/`](../../parp14/); see [`../../docs/PARP14_PROJECT.md`](../../docs/PARP14_PROJECT.md)
for the full pipeline protocol and [`../../CLAUDE.md`](../../CLAUDE.md) for
architectural context.
