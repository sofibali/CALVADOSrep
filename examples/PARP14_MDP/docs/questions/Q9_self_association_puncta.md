# Q9 — Which constructs self-associate into puncta with RNA and ADP-ribosylated substrate?

**Status:** OPEN — **all 10 active slab runs complete (2026-09-24): none show liquid-liquid phase separation.** Experimental
data ~months out (cell lines in construction as of Sept 2026).

---

## The question

PARP14 forms puncta in cells. Which domain-deletion constructs do so, and is it
driven by multivalent RNA binding, multivalent ADPr reading, or the catalytic
domains that create and remove the marks?

## Why it matters

This is the one question with a **large experimental counterpart**: a pooled
plasmid library of ~2042 constructs, cells sorted into DIFFUSE vs PUNCTATE
pools, each pool sequenced. The punctate-vs-diffuse log-ratio gives a
per-construct self-association score. Single-construct controls exist for
full-length plus catalytic-dead **ART**, **MD1** and **MD2** mutants.

Those controls map one-to-one onto a mechanistic decomposition:

| module | role | catalytic-dead control removes |
|---|---|---|
| RRMs + KH domains | RNA-binding valence | — |
| MD2, MD3, WWE | ADPr-reader valence | MD2-dead → mark *reading* |
| ART | writer | ART-dead → mark *creation* |
| MD1 | eraser | MD1-dead → mark *removal* |

Each predicts a different direction of change, so the screen can discriminate
between them rather than merely rank constructs.

## The hypothesis under test

Puncta are a multivalency/percolation phenomenon: PARP14 cross-links RNA through
its RNA-binding modules and cross-links ADP-ribosylated substrate through its
reader modules, while its own catalytic domains change how many marks exist to
be read. Cross-linking needs **both** arms, so the natural form is a *product*
of the two valences, not a sum.

## How it is being tested

### Model side — two pieces

1. **Sequence/architecture features for all 2047 combinations**
   (`sim_analysis/predict_puncta_propensity.py`). Emits *features*, not a
   hand-tuned score: valences, domain presence, length, FCR, NCPR, κ, SCD,
   aromatic/Arg content. Weights are **fit** against the measured enrichment via
   `--fit` with 5-fold CV — that is what tests the hypothesis rather than
   assuming it.

2. **Multi-chain slab simulations** (`prepare_slab.py`) → saturation
   concentration c_sat, the actual quantitative phase-separation observable.
   Every PARP14 simulation to date is single-chain, so nothing currently on disk
   can predict self-association. Two arms prepared and validated
   (`Sim.build_system()` passes):
   - **homotypic** — construct alone; does it self-associate at all?
   - **+RNA** — with polyU40; tests the RNA-multivalency leg.

### Panel

10 constructs spanning RNA-binding valence × ADPr-reader valence, using only
constructs that already have structures. Includes three matched ±ART pairs.
Coverage is RNA valence {0, 2, 3, 9, 12} — **gap at 4–8**, and nothing at
ADPr-reader 0. Filling either needs new AF3 predictions.

## First result — the slabs dissolve

Four of the twenty slab runs have finished the full 2e8 steps (2 µs):
`fl` and `core_full_go`, both arms. **None of them forms a stable condensate.**

`slab_eq` compresses every chain into a thin slab and then removes that pull.
If the construct self-associates the slab persists. Instead it disperses,
monotonically, for the whole 2 µs — measured as the z-range holding 90% of the
beads:

| run | chains | box z (nm) | after `slab_eq` | production start | production end |
|---|---|---|---|---|---|
| `homotypic/core_full_go` | 47 | 455 | 11 nm | 80 | **193** |
| `homotypic/fl` | 30 | 520 | 6 nm | 65 | **130** |
| `rna/core_full_go` | 47 | 455 | 11 nm | 57 | **108** |
| `rna/fl` | 30 | 520 | 6 nm | 46 | **78** |

So **c_sat is not defined for these runs** — there is no dense/dilute
coexistence to measure. `SlabAnalysis`'s tanh interface fit returns NaN and
nonsensical cutoffs, because it is fitting interfaces that do not exist.
`slab/analyze_slab.py` now checks slab stability first and reports
`NO PHASE SEP` rather than quoting a meaningless number.

**This is a result, not a failure.** Read literally: at 293 K, 150 mM ionic
strength, in CALVADOS3, neither full-length PARP14 nor the KH7a-to-ART core
self-associates strongly enough to hold a condensate together at the simulated
concentration.

### `md_full` holds together — but it is arrested, not a condensate

`md_full` (MD1+MD2+MD3, 83 chains, the only high-concentration run at
**406 mM residues**) is the one construct whose slab does **not** disperse. It
sits at 14.2 nm for the entire 2 µs.

That is not evidence of phase separation. A slab that never disperses looks the
same whether it is a liquid condensate or a jammed solid, and this one is
jammed:

| | `md_full` | the runs that dispersed |
|---|---|---|
| slab width over 2 µs | 14.2 ± **0.00** nm | 107–134 ± 33–46 nm |
| chain z-order correlation, first vs last frame | **0.982** | 0.21–0.37 |
| per-chain displacement over 2 µs | 0.55 nm | — |
| dense-phase density | **~1020 mg/mL** | — |
| beads in the dilute phase | **exactly 0** | — |

A liquid condensate fluctuates in width and lets chains diffuse past one
another. Here the width does not vary *at all*, no chain ever changes places,
and the density is ~3× a physiological condensate. `slab_eq` compressed 83
chains into a 25×25 nm cross-section and the result vitrified. **c_sat is not
defined for it**, any more than for the runs that dissolved.

Whether this is "so cohesive it glasses" or "an artifact of over-compression in
a small box" cannot be separated from one run. Either way it is not LLPS.
`analyze_slab.py` now tests width fluctuation and chain-order scrambling and
labels such runs `TRAPPED` rather than reporting a c_sat.

### `kh1_art_full` behaves like the rest

Both arms disperse (170 nm homotypic, 146 nm +RNA from a 6 nm slab), with RNA
again retarding dispersal — 0.86×, the same direction as every other pair.

### Matched pair C: ART and RNA are not independent

`fl_wwe_full_go` (−ART) finished on both arms, completing the first ±ART pair.
Dispersal width at the end of production (smaller = more cohesive):

| | +ART (`fl`) | −ART (`fl_wwe_full_go`) | effect of removing ART |
|---|---|---|---|
| homotypic | 130 nm | 176 nm | 1.35× more dispersed |
| +RNA | **78 nm** | 182 nm | **2.33× more dispersed** |
| effect of adding RNA | **1.67× less dispersed** | 1.03× — none | |

Removing ART makes the protein disperse faster in both arms. But the
interesting part is the interaction: **RNA retards dispersal only when ART is
present.** With ART, adding polyU tightens the assembly 1.67×; without ART,
RNA does nothing at all (176 → 182 nm, within noise).

Neither condenses — every one of the six still dissolves — so this is a
difference in *cohesion*, not in phase behaviour.

**Composition suggests why.** There is no ADP-ribosylation in CALVADOS; ART here
is simply 199 more residues of a particular composition. Those residues are
unusual for this protein:

| region | length | %K+R | net charge | %aromatic |
|---|---|---|---|---|
| ART 1603–1801 | 199 | 10.6 | **+3** | **13.1** |
| everything except ART | 1602 | 12.9 | −8 | 7.3 |

ART is *not* more Arg/Lys-rich (0.82× the rest), so this is not simple
electrostatic RNA binding. It is the protein's only **net-positive** module,
sitting in an otherwise acidic chain, and it is **1.8× richer in aromatics** —
both of which raise cohesion under CALVADOS's Ashbaugh–Hatch λ (high for F/W/Y)
and favour association with RNA.

**Caveat:** the pair is not perfectly matched. `fl_wwe_full_go` runs at 99 mM
residues vs `fl`'s 108 mM (9% lower) with 31 vs 30 chains, and lower
concentration disperses faster on its own. A 9% offset is unlikely to produce a
35% difference, and cannot explain the 2.33× in the RNA arm at all — but the
homotypic comparison in particular is not clean.

### The one suggestive signal

**The +RNA arm dissolves consistently more slowly than homotypic** — 108 vs
193 nm for `core_full_go`, 78 vs 130 nm for `fl`, a ~1.7× narrower distribution
in both matched pairs. That is the direction the RNA-multivalency leg of the
hypothesis predicts: polyU cross-links chains and retards dispersal, just not
enough to stabilise a phase. Two pairs is not a trend, but it is the first
multi-chain evidence in the project and the remaining 16 runs test it directly.

### Caveats that matter before reading too much into this

- **The slabs are still expanding at the end of production.** 2 µs is enough to
  show the direction unambiguously but not to reach the equilibrium dispersed
  state, so no *upper bound* on c_sat can be quoted either.
- **Concentration is a design choice.** Chain counts were set to hit ~50k beads
  per run, not to hit a physiological concentration. A construct that does not
  condense here might at higher concentration.
- **Force-field scope.** CALVADOS3 is parameterised for single-chain dimensions
  and IDP phase behaviour; the folded-domain interactions that might drive
  PARP14 assembly are represented only through the same residue-level
  Ashbaugh-Hatch/Debye-Hückel terms, with domains held rigid by restraints.
- The RNA arm additionally mixes force fields (CALVADOS2 RNA beads with
  CALVADOS3 protein), as noted below.

## What is blocking

- **Experimental data**: sort-seq counts not yet available. Ingestion is already
  scaffolded — `--fit` accepts `punctate_reads`+`diffuse_reads` or a precomputed
  `enrichment`, joined on `barcode` / `combo_key` / `combo`.
- **Compute**: slab work is **GPU-only** — CPU throughput on pollux is 5.76e5
  bead-steps/s, which puts one 100-chain slab at ~2 years. **Now measured on
  lyra** (4× L40S, no SLURM), 2026-09-17: ~5,000 steps/s on an idle card, i.e.
  **~11 h per run and ~9.5 GPU-days for the full panel × both arms** — not the
  ~79 GPU-days previously estimated from an assumed 50× speedup (the real
  speedup is ~350×). All 20 configs are verified to build and run on CUDA.
  Contention is the binding constraint, not throughput: sharing a card with
  another job costs 3.4–6×, so wall-clock ranges from ~2.5 days (4 idle GPUs) to
  weeks. See `slab/README.md` for the launch policy and per-GPU split.

## Known limits of the model side

- **No ADPr-substrate arm is buildable as a config change.** `PTMProtein` has no
  example and no test, sets `c_termini` to the last PTM bead
  (`calvados/components.py:737`), and never reads from PDB — so it is
  incompatible with `restraint: True`, which every multi-domain PARP14 construct
  needs. No ADP-ribose bead parameters exist anywhere. Viable route: put the
  ADPr beads on a short **unrestrained substrate peptide**.
- RNA is a generic 2-bead/nucleotide model with **no base identity**, so polyU is
  the only meaningful choice; its parameters were fit alongside CALVADOS2 while
  the proteins here are CALVADOS3.
- Construct length correlates strongly with RNA valence across the panel, so the
  two will be hard to separate in the regression.

## Reproduce

```bash
python prepare_slab.py --arm both --panel full   # already done; inputs are on disk
# on lyra (no SLURM -- submit_slab.slurm does not apply there):
cd slab
GPU=0 nohup ./run_slab_queue.sh        homotypic > logs/queue_homotypic.log 2>&1 &
GPU=1 nohup ./run_slab_opportunistic.sh rna      > logs/opp_rna.log         2>&1 &
python monitor_slab.py --watch

cd sim_analysis
python predict_puncta_propensity.py --prior --workers 48
python predict_puncta_propensity.py --fit <counts.csv>    # when data exists
```

## Outputs

- `data/puncta_features.csv` — 2047 constructs × 39 features (built, validated)
- `slab/{homotypic,rna}/` — prepared inputs; `run_slab_queue.sh` (lyra),
  `submit_slab.slurm` (a SLURM cluster, if one is ever used)
