# Q9 — Which constructs self-associate into puncta with RNA and ADP-ribosylated substrate?

**Status:** OPEN — **first 4 slab runs complete (2026-09-22): none of them phase-separate.** Experimental
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
