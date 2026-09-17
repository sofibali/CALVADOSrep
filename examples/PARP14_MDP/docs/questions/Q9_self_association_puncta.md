# Q9 — Which constructs self-associate into puncta with RNA and ADP-ribosylated substrate?

**Status:** OPEN — simulations prepared and validated, not yet run. Experimental
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

## What is blocking

- **Experimental data**: sort-seq counts not yet available. Ingestion is already
  scaffolded — `--fit` accepts `punctate_reads`+`diffuse_reads` or a precomputed
  `enrichment`, joined on `barcode` / `combo_key` / `combo`.
- **Compute**: slab work is **GPU-only**. Measured CPU throughput on pollux is
  5.76e5 bead-steps/s, which puts one 100-chain slab at ~2 years. The runs target
  the SLURM GPU server that shares this filesystem. Estimated ~4 GPU-days per
  run, ~79 GPU-days for the full panel × both arms — at an **assumed** 50× GPU
  speedup that must be replaced with `prepare_slab.py --benchmark` before
  committing.

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
python prepare_slab.py --benchmark              # calibrate GPU first
python prepare_slab.py --arm both --panel core
# on the GPU server:
sbatch --array=0-4 slab/submit_slab.slurm slab/homotypic

cd sim_analysis
python predict_puncta_propensity.py --prior --workers 48
python predict_puncta_propensity.py --fit <counts.csv>    # when data exists
```

## Outputs

- `data/puncta_features.csv` — 2047 constructs × 39 features (built, validated)
- `slab/{homotypic,rna}/` — prepared inputs; `slab/submit_slab.slurm`
