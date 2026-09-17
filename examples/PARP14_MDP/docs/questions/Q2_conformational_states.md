# Q2 — Does PARP14's macrodomain module occupy discrete conformational states?

**Status:** ANSWERED — and the answer is an honest **negative**.

> Full derivation, theory and diagnostics: [`../../TICAresultDescreption.md`](../../TICAresultDescreption.md).
> This file is the short version; that file is the primary source.

---

## The question

Two questions, in order:

1. **Which coordinate** captures the slow motions of the MD1–MD2–MD3 module?
2. Given that coordinate, are there **discrete metastable states** — basins the
   protein sits in and hops between — or a continuum?

## Why it matters

The whole downstream design pipeline depends on the answer. If there are real
metastable states, "representative structures" means kinetically distinct
conformers worth targeting individually ([[Q8]]). If there is not, extracted
clusters are just samples along a continuum, and claiming they are states would
be wrong.

## How it was tested

On `md_full` (FL 790–1388 incl. the 1194–1206 linker), 25 replicates extended to
100 ns (2.5 µs total):

- **Feature selection** by cross-validated VAMP-2 at lags 0.8 / 2.0 / 4.0 ns
- **Metastability** by the implied-timescale (ITS) test — a valid Markov model
  requires timescales that *plateau* with lag time

## The answer

### 1. The featurization question has a clean answer: `interface_ca`

| feature | 0.8 ns | 2.0 ns | 4.0 ns | reading |
|---|---|---|---|---|
| **interface_ca** | 2.18 | 1.74 | 1.56 | high **and stable** → real slow modes |
| ca_stride25 | 4.12 | 1.71 | 0.68 | high then collapses → mostly fast |
| linker_ca | 0.41 | 0.007 | 0.006 | fast noise |
| com / orient | ≤0.24 | ~0 | ~0 | uninformative |

Domain-edge contacts win unambiguously, and the choice is robust to trajectory
length.

### 2. The metastable-states question has a negative answer

Implied timescales for `interface_ca` rise almost perfectly **linearly** with
lag (IC1 ≈ 8.3 × τ; λ₁ ≈ 0.89 at *every* lag from 0.8 to 32 ns). No plateau.
Extending 5 ns → 100 ns did not produce one; it just pushed the apparent
timescale up in lock-step.

**Physically**, combined with the finding that the macrodomains have no
persistent inter-domain interface ([[Q5]]: contact frequency ≈ 0.005, COM–COM
fluctuating 0.4–0.8 nm), this says MD1–MD2–MD3 behaves **less like a few
metastable states separated by barriers and more like a flexible, diffusive
multidomain chain** — macrodomains on flexible linkers, closer to a
semi-flexible polymer with a continuum of slow modes than to a two- or
three-state switch.

The MSM framework *assumes* metastable basins. When the landscape is shallow,
the diagnostics correctly refuse to certify states; no amount of re-tuning τ or
*k* manufactures basins that are not there.

## Consequence for downstream work

Extracted clusters are **structural representatives sampled along a continuum**,
not kinetically distinct states. That is still fine for the binder-design goal of
*enforcing a new interface* ([[Q8]]) — deliberately picking diverse poses to
clamp needs **geometric diversity**, not kinetic metastability. The claim that
these are long-lived states the protein "sits in" is not made.

## Outputs

- `data/its_results_*`, `data/vamp_*`, `data/cluster_*`, `data/msm_*`
- `figures/05_clustering/`
- `tica_pipeline/` (staged scripts + its own README)
