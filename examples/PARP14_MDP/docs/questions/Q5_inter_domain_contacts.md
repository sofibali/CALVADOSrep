# Q5 — Do PARP14's domains contact each other, or are they independent beads on a string?

**Status:** ANSWERED. Includes a **retracted earlier claim** — read the correction.

---

## The question

Does PARP14 have a compact, defined domain architecture with real inter-domain
interfaces, or is it a set of folded domains on flexible linkers that rarely
touch?

## Why it matters

This determines what kind of molecule PARP14 is, and therefore what can be done
with it. A protein with stable interfaces has a defined surface to target; a
diffusive beads-on-a-string protein does not, and any binder campaign must
*create* an interface rather than stabilise one ([[Q8]]). It also decides whether
a Markov state model is even the right tool ([[Q2]]).

## The answer

**The macrodomains are largely non-associating.** From `bindcraft_md/README.md:17-33`:

| pair | COM–COM distance (nm) | contact frequency |
|---|---|---|
| MD1–MD2 | 3.16 ± 0.39 | ~0.006 |
| MD2–MD3 | 3.65 ± 0.68 | ~0.003 |
| MD1–MD3 | 3.81 ± 0.79 | ~0.004 |
| *intra-domain reference* | — | *~0.12–0.13* |

Inter-domain contact frequency is **~20–40× lower than intra-domain**. The three
macrodomains behave as dynamic, largely non-associating units.

**The most persistent interface in the protein is KHb-KH8–WWE**, not anything
among the macrodomains (`analysis_summary.txt:121-126`, ~0.01 contact frequency).

**Truncation pulls the triad apart.** From `analysis_summary.txt:99-105`: in the
macrodomain-only construct the three macrodomains sit in a tight cluster with all
pairwise distances ~3–4 nm; in full-length all distances increase by 50–100 %,
with MD1–MD3 going from 3.8 to 7.3 nm. WWE–ART is the tightest pair in `core`
(2.9 nm) — consistent with the ART/WWE adjacency seen independently in [[Q4]].

MD1–MD2 contacts are ~5× more frequent in the smaller constructs (0.005–0.006)
than in full-length (0.001).

## ⚠ Correction — an earlier claim was wrong

`bindcraft_md/README.md:30-33` retracts it explicitly:

> an earlier draft reported "MD1–MD2 permanently in contact (0.38 nm)". That was
> wrong — it came from the `fl_optimized` set and measured the `MD1L1` block
> *including the flexible L1 linker*.

**The lesson generalises:** `MD1L1` is MD1 *plus* the linker to MD2. Any
distance or contact computed on the `MD1L1` block is contaminated by that
linker and will appear far more associated than MD1 really is. Use trimmed
domain cores for contact work.

## How it was tested

- `sim_analysis/figure_md_distances.py` — inter-domain COM–COM (or `--metric min`
  CA–CA) distance violins, replicates pooled, 1.0 nm contact cutoff
- `sim_analysis/compare_full_sims_*.py` — cross-set contact-fraction bars with
  pairwise FDR correction
- Statistics conventions are documented at `README.md:200-222`. Note the **N=5
  floor** on the Mann–Whitney p-value: with 5 replicates per group the smallest
  attainable p is 2/C(10,5) = 0.0079, so no comparison at that depth can be
  significant past correction at α = 0.001.

## Caveats

- `analysis_summary.txt` describes the older 6 constructs at 25 × 20 ns, while
  the bindcraft/TICA numbers come from the longer sets. The two **agree on
  direction**, but do not mix the exact numbers.
- Several recent dated dirs under `figures/04_md_distances/` are empty; the most
  recent non-empty content is `2026-08-14/`.

## Outputs

- `data/md_distances_min.npz`, `data/*_dist_*.npy`
- `figures/04_md_distances/2026-08-14/`
- Prose: `analysis_summary.txt:99-126, 257-267`; `bindcraft_md/README.md:17-33`
