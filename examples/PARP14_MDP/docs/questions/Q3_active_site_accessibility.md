# Q3 — Are PARP14's catalytic sites sterically reachable, and does that depend on domain context?

**Status:** ANSWERED (2026-09-16) · 19 constructs · 215 replicates

---

## The question

A catalytic pocket that is buried by the protein's own chain cannot engage
substrate no matter how good its chemistry. Across PARP14's five ligand-binding
/ catalytic sites (MD1, MD2, MD3, WWE, ART), how much approach space does each
one actually have — and does deleting domains change it?

## Why it matters

PARP14 is both a writer (ART) and an eraser (MD1) of ADP-ribosylation. If the
domain architecture systematically favours one pocket over the other, that is a
structural bias built into the protein, independent of catalytic rate constants.
It also sets which pockets are plausible drug/binder targets ([[Q8]]).

## How it was tested

Solid-angle accessibility (SAA), `sim_analysis/analyze_accessibility.py`:

- from each pocket's centre of mass, cast **200 rays** on a Fibonacci sphere
- a ray is **blocked** if any non-pocket bead lies within **0.5 nm** of its path
  out to **5 nm**
- SAA = unblocked fraction. 0 = fully occluded, 1 = fully open
- beads within **±10 residues** of the pocket are excluded, so a domain does not
  count as occluding itself

Sampling: `--target-frames 2000`, which picks the stride per trajectory. The
sets span 3–4 k-frame originals and 100 k-frame extensions, so a single fixed
stride cannot serve both.

## The answer

**A strictly invariant ordering, with 0 / 19 violations:**

```
MD1  <  MD2  <  ART  <  MD3  <  WWE
0.355   0.448   0.507   0.576   0.765      (mean SAA)
```

Whichever domains are deleted, the sites keep this rank order. The ordering is a
property of the fold, not of any particular truncation.

Two consequences:

1. **MD1 (the eraser) is always the most buried catalytic site.** See [[Q4]].
2. **WWE is always the most exposed** (mean 0.765, up to 0.891) — consistent
   with its role as a PAR reader that must engage a polymer in trans.

Absolute values shift a lot with truncation (MD1 spans 0.231–0.453 across
constructs); only the *ordering* is invariant. The ranges overlap between sites —
the invariant is the **within-construct** ordering, not a global separation.

## Confidence and caveats

- This ranks **steric opportunity, not rate.** CALVADOS is one bead per residue
  with implicit solvent — no NAD⁺, no ADP-ribose, no transition state.
- Tight replicate spread on the long runs (±0.001 SAA) reflects replicates
  sharing a starting model. That is precision, not validation of the underlying
  AF3/AF2 structure.
- Two sampling regimes exist and should not be mixed (see [[Q0]] conventions):
  the 25-replicate original sets sample *structural* diversity across distinct
  AF3 seeds (±0.04); the 5-replicate extensions sample *conformational* time
  from similar starts (±0.001).
- Analysing every frame rather than 2000 per replicate shifts the mean by
  ≤0.002 SAA — checked directly on `fl`.

## Reproduce

```bash
cd sim_analysis
./backfill_accessibility.sh                  # fills data/accessibility_stats.npz
python predict_enzyme_dominance.py           # tables + figures
```

## Outputs

- `data/accessibility_stats.npz` — per-replicate SAA / cone angle / shell density
- `data/enzyme_dominance_saa.csv` — per-construct per-site means
- `figures/03_accessibility/latest/enzyme_dominance/site_accessibility_heatmap.svg`
