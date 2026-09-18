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
  out to **3 nm**
- SAA = unblocked fraction. 0 = fully occluded, 1 = fully open
- beads within **±10 residues** of the pocket are excluded, so a domain does not
  count as occluding itself

Sampling: `--target-frames 2000`, which picks the stride per trajectory. The
sets span 3–4 k-frame originals and 100 k-frame extensions, so a single fixed
stride cannot serve both.

### The 3 nm probe length is a SELECTED cutoff

It is chosen, not derived — set to match the DSS CA–CA crosslink ceiling used in
[[Q6]], so "reachable" means the same distance in both analyses. A sweep over
1–10 nm (`sweep_maxdist.sh`) shows the choice does not drive anything:

| probe | MD1 | MD2 | ART | MD3 | WWE | ordering |
|---|---|---|---|---|---|---|
| 1 nm | 0.737 | 0.828 | 0.912 | **1.000** | **1.000** | saturated, unusable |
| 2 nm | 0.452 | 0.594 | 0.606 | 0.733 | 0.952 | ✓ |
| **3 nm** | **0.405** | **0.497** | **0.569** | **0.640** | **0.854** | ✓ |
| 5 nm | 0.369 | 0.446 | 0.549 | 0.580 | 0.781 | ✓ |
| 10 nm | 0.345 | 0.416 | 0.532 | 0.541 | 0.745 | ✓ |

(`fl` values.) Three things follow: below 2 nm the metric is **degenerate**
(MD3 and WWE both hit 1.000); the gradient collapses from −0.285 to −0.0025 SAA
per nm, so it is asymptotic and local; and **the ordering is identical at every
distance from 2 nm up**. No conclusion here depends on the choice.

Matching the crosslink ceiling also has an empirical payoff — the correlation
between crosslink count and exposure in [[Q10]] is stronger at 3 nm than at 5 nm
for all five sites.

## The answer

**A strictly invariant ordering, with 0 / 19 violations:**

```
MD1  <  MD2  <  ART  <  MD3  <  WWE
0.384   0.485   0.529   0.623   0.823      (mean SAA, 3 nm probe)
```

Whichever domains are deleted, the sites keep this rank order. The ordering is a
property of the fold, not of any particular truncation.

Two consequences:

1. **MD1 (the eraser) is always the most buried catalytic site.** See [[Q4]].
2. **WWE is always the most exposed** (mean 0.823, up to 0.924) — consistent
   with its role as a PAR reader that must engage a polymer in trans.

Absolute values shift a lot with truncation (MD1 spans 0.261–0.453 across
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
