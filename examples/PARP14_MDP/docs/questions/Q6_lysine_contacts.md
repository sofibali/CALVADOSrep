# Q6 — Which lysine pairs are close enough to form a DSS crosslink?

**Status:** ANSWERED for full-length (2026-09-16). Method corrected and rerun; other sets still to do.

---

## The question

Which lysine pairs sit within DSS crosslinking reach across the conformational
ensemble, and which of those are actually *formable* — i.e. reachable by a
crosslinker that has to travel around the protein rather than through it?

## Why it matters

This is the model side of an XL-MS experiment. Predicted crosslinks are
structural restraints that can be tested at the bench, and they are one of the
few PARP14 observables where a CG ensemble makes a directly falsifiable
prediction. Long-range crosslinks would also speak to the inter-domain geometry
question ([[Q5]]).

## The cutoff is chemistry, not a loose contact radius

This was previously mis-documented and is worth stating precisely:

```
DSS spacer arm                ~15 Å
+ 2 × Lys side-chain reach
--------------------------------
= ~30 Å = 3.0 nm  CA–CA ceiling
```

**3.0 nm is correct and deliberate** — it is the standard DSS CA–CA ceiling in
the XL-MS field. A CG model with one bead per residue has no side chains to
resolve, so a permissive CA–CA ceiling is the right level of description.

The script's docstring used to claim 1.6 nm / 0.8 nm. **The docstring was wrong;
the code was always right.** Fixed 2026-09-16.

## Two corrections the raw cutoff needs

### 1. Sequence separation — the threshold is derived, not chosen

At a 3.0 nm ceiling, residues close in sequence are within reach from chain
connectivity alone: they satisfy the cutoff in *every* frame regardless of
structure, so they dominate the persistence ranking while carrying zero
tertiary-structure information.

The cutoff for "too close" follows from chain geometry. A CG chain has a uniform
**0.38 nm** CA–CA bond, so residues n apart cannot exceed `n × 0.38` nm even
fully extended:

```
ceil(3.0 / 0.38) = 8     first separation that can exceed the ceiling at all
                         ... but only at 8 × 0.38 = 3.04 nm, a 1% margin that a
                         flexible chain effectively never reaches
+ 1              = 9     first separation that is informative in practice
```

Measured against the unfiltered `fl` data, the cliff falls exactly there:

| \|Δresid\| | fraction at persistence 1.0 |
|---|---|
| 1–8 | **100 %** (forced by connectivity) |
| 9 | 5 % |
| 10+ | **0 %**, mean falls smoothly 0.97 → 0.51 by 18 |

`MIN_SEQ_SEP` is therefore computed as
`ceil(LYS_LYS_CUTOFF / CG_BOND_LENGTH) + 1` = **9**, and tracks the cutoff if it
ever changes (2 nm → 7, 4 nm → 12).

> **Do not raise this casually.** An earlier run used 20, picked post-hoc to sit
> above an artefact band that had been *observed* rather than derived. It
> discarded the 10–19 separations, which hold 245 `fl` pairs with mean
> persistence 0.71 spanning 0.36–0.98 — real structural signal, thrown away.
> Raise it only to deliberately restrict to long-range restraints.

### 2. Euclidean distance overestimates reachability

A straight line between two lysines can pass straight through the protein core,
which a crosslinker cannot do. The physically meaningful quantity is the
**solvent-accessible surface distance (SASD)** — the shortest path that stays
*outside* the protein, which is what Xwalk computes.

`--sasd` now implements this (`compute_sasd`):

- voxelize a local box around the two lysines (0.15 nm grid)
- mark a voxel blocked if within 0.30 nm of any bead other than the endpoints
- Dijkstra with 26-connectivity over the free voxels
- return the geodesic length, or ∞ if no path within the DSS ceiling exists

Verified behaviour: free space → equals Euclidean; a wall between the endpoints →
path detours (2.571 nm vs 2.000 Euclidean); a fully enclosed endpoint →
unreachable. ~0.15 s per pair.

**Pairs that pass the Euclidean test but fail SASD are buried false positives**,
and for a globular multidomain protein they are common.

Approximation to state plainly: a CG bead is a whole residue, so this "surface"
is coarser than Xwalk's all-atom one. Treat SASD here as a **reachability
filter**, not a quantitative distance.

## First corrected run — full-length, 2026-09-16

```bash
python analyze_lys_contacts.py --set fl --sasd --sasd-top 150   # min-seq-sep 9 by default
```

25 replicates. Of the 150 most persistent pairs re-scored by SASD:

| | count |
|---|---|
| solvent path within the 3.0 nm DSS ceiling | **135** |
| **buried false positives removed** | **15** |

Sequence separation of the re-scored set: min 9 (the derived filter), median 47,
max 188.

### Top predictions — persistent *and* solvent-reachable

| Lys i | Lys j | seq sep | Euclid persistence | SASD min (nm) |
|---|---|---|---|---|
| 942 | 979 | 37 | 1.000 | 0.724 |
| 839 | 871 | 32 | 1.000 | 0.944 |
| 21 | 80 | 59 | 1.000 | 1.024 |
| **1223** | **1364** | **141** | 1.000 | 1.024 |
| 26 | 80 | 54 | 1.000 | 1.094 |
| 1125 | 1162 | 37 | 1.000 | 1.120 |
| **1223** | **1365** | **142** | 1.000 | 1.182 |
| 1106 | 1154 | 48 | 1.000 | 1.204 |

The K1223-K1364/K1365 pairs are the most interesting: 141-142 residues apart,
permanently within reach, with a real solvent route. Those are genuine
long-range structural restraints and the most directly testable predictions here.

### What SASD removed

Pairs that pass the Euclidean test but have no solvent route — the straight line
runs through a domain core:

| Lys i | Lys j | seq sep | Euclid persistence |
|---|---|---|---|
| 828 | 912 | 84 | 1.000 |
| 871 | 955 | 84 | 1.000 |
| 1020 | 1078 | 58 | 1.000 |
| 1020 | 1154 | 134 | 1.000 |
| 1078 | 1141 | 63 | 1.000 |
| 1141 | 1158 | 17 | 1.000 |

All sit within a single macrodomain, which is exactly where a straight-line
distance is most misleading. **Every one would have been reported as a confident
crosslink by the Euclidean pass alone.**

`reachable = 1` means a solvent path within the DSS ceiling exists in at least
one sampled frame — a crosslink only needs one accessible conformer, which is
why the minimum over frames is the right statistic, not the mean.

Output columns: `resid_i, resid_j, seq_sep, euclid_persistence, sasd_mean_nm,
sasd_min_nm, n_frames_reachable, reachable`.

### Superseded

The 10,983-row `data/lys_contacts/all_sets_lys_acidic_summary.csv` predates both
corrections — no sequence-separation filter, no SASD. **Do not rank from it.**

## To close this question

1. Run the remaining sets — only `fl` has been regenerated:
   `python analyze_lys_contacts.py --set <name> --sasd` (the derived
   `--min-seq-sep 9` is the default; do not override it without reason).
2. Cross-check the survivors against RSA ([[Q7]]) — a lysine buried by RSA
   should not be appearing as crosslinkable.
3. Take the long-range pairs (K1223-K1364/K1365 at 141-142 separation first) to
   the bench as testable restraints.
4. If XL-MS data exists, compare satisfied vs violated predictions — a direct
   test of whether the CG ensemble has the inter-domain geometry right.

## Outputs

- `data/lys_contacts/all_sets_lys_acidic_summary.csv` — **superseded**, see above
- `data/{set}_lys_sasd.csv` — SASD-filtered pairs (new)
- `figures/07_lysine_contacts/2026-08-12/`. Note `latest -> 2026-09-04` points at
  an **empty** directory.
