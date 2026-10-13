---
name: slab-coexistence-qc
description: "Validate a CALVADOS slab / direct-coexistence run BEFORE interpreting it — classify the phase (dispersed / arrested / liquid), prove the inter-chain contacts actually exchange, and check c_sat, PBC and restraint counts. Use before reading any c_sat number, inter-chain contact map, or ±RNA / ±ligand comparison from a multi-chain run; also when a slab result looks surprisingly strong, surprisingly reproducible, or when a contact map shows one dominant domain pair."
---

# Slab coexistence QC

A multi-chain slab run produces numbers whether or not they mean anything. Every
check here exists because skipping it produced a wrong answer at least once in
the PARP14 campaign. Run them **before** interpreting, not after.

Order matters: 1 and 2 can invalidate everything downstream.

---

## 1. Do the inter-chain contacts exchange? (the one that matters most)

**Without this check, a dispersing run's contact map is an initial condition
dressed up as an interaction.**

`slab_eq` compresses every chain into a ~6-7 nm slab before production. If the
construct then disperses, chains that were threaded through one another at the
end of equilibration *cannot* separate — the restrained domains are rigid and
cannot pass through each other. Those leftovers never exchange, never decay, and
dominate the raw contact map.

Real case (PARP14 run 2, `fl` homotypic): inter-chain MD2–ART contacts were
**6,796 at frame 0 and 6,295 at frame 1,999** — flat over 2 µs — while every
*other* inter-chain contact decayed ~10×. 85% of all inter-chain contact, 80% of
it literally the same bead pairs, confined to 12 *consecutive-index* chain pairs
(neighbours in the initial slab placement). Read naively it was a strong,
seed-reproducible MD2–ART interaction. It was packing.

```bash
python slab/analyze_contact_exchange.py --all --workers 12
```

Metric: `retention(dt) = |pairs(t) ∩ pairs(t+dt)| / |pairs(t)|` at dt ≈ 10 ns.

| retention @10 ns | verdict |
|---|---|
| **< ~0.25** | liquid — contacts re-partner; the contact map is meaningful |
| **> ~0.4** | frozen — the map is an initial condition; do not interpret it |

PARP14 run 2: `md_full` 0.07–0.15 (liquid) vs `fl` / `core` / `kh1_art_full`
0.42–0.75 (frozen). Only `md_full`'s map was reported.

**Tell-tales that you are looking at a kinetic trap**, even before running the
check: one domain pair sitting 50–100× above all others; a standard deviation
across independent seeds of <1% (real association is noisy); a contact count
that is flat from frame 0; contacts concentrated in *consecutive chain indices*.

## 2. Dispersed, arrested, or liquid?

A slab that never disperses looks identical whether it is a condensate or a
jammed solid. `slab/analyze_slab.py` classifies it:

```bash
python slab/analyze_slab.py <run_dir>          # or --all
```

- **dispersing** — slab width grows to a large fraction of the box. `c_sat` does
  not exist to be measured; an empty `c_dilute` column is the *result*, not a
  failure. The log says `NO PHASE SEP` explicitly — check for that string rather
  than inferring from a blank cell.
- **arrested** — width sd ≈ 0.00, z-order correlation ≈ 1, very high density.
  `slab_eq` compressed a solid into existence. PARP14 `md_full` at 406 mM:
  sd 0.00, order 0.98, ~1020 mg/mL — it would otherwise have been reported as
  the campaign's first condensate.
- **liquid** — width plateaus but fluctuates (sd ~1 nm), chains rearrange past
  one another (order correlation ~0.5), and check 1 shows contacts exchanging.

**Both an arrested jam and a true condensate hold their slab. Only check 1
separates them.**

## 3. Did the c_sat fit converge?

Use the repo's `calvados.analysis`, not the installed copy: only this checkout
has the fix moving `fit_profile`'s "NOT CONVERGED" check *ahead* of its
`return`. In the installed copy that check is dead code after a return, so a
failed tanh interface fit never warns and `c_sat` comes back silently wrong.
`analyze_slab.py` forces the repo copy onto `sys.path` and prints which one it
loaded — confirm that line.

## 4. Do restraint counts scale with the number of chains?

CALVADOS custom restraints carry an explicit **chain-copy index** and are **not**
broadcast across copies: `Sim.map_custom_restraints` resolves each line to one
bead pair via `start_bead + (copy-1)*nbeads + (bead-1)`. A file written for a
monomer (`nmol=1`) restrains **chain 1 only** in a 47-chain slab.

In the build log, `Number of custom restraints:` must equal
**npairs × nmol**, not npairs. PARP14: 141 × 47 = 6,627.

## 5. Periodic boundaries

Contacts must be computed with minimum-image. Wrap coordinates into `[0, L)` and
pass `boxsize=` to `scipy.spatial.cKDTree`; otherwise every contact across an
x/y boundary is dropped. If you compute a distance yourself for a sanity check,
apply MIC — a raw `norm(a - b)` on wrapped coordinates reads boundary-crossing
pairs as tens of nm apart.

OpenMM's DCD writer wraps whole *molecules*, so chains are never split and
analyses need no unwrapping; only visualisation does
(`sim_analysis/unwrap_traj.py`).

## 6. Structural zeros in a ±ligand panel

The control arm of a ±RNA (or ±ligand) panel contains no ligand, so its
"0 contacts" rows are **structural zeros, not measurements**. Averaging both
arms halves every per-residue number and inflates its SD to ~100% of the mean —
which reads convincingly as seed-to-seed noise. Filter to the arm that has the
ligand, or emit no rows for runs without it.

---

## Reporting

State the phase classification and the retention number alongside any contact
map or c_sat you report. If retention is high, say the map is not interpretable
rather than quietly omitting it — the frozen contacts are a real property of the
protocol and worth recording.
