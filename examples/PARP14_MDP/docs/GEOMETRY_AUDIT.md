# PBC handling and restrained-domain geometry — 2026-09-24

Two things checked after the question "was PBC centering used for the monomer
analyses?". Answers: **no, and it is not needed** — but the check turned up a
real geometric subtlety in how restrained domains are represented.

---

## 1. Periodic boundaries: analyses are fine, visualisation is not

CALVADOS writes trajectories through OpenMM's `DCDReporter` with
`enforcePeriodicBox` left at its default, which for a periodic system wraps
whole **molecules** — each bonded molecule is translated as a unit, never split.

Verified directly: if a chain were broken across a boundary, some consecutive
CA–CA distance would jump to roughly box size.

| set | box | max consecutive-bead distance | bonds > 1 nm |
|---|---|---|---|
| `fl` | 300 nm | 0.89 nm | 0 |
| `fl_optimized` | 300 nm | 0.91 nm | 0 |
| `norrm` | 250 nm | 0.44 nm | 0 |
| `core` | 120 nm | 0.44 nm | 0 |
| `md_full` | 80 nm | 0.45 nm | 0 |

Never close to box size in any frame. Atom ordering was independently confirmed
by checking 400 known restraint target distances against the trajectory (mean
error −0.001 nm, rms 0.013 nm).

**So Rg, inter-domain distances, contact maps, RMSF, accessibility and the
active-site metrics need no unwrapping.** `mdtraj.compute_distances` also
defaults to `periodic=True`, which gives the same answer for a whole chain in a
large box.

### What does break: watching the trajectory

When the molecule's centre of mass crosses a boundary the whole chain jumps one
box length. Rate depends on box size and length of run:

| set | box | frames | wrap events |
|---|---|---|---|
| `md_full` | 80 nm | 100k | 191 |
| `md1_md2` | 70 nm | 100k | 206 |
| `md2_md3` | 70 nm | 100k | 187 |
| `mka_full` | 100 nm | 100k | 45 |
| `core_full_go` | 120 nm | 100k | 29 |
| `fl`, `fl_optimized`, `norrm`, `noart`, `core` | 120–300 nm | 3–10k | **0** |

It is the **1 µs extension runs in small boxes** that wrap; the original runs in
large boxes never do. Maximum jump equals the box length exactly, confirming
clean whole-molecule wraps.

Fix: **`sim_analysis/unwrap_traj.py`**

```bash
python sim_analysis/unwrap_traj.py <run_dir>            # writes {sysname}_viz.dcd
python sim_analysis/unwrap_traj.py --all --check        # survey without writing
python sim_analysis/unwrap_traj.py <run_dir> --mode unwrap
```

`center` (default) removes the COM every frame so the molecule sits still — best for
viewing. `unwrap` follows the COM continuously instead, preserving real
diffusion; use that if you ever compute an MSD or diffusion coefficient, which
**would** be wrong on the raw trajectory.

---

## 2. Restrained domains: the bead is a residue centre of mass, the bond is a CA–CA length

Backbone bonds inside restrained domains sit ~29% longer than the bond
parameter says they should:

| | n bonds | simulated | input PDB (CA) |
|---|---|---|---|
| inside a restrained domain | 1048 | **0.4962 nm** | 0.3855 |
| outside any restraint | 752 | 0.3856 nm | 0.3874 |

The bond is `r0 = 0.38 nm, kb = 8033 kJ/mol/nm²`, so thermal spread should be
±0.017 nm. Outside restraints the bond sits exactly where it should. Inside, it
is stretched by 0.11 nm — about 54 kJ/mol of permanent strain per bond.

### Why

`use_com: true` (CALVADOS's shipped default, `data/default_component.yaml:17`)
makes `geometry_from_pdb` place each bead at its **residue centre of mass**, and
restraint target distances are generated in that same representation. Confirmed
exactly — restraint targets match COM distances with **rms error 0.0000 nm**,
versus 0.106 nm against CA:

| reference used to reproduce restraint targets | rms error |
|---|---|
| CA positions | 0.106 nm |
| CB (CA for Gly) | 0.091 nm |
| **residue centre of mass** | **0.0000 nm** |

Consecutive-residue spacing differs between the two representations:

| representation | consecutive spacing |
|---|---|
| CA | 0.386 nm |
| **residue COM** | **0.499 nm** |

So the restraint network encodes ~0.50 nm spacing while the bond parameter
encodes the CA–CA value of 0.38 nm. Inside a domain the restraints (≈5 per bead
at k=700) collectively win over the single bond (k=8033) and the bond settles at
0.496 — essentially the COM spacing. **220 of the 1001 `i,i+2` restraint targets
exceed 0.76 nm, which is geometrically impossible with 0.38 nm bonds**, so the
stretch is forced, not incidental.

This is upstream CALVADOS behaviour with its default settings, not something
introduced here.

### How much does it actually matter

Much less than the 29% bond figure suggests, because Rg is dominated by overall
shape rather than local spacing:

| domain | CA ref Rg | COM ref Rg | simulated Rg | sim / CA |
|---|---|---|---|---|
| RRM1 | 1.168 | 1.196 | 1.190 | 1.018 |
| RRM2 | 1.093 | 1.113 | 1.104 | 1.011 |
| RRM3 | 1.082 | 1.113 | 1.107 | 1.023 |
| MD1 | 1.505 | 1.520 | 1.513 | **1.005** |
| MD2 | 1.521 | 1.535 | 1.527 | **1.004** |
| MD3 | 1.469 | 1.487 | 1.477 | **1.006** |
| WWE | 1.209 | 1.243 | 1.235 | 1.022 |
| ART | 1.665 | 1.685 | 1.674 | **1.005** |

**Domain Rg is inflated only 0.4–2.3%**, and the simulated value tracks the COM
reference almost exactly — the restraints are doing their job faithfully; they
are simply faithful to a COM structure.

### The number to carry forward

Comparing a CALVADOS model against a **CA** crystal structure has a built-in
floor that has nothing to do with simulation quality:

| domain | per-residue \|CA − COM\| | Kabsch RMSD |
|---|---|---|
| MD1 | 0.116 nm | 0.129 nm |
| MD2 | 0.115 nm | 0.126 nm |
| MD3 | 0.114 nm | 0.126 nm |
| ART | 0.127 nm | 0.142 nm |
| WWE | 0.129 nm | 0.141 nm |
| whole chain | 0.122 nm | **0.135 nm** |

**≈1.3 Å of RMSD is representational.** This matters for the Phase 2.5 restraint
optimisation, which scored trim values by RMSD and Rg against crystal
references (3VFQ, 3Q6Z, 3GOY, 1X4R): roughly 1.3 Å of any such RMSD is the
CA-vs-COM difference, and ~0.5–2% of any Rg excess is the same effect. Neither
invalidates the ranking — both are near-constant offsets across trim values —
but they should not be read as restraint error.

### Fixed by converting the reference, not the model — 2026-09-25

`sim_analysis/ref_geometry.py` puts an experimental structure into the same
representation as a CALVADOS bead model. The trajectories are untouched; only
the yardstick changes.

```bash
python sim_analysis/ref_geometry.py --report    # the floor, per reference
```

| reference | residues | CA-vs-COM RMSD | Rg CA | Rg COM | ΔRg |
|---|---|---|---|---|---|
| MD1 (3Q6Z) | 185 | 0.128 nm | 1.509 | 1.521 | +0.8% |
| MD2 (3VFQ) | 187 | 0.123 nm | 1.533 | 1.548 | +0.9% |
| WWE (3GOY) | 69 | 0.137 nm | 1.514 | 1.518 | +0.3% |
| ART (3GOY) | 110 | 0.133 nm | 1.491 | 1.517 | +1.8% |
| WWE (1X4R) | 99 | 0.130 nm | 1.834 | 1.851 | +0.9% |

`archive_analysis/compare_restraints.py` now loads its references this way.

Two details that matter for matching CALVADOS exactly:

* **Heavy atoms only.** The AF2/AF3 structures the simulations were built from
  carry no hydrogens, so their bead is a heavy-atom COM. `1x4r` is an NMR entry
  and is ~50% hydrogen; including those would shift that one reference's COM
  relative to every other and re-introduce an offset for the WWE comparison
  alone.
* **One chain.** Found while doing this: **3GOY has four chains (A–D), all
  numbered 1532–1720.** Selecting on residue number alone — which the original
  `load_xtal_ca` did — returned all four copies stacked: 441 "residues" for a
  118-residue domain, Rg 3.2 nm instead of 1.5. So the WWE and ART
  RMSD-vs-crystal numbers from Phase 2.5 were computed against a four-copy
  reference and should not be trusted. MD1 (3Q6Z) and MD2 (3VFQ) are
  single-chain and were unaffected.

`compare_restraints.py` also had every data path resolving one level too deep
after it was moved into `archive_analysis/`, so it could not have run there at
all. Root resolution now walks up to the project root.

Setting `use_com: false` would instead place beads at CA and build restraints
from CA distances, matching the 0.38 nm bond — but that changes the model and
would require re-running everything, so it was not done.

### A second-order oddity worth knowing

Because the bond wins in unrestrained regions and the restraints win inside
domains, the same chain carries **two different effective bead spacings** —
0.386 nm in linkers and IDRs, 0.496 nm inside folded domains. Domains are
compact so this barely affects their Rg, but it does mean the model's contour
length is not uniform along the chain.
