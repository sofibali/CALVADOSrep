# PARP14 Macrodomain BindCraft Campaign

De novo binder design (BindCraft) against the PARP14 macrodomain module (MD1–MD3),
informed by the CALVADOS CG-MD simulations in `examples/PARP14_MDP/`.

## Goal

Two binder concepts:
- **Block a single macrodomain** — bind the ADP-ribose pocket of MD1, MD2, or MD3
  to compete with substrate/ADP-ribose binding.
- **Clamp two macrodomains into a chosen arrangement** — *enforce a new interface*:
  a binder that bridges the inter-domain cleft of a selected conformational state,
  holding the two domains in a juxtaposition they do not stably adopt on their own.

## What the simulations told us (and how it shaped the design)

The `md` set has **no stable native inter-domain interface**. Folded-core
(restraint-core) results from `analyze_all.py`:

| Pair | COM–COM (nm) | inter-domain contact freq |
|------|--------------|---------------------------|
| MD1–MD2 | 3.16 ± 0.39 | ~0.006 |
| MD2–MD3 | 3.65 ± 0.68 | ~0.003 |
| MD1–MD3 | 3.81 ± 0.79 | ~0.004 |

(intra-domain contact freq ≈ 0.12–0.13). The three macrodomains behave as dynamic,
largely non-associating units. **This is why the design goal is to *enforce* a new
interface, not stabilise an existing one.**

> ⚠️ Correction: an earlier draft reported "MD1–MD2 permanently in contact (0.38 nm)".
> That was wrong — it came from the `fl_optimized` set and measured the `MD1L1`
> block *including the flexible L1 linker* (bonded to MD2, so always touching). The
> folded MD1↔MD2 domain contact is dynamic (table above).

The clamp poses come from your **TICA/MSM clustering** (`md_ca25_tica`, 2026-05-20),
states 3 (80%), 5 (7%), 1 (6%). Each state is a *distinct* juxtaposition; at folded-
core resolution each makes genuine contact (core min CA–CA ≈ 0.5–0.6 nm) with its own
contact patch. We enforce three states (state 1's MD2–MD3 is open ⇒ skipped).

Two structural caveats handled in target prep:
- **CG-backmapped multi-domain poses clash at atomic detail** (one-bead-per-residue
  artifact). We rigid-body **declash** them (below).
- The construct **deletes FL 1194–1206** (the MD2–MD3 linker) — MD2 and MD3 are
  joined directly. MD1L1+MD2 are contiguous (FL 790–1193, L1 linker retained).

## Target construction

Construct = MD1L1(FL 790–1004) + MD2(FL 1005–1193) + MD3(FL 1207–1388), local
numbering **1–586** (matches `backmapped_states/*.pdb` and the AF3 model).
FL→local: MD1/MD2 = FL−789; MD3 = FL−802.

Two structural sources per design goal (`prepare_targets.py`):
- **`*_af3`** — coordinates from the clash-free AF3 `md1l1_md2_md3` model.
- **`*_sim`** — the CALVADOS conformer for that goal: AF3 all-atom domains are
  Kabsch-fit (on each domain's structured core; per-domain CA-RMSD 0.1–0.3 Å) onto
  the backmapped CA frame, then **rigid-body declashed** — interpenetrating domains
  are translated apart along their COM axis to van-der-Waals contact (min-dist
  ≈ 2.6 Å) while preserving orientation.

Targets are truncated to the design-relevant domains to keep AF2 tractable. Residue
numbers are preserved (not renumbered) so hotspot numbers are consistent.

**11 targets** (`targets_summary.csv`):

| Target | goal | source | residues kept | # hotspots |
|--------|------|--------|---------------|-----------|
| `md1_block_af3` / `_sim` | block MD1 pocket | AF3 / sim | 1–404 (MD1+MD2) | 15 |
| `md2_block_af3` / `_sim` | block MD2 pocket | AF3 / sim | 1–586 (all 3)   | 23 |
| `md3_block_af3` / `_sim` | block MD3 pocket | AF3 / sim | 216–586 (MD2+MD3) | 23 |
| `clamp_md1md2_state3/5/1` | enforce MD1–MD2 @ TICA state | sim-state | 1–404 | 14 |
| `clamp_md2md3_state3/5` | enforce MD2–MD3 @ TICA state | sim-state | 216–586 | 12–14 |

**Block hotspots** = catalytic residues (`parp14/input/active_sites.yaml`) + pocket
residues that are surface-exposed (RSA ≥ 0.15, `data/sasa_face_per_residue.csv`) on
the active-site face. **Clamp hotspots** = the contact patch of that TICA state — the
closest core residues on *both* domains across the cleft (the surfaces the binder
must bridge to enforce the juxtaposition). Clamps are sim-pose-only.

> Caveat: a clamp target presents the two domains already in the target pose; BindCraft
> designs a binder to that fixed pose. Whether the binder actually *holds* the domains
> there (vs just binding the cleft) is an experimental question — these are candidates.

Full per-target list: `targets_summary.csv`. Settings JSONs: `settings/`.

## Running

> **Status (2026-09-17): BindCraft is installed and the campaign has run.**
> The repo is at `~/BindCraft` (not `~/Projects/BindCraft`, which does not
> exist), and `bindcraft_job.slurm` points at
> `BC_ENV=/home/sbali/miniconda3/envs/BindCraft`.
>
> **Results so far — the default protocol produces nothing.** `md1_block_af3`,
> `md1_block_sim`, `clamp_md1md2_state3`, `sweep_af3_diag` and
> `sweep_af3_mpnnorig` all have **0 accepted designs**, failing at the
> interface-confidence filters (`i_pAE`, `i_pLDDT`) rather than at clashes.
>
> **What works is `predict_initial_guess: True`:**
> `sweep_af3_guess` → **18 accepted**, `guess_clamp_md1md2_state3` → **2**.
> It is the only setting present in the winning arm and absent from every
> failing one.
>
> ⚠ That flag "introduces bias by providing binder atom positions as a starting
> point" (BindCraft README). The improved `i_pAE`/`i_pTM` are therefore **not
> independent validation** — AF2 was given the answer. Re-predict without the
> initial guess before treating any of these as real hits.
>
> Full write-up, including the config diff and the accepted-design metrics:
> [`../docs/questions/Q8_binder_design.md`](../docs/questions/Q8_binder_design.md).

Submit (GPU partition, one job per target):
```bash
cd ~/CALVADOS/examples/PARP14_MDP/bindcraft_md
./submit_bindcraft.sh blocks        # the 6 well-defined pocket blockers (recommended pilot)
./submit_bindcraft.sh all           # + the exploratory MD2–MD3 lock
./submit_bindcraft.sh md1_block_af3 # a single target
```
Designs are written to `designs/<target>/`; final accepted binders in
`final_design_stats.csv` there. Pilot = 8 final designs/target (edit
`number_of_final_designs` in the settings JSON to scale up).

## Files
- `prepare_targets.py` — builds targets + hotspots + settings (run in CALVADOS env)
- `targets/` — 7 BindCraft-ready target PDBs
- `settings/` — 7 BindCraft settings JSONs
- `targets_summary.csv` — target/hotspot table
- `bindcraft_job.slurm` — single GPU BindCraft job
- `submit_bindcraft.sh` — batch submitter
- `designs/`, `logs/` — outputs
