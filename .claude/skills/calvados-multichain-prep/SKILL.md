---
name: calvados-multichain-prep
description: "Build a multi-chain CALVADOS system (slab, coexistence, or any nmol > 1 panel) without the silent failures — per-copy custom restraints, explicit RNA force-field parameters, concentration matching, checkpoint cadence for preemptible GPUs, and paired seeds. Use when writing or reviewing a prepare_*.py that sets nmol > 1, adding an RNA or ligand arm, or designing a ±condition simulation panel."
---

# Preparing a multi-chain CALVADOS system

Every item here is something that **builds, runs, and produces plausible
numbers** while being wrong. None of them show up in a smoke test — only in the
physics. Three of the six came from real bugs in the PARP14 slab run-2 panel.

---

## 1. Custom restraints must be emitted per chain copy

The file format is:

```
<name> <copy> <bead> | <name> <copy> <bead> | <r> <k>
```

`Sim.map_custom_restraints` resolves **each line to exactly one bead pair** via
`start_bead + (copy-1)*nbeads + (bead-1)`. It does **not** broadcast a line
across copies.

Every monomer run has `nmol=1`, so `copy 1` is always correct there — which is
exactly why this is easy to miss. In a 47-chain slab, emitting only `copy 1`
bridges one chain and leaves 46 unrestrained.

```python
return [f'{sysname} {c} {li} | {sysname} {c} {lj} | {d} {k}\n'
        for c in range(1, nmol + 1) for li, lj, d, k in pairs]
```

**Verify in the build log**: `Number of custom restraints:` must be
`npairs × nmol`.

## 2. RNA force-field parameters must be passed explicitly

`calvados/data/default_component.yaml` ships much stiffer RNA than the CALVADOS
RNA paper values used in the PARP14 campaign:

| | campaign (`prepare_slab.py` `RNA_OPTS`) | package default |
|---|---|---|
| `rna_kb1` / `rna_kb2` | 1400 / 2200 | 8033 / 8033 |
| `rna_ka` | 4.20 | 7.24 |
| `rna_nb_cutoff` | 2.0 | 0.6 |

Omitting `RNA_OPTS` falls back to the defaults **without warning** and
substitutes a different polymer, making arms incomparable across campaigns.
Pass them into `Components(...)` for any run with an RNA component, and diff the
emitted `components.yaml` against a known-good run.

Related: CALVADOS RNA is 2 beads per nucleotide (`rP` + `rN`) with **no base
identity**, so polyU is the only meaningful sequence, and `n_atoms` exceeds
`n_residues` by the nucleotide count — counting RNA "residues" as beads
undercounts its contact surface by half.

## 3. Match concentration, not box volume

Across a construct panel, hold **residue molarity** constant, not box size or
chain count. Concentration is what decides condensation, so a panel at varying
concentration confounds construct identity with density. PARP14 run 1 spanned
99–863 mM and could not separate the two; run 2 fixed every cell at 100 mM.

Compute chains from the target concentration and the box, then report the
achieved value per cell (`--report`) before generating anything.

## 4. `checkpoint_interval = wfreq` for anything preemptible

Stock `sim.py` hard-codes `nbatches = 10`, so a run checkpoints once per tenth
of its step count. On a shared or preemptible GPU that is the dominant loss. The
local patch derives `nbatches` from an optional `checkpoint_interval` key
(absent/0 keeps original behaviour); set it to `wfreq` so a kill costs ≤1 frame.

**This applies to the production loop only.** Equilibration must keep writing no
checkpoint — an OpenMM checkpoint is bound to its System and the equilibration
System still carries the slab-centering force, so a checkpoint taken there
cannot load into the production System. The consequence is that equilibration is
**unpreemptible**; plan an uninterrupted slot for it
(see `gpu-sharing-policy-lyra` memory).

## 5. Pair seeds across arms

In a ±condition panel, use the **same seed set in both arms** so the comparison
is matched rather than between independent draws. With n=3 and paired seeds, a
difference under one pooled SD is a real null; with unpaired n=1 it is nothing.

Run 1's "RNA retards dispersal 1.67×" came from one unpaired cell at mismatched
concentration (99 vs 108 mM) and did not survive replication.

## 6. Sanity checks before launching the panel

- `--report` first: constructs, spans, chain counts, box sizes, achieved
  concentration, restraint counts. Read it before generating.
- Build **and run** one config of each type on CUDA (a 200-step smoke test), not
  just `Sim.build_system()`. Confirm in the log: restraint counts, chain counts,
  presence/absence of the ligand component, `Running on CUDA`.
- `components.yaml` carries **absolute paths** — the tree is not portable. Same
  filesystem, different GPU is fine; a different machine means re-running the
  prepare script there.
- Use the env whose OpenMM actually has CUDA (`envs/calvados`, not
  `envs/CALVADOS` — see the `calvados-conda-envs` memory).

## 7. Slicing a construct out of a full-length model

Slice at the **text level** and renumber residues 1..N. `AtomGroup.copy()` in
MDAnalysis does not produce an independent Universe, so renumbering `res.resid`
invalidates a `resid`-based selection and silently yields an empty group
(`IndexError: Cannot write an AtomGroup with 0 atoms`).

Verify a slice reproduces the parent: Kabsch RMSD of the slice against the
corresponding region of the source model should be 0.000 nm.

---

After preparing, before interpreting anything, run the `slab-coexistence-qc`
skill.
