# Q8 — Can we design binders that block a macrodomain pocket, or clamp two macrodomains together?

**Status:** ONGOING — the sweep has produced accepted designs. The
`predict_initial_guess` protocol is what unlocked it; read the caveat on what
that metric improvement does and does not mean.

---

## The question

Two design concepts (`bindcraft_md/README.md:6-13`):

1. **Block** — bind the ADP-ribose pocket of MD1 / MD2 / MD3 to compete with
   substrate.
2. **Clamp** — bind across two macrodomains to *enforce a new interface*,
   "holding the two domains in a juxtaposition they do not stably adopt on their
   own."

## Why "enforce", not "stabilise"

Directly from [[Q5]]: the macrodomains are largely non-associating (inter-domain
contact frequency ~0.003–0.006 vs ~0.12–0.13 intra-domain). There is no stable
interface to stabilise. `bindcraft_md/README.md:26-28` states it: "**This is why
the design goal is to *enforce* a new interface, not stabilise an existing one.**"

[[Q2]] justifies how the clamp poses were chosen: because the module is diffusive
rather than metastable, the TICA clusters are "structural representatives sampled
along a continuum, not kinetically-distinct metastable states" — so picking
**diverse poses** is right, since clamping needs geometric diversity, not kinetic
metastability.

Hotspots come from [[Q7]]: catalytic residues + pocket residues that are
surface-exposed (RSA ≥ 0.15) on the active-site face.

---

## What is to be tested

### Production targets — `bindcraft_md/settings/` (11)

| settings file | concept | binder length | hotspots |
|---|---|---|---|
| `md1_block_af3.json` | block MD1 pocket, AF3 pose | 70–130 | 15 |
| `md1_block_sim.json` | block MD1 pocket, MD pose | 70–130 | 15 |
| `md2_block_af3.json` | block MD2 pocket, AF3 pose | 70–130 | 23 |
| `md2_block_sim.json` | block MD2 pocket, MD pose | 70–130 | 23 |
| `md3_block_af3.json` | block MD3 pocket, AF3 pose | 70–130 | 23 |
| `md3_block_sim.json` | block MD3 pocket, MD pose | 70–130 | 23 |
| `clamp_md1md2_state1.json` | clamp MD1–MD2, TICA state 1 | 90–150 | 14 |
| `clamp_md1md2_state3.json` | clamp MD1–MD2, TICA state 3 (80 %) | 90–150 | 14 |
| `clamp_md1md2_state5.json` | clamp MD1–MD2, TICA state 5 (7 %) | 90–150 | 14 |
| `clamp_md2md3_state3.json` | clamp MD2–MD3, TICA state 3 | 90–150 | 12 |
| `clamp_md2md3_state5.json` | clamp MD2–MD3, TICA state 5 | 90–150 | 14 |

All target chain A, 8 final designs each. The `_af3` / `_sim` split is a
controlled comparison: **same pocket, same hotspots, different starting
conformation** — AF3 prediction vs an MD-sampled pose.

### Diagnostic sweep — `bindcraft_md/sweep/` (3 arms)

All three use the `md1_block_af3` hotspot set, so the *target is held fixed* and
only the protocol varies:

| arm | settings | advanced | varies |
|---|---|---|---|
| `af3_diag` | `sweep/settings/af3_diag.json` | `diag.json` | uncapped trajectories, diagnostic filters |
| `af3_guess` | `af3_guess.json` | `guess_cap100.json` | `max_trajectories: 50` |
| `af3_mpnnorig` | `af3_mpnnorig.json` | `mpnnorig_cap100.json` | original MPNN settings, `max_trajectories: 50` |

Shared advanced settings across all arms: `4stage` algorithm,
`use_multimer_design: true`, `sample_models: true`,
`num_recycles_design: 1` / `validation: 3`, `optimise_beta: true`,
soft/temporary/hard iterations 75 / 45 / 5.

`sweep/filters/diagnostic_filters.json` holds 218 filter entries — the relaxed
set used to see what *would* pass. Key interface thresholds there:
`Average_pLDDT ≥ 0.8`, `Average_i_pTM ≥ 0.5`, `Average_i_pAE ≤ 0.35`,
`Average_ShapeComplementarity ≥ 0.6`, `Average_dG ≤ 0`.

---

## Result — see `bindcraft_md/GPU_SERVER_SETUP.md` (primary source)

That file is the authoritative write-up and is more rigorous than the summary
below: it contains a **paired replay** (`analysis/replay_mpnn.py`, 396 pairs —
identical backbone, sequence and model, only the validation protocol differs)
that isolates the flag properly. Read it first. Summary:

### The baseline campaigns all returned zero

| target | ran | attempts | scored | accepted |
|---|---|---|---|---|
| `md1_block_af3` | 14 days | 617 | 10 | **0** |
| `clamp_md1md2_state3` | 12 days | 605 | 0 | **0** |
| `md1_block_sim` | 7 days | 474 | 0 | **0** |

### The wall was `i_pAE` at AF2 re-prediction, and it was not marginal

Median i_pAE **0.878** against a 0.35 threshold; **95.5% of failures exceed
0.70**. Meanwhile pLDDT passes 69% — the binders *fold* fine, they do not
*dock*. With `mpnn_fix_interface: True`, MPNN preserves interface residues
verbatim and AF2 still cannot recover the binding mode.

→ **Relaxing the i_pAE threshold would not have been defensible.** These are not
near-misses. That rules out "the filters are too strict", which was one of the
two hypotheses the sweep was built to test.

### `predict_initial_guess` is the fix

Paired replay, same trajectories, only the validation protocol differing:

| | guess OFF | guess ON |
|---|---|---|
| i_pAE median | 0.877 | **0.508** |
| pass all three AF2 gates | 1/396 | **119/396** |

118 rescued, 0 lost, 96.5% of pairs improved. Full-pipeline confirmation on the
same target with the same filters: runtime **14 days → 17 hours**, attempts
617 → 21, accepted **0 → 9**, from 5 distinct trajectories rather than one lucky
backbone.

### Why this is not just a rubber stamp

The obvious worry is that seeding AF2 with binder atom positions
("introduce bias", per BindCraft's README) simply makes the metric easier. Three
things argue against dismissing the result on that basis, all from the primary
source:

- **It still discriminates.** 70% of sequences still fail and the median i_pAE
  stays above threshold. The wall moved (99.4% → 58% rejection) rather than
  vanished.
- **Independent metrics improved too.** The accepted designs clear Rosetta gates
  that killed every baseline design: 1.0–4.0 buried unsatisfied H-bonds (limit 4;
  baseline sat at 5–9), ~2× the interface H-bonds, dG −40 to −70.
- **`Binder_RMSD` 0.83–1.63 Å was computed with no initial guess**, so that one
  cannot be inflated by the flag.

The leading rejection is now Rosetta shape complementarity, not AF2 confidence.

*(An earlier version of this file warned that i_pAE/i_pTM from these runs were
unvalidated. That was too strong — the guess-off Binder_RMSD and the Rosetta
gates are exactly the independent check it asked for.)*

## What would close this question

Per `GPU_SERVER_SETUP.md`, `sweep/run_queue.sh` is **in flight** — all 9
remaining targets with the flag, serial on one GPU, stopping at 8 accepted or
`max_trajectories=50`. Progress: `logs/queue_master.log`.

Open after that:

1. **Does the fix generalise past MD1?** The clamp is the cleanest test — its
   geometry was already right (94% of relaxed trajectories contact both domains)
   and it died at exactly this gate. First evidence is positive:
   `guess_clamp_md1md2_state3` has 2 accepted.
2. **Receptor conformer matters and is unresolved.** Same pocket and hotspots,
   AF3 vs CALVADOS-backmapped: clash rate 10% → 35%, usable relaxed yield
   58% → 30%. The backmapped surface is materially harder to design against.
3. **`md2_block_*` OOMs on a 46 GB L40S** (586-res target → 683-res complex →
   31.4 GiB single alloc). Either crop the target — all 23 MD2 hotspots lie in
   216–404, so MD3 is dead weight — or run with `TF_FORCE_UNIFIED_MEMORY=1
   XLA_PYTHON_CLIENT_MEM_FRACTION=4.0`.
4. **Shape complementarity is the new bottleneck**, having replaced AF2
   confidence.

## Note on stale docs

These are now badly out of date and should be corrected in place:
`bindcraft_md/README.md:87-88` still says "BindCraft is **NOT** yet installed"
and `GPU_SERVER_SETUP.md` (2026-08-14) says "`designs/` and `logs/` are empty —
nothing has been run yet". Both predate the runs above. The `sweep/` directory is
referenced in no README.

## Honest caveats already recorded

`bindcraft_md/README.md:41-58` — clamp targets are CG-backmapped poses that clash
and are rigid-body declashed. Lines 79-81: "Whether the binder actually *holds*
the domains there (vs just binding the cleft) is an experimental question — these
are candidates."

## Outputs

- `bindcraft_md/settings/`, `bindcraft_md/sweep/`
- `bindcraft_md/designs/*/failure_csv.csv`, `*/Trajectory/`, `*/Rejected/`
- `bindcraft_md/analyze_campaign.py` (aggregates; tolerates 0 accepted)
