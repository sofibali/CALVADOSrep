# Q8 — Can we design binders that block a macrodomain pocket, or clamp two macrodomains together?

**Status:** OPEN — campaign ongoing. First pass returned 0 accepted designs;
a parameter sweep is in progress to determine whether that is target difficulty
or filter strictness.

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

## First-pass result (to be confirmed or overturned by the sweep)

Three targets have substantial output — `md1_block_af3` (5128 files),
`clamp_md1md2_state3` (4690), `md1_block_sim` (3348). Current state:

- every `designs/*/final_design_stats.csv` is **header-only**
- **zero** PDBs under any `designs/*/Accepted/` (the entries that look like
  results are empty BindCraft scaffold dirs: `Animation`, `Pickle`, `Plots`,
  `Ranked`)

Where they died (`failure_csv.csv`) — rejection concentrates at the
**interface-confidence** filters, not at clash or sequence stages:

| target | Trajectory_Clashes | pAE | i_pAE | i_pLDDT |
|---|---|---|---|---|
| `md1_block_af3` | 46 | 1283 | 5243 | 5370 |
| `clamp_md1md2_state3` | 110 | 1642 | 4319 | 4400 |

`i_pAE` and `i_pLDDT` dominate — AlphaFold is not confident about the
*interface* of the designed complexes. That is the signature of a hard target
rather than a misconfigured run, which is exactly what the sweep is set up to
discriminate.

Sweep runs so far are incomplete: `logs/sweep_af3_diag.log` ends in a
`KeyboardInterrupt` during `binder_hallucination`; `logs/sweep_chain.log` holds
one line, "waiting for replay PID 4003519 to finish".

## What would close this question

- Does any arm of the sweep produce accepted designs under `diagnostic_filters`?
  If yes → the production filters are too strict for this target class. If no →
  the targets themselves are not designable by this protocol.
- Does the `_af3` vs `_sim` split matter? That tells you whether starting
  conformation is the bottleneck.
- Do clamps behave differently from blocks? Clamps ask for a larger, composite
  interface and may fail for a different reason.

## Note on stale docs

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
