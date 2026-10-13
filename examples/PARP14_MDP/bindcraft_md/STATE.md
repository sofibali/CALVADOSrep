# Campaign state — 2026-10-13

Nothing is running. No BindCraft process is alive and the GPU is free.
This file is the single place to look when picking the project back up.

## Where it ended

BindCraft 1, with `predict_initial_guess` enabled, produced **14 accepted binders**:

| run | target | apolar % of patch | relaxed | scored | accepted |
|---|---|---|---|---|---|
| `sweep_af3_guess` | MD1 block | 43.8 | 12 | 51 | **9** |
| `guess_clamp_md1md2_state3` | MD1–MD2 clamp, TICA state3 | 42.9 | 50 | 197 | **5** |
| `guess_md2_block_af3` | MD2 block | 38.1 | 50 | 16 | 0 |
| `guess_clamp_md1md2_state1` | MD1–MD2 clamp, TICA state1 | 28.6 | 28 | 49 | 0 |
| `guess_md3_block_af3` | MD3 block | 28.6 | 50 | 0 | 0 |

Accepted structures and the four stats CSVs per run are committed under `results/`.
The full write-up is the artifact (see `analysis/` below).

Two results that drive what to do next:

1. **`predict_initial_guess` was the whole fix.** The MD1 baseline ran 617
   trajectories over 14 days for 0 accepted; the same target with the flag gave 9
   accepted in 17 hours. Trajectory yield was *identical* (44%) in both arms — the
   flag changed nothing about hallucination, only what survived AF2 re-prediction.
2. **Apolar fraction of the hotspot patch predicts designability as a threshold**,
   not a gradient: everything >=42.9% produced binders, everything <=38.1% produced
   none. Pre-screen a site with `analysis/skill/scripts/pocket_chem.py` before
   spending GPU-days on it.

## Deliberately left undone

- **The MD1-focused BC2 arm is prepped but not launched** — paused on request.
  Configs in `bc2/campaigns/`. Goal: improve i_pAE on the 9 MD1 binders, which sit
  at 0.21–0.34 with 3 of 9 within 0.02 of the 0.35 cutoff.
  For that arm, **`ShapeComplementarity` is the gate to attack** (63% failure in
  `sweep_af3_guess`), not unsatisfied H-bonds (39%).
- **The BC2 selective configs still contain a `forced_targeting` key that breaks
  them.** `epitope_targeting.py:55` raises because the flag applies to every target
  including the detargets, which intentionally carry no hotspots. Left in place and
  documented rather than silently changed — decide the intent before editing.
- Remaining BC1 queue entries were dropped; they are left ready in `sweep/queue/`.

## Rebuilding the report

The pipeline is repo-relative on purpose: an earlier version of these figures was
lost when a `/tmp` scratchpad was cleared, so everything regenerates from committed
inputs.

```bash
python analysis/make_figdata.py     # recompute derived numbers from designs/, targets/, bc2/
python analysis/make_figures.py     # PyMOL renders -> analysis/assets/
python analysis/crop_figures.py     # trim/downscale to the embedded names
python analysis/build_report.py     # -> analysis/parp14_geometry_report.html
```

Then publish that HTML to https://claude.ai/artifact/DcdGTzfwaiib2iaTSDmGC5.
**Omit the `capabilities` field when republishing** — the page declares `db` and
`user` so notes can be added under any heading in the browser, and passing an empty
value would drop them.

## Reading campaign output

Use `analysis/skill/` (installed as the `bindcraft-triage` skill). Several BindCraft
counters do not mean what their names suggest — `failure_csv.csv` counts
**per-sequence** events, `LowConfidence/` conflates two different failures, the
4-stage early aborts never write a file, and `InterfaceResidues` is keyed by binder
residue rather than target. `SKILL.md` has the full list.

```bash
python analysis/skill/scripts/triage.py designs/sweep_af3_guess
python analysis/skill/scripts/pocket_chem.py settings/md1_block_af3.json
```

## Environment note

`/home/sbali/BindCraft` is a clone of **upstream** `martinpacesa/BindCraft`, not a
fork — nothing there is pushable, and `params/` is 5.3 GB of AF2 weights that must
stay untracked. The local `functions/pyrosetta_utils.py` patch is superseded by a
more robust upstream fix, so pulling is safe.
