# PARP14 BindCraft — status, findings, and how to restart

**Last updated:** 2026-09-21, on `lyra.fraserlab.com` (4× NVIDIA L40S, 46 GB each, **no SLURM**).
Supersedes the 2026-08-14 handoff, which described the pre-install state.

## Install: DONE

`BindCraft` conda env at **`/home/sbali/miniconda3/envs/BindCraft`**
(*not* `~/.conda/envs` as the old handoff assumed — `bindcraft_job.slurm` was corrected).
AF2 weights present in `/home/sbali/BindCraft/params/` (16 files).
Verified: jax 0.6.0 sees all 4 GPUs; `colabdesign` and `pyrosetta` import cleanly.

```bash
source /home/sbali/miniconda3/bin/activate BindCraft
export LD_LIBRARY_PATH=/home/sbali/miniconda3/envs/BindCraft/lib:${LD_LIBRARY_PATH}
cd /home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md
CUDA_VISIBLE_DEVICES=<n> python -u /home/sbali/BindCraft/bindcraft.py \
  --settings settings/<target>.json \
  --filters  /home/sbali/BindCraft/settings_filters/default_filters.json \
  --advanced /home/sbali/BindCraft/settings_advanced/default_4stage_multimer.json
```

No SLURM on this host — `submit_bindcraft.sh` / `sbatch` do not work here. Launch directly,
pin with `CUDA_VISIBLE_DEVICES`, and check `nvidia-smi` first (shared machine).

### REQUIRED PATCH to the BindCraft clone (not upstreamable from here)

`/home/sbali/BindCraft` has remote `martinpacesa/BindCraft` (upstream, no write access), so this
fix lives only in the local working tree. **Re-apply it after any fresh clone or `git checkout`**,
or every run dies at interface scoring:

`functions/pyrosetta_utils.py` — PyRosetta 2026.x changed `InterfaceAnalyzerMover.set_interface()`
to take a `DockingPartners` object instead of a string:

```python
from pyrosetta.rosetta.core.pose import DockingPartners       # add import
iam.set_interface(DockingPartners.docking_partners_from_string("A_B"))   # was: iam.set_interface("A_B")
```

## Baseline campaigns (no guess flag -- all ZERO accepted)

| Target | Concept | Ran | Attempts | Relaxed | Scored | Accepted |
|---|---|---|---|---|---|---|
| `md1_block_af3` | block MD1 pocket, AF3 receptor | Aug 14–28 (14 d) | 617 | 270 | 10 | **0** |
| `clamp_md1md2_state3` | clamp MD1–MD2, TICA state 3 | Aug 28–Sep 9 (12 d) | 605 | 220 | 0 | **0** |
| `md1_block_sim` | block MD1 pocket, CALVADOS conformer | Sep 9–16 (7 d) | 474 | 91 | 0 | **0** |

## Key findings (evidence in `analysis/`)

1. **Hotspots are fine.** All 15 MD1 hotspots are solvent-exposed (RSA 0.16–0.79) and form one
   compact patch (24.7 Å span). Not the bottleneck.
2. **Clamp geometry works.** 206 of 220 relaxed clamp trajectories (94%) contact *both* MD1 and
   MD2 within 4 Å. The two-domain ask is being satisfied — loosening hotspots would not help.
3. **Receptor conformer matters a lot.** Same pocket/hotspots, AF3 vs CALVADOS-backmapped:
   clash rate 10% → 35%, usable relaxed yield 58% → 30%. The back-mapped surface is materially
   harder to design against (as the target-prep notes predicted).
4. **THE WALL: `i_pAE` at AF2 re-prediction.** Rejects ~100% of MPNN sequences in all three runs
   (5,370/5,400 · 1,820/1,820 · 4,400/4,400).
5. **The failures are catastrophic, not marginal** (`analysis/replay_mpnn.py`, 380 paired
   predictions): median i_pAE **0.878** vs threshold 0.35; **95.5% of failures exceed 0.70**.
   Meanwhile pLDDT passes 69% — binders *fold* fine but do not *dock*. Note `mpnn_fix_interface`
   is True, so MPNN preserves interface residues verbatim and AF2 still cannot recover the
   binding mode. The hallucinated pose looks adversarial to AF2 rather than real.
   → **Relaxing the i_pAE threshold is NOT defensible.** These are not near-misses.

## THE FIX: `predict_initial_guess`

A paired replay over already-generated trajectories (`analysis/replay_mpnn.py`, 396 pairs -- identical
backbone, sequence and model, only the validation protocol differs) showed the gate is passable:

| | guess OFF | guess ON |
|---|---|---|
| i_pAE median | 0.877 | **0.508** |
| pass all three AF2 gates | 1/396 | **119/396** |

118 rescued, 0 lost; 96.5% of pairs improved. Not a rubber stamp -- 70% still fail and the median stays
above threshold, so it still discriminates; it just stops demanding de-novo docking.

**Full-pipeline confirmation (`designs/sweep_af3_guess/`), same target, same filters, one flag:**

| | baseline | with guess |
|---|---|---|
| Runtime | 14 days | **17 hours** |
| Trajectory attempts | 617 | 21 |
| Reached Rosetta scoring | 10 | 51 |
| **Accepted** | **0** | **9** |

Accepted designs came from **5 distinct trajectories** (not one lucky backbone) and clear the Rosetta gates
that killed every baseline design: 1.0-4.0 buried unsat H-bonds (limit 4; baseline sat at 5-9), ~2x the
interface H-bonds, dG -40 to -70, and Binder_RMSD 0.83-1.63 A -- the last computed with no initial guess,
so the flag cannot inflate it. The wall moved rather than vanished: i_pAE now rejects 58% of sequences
instead of 99.4%, and the leading rejection is Rosetta shape complementarity.

Run it via the stock preset `settings_advanced/default_4stage_multimer_hardtarget.json`
(identical to default except `predict_initial_guess: true`).

## Queue results so far

`sweep/run_queue.sh` (detached, survives logout) runs all 9 remaining targets with the flag, serial on one
GPU, each stopping at 8 accepted or `max_trajectories=50`. Progress: `logs/queue_master.log` and
`logs/queue_progress.log` (6-hourly, written by `sweep/progress_logger.sh`).

| target | relaxed | scored | accepted | status |
|---|---|---|---|---|
| `clamp_md1md2_state3` | 50 | 197 | **5** | complete (hit 50-traj cap, 3.2 d) |
| `clamp_md1md2_state1` | 9 | 0 | 0 | running, externally capped at 25 traj |
| 7 more | | | | queued |

**state3 -- the gate moved to interface chemistry.** i_pAE is beaten (197 designs reached Rosetta), but of
those, 68% fail buried unsatisfied H-bonds (median 5.0 vs limit 4) and 68% fail shape complementarity
(median 0.58 vs 0.60). dG and interface H-bond count fail 0% (medians -58.8 and 10.5) -- energetics are
fine, polar burial is the limiter. This is the same metric that killed the original md1_block_af3 campaign,
reappearing one stage later. A clamp has more polar surface to satisfy than a pocket plug, so it may be
intrinsic to the concept. Accepted clamp designs verified to bridge BOTH domains (7 MD1 + 18 MD2 contacts).

**state1 -- the flag cannot rescue a bad starting interface.** Its trajectories fold as well as state3's
(pLDDT 0.80) but start with median i_pAE **0.540** vs state3's 0.345, and 9 relaxed trajectories produced
0 scored designs. predict_initial_guess helps AF2 re-find a pose the sequence already supports; it cannot
manufacture an interface that was never good. Consistent with the TICA clustering -- state3 was the dominant
~80% basin, state1 a ~6% minor one -- but 9 trajectories is far too few to conclude designability tracks
state population. Capped at 25 trajectories (`sweep/cap_state1.sh`) to bound the cost of finding out; the
running process had already loaded max_trajectories=50, so the cap is applied externally by SIGINT, after
which run_queue.sh advances on its own.

## Open questions

- Does the fix generalise past MD1? The clamp is the cleanest test: its geometry was already correct
  (94% contact both domains) and it died at exactly this gate.
- Arm C (`mpnn_weights=original`) was dropped: the replay showed folding was never the bottleneck
  (pLDDT passed 70% even with the guess off), so scaffold composition is unlikely to be the lever.
- `md2_block_*` **OOMs on a 46 GB L40S** (586-res target → 683-res complex → 31.4 GiB single alloc).
  Either crop the target (all 23 md2 hotspots lie in 216–404, so MD3 is dead weight) or run with
  `TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0`.

## Layout

`settings/` targets · `sweep/` A/B arm configs + chain runner · `analysis/` diagnostic code and
results (replay CSV, geometry scans, report builder) · `designs/`, `logs/`, `figures/` are
gitignored (3.4 GB of run output).
