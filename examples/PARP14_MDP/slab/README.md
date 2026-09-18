# PARP14 slab / assembly runs

Everything here is prepared, validated, and **costed on real hardware**. No
production run has been started yet.

Status as of 2026-09-17, on `lyra.fraserlab.com`:

- all 20 configs build and reach `STARTING SIMULATION` on CUDA with the correct
  bead counts
- throughput measured on an idle L40S — **~11 h per run, ~9.5 GPU-days for the
  whole panel**, 8× cheaper than the assumption this file used to carry
- the launcher, resume arithmetic and completion guard are tested end-to-end
- the env the old scripts pointed at **cannot run these at all** (no CUDA
  platform) — fixed, see *Things that will bite*

```
/home/sbali/CALVADOS/examples/PARP14_MDP/slab
```

> ### ⚠ lyra has no SLURM — use `run_slab_queue.sh`, not `submit_slab.slurm`
>
> `bindcraft_md/GPU_SERVER_SETUP.md` records that `lyra.fraserlab.com`
> (4× NVIDIA L40S, 46 GB) has **no SLURM**. `submit_slab.slurm` was written
> before that was known and only applies if you have a SLURM cluster.
>
> On lyra, launch directly, pinned and detached — the same pattern as
> `bindcraft_md/sweep/run_queue.sh`:
>
> ```bash
> cd /home/sbali/CALVADOS/examples/PARP14_MDP/slab
>
> # GPU 0 -- the campaign's own GPU: long runs, start to finish
> GPU=0 nohup ./run_slab_queue.sh homotypic > logs/queue_homotypic.log 2>&1 &
>
> # GPU 1 and 2 -- opportunistic: run only while empty, yield to anyone else
> GPU=1 nohup ./run_slab_opportunistic.sh rna > logs/opp_rna.log 2>&1 &
>
> python monitor_slab.py --watch              # progress, ns/day, ETA, GPU load
> ```
>
> The benchmark step is already done — see *Cost on lyra*. Throughput was
> measured at ~5,000 steps/s on an idle L40S (~11 h per run, ~9.5 GPU-days for
> all 20), which is **8× cheaper than this README used to assume**.
>
> See *Sharing lyra* for which script belongs on which GPU, and why.

---

## What this is

Multi-chain **direct-coexistence (slab)** simulations. Every PARP14 simulation
so far has been single-chain, so nothing on disk can say whether a construct
self-associates. These give the actual quantitative observable — the saturation
concentration **c_sat** — via `calvados.analysis.SlabAnalysis`.

Two arms per construct:

| arm | contents | question |
|---|---|---|
| `homotypic/` | N copies of one construct | does it self-associate at all? |
| `rna/` | the same, plus polyU40 RNA | does RNA multivalency drive it? |

Plus `benchmark/` — the short calibration run. Its job is done: throughput has
now been measured on lyra (see *Cost on lyra*), so the panel no longer rests on
an assumed speedup.

## Why this runs on lyra and not pollux

`pollux.fraserlab.com` has **no GPU** — `nvidia-smi` fails and OpenMM sees only
`Reference` and `CPU`. Measured CPU throughput from this project's own `md_full`
runs is **5.76e5 bead-steps/s**, which puts a single 100-chain slab at roughly
**two years**. Every config here is `platform: CUDA`.

`lyra.fraserlab.com` has 4× NVIDIA L40S and a working CUDA build of OpenMM (in
`envs/calvados` — **not** `envs/CALVADOS`, see *Things that will bite*). Measured
there: ~350× pollux's CPU rate, i.e. ~11 h per run instead of ~2 years. lyra has
no SLURM, so `run_slab_queue.sh` replaces `submit_slab.slurm`.

---

## Run order

### 1. Benchmark — now measured on lyra (2026-09-17)

The 50×-CPU assumption is **superseded**. Two constructs were timed on an idle
L40S with this tree's own configs:

| construct | beads | steps/s | ns/day | bead-steps/s |
|---|---|---|---|---|
| `md2_md3` | 38,500 | 5,240 | 4,550 | 2.0e8 |
| `fl` | 54,030 | 6,360 | 5,490 | 3.4e8 |

That is **~350× CPU**, not 50×, so the panel is far cheaper than the README
previously claimed — see *Cost* below.

> **The formula this section used to give was wrong.** It read
> `bead_steps_per_s = 38500 * 1e6 / wall_seconds`, but with `slab_eq: true` the
> job runs `steps_eq` (5e6) *and then* `steps` (1e6) — 6e6 steps in that wall
> time, so it understated throughput ~6×. Take the rate from the production
> reporter log instead, which starts counting only after equilibration:
>
> ```bash
> # Step / Speed(ns/day) / Elapsed(s), production steps only
> cat benchmark/md2_md3/parp14_md2_md3.log
> ```
>
> `monitor_slab.py` already does this (it differences Step/Elapsed rather than
> trusting the running-average Speed column).

To re-cost the panel from your own number:

```bash
cd ..
python prepare_slab.py --report --panel full --gpu-rate <bead_steps_per_s>
```

Note that `--report` costs production only; each run also pays a one-off 5e6-step
equilibration (+2.5%).

### 2. Launch an arm

```bash
GPU=0 nohup ./run_slab_queue.sh homotypic > logs/queue_homotypic.log 2>&1 &
GPU=1 nohup ./run_slab_queue.sh rna       > logs/queue_rna.log       2>&1 &
```

Serial within an arm, one pinned GPU each, detached so it survives logout.
A failed construct is logged and the queue continues to the next.

*(If you do have a SLURM cluster, `submit_slab.slurm` runs the same tree as an
array job: `sbatch --array=0-9 submit_slab.slurm homotypic`.)*

### 3. Analyse

Per completed run, in its own directory:

```python
from calvados.analysis import SlabAnalysis
slab = SlabAnalysis(name=<sysname>, input_path='.', output_path='.',
                    input_pdb='top.pdb', input_dcd=None, centered_dcd='traj.dcd',
                    ref_chains=(0, nchain-1), ref_name='X', verbose=True)
slab.center(start=0, step=1, center_target='all')
slab.calc_profiles()
slab.calc_concentrations()      # writes {name}_ps_results.csv
slab.plot_density_profiles()
```

**c_sat is the `c_dilute` column** of `{name}_ps_results.csv`, in mM.

---

## Cost on lyra

Measured 2026-09-17 on an idle L40S, using this tree's own configs (not a
scaling argument). Per run: 2e8 production steps + 5e6 equilibration steps.

| | rate | per run | all 20 runs |
|---|---|---|---|
| **on one idle L40S** | ~4,970 steps/s | **~11.4 h** | **~9.5 GPU-days** |
| assumption this README used to make | 50× CPU | ~4 d | ~79 GPU-days |

So the panel is roughly **8× cheaper than previously documented** — the real
speedup over pollux's CPU is ~350×, not 50×. With all four GPUs to itself the
whole campaign is **~2.5 days wall-clock**; one GPU, both arms serial, is ~10
days.

These are **condensed-state** rates: each was timed after the full 5e6-step
`slab_eq` pull had actually formed the slab, not on the freshly-placed grid, so
they reflect what production really costs.

| construct | beads | spread (pre-`slab_eq`) | condensed (production) | penalty |
|---|---|---|---|---|
| `md2_md3` | 38,500 | 5,236 steps/s | 5,098 steps/s | 3% |
| `fl` | 54,030 | 6,360 steps/s | 4,851 steps/s | 24% |

Condensing costs less than expected — more pairs fall inside the 2.0/4.0 nm
cutoffs, but not enough to change the headline. Planning rate ~4,970 steps/s.

### Contention — the number that actually decides the schedule

lyra has 4 L40S and no scheduler, so runs land wherever you pin them and
compete with whatever else is on the card. This is not hypothetical; at the
time of the benchmark three of four GPUs were occupied by other work
(a 4-day `P2DFlow` training job, another user's `PLACER` jobs, and this
project's own BindCraft campaign).

Measured effect on the *same* construct:

| situation | ns/day | vs idle |
|---|---|---|
| idle GPU | 4,550 | 1.0× |
| sharing with one DL training job | 1,330 | **3.4× slower** |
| sharing during the busiest observed period | ~720 | **~6× slower** |

Contention is size-independent — `md2_md3` (38.5k beads) and `fl` (54k beads)
both dropped to the same ~830 steps/s when sharing, because they are
time-slicing the card rather than saturating it.

**Plan against the shared case, not the idle one.** A realistic bracket for the
full 20-run panel:

| | wall-clock |
|---|---|
| 4 idle GPUs | ~2.5 days |
| 2 GPUs, typical sharing | ~2 weeks |
| 1 shared GPU | ~1–2 months |

Check `nvidia-smi` before pinning, and prefer a card with no other compute
process. `monitor_slab.py` prints per-GPU utilisation next to run progress so a
throughput collapse is attributable rather than mysterious.

### Sharing lyra — which script on which GPU

Every run in this panel is ~11 h, so all of them count as long jobs. Confining
long jobs to one GPU keeps the rest of the machine usable by other people.

| GPU | what runs there | script |
|---|---|---|
| **0** | the campaign's own GPU. Long runs, start to finish, not interrupted. | `run_slab_queue.sh` |
| **1, 2** | opportunistic. Runs **only while the GPU is empty** and **gets off the moment anyone else appears**. | `run_slab_opportunistic.sh` |
| **3** | left alone | — |

**Seed first, then spread.** The opportunistic GPUs **refuse a construct that has
no `restart.chk`**, because equilibration cannot be preempted safely (see below).
So the campaign starts with a seeding pass on GPU 0 that takes every construct
through equilibration plus one checkpoint (~20 min each, ~6 h for all 20), after
which any GPU can pick up any construct:

```bash
# 1. seed: equilibrate every construct and leave a checkpoint behind
GPU=0 LEG_STEPS=2000000 nohup ./run_slab_queue.sh homotypic > logs/seed_homotypic.log 2>&1 &

# 2. then run: GPU 0 grinds through them, 1 and 2 help whenever they are free
GPU=0 nohup ./run_slab_queue.sh         homotypic > logs/queue_homotypic.log 2>&1 &
GPU=1 nohup ./run_slab_opportunistic.sh rna       > logs/opp_rna.log         2>&1 &
GPU=2 nohup ./run_slab_opportunistic.sh homotypic > logs/opp_homotypic.log   2>&1 &
```

**Why equilibration cannot be preempted.** `sim.py` runs the 5e6 `slab_eq` steps
as a single `simulation.step(steps_eq)` with only a DCD reporter — **no
checkpoint is written during it**. A preemption anywhere in those ~17 min throws
all of it away and starts over, so on a GPU reclaimed more often than that a
fresh construct would livelock and never reach production. Once `restart.chk`
exists the problem disappears entirely: `sim.py` forces `slab_eq` off on
checkpoint restart, so equilibration is done once and never repeated, and
everything after it checkpoints every `LEG_STEPS/10`.

`ALLOW_FRESH=1` overrides the refusal, but only makes sense on a GPU you expect
to stay quiet for 20 minutes.

You can start the opportunistic runners **at the same time as the seeding pass**.
They re-walk their arm until everything in it is finished, reporting
`nothing runnable yet (N still awaiting a seed on GPU 0)` and re-checking every
5 min, so they pick each construct up the moment its seed lands rather than
exiting early and leaving the GPU idle.

**If an opportunistic GPU is churning**, the runner reports `NO PROGRESS ... after
N yields` — that GPU is being reclaimed faster than the checkpoint interval and
is doing nothing useful. Lower `LEG_STEPS` or give the construct to GPU 0. This
is a real state, not a hypothetical: during testing GPU 2 was being reclaimed
within ~6 s.

**How yielding works.** The opportunistic runner works in short legs
(`LEG_STEPS`, default 2e7 ≈ 1 h idle). sim.py splits a leg into 10 checkpointed
batches, so a checkpoint lands every `LEG_STEPS/10` — ~7 min of work at the
default. A watcher polls the GPU every `POLL` seconds; as soon as a process that
is not ours appears it SIGTERMs the run, so at most one checkpoint interval is
lost. The next leg resumes from `restart.chk` and `slab_steps.py` keeps the
arithmetic exact, so the construct still finishes on precisely 2e8 steps.

Set `LEG_STEPS` smaller to yield faster at slightly more restart overhead; each
restart re-reads the PDB and rebuilds the restraint lists (~1 min for these
systems), so legs much below ~5e6 steps start to waste real time.

**A yielded run is not stranded on that GPU.** Resuming needs only
`restart.chk` and the construct directory, so the same construct can be picked
up on GPU 0 — or the other opportunistic GPU — simply by launching there. That
is what makes these runs movable as well as pausable.

**Collisions are prevented, not just discouraged.** Both scripts take an
exclusive `flock` on `<construct>/.slab.lock` for as long as they are working on
it, and skip a construct that is already locked. Two processes in one construct
directory would otherwise both write `restart.chk` and the DCD and corrupt each
other — which is the *likely* case when a queue and an opportunistic runner
share an arm, not an exotic one. You can still give each invocation a different
arm, or split one arm with `ONLY=`, but nothing breaks if they overlap.

### Packing two runs per GPU — untested, default to one

Memory imposes no limit at all — a run holds ~450 MiB of a 46 GB card, so
dozens would fit. And a single 50k-bead CG system almost certainly does not
saturate an L40S, so packing *might* raise aggregate throughput. That was
**not measured** — every GPU was busy with other users' work during the
benchmark window, so there was no fair test to run.

Default to one run per GPU. If you want the answer, the test is: run `md2_md3`
alone on an idle card, note ns/day from its production log, then start a second
copy on the same card and compare the *sum*. If the sum beats the single-run
rate, pack; if each simply halves, don't.

---

## The panel

10 constructs × 2 arms, chosen to span the two valences the multivalency
hypothesis turns on. Chain count is traded against construct size to keep every
run in the same cost bracket (~50k beads).

| construct | res | chains | box (nm) | RNA valence | ADPr readers | ART |
|---|---|---|---|---|---|---|
| `md2_md3` | 385 | 100 | 20×20×260 | 0 | 2 | – |
| `md_full` | 599 | 83 | 25×25×325 | 0 | 2 | – |
| `md3_wwe_full` | 396 | 100 | 20×20×260 | 2 | 2 | – |
| `mka_wwe_full` | 813 | 62 | 30×30×390 | 2 | 3 | – |
| `mka_full` | 1012 | 49 | 30×30×390 | 2 | 3 | ✓ |
| `core_full_go` | 1064 | 47 | 35×35×455 | 3 | 3 | ✓ |
| `kh1_wwe_full` | 1288 | 39 | 35×35×455 | 9 | 3 | – |
| `kh1_art_full` | 1487 | 34 | 40×40×520 | 9 | 3 | ✓ |
| `fl_wwe_full_go` | 1602 | 31 | 40×40×520 | 12 | 3 | – |
| `fl` | 1801 | 30 | 40×40×520 | 12 | 3 | ✓ |

Three matched ±ART pairs are included: `mka_wwe_full`↔`mka_full`,
`kh1_wwe_full`↔`kh1_art_full`, `fl_wwe_full_go`↔`fl`.

**Coverage gap worth knowing:** RNA valence {0, 2, 3, 9, 12} — the 4–8 band is
empty, and nothing sits at ADPr-reader 0 (RNA-binding domains with no readers
would be a clean negative control). Both need new AF3 structures.

---

## Each construct directory

```
config.yaml        CALVADOS simulation config (platform: CUDA, topol: slab)
components.yaml    system composition; absolute paths, so do not move the tree
run.py             entry point, called by run_slab_queue.sh (or the SLURM script)
slab_meta.yaml     construct, arm, nchain, n_rna, box, steps, beads
input/
  <sysname>.pdb                starting structure
  domains.yaml                 restraint domains, LOCAL copy so the run is
                               self-contained
  residues_CALVADOS3.csv       (homotypic) or
  residues_CALVADOS3_RNA.csv   (rna arm: C3 amino acids + RNA bead rows)
  rna.fasta                    (rna arm only) polyU40
```

`components.yaml` carries **absolute paths**. Moving or copying the tree
elsewhere will break it — re-run `prepare_slab.py` at the new location instead.

## Tooling in this directory

| file | what it does |
|---|---|
| `run_slab_queue.sh` | **the launcher for GPU 0.** Serial queue over one arm, one pinned GPU, detached, runs each construct start to finish. Refuses to start if the env has no CUDA platform. `ONLY=a,b` restricts it to a subset. |
| `run_slab_opportunistic.sh` | **the launcher for GPUs 1–2.** Same queue, but in short checkpointed legs, and it only runs while the GPU is empty — it stops the moment another user's process appears and resumes when they are done. `LEG_STEPS`, `POLL`, `ONLY` tune it. |
| `monitor_slab.py` | progress / ns/day / ETA per run, plus GPU load. Finds live runs in `/proc`, so it works regardless of how they were started and has no state file to go stale. `--watch` to refresh. |
| `slab_steps.py` | production-step bookkeeping. `status` reports done/target/remaining; `prepare` rewrites `config.yaml`'s `steps` to the remainder so a re-run of the queue resumes to exactly 2e8 instead of overshooting. Called by the queue script. |
| `submit_slab.slurm` | for an actual SLURM cluster. **Not usable on lyra**, and it still names the CUDA-less `envs/CALVADOS`. |

## Protocol

- `topol: slab`, `slab_eq: true`, `k_eq: 0.02`, `steps_eq: 5e6`
  (a linear pull toward z = Lz/2 that is removed before production)
- `steps: 2e8` (2 µs at dt = 0.01 ps), `wfreq: 1e5`
- 293.15 K, 0.15 M ionic, pH 7.5
- domain restraints on, from each construct's own `domains.yaml`

**Restarts continue rather than restart.** CALVADOS restarting from
`restart.chk` runs `steps` *additional* steps and appends to the DCD, so
re-launching an interrupted construct picks up where it stopped.

---

## Things that will bite

- **The conda env matters, and the obvious one is wrong.**
  `/home/sbali/miniconda3/envs/CALVADOS` has an OpenMM installed from a **pip
  wheel, which ships no CUDA plugin** — it exposes only Reference, CPU and
  OpenCL, so every `platform: CUDA` config here dies at startup with
  `OpenMMException: There is no registered Platform called "CUDA"`. Verified on
  lyra 2026-09-17. Use the lowercase **`/home/sbali/miniconda3/envs/calvados`**,
  which has CUDA and numpy 1.24, and whose `sim.py` is byte-identical to this
  checkout's. `run_slab_queue.sh` now points there and refuses to start if the
  CUDA platform is missing; `submit_slab.slurm` still names the broken env.
- **Do not set `slab_eq` and `ext_force` together.** `sim.py:47-48` overwrites
  the slab-centering force with the external one.
- **NumPy ≥ 2.0 breaks setup.** `build.build_xyzgrid` uses `np.product`, removed
  in NumPy 2. `envs/calvados` pins 1.24 — do not run this under a NumPy 2
  interpreter.
- **Re-running the queue would overshoot the target.** A checkpoint restart runs
  `steps` *additional* steps, so a second pass over a finished tree would push
  every construct to 4e8. `run_slab_queue.sh` now calls `slab_steps.py prepare`
  first, which reads the production steps already done and rewrites `steps` to
  the remainder (skipping constructs that are complete). `slab_meta.yaml` holds
  the authoritative 2e8 target; `config.yaml`'s `steps` is the working value.
- **Restarting skips equilibration by design.** If `restart.chk` exists,
  `sim.py:40-42` silently forces `slab_eq` off. Correct, but surprising.
- **Check the interface fit converged.** `SlabAnalysis.fit_profile` prints
  `NOT CONVERGED` when the left/right cutoff ratio is off. That check used to be
  dead code after a `return` — it is fixed in this checkout, but if you run
  against a different CALVADOS install, a failed fit will pass silently and
  c_sat will be wrong. Sanity-check the `cutoffs_*` columns either way.
- **No SLURM on lyra** — see the box at the top. `submit_slab.slurm` is kept for
  a cluster that has one; `run_slab_queue.sh` is what works on lyra.
- **OOM is not remotely a risk.** Measured on the running seed job: a 50k-bead
  slab run holds **~450 MiB** on a 46 GB L40S — essentially just the CUDA
  context. These are 50k coarse-grained beads, not an all-atom system, so the
  BindCraft unified-memory escape hatch is irrelevant here. **Memory is never
  the constraint; compute is** (see *Contention*).

  (An earlier draft of this file claimed ~7 GB per run. That was wrong — it
  misattributed a co-resident PLACER job's memory to ours. Check
  `nvidia-smi --query-compute-apps=pid,used_memory` per PID before believing a
  per-run memory figure on this machine.)
- **lyra is shared, and that dominates the schedule.** See *Contention* below.

## Not included

**There is no ADPr-substrate arm.** `PTMProtein` in CALVADOS has no example and
no test, sets `c_termini` to the last PTM bead (`components.py:737`), and never
reads from PDB — so it is incompatible with `restraint: true`, which every
multi-domain PARP14 construct needs. No ADP-ribose bead parameters exist either.
The viable route is putting ADPr beads on a short **unrestrained substrate
peptide**, which sidesteps the restraint bug. That is separate work.

Also note the RNA arm mixes force fields: the RNA bead parameters were fit
alongside CALVADOS2 while the protein parameters here are CALVADOS3, and the RNA
model has **no base identity** (polyU is the only meaningful choice). Any RNA-arm
result inherits that approximation.

---

Regenerate any of this with:

```bash
cd /home/sbali/CALVADOS/examples/PARP14_MDP
python prepare_slab.py --benchmark
python prepare_slab.py --arm both --panel full
```
