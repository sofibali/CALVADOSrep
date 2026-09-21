# PARP14 slab campaign — lyra session log, 2026-09-17/18

Everything done in this session, why, and what state the tree is in now.
The campaign was **started and then stopped** on lyra; it is being moved to
slower but less contended GPUs. Nothing here is lost — see *State of the tree*.

---

## TL;DR

- The launchers pointed at a conda env **with no CUDA platform**. Every run
  would have died at startup. Fixed.
- Throughput was **assumed** (50× CPU). It is now **measured**: ~4,970 steps/s
  on an idle L40S → **~11.4 h per run, ~9.5 GPU-days for all 20**, not the
  ~79 GPU-days this tree used to claim.
- All 20 configs verified to build and run on CUDA with correct bead counts.
- Built three tools: a polite opportunistic runner, a monitor, and exact
  remaining-step arithmetic. All tested end-to-end against real contention.
- **lyra is the wrong machine for this.** Not because it is slow — it is fast —
  but because it is shared with no scheduler, and sharing costs 3.4–6×. GPU 0
  in particular sits under a BindCraft sweep with ~2 weeks left.
- One construct (`homotypic/core_full_go`) is equilibrated with 1e6 production
  steps banked. It resumes rather than restarts.

---

## 1. What was wrong before anything could run

### 1.1 The conda env had no CUDA (this was the blocker)

`run_slab_queue.sh` and `submit_slab.slurm` both used
`/home/sbali/miniconda3/envs/CALVADOS`. Its OpenMM is a **pip wheel with no
CUDA plugin**:

```
$ envs/CALVADOS/bin/python -c 'import openmm; openmm.Platform.getPlatformByName("CUDA")'
OpenMMException: There is no registered Platform called "CUDA"
```

It exposes only Reference, CPU and OpenCL. Every config here is
`platform: CUDA`, so **all 20 runs would have failed in the first seconds**.

The working env is the lowercase **`/home/sbali/miniconda3/envs/calvados`**:
CUDA present, numpy 1.24 (needed — `build.build_xyzgrid` uses `np.product`,
removed in numpy 2), and its `calvados/sim.py` is **byte-identical** to this
checkout's. Its `components.py` is slightly newer (adds `comp_setup` and an
`init_restraint_force` stub); its `analysis.py` still has the old
`fit_profile` bug, so **analyse with the repo copy, run with the env copy**.

Fixed in `run_slab_queue.sh`, and it now refuses to start at all if the CUDA
platform is missing rather than failing 20 times in a row.

### 1.2 The documented benchmark formula was wrong

The README said:

```
bead_steps_per_s = 38500 * 1e6 / wall_seconds
```

With `slab_eq: true` the job runs `steps_eq` (5e6) **and then** `steps` (1e6),
so that wall time covers 6e6 steps, not 1e6 — understating throughput ~6×.
The rate must come from the production reporter log, which is only attached
*after* equilibration (`sim.py` builds a fresh `Simulation` once the centering
force is removed), so its Step column is production-only.

---

## 2. The measurements

All on NVIDIA L40S (46 GB), env `envs/calvados`, this tree's own configs.

### 2.1 Throughput, condensed state

Timed **after** the full 5e6-step `slab_eq` had actually formed the slab — not
on the freshly-placed grid, which would flatter the numbers.

| construct | beads | spread (pre-eq) | condensed (production) | penalty |
|---|---|---|---|---|
| `md2_md3` | 38,500 | 5,236 steps/s | 5,098 steps/s | 3% |
| `fl` | 54,030 | 6,360 steps/s | 4,851 steps/s | 24% |

**Planning rate ~4,970 steps/s.** Per run (2e8 production + 5e6 eq =
2.05e8 steps): **~11.4 h**. All 20 runs: **~9.5 GPU-days**.

That is ~350× pollux's CPU rate (5.76e5 bead-steps/s), not the assumed 50×.
The old estimate of ~79 GPU-days was **8× too pessimistic**.

### 2.2 Contention — the number that actually mattered

| situation | ns/day | vs idle |
|---|---|---|
| idle GPU | 4,550 | 1.0× |
| sharing with one DL training job | 1,330 | 3.4× slower |
| busiest observed | ~720 | ~6× slower |

Size-independent: `md2_md3` (38.5k beads) and `fl` (54k) both fell to
~830 steps/s when sharing. They are time-slicing the card, not saturating it.

Measured on GPU 0 under the live BindCraft sweep, the rate swung between
**389 and 3,000 steps/s** depending on BindCraft's phase, averaging
**~940 steps/s** — a 5.3× penalty. At that rate a single construct's full
production is ~60 h, and the whole panel on GPU 0 alone would be **~50 days**
instead of ~9.5.

### 2.3 GPU memory — a figure I got wrong and corrected

A slab run holds **~452 MiB**, verified per-PID on the live job. Memory is
never a constraint; dozens of runs would fit on one card.

I earlier wrote ~7 GB in the README. That was wrong: I read total per-GPU
memory while another user's PLACER job (~6.8 GB) was co-resident and attributed
it to us. **On this machine always check
`nvidia-smi --query-compute-apps=pid,used_memory` per PID**, never the
per-GPU total.

### 2.4 Config validation

All 20 built and reached `STARTING SIMULATION` on CUDA with bead counts
matching `slab_meta.yaml` — 10/10 homotypic, 10/10 rna. The rna arm (different
residue table, RNA components) was included; nothing failed.

---

## 3. What was built

| file | status | what it does |
|---|---|---|
| `run_slab_queue.sh` | **modified** | env fixed; CUDA preflight; per-construct `flock`; `ONLY=` subset; `LEG_STEPS=` for seeding; calls `slab_steps.py` so re-running never overshoots |
| `run_slab_opportunistic.sh` | **new** | polite runner for shared GPUs: runs only while the GPU is empty, yields on any foreign process, resumes after |
| `monitor_slab.py` | **new** | progress / ns-per-day / ETA per run + GPU load |
| `slab_steps.py` | **new** | remaining-step arithmetic; `status` and `prepare [--max-leg N]` |
| `.gitignore` | **new** | run artifacts (the `.xml` system dumps are ~40 MB each and were being committed) |
| `README.md` | **rewritten in parts** | cost, contention, sharing policy, corrected gotchas |
| `docs/questions/Q9_self_association_puncta.md` | **modified** | compute section replaced with measured numbers |

### 3.1 `slab_steps.py` — why it exists

CALVADOS restarting from `restart.chk` runs `steps` **ADDITIONAL** steps and
appends to the DCD. So re-running a queue over a finished tree would push every
construct to 4e8, and resuming a half-done one would overshoot by whatever it
had already done. Nothing in the tree guarded against this.

`prepare` reads production steps already done and rewrites `config.yaml`'s
`steps` to the remainder, exiting 3 when the construct is complete.
`slab_meta.yaml` holds the authoritative 2e8 target; `config.yaml`'s `steps` is
the working value.

Verified: 2,000 → resume → ran exactly 3,000 more → landed on exactly 5,000,
with the log confirming `Reading check point file` / `Appending trajectory`.

### 3.2 `run_slab_opportunistic.sh` — the yield mechanism

Work is done in legs of `LEG_STEPS`. `sim.py` splits a leg into 10 checkpointed
batches, so a checkpoint lands every `LEG_STEPS/10` — that is the most a
SIGTERM can destroy. A watcher polls the GPU; when a PID that is not ours
appears it stops the run, waits for the GPU to clear, and resumes from the
checkpoint.

Tested against real traffic, not simulated: it waited while another user's job
held GPU 1, started when free, and yielded within seconds when the next job
landed, then went back to waiting.

Three things had to be fixed during testing:

1. **The watched PID must be the CUDA process.** The run is launched as
   `( cd "$d" && CUDA_VISIBLE_DEVICES=$GPU exec python run.py ) &` — the `exec`
   is load-bearing. Without it `$!` is the subshell, nvidia-smi reports the
   python child, and the runner would see its own job as foreign and yield in
   an infinite loop.
2. **Equilibration cannot be preempted.** `sim.py` runs the 5e6 `slab_eq` steps
   as a single `simulation.step(steps_eq)` with **no checkpointing**. A
   preemption anywhere in those ~17 min throws all of it away. On a GPU
   reclaimed more often than that, a fresh construct livelocks forever. Caught
   because a test yield reported `progress=0/4000000`. The runner now
   **refuses a construct with no `restart.chk`** (override: `ALLOW_FRESH=1`).
3. **A single pass was not enough.** It walked the arm once and exited, so
   starting it before seeding finished left the GPU idle permanently. It now
   re-walks until the arm is complete, reporting
   `nothing runnable yet (N awaiting a seed)` every 5 min.

### 3.3 Locking

Both launchers take an exclusive `flock` on `<construct>/.slab.lock` and skip a
construct that is already locked. Without it, a queue and an opportunistic
runner sharing an arm would both write the same `restart.chk` and DCD and
corrupt each other — the likely case, not an exotic one. Verified: with the
lock held, the second launcher printed
`already running under another launcher -- skipping`.

### 3.4 The seeding pass

Because opportunistic GPUs refuse un-equilibrated constructs, the campaign has
to start with a seeding pass on an uninterrupted GPU:

```bash
GPU=0 LEG_STEPS=2000000 ./run_slab_queue.sh homotypic
```

This takes each construct through equilibration plus one checkpoint, then moves
on. After that `sim.py` forces `slab_eq` off on every restart, so equilibration
never repeats and all remaining work is freely preemptible.

---

## 4. What was actually run, and stopped

1. Seeding launched on GPU 0 (homotypic, rna chained behind it).
2. GPU 0 turned out to be under a BindCraft sweep — averaged ~940 steps/s, with
   stalls to 389. Seeding one construct there was ~2.1 h instead of ~23 min.
3. On your instruction (BindCraft is priority), GPU 0 work was **stopped and the
   GPU handed back**. `core_full_go` lost 42/50 equilibration frames — the right
   trade, since redoing all 50 on a free GPU takes ~17 min.
4. Seeding moved to GPU 2, and GPU 1 when free. `core_full_go` completed
   equilibration on GPU 2 and banked 1e6 production steps.
5. GPU 0 was chained behind the **whole** BindCraft sweep (PID 427444, the
   `sweep/run_queue.sh` driver) rather than the current target, since the sweep
   runs 9 targets sequentially.
6. **Everything of ours was then stopped** for the move. Only BindCraft (GPU 0)
   and czou's P2DFlow (GPU 3) remain.

### BindCraft, for planning

It is target **1 of 9** (`clamp_md1md2_state3`), 26.5 h in, 3 accepted /
15 trajectories against a cap of 50. At ~1.8 h/trajectory that target has
~2 days left and the **full sweep is realistically 2+ weeks**. GPU 0 is not a
near-term option.

---

## 5. State of the tree

Clean and ready to move. Verified after stopping:

- **All 20 `config.yaml` restored to `steps: 200000000`.** The seeding pass had
  patched two of them down to 2e6; they are back at the authoritative target.
- **Stale `.slab.lock` files removed.**
- **One construct has real, resumable progress:**

  | construct | equilibration | production steps | resumes? |
  |---|---|---|---|
  | `homotypic/core_full_go` | **done** | 1,000,000 / 200,000,000 | yes, skips eq |

  Everything else is untouched, at 0 steps.

- `homotypic/fl` has leftover `bonds_*.txt` / `restr_*.txt` from an interrupted
  system build. Harmless — regenerated every run.

To confirm state anywhere at any time:

```bash
python monitor_slab.py --arm both
python slab_steps.py status homotypic/core_full_go   # -> "1000000 200000000 199000000"
```

---

## 6. Moving to the other GPUs

The tree is **not portable as-is**: `components.yaml` carries absolute paths
(`/home/sbali/CALVADOS/examples/PARP14_MDP/slab/...`).

- **Same filesystem, different GPUs** → nothing to do. Just launch there.
- **Different machine/path** → re-run `python prepare_slab.py --arm both
  --panel full` at the new location. Do not hand-edit the paths. Note this
  regenerates the inputs, so copy `homotypic/core_full_go/restart.chk` (plus its
  `.dcd` and `parp14_core_full_go.log`) back afterwards if you want to keep
  those 1e6 steps.

What to check on the new machine:

```bash
# 1. the env has CUDA -- this is the thing that silently breaks everything
<env>/bin/python -c 'import openmm; print(openmm.Platform.getPlatformByName("CUDA"))'
# 2. numpy < 2
<env>/bin/python -c 'import numpy; print(numpy.__version__)'
```

If the env path differs from `/home/sbali/miniconda3/envs/calvados`, update
`CAL_ENV` at the top of `run_slab_queue.sh` and `run_slab_opportunistic.sh`.

**If the new GPUs are dedicated** (no other users), the opportunistic runner is
unnecessary — use `run_slab_queue.sh` on each GPU with `ONLY=` to split the
arms, skip the seeding pass, and let each construct run start to finish.

**Re-benchmark there.** Scale the measured L40S numbers by whatever the new card
does; the panel cost is directly proportional:

```
hours_per_run = 2.05e8 / steps_per_second / 3600
```

An idle card at half an L40S's speed still finishes the panel in ~19 GPU-days,
which beats a contended L40S at ~50. **Dedicated and slower is the right call**
— that is the whole lesson of this session.

---

## 7. Gotchas worth carrying forward

- `envs/CALVADOS` has no CUDA; use `envs/calvados`. Run with the env's
  `calvados`, analyse with the repo's (only the repo has the `fit_profile`
  convergence fix).
- A checkpoint restart runs `steps` **additional** steps. Always go through
  `slab_steps.py prepare`.
- `slab_eq` is skipped entirely on checkpoint restart (`sim.py:40-42`) — correct,
  but surprising, and it is what makes the seed-then-spread workflow work.
- Equilibration is **not** checkpointed. Never start a fresh construct on a GPU
  you cannot hold for ~20 minutes.
- Do not set `slab_eq` and `ext_force` together — `sim.py:47-48` clobbers
  `rcent`.
- Per-run GPU memory is ~450 MiB. If you see GB-scale numbers, you are reading
  someone else's process.
- `SlabAnalysis.fit_profile` prints `NOT CONVERGED` on a bad interface fit —
  but only in this checkout. Sanity-check the `cutoffs_*` columns regardless;
  c_sat is the `c_dilute` column of `{name}_ps_results.csv`, in mM.
- The `.xml` system dumps are ~40 MB per construct and were being committed to
  git. `.gitignore` added; the already-tracked ones need
  `git rm --cached` to stop following them.
