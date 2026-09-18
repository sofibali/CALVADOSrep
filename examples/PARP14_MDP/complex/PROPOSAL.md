# PARP14 / PARP9 / DTX3L — small-box association simulations

Prepared 2026-09-17. Everything here is ready to run; nothing has been launched.

---

## 1. Why these runs exist

The existing `p14fl_test/complexes/` runs **cannot measure binding**. All 155 of them
were re-analysed (500 ps sampling, minimum-image inter-chain distances):

- **0 of 51** system/variant/chain-pair combinations keeps all replicates bound for 100 ns.
- Chains start at the AF3 docked pose (0.3–0.7 nm at frame 0) and finish **13–55 nm apart**.
- Bound fraction 0.04–0.36; median time to final unbinding 9–62 ns.

Cause: `prepare_complex.py` sets `BOX_MARGIN = 80.0`, so boxes are 94–110 nm. One
molecule in (100 nm)³ is **1.7 µM** — the re-encounter time is far beyond 100 ns.
Those runs measure unbinding only.

Two further results shaped this design:

- The dominant `parp14:MD3 – dtx3l:linker` contact (13.1 contacts/frame) is
  **restraint-induced**: a difference map against the unrestrained control shows
  +25 contacts/frame under Gō-KH and nothing without it. Gō-KH is therefore excluded here.
- HyRes, under matched conditions from the same AF3 start, keeps PARP14–PARP9 in
  contact **99.98%** of the time. The two force fields disagree about whether these
  complexes exist. That disagreement is the open question these runs address.

## 2. Design

| Choice | Value | Why |
|---|---|---|
| Box | 40 nm cubic | Smallest that does not distort PARP14 (see below) |
| Concentration | 25.9 µM per chain | 15× the old 1.7 µM; encounter becomes frequent |
| Length | 500 ns × 20 replicates | 10 µs per set, 120 µs total |
| Timestep | 0.01 ps | Stock CALVADOS value — **no patch required** |
| Restraints | Harmonic k=700 + split-KH custom | Preserves folds (FNC 0.62–0.67) |
| Gō-KH | excluded | Manufactures the MD3–linker contact |

**Box sizing.** A chain must not see its own periodic image: `L > chain span + cutoff_yu`
(4.0 nm). Spans over the existing trajectories, 99th percentile: PARP14 32.0 nm,
DTX3L 26.8 nm, PARP9 20.8 nm. PARP14 sets the floor at 36 nm. 40 nm gives margin and
puts every system at the **same concentration**, which is required before bound
fractions can be compared across systems.

**Two arms.** Separated start (`topol: grid`, chains ~17 nm apart) and docked start
(AF3 pose, `restart: pdb`, 0.42–0.49 nm). If they converge on the same bound fraction
the result is thermodynamic; if not, it is hysteresis and 500 ns is too short. This is
the test the HyRes docked/separated pair failed.

**Pair counting.** A 1+1 heterodimer box and a 2× homodimer box each contain exactly
one potential pair, so bound fractions compare directly. The ternary box holds three
pairs at the same per-partner concentration, so ternary − dimer isolates the third chain.

## 3. What is prepared

| Set | Contents | Chains | Beads | Docked arm |
|---|---|---|---|---|
| `ternary` | PARP14 + PARP9 + DTX3L | 3 | 3,395 | yes |
| `p14_p9` | PARP14 + PARP9 | 2 | 2,655 | yes |
| `p14_dtx3l` | PARP14 + DTX3L | 2 | 2,541 | yes |
| `p9_dtx3l` | PARP9 + DTX3L | 2 | 1,594 | yes |
| `p14_homo` | PARP14 × 2 | 2 | 3,602 | **no AF3 model** |
| `p9_homo` | PARP9 × 2 | 2 | 1,708 | **no AF3 model** |
| `dtx3l_homo` | DTX3L × 2 | 2 | 1,480 | yes |

7 separated sets (140 replicates) + 5 docked sets (100) = **240 replicates, 120 µs**.

`p14_homo` and `p9_homo` have no docked control because no AF3 complex prediction
exists. If the separated arm under-samples binding, those two carry the bias
undetected. An AF3 PARP14 homodimer prediction would close this — worth doing, since
PARP14 self-association dominated the slab condensate.

## 4. Transfer and setup

```bash
# on the GPU host
tar xzf binding.tar.gz                 # -> binding/  (12 sets, 240 replicates, 97 MB)
export PYTHON_EXE=/path/to/calvados/env/bin/python
bash preflight.sh                      # GPU + CALVADOS + a real 2000-step run
```

`preflight.sh` must pass before launching. It checks the GPU, that OpenMM can create a
CUDA context, that CALVADOS imports, whether the build is patched or stock (either is
fine at dt 0.01), and that one replicate actually produces a trajectory.

All paths inside the runs are **relative**; symlinks are relative and survive the
tarball (verified: 0 broken after extraction). The only absolute path is the default
`PYTHON_EXE` in each `run.sh`, which the env var above overrides.

## 5. Launch commands

Run from the directory containing `binding/`. The numeric argument is replicates
packed per GPU; with 8 GPUs, `3` gives 24 slots so all 20 replicates of a set run at once.

| # | Command | Beads | Sim length | Est. wall |
|---|---|---|---|---|
| 1 | `bash binding/p9_dtx3l/run.sh 3` | 1,594 | 20 × 500 ns = 10 µs | ~40 min |
| 2 | `bash binding/p9_dtx3l_docked/run.sh 3` | 1,594 | 10 µs | ~40 min |
| 3 | `bash binding/dtx3l_homo/run.sh 3` | 1,480 | 10 µs | ~40 min |
| 4 | `bash binding/dtx3l_homo_docked/run.sh 3` | 1,480 | 10 µs | ~40 min |
| 5 | `bash binding/p9_homo/run.sh 3` | 1,708 | 10 µs | ~40 min |
| 6 | `bash binding/p14_dtx3l/run.sh 3` | 2,541 | 10 µs | ~40 min |
| 7 | `bash binding/p14_dtx3l_docked/run.sh 3` | 2,541 | 10 µs | ~40 min |
| 8 | `bash binding/p14_p9/run.sh 3` | 2,655 | 10 µs | ~40 min |
| 9 | `bash binding/p14_p9_docked/run.sh 3` | 2,655 | 10 µs | ~40 min |
| 10 | `bash binding/ternary/run.sh 3` | 3,395 | 10 µs | ~40 min |
| 11 | `bash binding/ternary_docked/run.sh 3` | 3,395 | 10 µs | ~40 min |
| 12 | `bash binding/p14_homo/run.sh 3` | 3,602 | 10 µs | ~40 min |
| — | `bash binding/run_all.sh 3` | all 12 | **120 µs** | **~8 h** |

Every replicate is 5×10⁷ steps. Timing basis: Vincent's own logs at the identical step
count, 1,180–2,384 s per replicate. The spread tracks GPU contention, not system size.
Total ≈ **105–125 GPU-hours**, 50–90 GB of trajectory.

`run.sh` skips any replicate that already has a `.dcd`, so re-running after an
interruption is safe.

**Restricting GPUs:**
```bash
GPU_LIST="0 1" bash binding/ternary/run.sh 2    # 2 GPUs, 2 each = 4 concurrent
GPU_LIST="5"   bash binding/p9_dtx3l/run.sh 1   # single GPU, serial
```

**Recommended order:** run rows 1–2 first and stop. If separated and docked disagree on
PARP9–DTX3L — the one interface that reproduced across molecular context — then 500 ns
is too short and the remaining sets need a longer budget before they are worth spending.

## 6. Box-size control

If bound fractions are thermodynamic they must shift predictably with concentration.
Generate a 50 nm (13.3 µM) arm and compare:

```bash
python prepare_binding.py --all --start both --box 50 --tag box50
bash binding/run_all.sh 3
```

`prepare_binding.py` is included. Other useful flags: `--ns` (length per replicate),
`--nreps`, `--box`, `--start {separated,docked,both}`.

## 7. Analysis

The three scripts used for the convergence work are included and apply unchanged to
the new runs — edit the `ROOT` path at the top of `analyze_convergence.py`:

```bash
python analyze_convergence.py    # per-replicate observables -> data/raw.npy
python report2.py                # bound fractions, dissociation kinetics, tables
python figmaps.py                # bound-window contact maps + difference maps
```

The primary readout is `dissociation_table.csv`: `bound_frac_mean` per set. The
comparisons that matter are separated vs docked (is it equilibrium?), ternary vs dimer
(does the third chain change the interface?), and hetero vs homo (does PARP14 prefer
itself?).

## 8. Known limitations

- 25.9 µM is far above physiological. This is an enhanced-sampling device to make
  encounter frequent, not a physiological concentration.
- CALVADOS is parameterised for disordered regions; folded-domain surfaces may be
  under-sticky. That is the likely source of the disagreement with HyRes, and these
  runs will quantify it but not resolve which force field is right.
- No ADP-ribose or ubiquitin in any model, so nothing here addresses catalysis.
- Starting structures inherit AF3 interfaces at ipTM 0.30–0.51. The separated arm is
  the mitigation: it does not use the predicted pose at all.
