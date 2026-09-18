#!/bin/bash
# PARP14 slab coexistence queue -- direct launch, no SLURM.
#
# lyra.fraserlab.com has 4x NVIDIA L40S (46 GB) and NO SLURM
# (bindcraft_md/GPU_SERVER_SETUP.md), so submit_slab.slurm does not apply there.
# This mirrors the pattern that actually works on that host:
# sweep/run_queue.sh -- serial on one pinned GPU, detached, survives logout.
#
#   GPU=0 nohup ./run_slab_queue.sh homotypic > logs/queue_homotypic.log 2>&1 &
#   GPU=1 nohup ./run_slab_queue.sh rna       > logs/queue_rna.log       2>&1 &
#
# Optionally wait for another PID to finish first (chaining onto a running job):
#   GPU=0 ./run_slab_queue.sh homotypic 12345
#
# ONLY=a,b restricts the queue to those constructs, so the 4 GPUs can be split
# across one arm rather than one GPU per arm:
#   GPU=0 ONLY=md2_md3,md_full,md3_wwe_full,mka_wwe_full,mka_full ./run_slab_queue.sh homotypic
#   GPU=1 ONLY=core_full_go,kh1_wwe_full,kh1_art_full,fl_wwe_full_go,fl ./run_slab_queue.sh homotypic
#
# Progress:  tail -f logs/queue_<arm>.log
#            python monitor_slab.py --watch
#            nvidia-smi          (shared machine -- check before pinning)
set +u
SLAB=$(cd "$(dirname "$0")" && pwd)
# The env whose OpenMM actually has a CUDA platform. NOT envs/CALVADOS: that
# one's openmm is a pip wheel with no CUDA plugin, so it exposes only
# Reference/CPU/OpenCL and every `platform: CUDA` run dies at startup with
#   OpenMMException: There is no registered Platform called "CUDA"
# (verified on lyra 2026-09-17). envs/calvados has CUDA + numpy 1.24, and its
# sim.py is byte-identical to this checkout's.
CAL_ENV=/home/sbali/miniconda3/envs/calvados
ARM=${1:?usage: [GPU=n] ./run_slab_queue.sh <homotypic|rna|benchmark> [wait_pid]}
WAIT_PID=${2:-0}
GPU=${GPU:-0}
mkdir -p "$SLAB/logs"

if ! "$CAL_ENV/bin/python" -c 'import openmm; openmm.Platform.getPlatformByName("CUDA")' 2>/dev/null; then
  echo "FATAL: $CAL_ENV has no CUDA platform in OpenMM -- every run would die at startup."
  exit 1
fi

ARM_DIR="$SLAB/$ARM"
[ -d "$ARM_DIR" ] || { echo "no such arm: $ARM_DIR"; exit 1; }

if [ "$WAIT_PID" != "0" ]; then
  echo "[slab $(date '+%F %T')] waiting for PID $WAIT_PID"
  while kill -0 "$WAIT_PID" 2>/dev/null; do sleep 120; done
fi

echo "[slab $(date '+%F %T')] arm=$ARM gpu=$GPU host=$(hostname)"
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv,noheader 2>/dev/null \
  || echo "  (nvidia-smi unavailable -- are you on the GPU host?)"

for d in $(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort); do
  name=$(basename "$d")
  if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then continue; fi
  beads=$(grep -oP 'beads:\s*\K[0-9]+' "$d/slab_meta.yaml" 2>/dev/null)
  echo ""
  echo "[slab $(date '+%F %T')] START $name (${beads:-?} beads)"
  # One launcher per construct directory. Two processes in the same directory
  # would both write restart.chk and the DCD and corrupt each other; this is
  # easy to hit when a queue and an opportunistic runner share an arm.
  exec 9>"$d/.slab.lock"
  if ! flock -n 9; then
    echo "  ^ already running under another launcher -- skipping"
    exec 9>&-
    continue
  fi
  # Set `steps` to what is actually LEFT. A checkpoint restart runs `steps`
  # ADDITIONAL steps, so without this a re-run of the queue would push
  # finished constructs past 2e8 and overshoot partial ones. rc=3 means the
  # construct already hit its target.
  # LEG_STEPS caps this leg instead of running the construct to completion.
  # Its main use is a SEEDING PASS: `LEG_STEPS=2000000` takes each construct
  # through equilibration and one checkpoint, then moves on, so the
  # opportunistic GPUs have something they are allowed to pick up (they refuse
  # un-equilibrated constructs -- see run_slab_opportunistic.sh).
  if [ -n "$LEG_STEPS" ]; then
    "$CAL_ENV/bin/python" "$SLAB/slab_steps.py" prepare "$d" --max-leg "$LEG_STEPS"
  else
    "$CAL_ENV/bin/python" "$SLAB/slab_steps.py" prepare "$d"
  fi
  if [ $? -eq 3 ]; then exec 9>&-; continue; fi
  start=$(date +%s)
  # A 46 GB L40S handles ~50k CG beads comfortably; the OOM seen with BindCraft
  # was a 683-residue all-atom AF2 complex, a different regime. If a large
  # construct does OOM, the same escape hatch applies:
  #   TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
  ( cd "$d" && CUDA_VISIBLE_DEVICES=$GPU "$CAL_ENV/bin/python" -u run.py ) \
      >> "$SLAB/logs/${ARM}_${name}.log" 2>&1
  rc=$?
  echo "[slab $(date '+%F %T')] END $name exit=$rc elapsed=$(( $(date +%s) - start ))s"
  if [ $rc -ne 0 ]; then
    echo "  ^ FAILED -- see logs/${ARM}_${name}.log. Continuing to the next construct."
  fi
  exec 9>&-    # release the per-construct lock
done

echo ""
echo "[slab $(date '+%F %T')] arm $ARM complete"
