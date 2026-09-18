#!/bin/bash
# PARP14 slab -- OPPORTUNISTIC runner for the shared GPUs (1 and 2).
#
# POLICY THIS IMPLEMENTS
# ----------------------
# lyra is shared and has no scheduler, so the campaign is split in two:
#
#   GPU 0        long, uninterrupted runs. `run_slab_queue.sh`, one construct
#                start to finish (~11 h each). This is the campaign's own GPU.
#   GPU 1, 2     opportunistic. THIS script. It only runs while the GPU is
#                otherwise empty, and it GETS OFF as soon as anyone else's
#                process appears -- their job should not have to queue behind
#                a multi-day MD campaign.
#   GPU 3        left alone.
#
# HOW YIELDING WORKS
# ------------------
# Work is done in short legs (default 2e7 steps, ~1 h on an idle L40S).
# sim.py splits a leg into 10 checkpointed batches, so a checkpoint lands
# every leg/10 -- ~7 min of work at the default. A watcher polls the GPU; the
# moment a process that is not ours shows up, it SIGTERMs the run. At most one
# checkpoint interval is lost, and the next leg resumes from `restart.chk`
# exactly where it stopped (slab_steps.py keeps the arithmetic honest, so the
# construct still lands on exactly its 2e8-step target).
#
# Because resuming only needs the checkpoint, a yielded run is not stuck on
# this GPU -- picking it up later on GPU 0, or on the other opportunistic GPU,
# is just another launch.
#
# USAGE
#   GPU=1 nohup ./run_slab_opportunistic.sh homotypic > logs/opp_homotypic.log 2>&1 &
#   GPU=2 nohup ./run_slab_opportunistic.sh rna       > logs/opp_rna.log       2>&1 &
#
#   LEG_STEPS=10000000   smaller legs -> yields faster, slightly more overhead
#   POLL=30              seconds between occupancy checks
#   ONLY=a,b             restrict to these constructs
#
# Stop it with: kill <pid of this script>   (the current leg checkpoints out)
set +u
SLAB=$(cd "$(dirname "$0")" && pwd)
CAL_ENV=/home/sbali/miniconda3/envs/calvados
ARM=${1:?usage: [GPU=n] ./run_slab_opportunistic.sh <homotypic|rna|benchmark>}
GPU=${GPU:-1}
LEG_STEPS=${LEG_STEPS:-20000000}
POLL=${POLL:-30}
mkdir -p "$SLAB/logs"

ARM_DIR="$SLAB/$ARM"
[ -d "$ARM_DIR" ] || { echo "no such arm: $ARM_DIR"; exit 1; }

if ! "$CAL_ENV/bin/python" -c 'import openmm; openmm.Platform.getPlatformByName("CUDA")' 2>/dev/null; then
  echo "FATAL: $CAL_ENV has no CUDA platform in OpenMM."
  exit 1
fi

# PIDs on $GPU that are not this script's descendants.
foreign_pids () {
  local uuid pids p
  uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader \
         | awk -F', *' -v g="$GPU" '$1==g {print $2}')
  [ -n "$uuid" ] || return 0
  pids=$(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader \
         | awk -F', *' -v u="$uuid" '$1==u {print $2}')
  for p in $pids; do
    # ours if it is the leg we launched
    [ -n "$RUN_PID" ] && [ "$p" = "$RUN_PID" ] && continue
    echo "$p"
  done
}

wait_for_empty () {
  local first=1
  while true; do
    local f; f=$(foreign_pids)
    [ -z "$f" ] && { [ $first -eq 0 ] && echo "[opp $(date '+%F %T')] GPU $GPU is free again"; return; }
    if [ $first -eq 1 ]; then
      echo "[opp $(date '+%F %T')] GPU $GPU busy (pids: $(echo $f | tr '\n' ' ')) -- waiting"
      first=0
    fi
    sleep "$POLL"
  done
}

echo "[opp $(date '+%F %T')] arm=$ARM gpu=$GPU leg=$LEG_STEPS host=$(hostname)"
echo "[opp $(date '+%F %T')] policy: run only while GPU $GPU is empty; yield on any foreign process"

# Outer pass loop. A single walk of the arm is not enough: constructs are
# seeded progressively on GPU 0, so on an early pass most of them have no
# restart.chk yet and are skipped. Without this loop the runner would walk the
# list once, find nothing it is allowed to touch, and exit -- leaving an idle
# GPU while work appears minutes later. Instead it keeps re-walking until every
# construct in the arm is finished, sleeping when a whole pass found nothing.
while true; do
did_work=0
for d in $(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort); do
  name=$(basename "$d")
  if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then continue; fi

  # One launcher per construct directory. The whole point of this script is
  # that it shares an arm with the GPU 0 queue, so the collision it guards
  # against -- two processes writing the same restart.chk and DCD -- is the
  # likely case, not the exotic one.
  exec 9>"$d/.slab.lock"
  if ! flock -n 9; then
    echo "[opp $(date '+%F %T')] $name already running under another launcher -- skipping"
    exec 9>&-
    continue
  fi

  # A fresh construct must NOT be started here. Equilibration is a single
  # uninterruptible `simulation.step(steps_eq)` with no checkpointing (sim.py),
  # so a preemption anywhere in those 5e6 steps throws away the whole ~17 min
  # and starts over. On a GPU that is interrupted more often than that the
  # construct would livelock, never reaching production. Seed it on GPU 0
  # instead: once `restart.chk` exists, equilibration is done and forever
  # skipped (sim.py forces slab_eq off on checkpoint restart), and everything
  # after that is checkpointed every LEG_STEPS/10 and safe to preempt.
  if [ ! -f "$d/restart.chk" ] && [ "$ALLOW_FRESH" != "1" ]; then
    echo "[opp $(date '+%F %T')] SKIP $name -- not equilibrated yet (no restart.chk)."
    echo "    Seed it on GPU 0 first:  GPU=0 ONLY=$name ./run_slab_queue.sh $ARM"
    echo "    (ALLOW_FRESH=1 overrides, but risks losing equilibration repeatedly)"
    exec 9>&-
    continue
  fi

  stalled=0
  while true; do
    # Exact remaining-step arithmetic, capped to one leg.
    RUN_PID=''
    "$CAL_ENV/bin/python" "$SLAB/slab_steps.py" prepare "$d" --max-leg "$LEG_STEPS"
    rc=$?
    [ $rc -eq 3 ] && break          # construct complete -> next construct
    [ $rc -ne 0 ] && { echo "  ^ slab_steps failed for $name"; break; }
    read before _ _ < <("$CAL_ENV/bin/python" "$SLAB/slab_steps.py" status "$d")

    wait_for_empty

    echo "[opp $(date '+%F %T')] START $name leg on GPU $GPU"
    start=$(date +%s)
    ( cd "$d" && CUDA_VISIBLE_DEVICES=$GPU exec "$CAL_ENV/bin/python" -u run.py ) \
        >> "$SLAB/logs/${ARM}_${name}.log" 2>&1 &
    RUN_PID=$!

    # Watch for anyone else landing on this GPU; yield immediately if so.
    yielded=0
    while kill -0 "$RUN_PID" 2>/dev/null; do
      sleep "$POLL"
      f=$(foreign_pids)
      if [ -n "$f" ]; then
        echo "[opp $(date '+%F %T')] YIELD -- foreign pid(s) $(echo $f | tr '\n' ' ') on GPU $GPU; stopping $name"
        kill -TERM "$RUN_PID" 2>/dev/null
        # give OpenMM a moment to unwind; the last checkpoint is already on disk
        for _ in $(seq 20); do kill -0 "$RUN_PID" 2>/dev/null || break; sleep 1; done
        kill -KILL "$RUN_PID" 2>/dev/null
        yielded=1
        break
      fi
    done
    wait "$RUN_PID" 2>/dev/null
    rc=$?
    RUN_PID=''
    read done target remaining < <("$CAL_ENV/bin/python" "$SLAB/slab_steps.py" status "$d")
    echo "[opp $(date '+%F %T')] END $name exit=$rc elapsed=$(( $(date +%s) - start ))s progress=${done}/${target}"

    if [ $yielded -eq 1 ]; then
      # If we keep getting preempted before a single checkpoint lands, this GPU
      # is churning faster than LEG_STEPS/10 and we are making no progress at
      # all. Say so rather than spinning silently forever.
      if [ "$done" = "$before" ]; then
        stalled=$((stalled + 1))
        if [ $stalled -ge 3 ]; then
          echo "[opp $(date '+%F %T')] NO PROGRESS on $name after $stalled yields --"
          echo "    GPU $GPU is being reclaimed faster than the $((LEG_STEPS / 10))-step"
          echo "    checkpoint interval. Lower LEG_STEPS, or run this construct on GPU 0."
          stalled=0
        fi
      else
        stalled=0
      fi
      wait_for_empty            # sit out until they are done, then continue
      continue
    fi
    if [ $rc -ne 0 ]; then
      echo "  ^ FAILED (not a yield) -- see logs/${ARM}_${name}.log. Next construct."
      break
    fi
    [ "$remaining" = "0" ] && break
  done
  did_work=1
  exec 9>&-    # release the per-construct lock
done

# Is anything in this arm still unfinished? If not, we are genuinely done.
pending=0
for d in $(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort); do
  [ -f "$d/slab_meta.yaml" ] || continue
  read _ _ rem < <("$CAL_ENV/bin/python" "$SLAB/slab_steps.py" status "$d")
  [ "$rem" != "0" ] && pending=$((pending + 1))
done
[ $pending -eq 0 ] && break

if [ $did_work -eq 0 ]; then
  echo "[opp $(date '+%F %T')] nothing runnable yet ($pending construct(s) still awaiting a seed on GPU 0) -- re-checking in 5 min"
  sleep 300
fi
done

echo "[opp $(date '+%F %T')] arm $ARM complete on GPU $GPU"
