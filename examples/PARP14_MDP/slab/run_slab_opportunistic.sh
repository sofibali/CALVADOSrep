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
ARM=${1:?usage: [GPU=n] ./run_slab_opportunistic.sh <arm> [wait_pid]}
WAIT_PID=${2:-0}
GPU=${GPU:-1}
LEG_STEPS=${LEG_STEPS:-10000000}  # 1e7: ~3% restart overhead, 10 dup frames/yield. See README "Yielding leaves overlapping frames".
POLL=${POLL:-30}
mkdir -p "$SLAB/logs"

ARM_DIR="$SLAB/$ARM"
[ -d "$ARM_DIR" ] || { echo "no such arm: $ARM_DIR"; exit 1; }

if ! "$CAL_ENV/bin/python" -c 'import openmm; openmm.Platform.getPlatformByName("CUDA")' 2>/dev/null; then
  echo "FATAL: $CAL_ENV has no CUDA platform in OpenMM."
  exit 1
fi

# PIDs on $GPU that are not this script's descendants.
# NOTE: this FAILS CLOSED. If nvidia-smi cannot be read we report a sentinel
# "unknown" rather than "nobody is there" -- the whole point of this script is
# deference, and a transient nvidia-smi failure must not be indistinguishable
# from an idle GPU or we would keep running straight through someone else's job.
foreign_pids () {
  local uuid pids p
  uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader 2>/dev/null \
         | awk -F', *' -v g="$GPU" '$1==g {print $2}')
  if [ -z "$uuid" ]; then echo "nvidia-smi-unreadable"; return 0; fi
  pids=$(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader 2>/dev/null) \
      || { echo "nvidia-smi-unreadable"; return 0; }
  pids=$(echo "$pids" | awk -F', *' -v u="$uuid" '$1==u {print $2}')
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

# Terminate the current leg cleanly if this script is killed. Without this the
# backgrounded run.py is orphaned: nothing polls the GPU for it any more (so it
# never yields, defeating the whole point), and it keeps the inherited fd-9
# flock, so every other launcher then skips that construct forever.
RUN_PID=''
cleanup () {
  trap - TERM INT EXIT
  if [ -n "$RUN_PID" ] && kill -0 "$RUN_PID" 2>/dev/null; then
    echo "[opp $(date '+%F %T')] stopping -- terminating current leg (pid $RUN_PID)"
    kill -TERM "$RUN_PID" 2>/dev/null
    for _ in $(seq 20); do kill -0 "$RUN_PID" 2>/dev/null || break; sleep 1; done
    kill -KILL "$RUN_PID" 2>/dev/null
  fi
  exec 9>&- 2>/dev/null
  exit 143
}
trap cleanup TERM INT

# Chain behind another launcher, so a follow-on arm (e.g. the FL
# concentration ladder) starts only once the main panel is finished.
if [ "$WAIT_PID" != "0" ]; then
  echo "[opp $(date '+%F %T')] waiting for PID $WAIT_PID before starting arm $ARM"
  while kill -0 "$WAIT_PID" 2>/dev/null; do sleep 120; done
  echo "[opp $(date '+%F %T')] PID $WAIT_PID finished -- starting"
fi

echo "[opp $(date '+%F %T')] arm=$ARM gpu=$GPU leg=$LEG_STEPS host=$(hostname)"
echo "[opp $(date '+%F %T')] policy: run only while GPU $GPU is empty; yield on any foreign process"

# Outer pass loop. A single walk of the arm is not enough: constructs are seeded
# progressively on GPU 0, so on an early pass most have no restart.chk yet and
# are skipped. Without this the runner would walk the list once, find nothing it
# may touch, and exit -- leaving an idle GPU while work appears minutes later.
while true; do
did_work=0
for d in $(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort); do
  name=$(basename "$d")
  if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then continue; fi
  [ -f "$d/slab_meta.yaml" ] || continue
  # A construct can be taken out of the queue without deleting or moving
  # it (moving would break components.yaml's absolute paths). The marker
  # file says why it was parked and how to un-park it.
  if [ -f "$d/.parked" ]; then
    echo "[opp $(date '+%F %T')] PARKED $name -- skipping (see $d/.parked)"
    continue
  fi    # not a construct directory

  # A fresh construct must NOT be started here. Equilibration is a single
  # uninterruptible `simulation.step(steps_eq)` with NO checkpointing (sim.py),
  # so a preemption anywhere in those 5e6 steps throws away the whole ~17 min
  # and starts over; on a GPU reclaimed more often than that the construct
  # livelocks and never reaches production. Seed it on GPU 0 instead: once
  # restart.chk exists, slab_eq is forced off on restart and everything after
  # is checkpointed every LEG_STEPS/10 and safe to preempt.
  if [ ! -f "$d/restart.chk" ] && [ "$ALLOW_FRESH" != "1" ]; then
    echo "[opp $(date '+%F %T')] SKIP $name -- not equilibrated yet (no restart.chk)."
    echo "    Seed it on GPU 0 first:  GPU=0 ONLY=$name ./run_slab_queue.sh $ARM"
    echo "    (ALLOW_FRESH=1 overrides, but risks losing equilibration repeatedly)"
    continue
  fi

  stalled=0
  while true; do
    # Wait for the GPU BEFORE taking the lock. Holding the per-construct lock
    # across an unbounded wait would pin the construct to this GPU: if another
    # user sat on GPU $GPU for two days, run_slab_queue.sh on GPU 0 would find
    # the lock held and skip the construct, which is the opposite of the
    # "a yielded run is movable" property this script is supposed to provide.
    wait_for_empty

    exec 9>"$d/.slab.lock"
    if ! flock -n 9; then
      echo "[opp $(date '+%F %T')] $name taken by another launcher -- moving on"
      exec 9>&-
      break
    fi

    # Exact remaining-step arithmetic, capped to one leg.
    RUN_PID=''
    "$CAL_ENV/bin/python" "$SLAB/slab_steps.py" prepare "$d" --max-leg "$LEG_STEPS"
    rc=$?
    if [ $rc -eq 3 ]; then exec 9>&-; break; fi          # complete
    if [ $rc -ne 0 ]; then
      echo "  ^ slab_steps.py failed (rc=$rc) for $name -- skipping rather than"
      echo "    running with an unpatched steps value and overshooting the target"
      exec 9>&-; break
    fi
    read before _ _ < <("$CAL_ENV/bin/python" "$SLAB/slab_steps.py" status "$d")

    # A previous yield may have left frames in the DCD that the checkpoint does
    # not back -- the resumed leg would re-simulate that window and append a
    # second trajectory through it. Truncate back to the checkpoint first so
    # this leg appends onto a clean boundary. No-op when there is nothing to
    # drop, and O(1) regardless of trajectory size.
    "$CAL_ENV/bin/python" "$SLAB/trim_dcd.py" "$d" --apply

    # Somebody may have landed while we were taking the lock.
    if [ -n "$(foreign_pids)" ]; then exec 9>&-; continue; fi

    echo "[opp $(date '+%F %T')] START $name leg on GPU $GPU"
    start=$(date +%s)
    ( cd "$d" && CUDA_VISIBLE_DEVICES=$GPU exec "$CAL_ENV/bin/python" -u run.py ) \
        >> "$SLAB/logs/${ARM}_${name}.log" 2>&1 &
    RUN_PID=$!

    yielded=0
    while kill -0 "$RUN_PID" 2>/dev/null; do
      sleep "$POLL"
      f=$(foreign_pids)
      [ -z "$f" ] && continue
      # The run may have finished during that sleep; killing a corpse and
      # calling it a yield would hide a real failure and retry it forever.
      kill -0 "$RUN_PID" 2>/dev/null || break
      echo "[opp $(date '+%F %T')] YIELD -- foreign pid(s) $(echo $f | tr '\n' ' ') on GPU $GPU; stopping $name"
      kill -TERM "$RUN_PID" 2>/dev/null
      for _ in $(seq 20); do kill -0 "$RUN_PID" 2>/dev/null || break; sleep 1; done
      kill -KILL "$RUN_PID" 2>/dev/null
      yielded=1
      break
    done
    wait "$RUN_PID" 2>/dev/null
    rc=$?
    RUN_PID=''
    did_work=1
    read done target remaining < <("$CAL_ENV/bin/python" "$SLAB/slab_steps.py" status "$d")
    echo "[opp $(date '+%F %T')] END $name exit=$rc elapsed=$(( $(date +%s) - start ))s progress=${done}/${target}"
    exec 9>&-    # release before any further waiting

    # rc=0 means it exited cleanly on its own, whatever the GPU looked like.
    if [ $yielded -eq 1 ] && [ $rc -ne 0 ]; then
      # Preempted before a checkpoint landed => no progress at all. If that
      # keeps happening this GPU is churning faster than the checkpoint
      # interval and we are burning setup cost for nothing.
      if [ "$done" = "$before" ]; then
        stalled=$((stalled + 1))
        if [ $stalled -ge 3 ]; then
          echo "[opp $(date '+%F %T')] NO PROGRESS on $name after $stalled yields --"
          echo "    GPU $GPU is reclaimed faster than the $((LEG_STEPS / 10))-step"
          echo "    checkpoint interval. Lower LEG_STEPS, or run this one on GPU 0."
          stalled=0
        fi
      else
        stalled=0
      fi
      continue
    fi
    if [ $rc -ne 0 ]; then
      echo "  ^ FAILED (not a yield) -- see logs/${ARM}_${name}.log. Next construct."
      break
    fi
    [ "$remaining" = "0" ] && break
  done
done

# Anything left unfinished in this arm?
pending=0
for d in $(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort); do
  [ -f "$d/slab_meta.yaml" ] || continue
  read _ _ rem < <("$CAL_ENV/bin/python" "$SLAB/slab_steps.py" status "$d")
  [ "$rem" != "0" ] && pending=$((pending + 1))
done
[ $pending -eq 0 ] && break

if [ $did_work -eq 0 ]; then
  echo "[opp $(date '+%F %T')] nothing runnable yet ($pending construct(s) awaiting a seed on GPU 0) -- re-checking in 5 min"
  sleep 300
fi
done

echo "[opp $(date '+%F %T')] arm $ARM complete on GPU $GPU"
