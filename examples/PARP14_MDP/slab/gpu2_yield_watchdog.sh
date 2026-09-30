#!/bin/bash
# Watch GPU2. If other users show up and stay, pause the continuous rnapanel
# queue AT A CHECKPOINT and hand the card over to the opportunistic runner.
#
#   nohup ./gpu2_yield_watchdog.sh > logs/gpu2_watchdog.log 2>&1 &
#
#   POLL=60      seconds between occupancy checks
#   PINGS=3      consecutive busy polls before switching (3 x 60s = 3 min)
#   MAX_WAIT=3600  give up waiting for a checkpoint after this and stop anyway
#
# WHY IT WAITS INSTEAD OF KILLING IMMEDIATELY
# -------------------------------------------
# sim.py writes restart.chk via saveCheckpoint() after each batch (sim.py:642),
# and NOT AT ALL during the 5e6-step slab equilibration. So the cost of a
# SIGTERM depends entirely on where the run is:
#
#   mid-production     <= one checkpoint_interval = 100k steps ~ 20 s. Trivial.
#   mid-equilibration  the whole equilibration so far, up to ~17 min, and the
#                      run then has to redo it from scratch.
#
# So this does not kill on detection. It arms, then waits for restart.chk to be
# written -- a fresh mtime if production is already under way, or the file
# appearing at all if the run is still equilibrating -- and terminates in the
# seconds right after. Losing ~0 either way, and never throwing away an
# equilibration that is nearly done.
#
# Deference cost of that wait is bounded: production checkpoints every ~20 s,
# so in the common case the handover takes well under a minute.
set +u
SLAB=$(cd "$(dirname "$0")" && pwd)
POLL=${POLL:-60}
PINGS=${PINGS:-3}
MAX_WAIT=${MAX_WAIT:-3600}
GPU=2
log () { echo "[watchdog $(date '+%F %T')] $*"; }

# PIDs on GPU2 that are not ours. Ours = descendants of the driver, which is
# every process the campaign launched on this card.
foreign_pids () {
  local uuid pids p root
  uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader 2>/dev/null \
         | awk -F', *' -v g="$GPU" '$1==g {print $2}')
  # Fail CLOSED, but toward yielding: if we cannot see the GPU we assume
  # someone is on it. Switching to opportunistic is the deferential direction,
  # so an unreadable nvidia-smi should push us that way, not paper over it.
  [ -z "$uuid" ] && { echo "nvidia-smi-unreadable"; return 0; }
  pids=$(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader 2>/dev/null) \
      || { echo "nvidia-smi-unreadable"; return 0; }
  for p in $(echo "$pids" | awk -F', *' -v u="$uuid" '$1==u {print $2}'); do
    root=$p
    while [ -n "$root" ] && [ "$root" != "1" ] && [ "$root" != "$DRIVER" ]; do
      root=$(ps -o ppid= -p "$root" 2>/dev/null | tr -d ' ')
    done
    [ "$root" = "$DRIVER" ] && continue     # ours
    echo "$p"
  done
}

DRIVER=$(pgrep -f 'bash .*run_rnapanel_gpu2.sh' | head -1)
[ -z "$DRIVER" ] && { log "no run_rnapanel_gpu2.sh running -- nothing to watch"; exit 1; }
log "watching GPU$GPU; driver pid $DRIVER; switch after $PINGS consecutive busy polls of ${POLL}s"

busy=0
while kill -0 "$DRIVER" 2>/dev/null; do
  f=$(foreign_pids)
  if [ -n "$f" ]; then
    busy=$((busy+1))
    log "GPU$GPU busy ($(echo $f | tr '\n' ' ')) -- ping $busy/$PINGS"
  else
    [ "$busy" -gt 0 ] && log "GPU$GPU clear again -- resetting ping count"
    busy=0
  fi
  [ "$busy" -ge "$PINGS" ] && break
  sleep "$POLL"
done

if ! kill -0 "$DRIVER" 2>/dev/null; then
  log "driver exited on its own -- panel finished or was stopped; not switching"
  exit 0
fi

log "SUSTAINED CONTENTION -- arming a checkpoint-aligned pause"

# The leg currently on the GPU, and the construct directory it is running in.
# Identify it the same way ownership is decided: the process ON GPU2 whose
# ancestry leads back to the driver. Matching on pgrep ordering instead would
# be fragile -- there are two nested run_slab_queue.sh shells.
RUNPID=''
_uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader 2>/dev/null \
        | awk -F', *' -v g="$GPU" '$1==g {print $2}')
for p in $(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader 2>/dev/null \
           | awk -F', *' -v u="$_uuid" '$1==u {print $2}'); do
  root=$p
  while [ -n "$root" ] && [ "$root" != "1" ] && [ "$root" != "$DRIVER" ]; do
    root=$(ps -o ppid= -p "$root" 2>/dev/null | tr -d ' ')
  done
  [ "$root" = "$DRIVER" ] && { RUNPID=$p; break; }
done
if [ -n "$RUNPID" ]; then
  DIR=$(readlink -f "/proc/$RUNPID/cwd" 2>/dev/null)
  CHK="$DIR/restart.chk"
  before=$(stat -c %Y "$CHK" 2>/dev/null || echo none)
  if [ "$before" = "none" ]; then
    log "leg $(basename "$DIR") (pid $RUNPID) is still EQUILIBRATING -- waiting for its first checkpoint so the equilibration is not thrown away"
  else
    log "leg $(basename "$DIR") (pid $RUNPID) is in production -- waiting for the next checkpoint (~20 s)"
  fi
  waited=0
  while kill -0 "$RUNPID" 2>/dev/null && [ "$waited" -lt "$MAX_WAIT" ]; do
    now=$(stat -c %Y "$CHK" 2>/dev/null || echo none)
    [ "$now" != "$before" ] && [ "$now" != "none" ] && { log "checkpoint written -- pausing now (0 steps lost)"; break; }
    sleep 5; waited=$((waited+5))
  done
  [ "$waited" -ge "$MAX_WAIT" ] && log "WARNING: no checkpoint after ${MAX_WAIT}s -- stopping anyway, this leg loses its progress since the last one"
fi

# Stop the shells first so the queue cannot start the next construct, then the
# leg itself. Killing the queue shell releases its fd-9 flock on the construct.
log "stopping continuous queue (driver $DRIVER)"
kill -TERM "$DRIVER" 2>/dev/null
pkill -TERM -f 'run_slab_queue.sh rnapanel' 2>/dev/null
[ -n "$RUNPID" ] && kill -TERM "$RUNPID" 2>/dev/null
for _ in $(seq 30); do kill -0 "$RUNPID" 2>/dev/null || break; sleep 1; done
kill -KILL "$RUNPID" 2>/dev/null
sleep 5

log "starting opportunistic runner on GPU$GPU (yields on any foreign process)"
cd "$SLAB" || exit 1
GPU=2 setsid nohup ./run_slab_opportunistic.sh rnapanel \
    > logs/opp_rnapanel_gpu2.log 2>&1 < /dev/null &
disown
log "handed over. opportunistic pid $!; log: logs/opp_rnapanel_gpu2.log"
log "to go back to continuous: kill the opportunistic runner, then ./run_rnapanel_gpu2.sh"
