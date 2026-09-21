#!/bin/bash
# One-shot watcher: stop guess_clamp_md1md2_state1 at 25 relaxed trajectories.
# The running bindcraft process loaded max_trajectories=50 at startup, so the cap
# has to be applied externally. SIGINT lets it exit between steps; run_queue.sh
# then advances to the next target on its own.
set +u
CD=/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md
D=$CD/designs/guess_clamp_md1md2_state1
LIMIT=25
LOG=$CD/logs/queue_progress.log
while true; do
  PID=$(pgrep -f "bindcraft[.]py --settings .*clamp_md1md2_state1" | head -1)
  [ -z "$PID" ] && { echo "$(date '+%F %T') cap_state1: target no longer running, watcher exiting" >> "$LOG"; exit 0; }
  N=$(ls "$D"/Trajectory/Relaxed/*.pdb 2>/dev/null | wc -l)
  if [ "$N" -ge "$LIMIT" ]; then
    echo "$(date '+%F %T') cap_state1: reached $N relaxed >= $LIMIT, sending SIGINT to $PID" >> "$LOG"
    kill -SIGINT "$PID"
    exit 0
  fi
  sleep 600
done
