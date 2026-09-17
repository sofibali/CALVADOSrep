#!/bin/bash
# Detached 6-hourly progress logger for the predict_initial_guess queue.
# Survives session end; appends one timestamped line per check to logs/queue_progress.log
set +u
CD=/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md
LOG=$CD/logs/queue_progress.log
while true; do
  CUR=$(pgrep -af "bindcraft[.]py --settings" | sed 's|.*/queue/||; s|\.json.*||' | head -1)
  LINE=""
  for d in "$CD"/designs/guess_*/; do
    [ -d "$d" ] || continue
    t=$(basename "$d" | sed 's|^guess_||')
    a=$(ls "$d"Accepted/*.pdb 2>/dev/null | wc -l)
    r=$(ls "$d"Trajectory/Relaxed 2>/dev/null | wc -l)
    s=$(( $(wc -l < "$d"mpnn_design_stats.csv 2>/dev/null || echo 1) - 1 ))
    LINE="$LINE ${t}:${a}acc/${s}scored/${r}traj"
  done
  if [ -z "$CUR" ]; then
    echo "$(date '+%F %T') IDLE --$LINE" >> "$LOG"
    pgrep -f "[r]un_queue.sh" >/dev/null || { echo "$(date '+%F %T') QUEUE PROCESS GONE - logger exiting" >> "$LOG"; exit 0; }
  else
    echo "$(date '+%F %T') running=$CUR --$LINE" >> "$LOG"
  fi
  sleep 21600
done
