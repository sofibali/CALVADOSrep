#!/bin/bash
# Run-2 RNA panel on GPU2: seeding pass, then continuous production.
#
#   nohup ./run_rnapanel_gpu2.sh > logs/rnapanel_gpu2.log 2>&1 &
#
# Pass 1 (LEG_STEPS=200000) takes all 24 unparked runs through the 5e6-step
# slab equilibration and two checkpoints, then moves on. After it, every run
# has a restart.chk, so any of them can be picked up later by an opportunistic
# runner on another card -- the opportunistic launcher refuses constructs that
# have not equilibrated. ~20 min/run, so ~8 h for the pass.
#
# Pass 2 runs each construct to its 2e8-step target. slab_steps.py sets `steps`
# to what is actually LEFT on every entry, so re-running this script is safe
# and finished constructs are skipped (rc=3).
#
# The -WWE-ART tier is parked (.parked markers); remove those to include it.
set -u
SLAB=$(cd "$(dirname "$0")" && pwd)

echo "=== pass 1: seeding (equilibration + first checkpoints) ==="
GPU=2 LEG_STEPS=200000 "$SLAB/run_slab_queue.sh" rnapanel

echo ""
echo "=== pass 2: production to 2e8 steps ==="
while true; do
  GPU=2 "$SLAB/run_slab_queue.sh" rnapanel
  # Re-sweep: a run that crashed mid-leg is retried on the next pass rather
  # than abandoning the rest of the panel. Exit once nothing has work left.
  left=0
  for d in "$SLAB"/rnapanel/*/; do
    [ -f "$d/slab_meta.yaml" ] || continue
    [ -f "$d/.parked" ] && continue
    read -r done target _ < <("/home/sbali/miniconda3/envs/calvados/bin/python" \
        "$SLAB/slab_steps.py" status "$d" 2>/dev/null) || continue
    [ "${done:-0}" -lt "${target:-1}" ] && left=$((left+1))
  done
  echo "[rnapanel $(date '+%F %T')] runs still short of target: $left"
  [ "$left" -eq 0 ] && break
  sleep 60
done
echo "[rnapanel $(date '+%F %T')] PANEL COMPLETE"
