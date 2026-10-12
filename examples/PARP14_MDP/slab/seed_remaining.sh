#!/bin/bash
# Seed the rnapanel constructs that never got an equilibration, in parallel
# across all four GPUs, WITHOUT yielding part-way.
#
#   nohup ./seed_remaining.sh > logs/seed_remaining.log 2>&1 &
#
# WHY THIS DOES NOT YIELD
# -----------------------
# Everything else in this campaign defers to other users mid-run. This step
# deliberately does not, because it cannot:
#
#   * Production is chunked and resumes from restart.chk, so it survives being
#     preempted -- that is how the opportunistic runner finished 13 runs.
#   * Equilibration is 5e6 steps (~20 min) and is NOT resumable. An OpenMM
#     checkpoint is bound to its System, and the equilibration System still
#     carries the slab-centering force, so a checkpoint taken there cannot be
#     loaded into the production System (sim.py:632).
#
# Measured contention on GPU2 over the week of 2026-09-30..10-07 was ~2 min
# busy / ~5 min free, 252 cycles. No free window is ever 20 min, so a yielding
# seeder would restart the same equilibration forever and never finish. These
# 11 constructs are the whole md_full arm (the positive control -- the only
# construct that condensed in run 1) plus 5 of 6 kh1_art_full, so they are not
# droppable.
#
# Mitigation: each job is ~20 min and ~450 MB, against 23-39 GB free per card.
# It shares compute for a short window rather than occupying a card for days.
set -u
SLAB=$(cd "$(dirname "$0")" && pwd)
CAL_ENV=/home/sbali/miniconda3/envs/calvados
LEG=${LEG:-200000}
mkdir -p "$SLAB/logs"

# Whatever still has no restart.chk, split round-robin over the 4 cards.
mapfile -t TODO < <(for d in "$SLAB"/rnapanel/*/; do
  [ -f "$d/.parked" ] && continue
  [ -f "$d/slab_meta.yaml" ] || continue
  [ -f "$d/restart.chk" ] || basename "$d"
done)
if [ "${#TODO[@]}" -eq 0 ]; then
  echo "[seed $(date '+%F %T')] nothing to seed -- every unparked construct already has restart.chk"
  exit 0
fi
echo "[seed $(date '+%F %T')] seeding ${#TODO[@]} construct(s) across GPUs 0-3, leg=$LEG, no yielding"
printf '  %s\n' "${TODO[@]}"

declare -A GROUP
for i in "${!TODO[@]}"; do
  g=$(( i % 4 ))
  GROUP[$g]="${GROUP[$g]:+${GROUP[$g]},}${TODO[$i]}"
done

pids=()
for g in 0 1 2 3; do
  [ -z "${GROUP[$g]:-}" ] && continue
  echo "[seed $(date '+%F %T')] GPU$g <- ${GROUP[$g]}"
  GPU=$g ONLY="${GROUP[$g]}" LEG_STEPS=$LEG \
    "$SLAB/run_slab_queue.sh" rnapanel > "$SLAB/logs/seed_gpu${g}.log" 2>&1 &
  pids+=($!)
done

rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "[seed $(date '+%F %T')] all seeding queues finished (rc=$rc)"

left=0
for d in "$SLAB"/rnapanel/*/; do
  [ -f "$d/.parked" ] && continue
  [ -f "$d/slab_meta.yaml" ] || continue
  [ -f "$d/restart.chk" ] || { echo "  STILL UNSEEDED: $(basename "$d")"; left=$((left+1)); }
done
echo "[seed $(date '+%F %T')] unseeded remaining: $left"

# Hand the finished work to opportunistic runners. GPU2 already has one from
# the 2026-09-30 handover; add 1 and 3 so the remaining ~128 GPU-h of
# production is not funnelled through a single contended card.
for g in 1 3; do
  if pgrep -f "run_slab_opportunistic.sh rnapanel" >/dev/null && \
     pgrep -af "run_slab_opportunistic.sh rnapanel" | grep -q "gpu${g}marker"; then continue; fi
  if [ -f "$SLAB/logs/opp_rnapanel_gpu${g}.pid" ] && \
     kill -0 "$(cat "$SLAB/logs/opp_rnapanel_gpu${g}.pid")" 2>/dev/null; then
    echo "[seed $(date '+%F %T')] opportunistic runner already on GPU$g"
    continue
  fi
  GPU=$g setsid nohup "$SLAB/run_slab_opportunistic.sh" rnapanel \
      > "$SLAB/logs/opp_rnapanel_gpu${g}.log" 2>&1 < /dev/null &
  echo $! > "$SLAB/logs/opp_rnapanel_gpu${g}.pid"
  echo "[seed $(date '+%F %T')] started opportunistic runner on GPU$g (pid $(cat "$SLAB/logs/opp_rnapanel_gpu${g}.pid"))"
done
echo "[seed $(date '+%F %T')] done. GPU2's existing opportunistic runner picks up the rest."
