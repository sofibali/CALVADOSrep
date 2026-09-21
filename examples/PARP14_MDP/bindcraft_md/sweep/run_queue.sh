#!/bin/bash
# PARP14 BindCraft queue -- remaining targets with predict_initial_guess=True.
# Serial, one GPU. Each target stops at 8 accepted designs or max_trajectories=50.
#
# Order set 2026-09-21: AF3-receptor blockers first, MD3 before MD2 (MD3 is the
# smaller target and needs no unified memory, so it reads out sooner), then the
# two sim-conformer blockers, then the MD2-MD3 clamp
# states, then the remaining MD1-MD2 clamp. Block-target hotspots were
# regenerated from the corrected UniProt-referenced active sites
# (docs/NUMBERING_AUDIT.md). clamp_md1md2_state3 (complete) and _state1
# (running when this order was set) are not repeated here.
set +u
BC=/home/sbali/BindCraft
CD=/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md
SW=$CD/sweep
GPU=${GPU:-0}
WAIT_PID=${1:-0}

source /home/sbali/miniconda3/bin/activate BindCraft
export LD_LIBRARY_PATH=/home/sbali/miniconda3/envs/BindCraft/lib:${LD_LIBRARY_PATH}
cd "$CD"

if [ "$WAIT_PID" != "0" ]; then
  echo "[queue $(date '+%F %T')] waiting for PID $WAIT_PID"
  while kill -0 "$WAIT_PID" 2>/dev/null; do sleep 120; done
fi

run () {
  local t=$1 unified=${2:-0}
  local dir="$CD/designs/guess_${t}"
  echo "[queue $(date '+%F %T')] START $t (unified_mem=$unified)"
  if [ "$unified" = "1" ]; then
    # 586-res target -> ~716-res complex OOMs a 46GB L40S; spill to host RAM
    export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
  else
    unset TF_FORCE_UNIFIED_MEMORY XLA_PYTHON_CLIENT_MEM_FRACTION
  fi
  CUDA_VISIBLE_DEVICES=$GPU python -u "$BC/bindcraft.py" \
    --settings "$SW/queue/${t}.json" \
    --filters  "$BC/settings_filters/default_filters.json" \
    --advanced "$SW/advanced/guess_cap50.json" \
    >> "$CD/logs/queue_${t}.log" 2>&1
  echo "[queue $(date '+%F %T')] END $t exit=$? accepted=$(ls $dir/Accepted/*.pdb 2>/dev/null|wc -l) relaxed=$(ls $dir/Trajectory/Relaxed 2>/dev/null|wc -l)"
}

run md3_block_af3
run md2_block_af3  1
run md3_block_sim
run md2_block_sim  1
run clamp_md2md3_state3
run clamp_md2md3_state5
run clamp_md1md2_state5
echo "[queue $(date '+%F %T')] QUEUE COMPLETE"
