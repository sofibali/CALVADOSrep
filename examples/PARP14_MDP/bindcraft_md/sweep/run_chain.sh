#!/bin/bash
# Serial sweep on GPU 0: wait for the replay to finish, then Arm A, then Arm C.
# Each arm self-terminates at max_trajectories=100 (relaxed trajectories).
set +u   # conda activate scripts reference unbound vars
BC=/home/sbali/BindCraft
CD=/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md
SW=$CD/sweep
LOG=$CD/logs
REPLAY_PID=${1:-0}

source /home/sbali/miniconda3/bin/activate BindCraft
export LD_LIBRARY_PATH=/home/sbali/miniconda3/envs/BindCraft/lib:${LD_LIBRARY_PATH:-}
cd "$CD"

echo "[chain $(date '+%F %T')] waiting for replay PID $REPLAY_PID to finish"
while kill -0 "$REPLAY_PID" 2>/dev/null; do sleep 120; done
echo "[chain $(date '+%F %T')] replay done"

run_arm () {
  local name=$1 adv=$2
  echo "[chain $(date '+%F %T')] START $name"
  CUDA_VISIBLE_DEVICES=0 python -u "$BC/bindcraft.py" \
    --settings "$SW/settings/af3_${name}.json" \
    --filters  "$BC/settings_filters/default_filters.json" \
    --advanced "$adv" >> "$LOG/sweep_af3_${name}.log" 2>&1
  echo "[chain $(date '+%F %T')] END $name (exit $?) relaxed=$(ls $CD/designs/sweep_af3_${name}/Trajectory/Relaxed 2>/dev/null | wc -l)"
}

run_arm mpnnorig "$SW/advanced/mpnnorig_cap100.json"
run_arm guess    "$SW/advanced/guess_cap100.json"
echo "[chain $(date '+%F %T')] CHAIN COMPLETE"
