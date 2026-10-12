#!/bin/bash
# BC2 PARP14 campaigns, serial on one GPU.
# Order is experimental logic, not preference:
#   1 md1_baseline   - does BC2 work at all? head-to-head vs BC1 (9 accepted / 21 traj)
#   2 md1_selective  - does detargeting work on the target we KNOW is designable?
#   3 md3_selective  - the hard case: BC1 gave 0 scored from 50 trajectories
#   4 md2_selective  - BC1 gave 16 scored, 0 accepted
set +u
cd /home/sbali/BindCraft2
source .venv/bin/activate
GPU=${GPU:-0}
run () {
  local c=$1
  echo "[bc2 $(date '+%F %T')] START $c"
  CUDA_VISIBLE_DEVICES=$GPU bindcraft design "parp14_campaigns/${c}.json" \
    >> "logs/${c}.log" 2>&1
  echo "[bc2 $(date '+%F %T')] END $c exit=$?"
}
run md1_baseline
run md1_selective
run md3_selective
run md2_selective
echo "[bc2 $(date '+%F %T')] BC2 CHAIN COMPLETE"
