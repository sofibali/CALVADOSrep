#!/bin/bash
# Wait until a GPU is genuinely free (no compute processes from anyone), then run Arm 1.
# Polls every 10 min. Logs each check so the wait is auditable.
set +u
cd /home/sbali/BindCraft2
LOG=logs/arm1_waiter.log
CAMPAIGN=md1_focused

free_gpu () {
  # print the index of the first GPU with zero compute processes, else nothing
  for i in $(nvidia-smi --query-gpu=index --format=csv,noheader); do
    uuid=$(nvidia-smi --query-gpu=uuid --format=csv,noheader -i "$i")
    n=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader | grep -c "$uuid")
    if [ "$n" -eq 0 ]; then echo "$i"; return; fi
  done
}

echo "[waiter $(date '+%F %T')] waiting for a free GPU before starting $CAMPAIGN" >> "$LOG"
while true; do
  G=$(free_gpu)
  if [ -n "$G" ]; then
    echo "[waiter $(date '+%F %T')] GPU $G is free -> starting $CAMPAIGN" >> "$LOG"
    source .venv/bin/activate
    CUDA_VISIBLE_DEVICES=$G bindcraft design "parp14_campaigns/${CAMPAIGN}.json" \
      >> "logs/${CAMPAIGN}.log" 2>&1
    echo "[waiter $(date '+%F %T')] $CAMPAIGN finished exit=$?" >> "$LOG"
    exit 0
  fi
  echo "[waiter $(date '+%F %T')] all GPUs busy: $(nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader | tr '\n' ' | ')" >> "$LOG"
  sleep 600
done
