#!/bin/bash
# One-glance campaign status.
cd "$(dirname "${BASH_SOURCE[0]}")"
echo "=== daemons ==="
pgrep -f run_campaign.sh     >/dev/null && echo "  campaign driver  RUNNING" || echo "  campaign driver  STOPPED"
pgrep -f process_completed.sh >/dev/null && echo "  analysis watcher RUNNING" || echo "  analysis watcher STOPPED"
echo
echo "=== sets (5 replicates each) ==="
tot=0
for d in binding/*/; do
    s=$(basename "$d"); n=$(find "$d" -name '*.dcd' | wc -l); tot=$((tot+n))
    m=""; [ -e "logs/$s.complete" ] && m="run-done"; [ -e "logs/$s.analyzed" ] && m="ANALYSED"
    printf "  %-20s %2d/5  %s\n" "$s" "$n" "$m"
done
echo "  ------------------------------------"
printf "  %-20s %3d/60 replicates started\n" TOTAL "$tot"
echo
echo "=== live replicates ==="
for p in $(pgrep -f "bin/python run.py"); do
    gpu=$(tr '\0' '\n' < /proc/$p/environ 2>/dev/null | grep CUDA_VISIBLE_DEVICES | cut -d= -f2)
    cwd=$(readlink /proc/$p/cwd 2>/dev/null | sed 's|.*/binding/||')
    [ -n "$cwd" ] && [ -n "$gpu" ] && printf "  GPU %s  %s\n" "$gpu" "$cwd"
done
echo
echo "=== GPUs ==="
nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader | sed 's/^/  /'
[ -f analysis/dissociation_table.csv ] && { echo; echo "=== latest results ==="; column -s, -t analysis/dissociation_table.csv | sed 's/^/  /'; }
