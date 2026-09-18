#!/bin/bash
# PARP14/PARP9/DTX3L binding campaign driver — castor.fraserlab.com
# 240 replicates, 12 sets, 4 GPUs, 1 replicate per GPU.
#
# PER_GPU=1 is deliberate and measured: castor has CUDA MPS off and compute
# mode Default, so packing 2-4 sims per GPU costs ~55% of aggregate throughput
# (ternary: 1.98e4 ns/day at 1/GPU vs 8.1e3 aggregate at 2/GPU).
#
# Set order: the p9_dtx3l gate pair first (separated vs docked on the one
# interface that reproduced across molecular context), then ascending bead
# count so cheap results land early.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHON_EXE=/home/sbali/miniconda3/envs/calvados/bin/python
export GPU_LIST="0 1 2 3"
PER_GPU=1

SETS="p9_dtx3l p9_dtx3l_docked dtx3l_homo dtx3l_homo_docked p9_homo \
p14_dtx3l p14_dtx3l_docked p14_p9 p14_p9_docked ternary ternary_docked p14_homo"

mkdir -p logs
for s in $SETS; do
    echo "=== $s : START $(date -Is) ==="
    bash binding/$s/run.sh "$PER_GPU" 2>&1 | sed "s/^/[$s] /"
    ndcd=$(find binding/$s -name '*.dcd' | wc -l)
    echo "=== $s : DONE $(date -Is)  ($ndcd/20 trajectories) ==="
    touch logs/${s}.complete
done
echo "=== CAMPAIGN COMPLETE $(date -Is) ==="
