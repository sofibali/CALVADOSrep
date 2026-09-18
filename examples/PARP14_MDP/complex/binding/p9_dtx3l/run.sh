#!/bin/bash
# Launch this system's rep-* replicates across GPUs (round-robin).
#   bash run.sh [PER_GPU]          GPU_LIST="0 1" bash run.sh 2
set -e
CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHON_EXE="${PYTHON_EXE:-/home/vle/.conda/envs/calvados/bin/python}"
read -r -a GPUS <<< "${GPU_LIST:-0 1 2 3 4 5 6 7}"
NGPU=${#GPUS[@]}; PER_GPU=${1:-3}; PARALLEL=$(( NGPU * PER_GPU ))
cd "$CWD"
i=0
for d in $(ls -d rep-* | sort -V); do echo "${GPUS[$((i % NGPU))]} $d"; i=$((i + 1)); done | xargs -P "$PARALLEL" -L1 bash -c '
    gpu="$0"; d="$1"
    cd "$d" || exit 1
    if compgen -G "*.dcd" > /dev/null; then echo "  $d: has trajectory, skipping"; exit 0; fi
    echo "  starting $d on GPU $gpu"
    CUDA_VISIBLE_DEVICES="$gpu" "$PYTHON_EXE" run.py > run.log 2>&1 || echo "  $d: FAILED"
    echo "  $d: done"
'
