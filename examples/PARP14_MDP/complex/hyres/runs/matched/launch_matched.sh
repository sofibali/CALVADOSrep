#!/bin/bash
# 2 arms x 3 replicates, one per GPU. ~3.5 h each, so ~3.5 h total on 6 GPUs.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
P=/home/sbali/miniconda3/envs/hyres/bin/python
i=0
for arm in separated docked; do
  for rep in 1 2 3; do
    d=${arm}_rep${rep}; mkdir -p $d
    cp complex.psf top_hyres_GPU.inp param_hyres_GPU.inp run_matched.py _ff_body.py $d/
    cp start_${arm}.pdb $d/start.pdb
    cp domain_restraints_${arm}.txt $d/restraints.txt
    ( cd $d && CUDA_VISIBLE_DEVICES=$i nohup $P run_matched.py start.pdb complex.psf \
        restraints.txt 0 $((1000+rep)) > run.log 2>&1 ) &
    echo "  launched $d on GPU $i (seed $((1000+rep)))"
    i=$((i+1))
  done
done
wait
echo "=== all matched HyRes runs finished $(date -Is) ==="
