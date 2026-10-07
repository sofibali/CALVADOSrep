# Resuming the BC1 queue (stopped 2026-09-30)

Stopped deliberately to move to BindCraft2, not because anything failed. Everything
is intact and resumable — BindCraft skips trajectories whose PDB already exists, and
`check_accepted_designs` re-reads from disk each loop.

## Completed (do not re-run)

| target | relaxed | scored | accepted |
|---|---|---|---|
| clamp_md1md2_state3 | 50 | 197 | **5** |
| clamp_md1md2_state1 | 28 | 49 | 0 |
| md3_block_af3 | 50 | 0 | 0 |
| md2_block_af3 | 50 | 16 | 0 |
| md3_block_sim | 12 (partial) | 0 | 0 |

## Remaining, in `sweep/run_queue.sh` order

md3_block_sim (partial, 12 relaxed) · md2_block_sim · clamp_md2md3_state3 ·
clamp_md2md3_state5 · clamp_md1md2_state5

## To resume everything

    cd /home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md
    setsid nohup ./sweep/run_queue.sh 0 > logs/queue_master.log 2>&1 < /dev/null &
    setsid nohup ./sweep/progress_logger.sh > /dev/null 2>&1 < /dev/null &

`run_queue.sh` still lists md3_block_af3 and md2_block_af3 first; both are complete
and will exit immediately at their 50-trajectory cap, so it is safe to leave them.

## To resume a single target

    source /home/sbali/miniconda3/bin/activate BindCraft
    export LD_LIBRARY_PATH=/home/sbali/miniconda3/envs/BindCraft/lib:${LD_LIBRARY_PATH}
    CUDA_VISIBLE_DEVICES=<n> python -u /home/sbali/BindCraft/bindcraft.py \
      --settings sweep/queue/<target>.json \
      --filters  /home/sbali/BindCraft/settings_filters/default_filters.json \
      --advanced sweep/advanced/guess_cap50.json

## Why they were dropped

The two `_sim` blockers were judged unlikely to inform: MD3-AF3 gave 0 scored and
MD2-AF3 gave 0 accepted, and the CG-backmapped sim conformer was strictly worse than
AF3 on MD1 (35% vs 10% clash rate, 30% vs 58% relaxed yield). The two clamp_md2md3
targets DO still carry information — whether clamp designability tracks TICA state
population in a second cleft — and are worth resuming if BC2 does not supersede them.
