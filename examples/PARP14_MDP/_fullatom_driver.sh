#!/bin/bash
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
cd /home/sbali/CALVADOS/examples/PARP14_MDP
run(){ echo "##### $1"; $PY fullatom_minimize_states.py --states-dir "$1" --reference "$2" --fdomains "$3" 2>&1 | grep -vE "Warning|warn"; }
MDREF=md_full/input/ref_allatom.pdb; MDDOM=md_full/input/domains.yaml
FLREF=input/parp14.pdb; FLDOM=input/domains.yaml
run representative_frames/2026-06-25/md_full_poseiface_tica  $MDREF $MDDOM
run representative_frames/2026-06-24/md_full_ca25_tica       $MDREF $MDDOM
run representative_frames/2026-06-24/md_full_iface10_tica    $MDREF $MDDOM
run representative_frames/2026-06-25/fl_optimized_ca25_tica  $FLREF $FLDOM
run representative_frames/2026-06-25/fl_optimized_iface10_tica $FLREF $FLDOM
echo "ALL DONE"
