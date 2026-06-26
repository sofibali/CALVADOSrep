#!/bin/bash
# Submit PARP14 macrodomain BindCraft targets to the GPU queue.
#
#   ./submit_bindcraft.sh            # submit all 7 targets
#   ./submit_bindcraft.sh blocks     # the 6 single-MD pocket blockers only
#   ./submit_bindcraft.sh md1_block_af3   # one named target
#
# Each target -> one GPU job (see bindcraft_job.slurm). Pilot = 8 final designs/target.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SET="${HERE}/settings"
mkdir -p "${HERE}/logs"

case "${1:-all}" in
  all)    list=$(ls "${SET}"/*.json) ;;
  blocks) list=$(ls "${SET}"/*_block_*.json) ;;
  clamps) list=$(ls "${SET}"/clamp_*.json) ;;
  *)      list="${SET}/$1.json" ;;
esac

for s in $list; do
  name=$(basename "$s" .json)
  echo "Submitting ${name}"
  sbatch --job-name="bc_${name}" "${HERE}/bindcraft_job.slurm" -s "$s"
done
echo "Done. Monitor with: squeue -u $USER   |   logs in ${HERE}/logs/"
