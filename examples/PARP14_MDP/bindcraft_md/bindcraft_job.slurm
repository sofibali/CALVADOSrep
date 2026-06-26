#!/bin/bash
#SBATCH --job-name=bc_parp14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=42gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/bc_%x_%j.log
#
# One BindCraft run for a single PARP14 macrodomain target.
# Usage (via submit_bindcraft.sh, or directly):
#   sbatch bindcraft_job.slurm -s settings/md1_block_af3.json
#
# IMPORTANT: set BC_ENV to the env created by install_bindcraft.sh.
# On this system the env is in ~/.conda/envs (NOT $CONDA_BASE/envs), so the
# repo's stock bindcraft.slurm path is wrong -- we activate explicitly here.

BC_REPO=/home/sbali/Projects/BindCraft
BC_ENV=/home/sbali/.conda/envs/BindCraft          # <-- adjust if install put it elsewhere
ADV=${BC_REPO}/settings_advanced/default_4stage_multimer.json
FIL=${BC_REPO}/settings_filters/default_filters.json

# --- parse -s/--settings (also accept -f/-a overrides) ---
SETTINGS=""
TEMP=$(getopt -o s:f:a: --long settings:,filters:,advanced: -n 'bindcraft_job' -- "$@")
eval set -- "$TEMP"
while true; do
  case "$1" in
    -s|--settings) SETTINGS="$2"; shift 2 ;;
    -f|--filters)  FIL="$2"; shift 2 ;;
    -a|--advanced) ADV="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done
[ -z "$SETTINGS" ] && { echo "Error: -s <settings.json> required"; exit 1; }

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate ${BC_ENV}
export LD_LIBRARY_PATH=${BC_ENV}/lib:${LD_LIBRARY_PATH}

echo "Running BindCraft: settings=${SETTINGS}"
echo "  advanced=${ADV}"
echo "  filters =${FIL}"
python -u "${BC_REPO}/bindcraft.py" --settings "${SETTINGS}" --filters "${FIL}" --advanced "${ADV}"
