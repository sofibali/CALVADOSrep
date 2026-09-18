#!/bin/bash
#SBATCH --job-name=parp14_slab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/slab_%x_%A_%a.log
#
# One PARP14 slab coexistence run per array task.
#
#   sbatch --array=0-5 submit_slab.slurm slab/homotypic
#
# ARM_DIR (positional) holds one subdirectory per construct; the array index
# selects one. Partition/qos/gres mirror bindcraft_job.slurm, which is the only
# verified-working GPU submission in this repo -- adjust if the GPU server uses
# different names.
#
# CALVADOS restarts from restart.chk run `steps` ADDITIONAL steps and append to
# the DCD, so re-submitting a timed-out task continues rather than restarting.

set -u
ARM_DIR=${1:?usage: sbatch --array=0-N submit_slab.slurm <arm_dir>}
CAL_ENV=/home/sbali/miniconda3/envs/CALVADOS

mapfile -t DIRS < <(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
D=${DIRS[$SLURM_ARRAY_TASK_ID]}
echo "host=$(hostname)  task=$SLURM_ARRAY_TASK_ID  dir=$D"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

cd "$D" || exit 1
srun "$CAL_ENV/bin/python" -u run.py
