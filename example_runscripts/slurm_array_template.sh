#!/bin/bash
#============================================================================
# SLURM Array Job Template for CALVADOS Simulations
#============================================================================
# Usage:
#   sbatch slurm_array_template.sh
#
# For PARP14 library (75,000 jobs = 3000 compositions × 25 structures):
#   sbatch --array=1-75000%500 slurm_array_template.sh
#
# For testing (first 10 jobs):
#   sbatch --array=1-10 slurm_array_template.sh
#============================================================================

#SBATCH --job-name=calvados
#SBATCH --array=1-100%50          # 100 jobs, max 50 concurrent
#SBATCH --time=02:00:00           # 2 hours per job
#SBATCH --gres=gpu:1              # 1 GPU per job
#SBATCH --mem=8G                  # 8 GB RAM
#SBATCH --cpus-per-task=4         # 4 CPU cores
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

# Optional: Email notifications
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=your.email@example.com

#============================================================================
# Configuration - MODIFY THESE
#============================================================================

# Base directory for simulations
BASE_DIR="${HOME}/PARP14_library"

# Conda environment name
CONDA_ENV="calvados"

# Number of replicates per composition (for PARP14 library mode)
N_REPLICATES=25

# Simulation mode: "single" or "library"
MODE="single"  # Change to "library" for PARP14-style batch

#============================================================================
# Job Decoding
#============================================================================

# For single mode: SLURM_ARRAY_TASK_ID = simulation number
# For library mode: Decode composition and replicate from task ID

if [ "$MODE" == "library" ]; then
    # PARP14 library mode: decode composition and replicate
    COMP_ID=$(( (SLURM_ARRAY_TASK_ID - 1) / N_REPLICATES + 1 ))
    REP_ID=$(( (SLURM_ARRAY_TASK_ID - 1) % N_REPLICATES + 1 ))

    COMP_DIR=$(printf "comp_%04d" $COMP_ID)
    REP_DIR=$(printf "rep_%02d" $REP_ID)
    SIM_PATH="${BASE_DIR}/simulations/${COMP_DIR}/${REP_DIR}"
else
    # Single mode: use task ID directly
    SIM_ID=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
    SIM_PATH="${BASE_DIR}/simulations/sim_${SIM_ID}"
fi

#============================================================================
# Environment Setup
#============================================================================

echo "=============================================="
echo "CALVADOS Simulation Job"
echo "=============================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Simulation path: ${SIM_PATH}"
echo "=============================================="

# Load conda
source ~/.bashrc
conda activate ${CONDA_ENV}

# Verify GPU is available
if command -v nvidia-smi &> /dev/null; then
    echo "GPU Info:"
    nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv
fi

#============================================================================
# Run Simulation
#============================================================================

# Change to simulation directory
if [ ! -d "${SIM_PATH}" ]; then
    echo "Error: Simulation directory not found: ${SIM_PATH}"
    exit 1
fi

cd ${SIM_PATH}

# Check for run.py
if [ ! -f "run.py" ]; then
    echo "Error: run.py not found in ${SIM_PATH}"
    exit 1
fi

# Run simulation
echo ""
echo "Starting simulation at $(date)"
echo "----------------------------------------------"

python run.py --path .

EXIT_CODE=$?

echo "----------------------------------------------"
echo "Simulation finished at $(date)"
echo "Exit code: ${EXIT_CODE}"

#============================================================================
# Post-processing (optional)
#============================================================================

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Simulation completed successfully"

    # Optional: Run analysis
    # python ${BASE_DIR}/scripts/analyze_trajectory.py --top top.pdb --traj *.dcd

    # Optional: Compress trajectory
    # gzip *.dcd

    # Optional: Clean up
    # rm -f *.xml *.log
else
    echo "Simulation failed with exit code ${EXIT_CODE}"
fi

echo "=============================================="
echo "Job completed"
echo "=============================================="

exit ${EXIT_CODE}
