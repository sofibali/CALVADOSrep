#!/bin/bash
#
# Run SASA analysis in background with progress logging
#
# Usage: bash run_sasa_analysis.sh

LOG_FILE="sasa_analysis_$(date +%Y%m%d_%H%M%S).log"

echo "Starting SASA analysis..."
echo "Output log: $LOG_FILE"
echo ""
echo "To watch progress in real-time:"
echo "  tail -f $LOG_FILE"
echo ""

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate calvados

# Run in background, redirect output to log file
nohup python 06_sasa_analysis.py \
    --input /home/sbali/parp14/alphafold_outputs \
    --domains /home/sbali/parp14/Domain_boundries.csv \
    --output /home/sbali/parp14/sasa_per_residue.csv \
    > "$LOG_FILE" 2>&1 &

PID=$!

echo "Process started with PID: $PID"
echo "Command to check if still running:"
echo "  ps -p $PID"
echo ""
echo "Command to stop if needed:"
echo "  kill $PID"
