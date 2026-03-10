#!/bin/bash
#
# Run pLDDT analysis in background with progress logging
# activate calvados enviornment before running 
# Usage: bash run_pLDDT_analysis.sh

LOG_FILE="pLDDT_analysis_$(date +%Y%m%d_%H%M%S).log"

echo "Starting pLDDT analysis..."
echo "Output log: $LOG_FILE"
echo ""
echo "To watch progress in real-time:"
echo "  tail -f $LOG_FILE"
echo ""

# Run in background, redirect output to log file
nohup python 05_pLDDT_analysis.py \
    --input /home/sbali/parp14/alphafold_outputs \
    --output-residue /home/sbali/parp14/pipeline/pLDDT_per_residue.csv \
    --output-structure /home/sbali/parp14/pipeline/confidence_scores.csv \
    > "$LOG_FILE" 2>&1 &

PID=$!

echo "Process started with PID: $PID"
echo "Command to check if still running:"
echo "  ps -p $PID"
echo ""
echo "Command to stop if needed:"
echo "  kill $PID"
