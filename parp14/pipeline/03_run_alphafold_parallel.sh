#!/bin/bash
#
# Parallel AlphaFold3 Batch Processing Script
# Runs multiple predictions in parallel on a single GPU to maximize compute utilization
#
# Usage:
#   bash 03_run_alphafold_parallel.sh
#   or with nohup:
#   nohup bash 03_run_alphafold_parallel.sh > alphafold_parallel.log 2>&1 &
#
# Configuration:
#   - Edit PARALLEL_JOBS to control how many predictions run simultaneously
#   - Typical values: 2-4 for a single GPU depending on VRAM and model size
#   - Monitor GPU usage with: watch -n 1 nvidia-smi

set -e

# Configuration
LOG_FILE="alphafold_parallel_$(date +%Y%m%d_%H%M%S).log"
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"
OUTPUT_DIR="/home/sbali/parp14/alphafold_outputs_missing"
JSON_DIR="/home/sbali/parp14/pipeline/alphafold_inputs_missing"
MISSING_LIST="/home/sbali/parp14/missing_predictions.txt"

# Parallel job control
PARALLEL_JOBS=3  # Number of simultaneous predictions (adjust based on GPU memory)
export CUDA_VISIBLE_DEVICES=0
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas

echo "===================================================================" | tee -a "$LOG_FILE"
echo "AlphaFold3 Parallel Batch Processing Script" | tee -a "$LOG_FILE"
echo "Started at: $(date)" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "Parallel jobs: $PARALLEL_JOBS" | tee -a "$LOG_FILE"
echo "GPU: $CUDA_VISIBLE_DEVICES" | tee -a "$LOG_FILE"
echo "===================================================================" | tee -a "$LOG_FILE"

echo "Sourcing SBGrid environment..." | tee -a "$LOG_FILE"
source /programs/sbgrid.shrc 2>&1 | grep -v "^[*]" | tee -a "$LOG_FILE" || true

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if missing predictions list exists
if [ ! -f "$MISSING_LIST" ]; then
    echo "Error: Missing predictions list not found: $MISSING_LIST" | tee -a "$LOG_FILE"
    echo "Run check_missing.py first to generate the list." | tee -a "$LOG_FILE"
    exit 1
fi

# Count total missing
TOTAL_MISSING=$(wc -l < "$MISSING_LIST")
echo "Total predictions to run: $TOTAL_MISSING" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Function to run a single prediction
run_prediction() {
    local basename=$1
    local json_file="$JSON_DIR/${basename}.json"
    local output_subdir="$OUTPUT_DIR/${basename,,}"
    local job_log="${OUTPUT_DIR}/logs/${basename,,}.log"
    
    # Create logs directory
    mkdir -p "${OUTPUT_DIR}/logs"
    
    # Skip if already exists
    if [ -d "$output_subdir" ] && [ "$(ls -A "$output_subdir" 2>/dev/null)" ]; then
        echo "SKIPPED: $basename (already exists)" >> "$LOG_FILE"
        return 0
    fi
    
    # Check JSON file exists
    if [ ! -f "$json_file" ]; then
        echo "ERROR: JSON file not found: $json_file" >> "$LOG_FILE"
        return 1
    fi
    
    echo "STARTED: $basename at $(date)" >> "$LOG_FILE"
    
    # Run AlphaFold3
    if /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --json_path "$json_file" > "$job_log" 2>&1; then
        echo "COMPLETED: $basename at $(date)" >> "$LOG_FILE"
        return 0
    else
        echo "FAILED: $basename at $(date)" >> "$LOG_FILE"
        return 1
    fi
}

export -f run_prediction
export DB_DIR MODEL_DIR OUTPUT_DIR JSON_DIR LOG_FILE

# Run predictions in parallel using GNU parallel or xargs
if command -v parallel &> /dev/null; then
    # Use GNU parallel if available (preferred)
    echo "Using GNU parallel for job management..." | tee -a "$LOG_FILE"
    cat "$MISSING_LIST" | parallel -j $PARALLEL_JOBS --bar --joblog "${OUTPUT_DIR}/parallel_jobs.log" run_prediction {}
else
    # Fallback to xargs with manual job control
    echo "Using xargs for job management (install GNU parallel for better progress tracking)..." | tee -a "$LOG_FILE"
    echo "Install with: sudo apt-get install parallel  or  brew install parallel" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Process with xargs
    cat "$MISSING_LIST" | xargs -P $PARALLEL_JOBS -I {} bash -c 'run_prediction "$@"' _ {}
fi

# Summary
echo "" | tee -a "$LOG_FILE"
echo "=======================================================================" | tee -a "$LOG_FILE"
echo "Parallel Batch Processing Summary" | tee -a "$LOG_FILE"
echo "Finished at: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Count results
COMPLETED_COUNT=$(find "$OUTPUT_DIR" -maxdepth 1 -type d | wc -l)
COMPLETED_COUNT=$((COMPLETED_COUNT - 1))  # Subtract 1 for parent dir
REMAINING=$((TOTAL_MISSING - COMPLETED_COUNT))

echo "Results:" | tee -a "$LOG_FILE"
echo "  Total to process: $TOTAL_MISSING" | tee -a "$LOG_FILE"
echo "  Completed: $COMPLETED_COUNT" | tee -a "$LOG_FILE"
echo "  Remaining: $REMAINING" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Check for failures
if [ -d "${OUTPUT_DIR}/logs" ]; then
    FAILED_COUNT=$(grep -l "Error\|error\|FAILED" "${OUTPUT_DIR}/logs"/*.log 2>/dev/null | wc -l || echo 0)
    echo "  Failed jobs: $FAILED_COUNT" | tee -a "$LOG_FILE"
    if [ "$FAILED_COUNT" -gt 0 ]; then
        echo "  Check individual logs in: ${OUTPUT_DIR}/logs/" | tee -a "$LOG_FILE"
    fi
fi

echo "=======================================================================" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "Job logs: ${OUTPUT_DIR}/logs/" | tee -a "$LOG_FILE"
if command -v parallel &> /dev/null; then
    echo "Parallel job log: ${OUTPUT_DIR}/parallel_jobs.log" | tee -a "$LOG_FILE"
fi
echo "=======================================================================" | tee -a "$LOG_FILE"
