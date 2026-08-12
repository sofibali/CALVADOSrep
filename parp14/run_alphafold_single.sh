#!/bin/bash
set -e

source /programs/sbgrid.shrc
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
export CUDA_VISIBLE_DEVICES=0

# Usage: bash run_alphafold_single.sh <path/to/input.json> [output_dir]

JSON_PATH="$1"
OUTPUT_DIR="${2:-alphafold_outputs}"

if [ -z "$JSON_PATH" ]; then
    echo "Usage: $0 <path/to/input.json> [output_dir]"
    exit 1
fi

if [ ! -f "$JSON_PATH" ]; then
    echo "Error: JSON file not found: $JSON_PATH"
    exit 1
fi

# Configuration
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"

mkdir -p "$OUTPUT_DIR"

BASENAME=$(basename "$JSON_PATH" .json)
echo "Processing: $BASENAME"

if /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
    --db_dir "$DB_DIR" \
    --model_dir "$MODEL_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --json_path "$JSON_PATH"; then
    echo "  ✓ Completed"
else
    echo "  ✗ Failed"
    exit 1
fi
