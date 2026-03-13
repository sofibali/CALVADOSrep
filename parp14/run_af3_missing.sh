#!/bin/bash
#
# Run missing PARP14 AlphaFold3 predictions (1288 combinations)
#
# Two-phase approach for maximum GPU efficiency:
#   Phase 1: Data pipeline (MSA/template search) -- CPU only, highly parallel
#   Phase 2: Structure inference -- GPU, model loads ONCE via --input_dir
#
# Usage:
#   ./run_af3_missing.sh                    # default: 2 GPUs, 20 MSA parallel
#   ./run_af3_missing.sh --gpu0 0 --gpu1 1  # specify GPU IDs
#   ./run_af3_missing.sh --msa-parallel 16  # control CPU parallelism
#   ./run_af3_missing.sh --skip-msa         # skip phase 1 (already done)
#   ./run_af3_missing.sh --dry-run          # show plan without running
#
# Monitor:
#   tail -f run_af3_phase1_msa.log
#   tail -f run_af3_gpu0_inference.log
#   tail -f run_af3_gpu1_inference.log
#   find alphafold_outputs -name 'model.cif' | wc -l

set +e  # Continue even if individual predictions fail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# Defaults
# ============================================================================
GPU0=0
GPU1=1
MSA_PARALLEL=8
MSA_CPUS_PER_JOB=2         # CPU cores per MSA job (MSA_PARALLEL * MSA_CPUS_PER_JOB = total cores used)
DRY_RUN=0
SKIP_MSA=0
SKIP_INFERENCE=0
JAX_CACHE_DIR="$SCRIPT_DIR/.jax_cache"
NUM_DIFFUSION_SAMPLES=5

# Paths
JSON_DIR="$SCRIPT_DIR/alphafold_inputs_missing"
OUTPUT_DIR="$SCRIPT_DIR/alphafold_outputs"
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"

# ============================================================================
# Parse arguments
# ============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpu0)           GPU0="$2"; shift 2 ;;
        --gpu1)           GPU1="$2"; shift 2 ;;
        --msa-parallel)   MSA_PARALLEL="$2"; shift 2 ;;
        --msa-cpus)       MSA_CPUS_PER_JOB="$2"; shift 2 ;;
        --skip-msa)       SKIP_MSA=1; shift ;;
        --skip-inference) SKIP_INFERENCE=1; shift ;;
        --msa-only)       SKIP_INFERENCE=1; shift ;;
        --inference-only) SKIP_MSA=1; shift ;;
        --dry-run)        DRY_RUN=1; shift ;;
        --diffusion-samples) NUM_DIFFUSION_SAMPLES="$2"; shift 2 ;;
        *)                echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ============================================================================
# Environment
# ============================================================================
source /programs/sbgrid.shrc 2>/dev/null
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas

mkdir -p "$JAX_CACHE_DIR" "$OUTPUT_DIR"

# ============================================================================
# Count inputs, split between GPUs
# ============================================================================
ALL_JSONS=( "$JSON_DIR"/*.json )
TOTAL=${#ALL_JSONS[@]}

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No JSON files found in $JSON_DIR"
    exit 1
fi

# Count how many still need to be run (no model.cif in output)
NEED_RUN=()
ALREADY_DONE=0
for json_file in "${ALL_JSONS[@]}"; do
    name=$(basename "$json_file" .json)
    out_dir="$OUTPUT_DIR/${name,,}"
    if [ -n "$(find "$out_dir" -name 'model.cif' -print -quit 2>/dev/null)" ]; then
        ((ALREADY_DONE++))
    else
        NEED_RUN+=("$json_file")
    fi
done

N_TODO=${#NEED_RUN[@]}

# Split between 2 GPUs
HALF=$(( (N_TODO + 1) / 2 ))
GPU0_JSONS=("${NEED_RUN[@]:0:$HALF}")
GPU1_JSONS=("${NEED_RUN[@]:$HALF}")

# Time estimates (~5 min/prediction on A100 for avg ~930 residue protein)
INF_MIN=5
GPU0_HRS=$(( ${#GPU0_JSONS[@]} * INF_MIN / 60 ))
GPU1_HRS=$(( ${#GPU1_JSONS[@]} * INF_MIN / 60 ))
MSA_WALL_HRS=$(( N_TODO * 30 / MSA_PARALLEL / 60 ))  # ~30 min/seq MSA

echo "============================================================================"
echo "PARP14 AlphaFold3 -- Missing Predictions"
echo "============================================================================"
echo ""
echo "Total input JSONs:     $TOTAL"
echo "Already completed:     $ALREADY_DONE"
echo "Still needed:          $N_TODO"
echo "JAX cache:             $JAX_CACHE_DIR"
echo ""
echo "--- Phase 1: Data Pipeline ($MSA_PARALLEL jobs × $MSA_CPUS_PER_JOB cores = $(( MSA_PARALLEL * MSA_CPUS_PER_JOB )) total cores) ---"
if [ $SKIP_MSA -eq 1 ]; then
    echo "  SKIPPED (--skip-msa)"
else
    echo "  Est. wall time: ~${MSA_WALL_HRS} hours"
fi
echo ""
echo "--- Phase 2: Inference (2 GPUs, model loads once per GPU) ---"
echo "  GPU 0 ($GPU0): ${#GPU0_JSONS[@]} predictions (~${GPU0_HRS} hours)"
echo "  GPU 1 ($GPU1): ${#GPU1_JSONS[@]} predictions (~${GPU1_HRS} hours)"
echo ""
TOTAL_HRS=$(( MSA_WALL_HRS + (GPU0_HRS > GPU1_HRS ? GPU0_HRS : GPU1_HRS) ))
echo "  Est. total wall time: ~${TOTAL_HRS} hours (~$(( (TOTAL_HRS + 23) / 24 )) days)"
echo "============================================================================"
echo ""

if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Exiting."
    exit 0
fi

START_TIME=$(date +%s)

# ============================================================================
# PHASE 1: Data Pipeline (MSA/template search) -- CPU only, parallel
# ============================================================================
run_msa_for_json() {
    local json_file=$1
    local name=$(basename "$json_file" .json)
    local out_dir="$OUTPUT_DIR/${name,,}"

    # Skip if data pipeline or inference already done
    if [ -n "$(find "$out_dir" -name '*_data.json' -print -quit 2>/dev/null)" ]; then
        return 0
    fi
    if [ -n "$(find "$out_dir" -name 'model.cif' -print -quit 2>/dev/null)" ]; then
        return 0
    fi

    # Limit CPU cores: cap threads for hmmer/hhblits/jackhmmer inside AF3
    export OMP_NUM_THREADS=$MSA_CPUS_PER_JOB
    export OPENMM_CPU_THREADS=$MSA_CPUS_PER_JOB
    export MKL_NUM_THREADS=$MSA_CPUS_PER_JOB
    export NUMEXPR_MAX_THREADS=$MSA_CPUS_PER_JOB
    export TF_NUM_INTEROP_THREADS=$MSA_CPUS_PER_JOB
    export TF_NUM_INTRAOP_THREADS=$MSA_CPUS_PER_JOB
    export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false"

    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --json_path "$json_file" \
        --norun_inference \
        --nhmmer_n_cpu "$MSA_CPUS_PER_JOB" \
        --jackhmmer_n_cpu "$MSA_CPUS_PER_JOB" \
        &> "$OUTPUT_DIR/logs/msa_${name,,}.log"
}

if [ $SKIP_MSA -eq 0 ]; then
    echo "============================================================"
    echo "PHASE 1: Data Pipeline (MSA, $MSA_PARALLEL parallel)"
    echo "Started: $(date)"
    echo "============================================================"

    mkdir -p "$OUTPUT_DIR/logs"
    export -f run_msa_for_json
    export OUTPUT_DIR DB_DIR MODEL_DIR MSA_CPUS_PER_JOB

    ACTIVE=0
    SUBMITTED=0

    for json_file in "${NEED_RUN[@]}"; do
        run_msa_for_json "$json_file" &
        ((ACTIVE++))
        ((SUBMITTED++))

        if [ $ACTIVE -ge $MSA_PARALLEL ]; then
            wait -n 2>/dev/null || wait
            ((ACTIVE--))
        fi

        if [ $((SUBMITTED % 50)) -eq 0 ]; then
            echo "[$(date '+%H:%M:%S')] Phase 1: $SUBMITTED/$N_TODO submitted"
        fi
    done

    wait
    echo "[$(date '+%H:%M:%S')] Phase 1 complete: $SUBMITTED MSA jobs finished"
    echo ""
else
    echo "PHASE 1: SKIPPED"
    echo ""
fi

# ============================================================================
# PHASE 2: Inference -- model loads ONCE per GPU via --input_dir
# ============================================================================
prepare_batch() {
    # Create a temp dir with symlinks to _data.json or raw .json for one GPU's share
    local batch_dir=$1
    shift
    local json_files=("$@")
    local n_data=0
    local n_raw=0

    rm -rf "$batch_dir"
    mkdir -p "$batch_dir"

    for json_file in "${json_files[@]}"; do
        local name=$(basename "$json_file" .json)
        local out_dir="$OUTPUT_DIR/${name,,}"

        # Skip already completed
        if [ -n "$(find "$out_dir" -name 'model.cif' -print -quit 2>/dev/null)" ]; then
            continue
        fi

        # Prefer _data.json (MSA pre-computed)
        local data_json=$(find "$out_dir" -name '*_data.json' -print -quit 2>/dev/null)
        if [ -n "$data_json" ]; then
            ln -sf "$data_json" "$batch_dir/"
            ((n_data++))
        else
            ln -sf "$json_file" "$batch_dir/"
            ((n_raw++))
        fi
    done

    echo "$n_data pre-computed MSA, $n_raw need full pipeline"
}

run_gpu_inference() {
    local gpu_id=$1
    local label=$2
    local batch_dir=$3
    local log_file="$SCRIPT_DIR/run_af3_gpu${label}_inference.log"

    local n_inputs=$(ls "$batch_dir"/*.json 2>/dev/null | wc -l)
    if [ "$n_inputs" -eq 0 ]; then
        echo "GPU $label: nothing to process" | tee "$log_file"
        return 0
    fi

    echo "[$(date '+%H:%M:%S')] GPU $label ($gpu_id): $n_inputs inputs" | tee "$log_file"

    # Check if all are _data.json -> skip data pipeline
    local n_data=$(ls "$batch_dir"/*_data.json 2>/dev/null | wc -l)
    local pipeline_flag=""
    if [ "$n_data" -eq "$n_inputs" ]; then
        pipeline_flag="--norun_data_pipeline"
        echo "  All inputs have pre-computed MSAs" | tee -a "$log_file"
    fi

    CUDA_VISIBLE_DEVICES=$gpu_id \
    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --input_dir "$batch_dir" \
        --jax_compilation_cache_dir "$JAX_CACHE_DIR" \
        --flash_attention_implementation triton \
        --num_diffusion_samples $NUM_DIFFUSION_SAMPLES \
        $pipeline_flag \
        >> "$log_file" 2>&1

    echo "[$(date '+%H:%M:%S')] GPU $label: done (exit $?)" | tee -a "$log_file"
}

if [ $SKIP_INFERENCE -eq 0 ]; then
    echo "============================================================"
    echo "PHASE 2: Inference (2 GPUs, model loads once each)"
    echo "Started: $(date)"
    echo "============================================================"
    echo ""

    BATCH0="$SCRIPT_DIR/.af3_batch_gpu0"
    BATCH1="$SCRIPT_DIR/.af3_batch_gpu1"

    echo "Preparing GPU 0 batch..."
    GPU0_SUMMARY=$(prepare_batch "$BATCH0" "${GPU0_JSONS[@]}")
    echo "  GPU 0: $GPU0_SUMMARY"

    echo "Preparing GPU 1 batch..."
    GPU1_SUMMARY=$(prepare_batch "$BATCH1" "${GPU1_JSONS[@]}")
    echo "  GPU 1: $GPU1_SUMMARY"
    echo ""

    echo "Launching inference on both GPUs..."
    echo "  Monitor: tail -f $SCRIPT_DIR/run_af3_gpu0_inference.log"
    echo "           tail -f $SCRIPT_DIR/run_af3_gpu1_inference.log"
    echo ""

    run_gpu_inference "$GPU0" "0" "$BATCH0" &
    PID0=$!
    echo "  GPU 0: PID $PID0"

    run_gpu_inference "$GPU1" "1" "$BATCH1" &
    PID1=$!
    echo "  GPU 1: PID $PID1"

    echo ""
    echo "Waiting for both GPUs..."
    FAILURES=0
    wait $PID0 || ((FAILURES++))
    wait $PID1 || ((FAILURES++))

    rm -rf "$BATCH0" "$BATCH1"
else
    echo "PHASE 2: SKIPPED"
fi

# ============================================================================
# Summary
# ============================================================================
END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))

echo ""
echo "============================================================================"
echo "PARP14 AlphaFold3 Batch Complete"
echo "============================================================================"
echo "  Elapsed: $(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m"
echo ""

# Count completed
N_COMPLETE=$(find "$OUTPUT_DIR" -maxdepth 3 -name 'model.cif' -print 2>/dev/null | \
    sed 's|/seed-[0-9]*_sample-[0-9]*/model.cif||' | sort -u | wc -l)

echo "  Completed combinations: $N_COMPLETE / 2047 target"
echo "  Total model.cif files:  $(find "$OUTPUT_DIR" -name 'model.cif' | wc -l)"

if [ "${FAILURES:-0}" -gt 0 ]; then
    echo "  WARNING: $FAILURES GPU(s) had errors. Check inference logs."
fi
echo "============================================================================"
