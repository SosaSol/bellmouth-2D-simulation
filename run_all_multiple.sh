#!/bin/bash
set -euo pipefail # Exit on error, unset variable, or pipe failure
trap 'echo "[ERROR] Script interrupted"; exit 1' SIGINT # Handle Ctrl+C

# ---------------------- Usage ----------------------
# ./run_all.sh [--serial|--parallel N] [--overwrite]

# ---------------------- Default Parameters ----------------------
RUN_MODE="serial"
NP=12
OVERWRITE=false

# ---------------------- Parse Named Arguments ----------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel)
            RUN_MODE="parallel"
            if [[ -n "${2:-}" && "$2" =~ ^[0-9]+$ ]]; then
                NP="$2"
                shift 2
            else
                echo "[ERROR] --parallel requires a numeric argument (e.g., --parallel 8)"
                exit 1
            fi
            ;;
        --serial)
            RUN_MODE="serial"
            shift
            ;;
        --overwrite)
            OVERWRITE=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [--serial|--parallel N] [--overwrite]"
            echo "    --serial               Run in serial mode (default)"
            echo "    --parallel N           Run in parallel using N processors"
            echo "    --overwrite            Overwrite existing case folders"
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1"
            echo "Use --help for usage."
            exit 1
            ;;
    esac
done

# ---------------------- Summary & Validation ----------------------
echo "[INFO] Run mode: $RUN_MODE"
[[ "$RUN_MODE" == "parallel" ]] && echo "[INFO] Number of processors: $NP"
echo "[INFO] Overwrite existing cases: $OVERWRITE"

if [[ "$RUN_MODE" == "parallel" && "$NP" -le 1 ]]; then
    echo "[ERROR] Cannot run in parallel with $NP processor(s). Use at least 2."
    exit 1
fi

# ---------------------- Constants ----------------------
Kx=0.33;    Ky=0.33;
r=0.020;    t=0.005
L=0.339
xmin=0.0;   ymin=0.0
xmax=25.0;  ymax=25.0
nt=$NP

MW_START=2; MW_END=12
MB_START=2; MB_END=12

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
CASE_TEMPLATE="${SCRIPT_DIR}/ELL-case-template"
OUTPUT_DIR="${SCRIPT_DIR}/outputsTest"
LOG_DIR="${SCRIPT_DIR}/logsTest"

# ---------------------- Function: Run Single Case ----------------------
run_case() {
    local Mw=$1
    local Mb=$2
    local fname=$3
    local case_path=$4
    local log_file=$5

    {
        echo "=================================================="
        echo "[INFO] Case ${COUNTER}/${TOTAL_ITER}: $fname"
        echo "=================================================="

        # 1) Prepare case directory
        if [ -d "$case_path" ]; then
            echo "[INFO] Overwriting existing case"
            rm -rf "$case_path"
        else
            echo "[INFO] Creating new case from template"
        fi

        cp -r "${CASE_TEMPLATE}" "$case_path"
        chmod +x "${case_path}/Allrun"

        # 2) Generate mesh
        msh_file="${case_path}/constant/triSurface/${fname}.msh"
        if [ -f "$msh_file" ]; then
            echo "[SKIP] Mesh already exists for $fname"
        else
            save_dir="${case_path}/constant/triSurface"
            echo "[INFO] Generating mesh..."
            python3 mesh.py \
                --Mw "$Mw" --Mb "$Mb" \
                --Kx "$Kx" --Ky "$Ky" \
                --r "$r" --t "$t" --L "$L" \
                --xmin "$xmin" --ymin "$ymin" \
                --xmax "$xmax" --ymax "$ymax" \
                --nt "$nt" \
                --sd "$save_dir"
        fi

        # 3) Run simulation
        echo "[INFO] Running OpenFOAM simulation..."
        pushd "$case_path" > /dev/null
        if [[ "$RUN_MODE" == "parallel" ]]; then
            ./Allrun parallel "$NP"
        else
            ./Allrun
        fi
        popd > /dev/null

        echo "[DONE] Case finished: $fname"
        echo
    } 2>&1 | tee "$log_file"
}

# ---------------------- Beginning of Script ----------------------
# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Compute total iterations and initialize counter
TOTAL_ITER=$(( (MW_END - MW_START + 1) * (MB_END - MB_START + 1) ))
COUNTER=1

MAX_JOBS=3
JOBS=0

SECONDS=0
echo "##################################################"
echo "# Parametric sweep started: $(date)"
echo "# Total cases: $TOTAL_ITER"
echo "##################################################"
echo

# ---------------------- Main Loop ----------------------
for Mw in $(seq "$MW_START" "$MW_END"); do
    for Mb in $(seq "$MB_START" "$MB_END"); do
        fname=$(printf "ELL-%d-%d-%.0f-%.0f-%.0f-%.0f-%.0f" \
            "$Mw" "$Mb" \
            "$(awk "BEGIN{printf \"%.0f\", $Kx * 100}")" \
            "$(awk "BEGIN{printf \"%.0f\", $Ky * 100}")" \
            "$(awk "BEGIN{printf \"%.0f\", $r * 1000}")" \
            "$(awk "BEGIN{printf \"%.0f\", $L * 1000}")" \
            "$(awk "BEGIN{printf \"%.0f\", $t * 1000}")")

        case_path="${OUTPUT_DIR}/${fname}"
        log_file="${LOG_DIR}/${fname}.log"

        if [ -d "$case_path" ] && [ "$OVERWRITE" != "true" ]; then
            echo "[SKIP] Case folder: $fname already exists and OVERWRITE is false" | tee "$log_file"
            ((COUNTER++))
            continue
        fi

        # Run case in background
        run_case "$Mw" "$Mb" "$fname" "$case_path" "$log_file" &
        ((JOBS++))
        ((COUNTER++))

        # Limit concurrent jobs
        if (( JOBS >= MAX_JOBS )); then
            wait -n  # Wait for *any* job to finish
            ((JOBS--))
        fi
    done
done

# Wait for all remaining background jobs to finish
wait
done
# ---------------------- End Main Loop ----------------------

# ---------------------- End Timing ----------------------
ELAPSED_TIME=$SECONDS
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo
echo "##################################################"
echo "# All ${TOTAL_ITER} cases finished: $(date)"
echo "# Total cases run: $((COUNTER - 1))"
echo "# Time finished: $(date)"
echo "# Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "##################################################"