#!/bin/bash
set -euo pipefail # Exit on error, unset variable, or pipe failure
trap 'echo "[ERROR] Script interrupted"; exit 1' SIGINT # Handle Ctrl+C

# Usage:
#   ./run_all.sh [serial|parallel] [number of cores]
# Default is serial.

# ---------------------- Parse Run Mode ----------------------
RUN_MODE=${1:-serial}       # Default to serial if no argument is passed
NP=${2:-12}                 # Default to 12 processors if none provided

if [[ "$RUN_MODE" != "serial" && "$RUN_MODE" != "parallel" ]]; then
    echo "Invalid run mode: $RUN_MODE"
    echo "Usage: $0 [serial|parallel] [nCores]"
    exit 1
fi

echo "[INFO] Run mode: $RUN_MODE"
[[ "$RUN_MODE" == "parallel" ]] && echo "[INFO] Number of processors: $NP"

# Validate core count if running in parallel
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
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
LOG_DIR="${SCRIPT_DIR}/logs"

# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Compute total iterations and initialize counter
TOTAL_ITER=$(( (MW_END - MW_START + 1) * (MB_END - MB_START + 1) ))
COUNTER=1

SECONDS=0
echo "##################################################"
echo "# Parametric sweep started: $(date)"
echo "# Total cases: $TOTAL_ITER"
echo "##################################################"
echo

# ---------------------- Main Loop ----------------------
# Loop over Mw and Mb
for Mw in $(seq $MW_START $MW_END); do
    for Mb in $(seq $MB_START $MB_END); do
        fname=$(printf "ELL-%d-%d-%.0f-%.0f-%.0f-%.0f-%.0f" \
            "$Mw" "$Mb" \
            "$(awk "BEGIN{printf \"%.0f\", $Kx * 100}")" \
            "$(awk "BEGIN{printf \"%.0f\", $Ky * 100}")" \
            "$(awk "BEGIN{printf \"%.0f\", $r * 1000}")" \
            "$(awk "BEGIN{printf \"%.0f\", $L * 1000}")" \
            "$(awk "BEGIN{printf \"%.0f\", $t * 1000}")")

        case_path="${OUTPUT_DIR}/${fname}"
        log_file="${LOG_DIR}/${fname}.log"
        {
            echo "=================================================="
            echo "[INFO] Case ${COUNTER}/${TOTAL_ITER}: $fname"
            echo "=================================================="
        
            # 1) Copy template if it doesn't exist
            if [ ! -d "$case_path" ]; then
                echo "[INFO] Copying template case"
                cp -r "${CASE_TEMPLATE}" "$case_path"
                chmod +x "${case_path}/Allrun"
            else
                echo "[SKIP] Case folder exists: $case_path"
            fi

            # 2) Generate mesh if it doesn't exists
            msh_file="${case_path}/constant/triSurface/${fname}.msh"
            if [ -f "$msh_file" ]; then # mesh exists
                echo "[SKIP] Mesh already exists for $fname"
            else # mesh doesn't exist
                save_dir="${case_path}/constant/triSurface"
                echo "[INFO] Generating mesh..."
                python3 mesh.py \
                    --Mw $Mw --Mb $Mb \
                    --Kx $Kx --Ky $Ky \
                    --r $r --t $t --L $L \
                    --xmin $xmin --ymin $ymin \
                    --xmax $xmax --ymax $ymax \
                    --nt $nt \
                    --sd $save_dir
            fi

            # 3) Run Allrun
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
        ((COUNTER++))
    done
done

# ---------------------- End Timing ----------------------
ELAPSED_TIME=$SECONDS
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo
echo "##################################################"
echo "# All ${TOTAL_ITER} cases finished: $(date)"
echo "# Total cases run: "
echo "# Time finished: "
echo "# Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS_REMAINING}s"
echo "##################################################"

