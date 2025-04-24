#!/bin/bash

# Usage:
#   ./run_all.sh [serial|parallel]
# Default is serial.

# ---------------------- Parse Run Mode ----------------------
RUN_MODE=${1:-serial}
if [[ "$RUN_MODE" != "serial" && "$RUN_MODE" != "parallel" ]]; then
    echo "Invalid run mode: $RUN_MODE"
    echo "Usage: $0 [serial|parallel]"
    exit 1
fi
echo "[INFO] Run mode: $RUN_MODE"

# ---------------------- Constants ----------------------
Kx=0.33;    Ky=0.9;
r=0.010;    t=0.005
L=0.300
xmin=0.0;   ymin=0.0
xmax=25.0;  ymax=25.0
nt=15

MW_START=2; MW_END=12
MB_START=2; MB_END=12

CASE_TEMPLATE="ELL-case-template"
OUTPUT_DIR="./outputs"
LOG_DIR="./logs"

# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Compute total iterations and initialize counter
TOTAL_ITER=$(( (MW_END - MW_START + 1) * (MB_END - MB_START + 1) ))
COUNTER=1

SECONDS=0
echo "##################################################"
echo "# Starting parametric sweep"
echo "# Total cases to run: ${TOTAL_ITER}"
echo "# Time started: $(date)"
echo "##################################################"
echo

# Loop over Mw and Mb
for Mw in $(seq $MW_START $MW_END); do
    for Mb in $(seq $MB_START $MB_END); do
        fname=$(printf "ELL-%d-%d-%.0f-%.0f-%.0f-%.0f-%.0f" \
            "$Mw" "$Mb" \
            "$(echo "$Kx * 100" | bc -l)" \
            "$(echo "$Ky * 100" | bc -l)" \
            "$(echo "$r * 1000" | bc -l)" \
            "$(echo "$L * 1000" | bc -l)" \
            "$(echo "$t * 1000" | bc -l)" )

        case_path="${OUTPUT_DIR}/${fname}"
        log_file="${LOG_DIR}/${fname}.log"
        {
            echo "=================================================="
            echo "Running case ${COUNTER}/${TOTAL_ITER}: $fname"
            echo "=================================================="
        
            # 1) Copy template if it doesn't exist
            if [ ! -d "$case_path" ]; then
                echo "[INFO] Copying template case"
                cp -r "${CASE_TEMPLATE}" "$case_path"
                chmod +x "${case_path}/Allrun"
            else
                echo "[SKIP] Case folder already exists"
            fi

            # 2) Generate mesh if it doesn't exists
            msh_file="${case_path}/constant/triSurface/${fname}.msh"
            if ls $msh_file 1> /dev/null 2>&1; then
                echo "[SKIP] Mesh already exists for $fname"
            else
                save_dir="${case_path}/constant/triSurface"
                echo "[INFO] Generating mesh for $fname"
                python3 mesh.py \
                    --Mw $Mw --Mb $Mb \
                    --Kx $Kx --Ky $Ky \
                    --r $r --t $t --L $L \
                    --xmin $xmin --ymin $ymin \
                    --xmax $xmax --ymax $ymax \
                    --nt $nt \
                    --sd $save_dir
            fi

            # 3) Run Allrun (pass parallel flag if requested)
            echo "[INFO] Running OpenFOAM simulation ($RUN_MODE)"
            if [ "$RUN_MODE" = "parallel" ]; then
                (cd "$case_path" && ./Allrun parallel)
            else
                (cd "$case_path" && ./Allrun)
            fi

            echo "[DONE] Case finished: $fname"
            echo
        } 2>&1 | tee "$log_file"
        ((COUNTER++))
    done
done

ELAPSED_TIME=$SECONDS
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS_REMAINING=$((ELAPSED_TIME % 60))

echo
echo "##################################################"
echo "# All cases completed!"
echo "# Total cases run: ${TOTAL_ITER}"
echo "# Time finished: $(date)"
echo "# Elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS_REMAINING}s"
echo "##################################################"

