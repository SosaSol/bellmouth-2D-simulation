#!/bin/bash
set -euo pipefail
trap 'echo "[ERROR] Script interrupted"; exit 1' SIGINT

# ---------------------- Usage ----------------------
# ./run_all.sh [--parallel N] [--overwrite] [--advanced]

# ---------------------- Default Parameters ----------------------
RUN_MODE="serial"
NP=12
OVERWRITE=false
ADVANCED=false

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
        --overwrite)
            OVERWRITE=true
            shift
            ;;
        --advanced)
            ADVANCED=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [--parallel N] [--overwrite] [--advanced]"
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
echo "[INFO] Advanced simulation: $ADVANCED"

if [[ "$RUN_MODE" == "parallel" && "$NP" -le 1 ]]; then
    echo "[ERROR] Cannot run in parallel with $NP processor(s). Use at least 2."
    exit 1
fi

# ---------------------- Constants ----------------------
Kx=0.33; Ky=0.33
r=0.050

L1=0.16866  # m
L2=0.13029  # m
L3=0.03986  # m
L=$(awk "BEGIN {print $L1 + $L2 + $L3}")  # m

GAP2=0.04256  # m
t=$(awk "BEGIN {print $GAP2 / 2}")  # m

xmin=0.0; ymin=0.0
xmax=25.0; ymax=25.0

nt=$NP

MW_START=2; MW_END=12
MB_START=2; MB_END=12

# ---------------------- Check Dependencies ----------------------

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
CASE_TEMPLATE="${SCRIPT_DIR}/ELL-case-template"
if [ ! -d "$CASE_TEMPLATE" ]; then
    echo "[ERROR] Case template directory not found: $CASE_TEMPLATE"
    exit 1
fi

if [ "$ADVANCED" == true ]; then
    OUTPUT_DIR="${SCRIPT_DIR}/outputs/advanced"
    LOG_DIR="${SCRIPT_DIR}/logs/advanced"
    MESH_SCRIPT="mesh_advanced_with_bellmouth.py"
    MESH_SCRIPT_STRAIGHT="mesh_advanced_no_bellmouth.py"
    MESH_FLAGS="--Kx $Kx --Ky $Ky --r $r --xmin $xmin --ymin $ymin --xmax $xmax --ymax $ymax --nt $nt"
    MESH_FLAGS_STRAIGHT="--xmin $xmin --ymin $ymin --xmax $xmax --ymax $ymax --nt $nt"
else
    OUTPUT_DIR="${SCRIPT_DIR}/outputs/simple"
    LOG_DIR="${SCRIPT_DIR}/logs/simple"
    MESH_SCRIPT="mesh.py"
    MESH_SCRIPT_STRAIGHT="mesh_straight.py"
    MESH_FLAGS="--Kx $Kx --Ky $Ky --r $r --t $t --L $L --xmin $xmin --ymin $ymin --xmax $xmax --ymax $ymax --nt $nt"
    MESH_FLAGS_STRAIGHT="--t $t --L $L --xmin $xmin --ymin $ymin --xmax $xmax --ymax $ymax --nt $nt"
fi

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

TOTAL_ITER=$(( (MW_END - MW_START + 1) * (MB_END - MB_START + 1) + (MW_END - MW_START + 1)))
COUNTER=1
SECONDS=0

echo "##################################################"
echo "# Parametric sweep started: $(date)"
echo "# Total cases: $TOTAL_ITER"
echo "##################################################"
echo

# ---------------------- Main Loop ----------------------
for Mw in $(seq "$MW_END" -1 "$MW_START"); do

    fname=$(printf "IN-%d-%d-%d" \
        "$Mw" \
        "$(awk "BEGIN{print int($L * 1000)}")" \
        "$(awk "BEGIN{print int($t * 1000)}")")

    case_path="${OUTPUT_DIR}/${fname}"
    log_file="${LOG_DIR}/${fname}.log"

    if [ -d "$case_path" ] && [ "$OVERWRITE" != "true" ]; then
        echo "[SKIP] Case folder: $fname already exists and OVERWRITE is false" | tee "$log_file"
        ((COUNTER++))
        continue
    fi

    {
        echo "=================================================="
        echo "[INFO] Case ${COUNTER}/${TOTAL_ITER}: $fname"
        echo "[INFO] Start: $(date)"
        echo "=================================================="

        if [ -d "$case_path" ]; then
            echo "[INFO] Overwriting existing case"
            rm -rf "$case_path"
        else
            echo "[INFO] Creating new case from template"
        fi
        cp -r "$CASE_TEMPLATE" "$case_path"
        chmod +x "${case_path}/Allrun"

        msh_file="${case_path}/constant/triSurface/${fname}.msh"
        if [ -f "$msh_file" ]; then
            echo "[SKIP] Mesh already exists for $fname"
        else
            save_dir="${case_path}/constant/triSurface"
            echo "[INFO] Generating mesh using $MESH_SCRIPT..."
            python3 "$SCRIPT_DIR/$MESH_SCRIPT_STRAIGHT" --Mw $Mw $MESH_FLAGS_STRAIGHT --sd "$save_dir"
        fi

        echo "[INFO] Running OpenFOAM simulation..."
        pushd "$case_path" > /dev/null
        if [[ "$RUN_MODE" == "parallel" ]]; then
            ./Allrun parallel "$NP"
        else
            ./Allrun
        fi
        popd > /dev/null

        echo "[DONE] Case finished: $fname"
        echo "[INFO] End: $(date)"
        echo
    } 2>&1 | tee "$log_file"

    ((COUNTER++))

    for Mb in $(seq "$MB_END" -1 "$MB_START"); do

        fname=$(printf "ELL-%d-%d-%d-%d-%d-%d-%d" \
            "$Mw" "$Mb" \
            "$(awk "BEGIN{print int($Kx * 100)}")" \
            "$(awk "BEGIN{print int($Ky * 100)}")" \
            "$(awk "BEGIN{print int($r * 1000)}")" \
            "$(awk "BEGIN{print int($L * 1000)}")" \
            "$(awk "BEGIN{print int($t * 1000)}")")

        case_path="${OUTPUT_DIR}/${fname}"
        log_file="${LOG_DIR}/${fname}.log"

        if [ -d "$case_path" ] && [ "$OVERWRITE" != "true" ]; then
            echo "[SKIP] Case folder: $fname already exists and OVERWRITE is false" | tee "$log_file"
            ((COUNTER++))
            continue
        fi

        {
            echo "=================================================="
            echo "[INFO] Case ${COUNTER}/${TOTAL_ITER}: $fname"
            echo "[INFO] Start: $(date)"
            echo "=================================================="

            if [ -d "$case_path" ]; then
                echo "[INFO] Overwriting existing case"
                rm -rf "$case_path"
            else
                echo "[INFO] Creating new case from template"
            fi
            cp -r "$CASE_TEMPLATE" "$case_path"
            chmod +x "${case_path}/Allrun"

            msh_file="${case_path}/constant/triSurface/${fname}.msh"
            if [ -f "$msh_file" ]; then
                echo "[SKIP] Mesh already exists for $fname"
            else
                save_dir="${case_path}/constant/triSurface"
                echo "[INFO] Generating mesh using $MESH_SCRIPT..."
                python3 "$SCRIPT_DIR/$MESH_SCRIPT" --Mw $Mw --Mb $Mb $MESH_FLAGS --sd "$save_dir"
            fi

            echo "[INFO] Running OpenFOAM simulation..."
            pushd "$case_path" > /dev/null
            if [[ "$RUN_MODE" == "parallel" ]]; then
                ./Allrun parallel "$NP"
            else
                ./Allrun
            fi
            popd > /dev/null

            echo "[DONE] Case finished: $fname"
            echo "[INFO] End: $(date)"
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
echo "# Total cases run: $((COUNTER - 1))"
echo "# Time finished: $(date)"
echo "# Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "##################################################"
