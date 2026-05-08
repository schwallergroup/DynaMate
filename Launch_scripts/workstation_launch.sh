#!/bin/bash

echo "[INFO] Loading environment..."
source setup.sh

CSV_FILE="Launch/systems.csv"
CHECKPOINT=".csv_checkpoint"
FAILED_SYSTEMS="Launch/failed_systems.txt"

if [[ ! -f "$CSV_FILE" ]]; then
    echo "[ERROR] CSV file not found: $CSV_FILE"
    exit 1
fi

echo "[INFO] Starting sequential execution..."
echo "[INFO] Each system will be run 3 times."
echo

# Load starting line (default=2)
START_LINE=2
[[ -s "$CHECKPOINT" ]] && START_LINE=$(<"$CHECKPOINT")

# Iterate through CSV (skip header)
LINE_NUM=$((START_LINE - 1))
tail -n +$START_LINE "$CSV_FILE" | while IFS=',' read -r PROTEIN LIGAND MODEL; do
    ((LINE_NUM++))

    echo "[INFO] Processing CSV line $LINE_NUM"

    # Trim leading/trailing spaces
    PROTEIN=$(echo "$PROTEIN" | xargs)
    LIGAND=$(echo "$LIGAND" | xargs)
    MODEL=$(echo "$MODEL" | xargs)
    
    # Take only the last part of the model name after the last '/'
    MODEL_SHORT=${MODEL##*/}
    echo "Model='$MODEL'"

    # Handle missing ligand
    if [[ "$LIGAND" == "" || "$LIGAND" == "None" || "$LIGAND" == "none" ]]; then
        LIG_ARG="None"
    else
        LIG_ARG="$LIGAND"
    fi

    for RUN in 1 2 3; do
        
        LOGFILE="agent_logs/${PROTEIN}_${LIG_ARG}_${MODEL_SHORT}_run${RUN}.log"

        echo "------------------------------------------------------"
        echo "[INFO] Running: Protein=$PROTEIN | Ligand=$LIG_ARG | Model=$MODEL (Repeat $RUN/3)"
        echo "[INFO] Output → $LOGFILE"
        echo "------------------------------------------------------"

        # Run the Python script
        if ! python main.py \
            --pdb_id "$PROTEIN" \
            --ligand "$LIG_ARG" \
            --model "$MODEL" \
            --temp 300 \
            --duration 0.01 \
            > "$LOGFILE" 2>&1
        then
            echo "[WARN] DynaMate failed for $PROTEIN $LIG_ARG, run $RUN."  | tee -a "$FAILED_SYSTEMS"
        fi

        echo "[INFO] Finished run $RUN for $PROTEIN ($MODEL)."

    done

    # Save progress
    echo $((LINE_NUM + 1)) > "$CHECKPOINT"

done

echo "[INFO] All runs completed!"
