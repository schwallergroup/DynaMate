#!/bin/bash
#SBATCH --account=a131
#SBATCH --partition=debug
#SBATCH --job-name=dynamate
#SBATCH --output=$SCRATCH/logs/dynamate/dynamate_%j_%N.out
#SBATCH --error=$SCRATCH/logs/dynamate/dynamate_%j_%N.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --environment="/capstor/store/cscs/swissai/a131/bnaida/ce-images/dynamate/md_tools_v0.1.0.toml"

set -e
mkdir -p agent_logs

PYTHON=/capstor/store/cscs/swissai/a131/bnaida/conda-envs/dynamate/bin/python
cd /capstor/store/cscs/swissai/a131/$USER/projects/DynaMate

CSV_FILE="Launch/systems.csv"
FAILED_SYSTEMS="Launch/failed_systems.txt"

if [[ ! -f "$CSV_FILE" ]]; then
    echo "[ERROR] CSV file not found: $CSV_FILE"
    exit 1
fi

NODE_ID=$SLURM_NODEID
GPU_ID=$SLURM_LOCALID
GPUS_PER_NODE=4

TOTAL_ROWS=$(tail -n +2 "$CSV_FILE" | wc -l)
DATA_PER_NODE=$(( TOTAL_ROWS / SLURM_NNODES ))
DATA_PER_GPU=$(( DATA_PER_NODE / GPUS_PER_NODE ))

START_ROW=$(( NODE_ID * DATA_PER_NODE + GPU_ID * DATA_PER_GPU ))
END_ROW=$(( START_ROW + DATA_PER_GPU - 1 ))

LINE_NUM=$((START_ROW - 1))
tail -n +2 "$CSV_FILE" | tail -n +$START_ROW | head -n $DATA_PER_GPU \
| while IFS=',' read -r PROTEIN LIGAND MODEL; do
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
        if ! $PYTHON -u main.py \
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
