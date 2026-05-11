#!/bin/bash
# worker.sh — executed once per srun task
#
# SLURM variables set per-task:
#   $SLURM_NODEID    — node index (0..NODES-1)
#   $SLURM_LOCALID   — local GPU rank on this node (0..GPUS_PER_NODE-1)
#   $SLURM_NNODES    — total number of nodes
#
# Inherited from submit.sh:
#   $PYTHON, $WORK_DIR, $CSV_FILE, $FAILED_SYSTEMS, $GPUS_PER_NODE

cd $WORK_DIR

NODE_ID=$SLURM_NODEID
GPU_ID=$SLURM_LOCALID   # 0-indexed within the node

echo "[INFO] Node $NODE_ID | GPU $GPU_ID | sees CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

# --- Chunk the CSV across nodes and GPUs ---
TOTAL_ROWS=$(tail -n +2 "$CSV_FILE" | wc -l)
DATA_PER_NODE=$(( TOTAL_ROWS / SLURM_NNODES ))
DATA_PER_GPU=$(( DATA_PER_NODE / GPUS_PER_NODE ))

START_ROW=$(( NODE_ID * DATA_PER_NODE + GPU_ID * DATA_PER_GPU + 1 ))  # +1: tail is 1-indexed

echo "[INFO] Handling rows $START_ROW to $(( START_ROW + DATA_PER_GPU - 1 )) of $TOTAL_ROWS"

tail -n +2 "$CSV_FILE" | tail -n +"$START_ROW" | head -n "$DATA_PER_GPU" \
| while IFS=',' read -r PROTEIN LIGAND MODEL; do

    PROTEIN=$(echo "$PROTEIN" | xargs)
    LIGAND=$(echo "$LIGAND"   | xargs)
    MODEL=$(echo "$MODEL"     | xargs)
    MODEL_SHORT=${MODEL##*/}

    LIG_ARG="None"
    if [[ -n "$LIGAND" && "$LIGAND" != "None" && "$LIGAND" != "none" ]]; then
        LIG_ARG="$LIGAND"
    fi

    for RUN in 1 2 3; do
        LOGFILE="agent_logs/${PROTEIN}_${LIG_ARG}_${MODEL_SHORT}_run${RUN}.log"

        echo "------------------------------------------------------"
        echo "[INFO] Node=$NODE_ID GPU=$GPU_ID | $PROTEIN | $LIG_ARG | $MODEL_SHORT | Run $RUN/3"
        echo "[INFO] Output → $LOGFILE"
        echo "------------------------------------------------------"

        if ! $PYTHON -u main.py \
            --pdb_id  "$PROTEIN" \
            --ligand  "$LIG_ARG" \
            --model   "$MODEL" \
            --temp    300 \
            --duration 0.01 \
            > "$LOGFILE" 2>&1
        then
            echo "[WARN] Failed: $PROTEIN $LIG_ARG run $RUN" | tee -a "$FAILED_SYSTEMS"
        fi

        echo "[INFO] Finished run $RUN for $PROTEIN ($MODEL_SHORT)."
    done
done

echo "[INFO] Node $NODE_ID | GPU $GPU_ID | All rows complete."