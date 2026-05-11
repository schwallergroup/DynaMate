#!/bin/bash
#SBATCH --account=a131
#SBATCH --partition=debug
#SBATCH --job-name=dynamate_launch
#SBATCH --output=slurm_logs/launch_%j_%N.out
#SBATCH --error=slurm_logs/launch_%j_%N.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --time=00:10:00

set -e

mkdir -p agent_logs slurm_logs

# --- Exported config inherited by all worker tasks ---
export PYTHON=/capstor/store/cscs/swissai/a131/bnaida/conda-envs/dynamate/bin/python
export WORK_DIR=/capstor/store/cscs/swissai/a131/cassandra/projects/DynaMate
export CSV_FILE="Launch_scripts/systems.csv"
export FAILED_SYSTEMS="Launch_scripts/failed_systems.txt"
export GPUS_PER_NODE=4

cd $WORK_DIR

if [[ ! -f "$CSV_FILE" ]]; then
    echo "[ERROR] CSV file not found: $CSV_FILE"
    exit 1
fi

srun --ntasks=4 --ntasks-per-node=4 --gpus-per-task=1 \
    --environment="/capstor/store/cscs/swissai/a131/bnaida/ce-images/dynamate/md_tools_v0.1.0.toml" \
    bash /capstor/store/cscs/swissai/a131/cassandra/projects/DynaMate/Launch_scripts/worker.sh
