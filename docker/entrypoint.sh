#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate dynamate
source /opt/software/ambertools25/amber.sh

# First arg is a flag (e.g. --model) -> run the agent; else exec it (e.g. /bin/bash)
if [ "${1:0:1}" = '-' ]; then
    set -- python /app/main.py "$@"
fi
exec "$@"