#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate dynagent
source /opt/software/ambertools25/amber.sh
exec "$@"