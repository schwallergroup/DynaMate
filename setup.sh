#!/bin/bash

# load your environment
source $(conda info --base)/bin/activate dynamate

# Activate GROMACS environment
if [ -f /usr/local/gromacs/bin/GMXRC ]; then
    source /usr/local/gromacs/bin/GMXRC
else
    echo "Error: /usr/local/gromacs/bin/GMXRC not found!"
    return 1 2>/dev/null || exit 1
fi

echo "PYTHONPATH set to: $PYTHONPATH"
echo "Environment setup complete!"

