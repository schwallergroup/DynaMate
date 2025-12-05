#!/bin/bash
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 input_gro npt_cpt_file"
    exit 1
fi

GMX='gmx'
INPUT_GRO="$1"
NPT_CPT_FILE="$2"
LOG_FILE="$3"

> $LOG_FILE 

#------- PRODUCTION MD ------------
if ! ls md.gro 1> /dev/null 2>&1; then
    $GMX grompp -f md.mdp -c $INPUT_GRO -t $NPT_CPT_FILE -p topol.top -n index.ndx -o md.tpr >> $LOG_FILE 2>&1
    echo "y" | $GMX mdrun -v -deffnm md >> $LOG_FILE 2>&1
else
    echo "'md.gro' already exists. Skipping production MD."
fi