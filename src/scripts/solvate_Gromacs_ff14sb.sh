#!/bin/bash
# location input mdp files
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDP_FILES="$SCRIPT_DIR/../../sandbox/experiments/inputs_unconstrained_MD"

# get input file
if [ $# -ne 1 ]; then
    echo "Usage: $0 input.pdb"
    exit 1
fi

PDBFILE=$(basename "$1" .pdb)

gmx editconf -f "$PDBFILE" -o box.gro -c -d 1.2 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f "$MDP_FILES"/ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo "13" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral