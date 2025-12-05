#!/bin/bash
# Usage: ./run_tleap.sh input.pdb

if [ $# -ne 4 ]; then
    echo "Usage: $0 sandbox_dir complex_pdb ligand_name prepi_file"
    exit 1
fi

SANDBOX_DIR=$1 # PDBFILE already has sandbox path
PDBFILE=$2
LIGNAME=$3
PREPI_FILE=$4

# Create tleap input file
cat > leap.in << EOF
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p

# Map PDB atom names to template atom names
addPdbAtomMap { { "CH3"  "C" } { "HH31" "H1" } { "HH32" "H2" } { "HH33" "H3" } { "CL1" "Cl1" } { "CL2" "Cl2" } }

# Load ligand parameters
loadamberprep ${SANDBOX_DIR}/${PREPI_FILE}
loadamberparams ${SANDBOX_DIR}/${LIGNAME}.frcmod

# PDBFILE already has sandbox path
mol = loadpdb ${PDBFILE}
solvatebox mol TIP3PBOX 16
addions mol Cl- 0 # Neutralize system
addions mol Na+ 0 # Neutralize system

saveamberparm mol ${SANDBOX_DIR}/complex.prmtop ${SANDBOX_DIR}/complex.inpcrd
savepdb mol ${SANDBOX_DIR}/complex_tleap.pdb

quit
EOF

# Run tleap
tleap -f leap.in