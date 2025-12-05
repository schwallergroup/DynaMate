#!/bin/bash
# Usage: ./run_tleap.sh input.pdb

if [ $# -ne 3 ]; then
    echo "Usage: $0 sandbox_dir input.pdb pdb_id"
    exit 1
fi

SANDBOX_DIR=$1
PDBFILE=$2
PDB_ID=$3

# Create tleap input file
cat > leap.in << EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p

# Map PDB atom names to template atom names
addPdbAtomMap { { "CH3"  "C" } { "HH31" "H1" } { "HH32" "H2" } { "HH33" "H3" } }

mol = loadpdb ${SANDBOX_DIR}/${PDBFILE}
solvatebox mol TIP3PBOX 16
addions mol Cl- 0 # Neutralize system
addions mol Na+ 0 # Neutralize system

saveamberparm mol ${SANDBOX_DIR}/${PDB_ID}.prmtop ${SANDBOX_DIR}/${PDB_ID}.inpcrd
savepdb mol ${SANDBOX_DIR}/${PDB_ID}_tleap.pdb

quit
EOF

# Run tleap
tleap -f leap.in

# if output files not in directory, echo that tleap failed, otyherwise confirm success
ls ${SANDBOX_DIR}/${PDB_ID}.prmtop ${SANDBOX_DIR}/${PDB_ID}.inpcrd ${SANDBOX_DIR}/${PDB_ID}_tleap.pdb > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "tleap failed to generate output files."
else
    echo "${PDBFILE} processed. Generated ${PDB_ID}.prmtop, ${PDB_ID}.inpcrd, and ${PDB_ID}_tleap.pdb."
fi