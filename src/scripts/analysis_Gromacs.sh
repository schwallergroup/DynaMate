#!/bin/bash
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 input_xtc [ligand_name]"
    exit 1
fi

GMX='gmx'
INPUT_XTC="$1"
LOG_FILE="$2"
> $LOG_FILE 
FILENAME="${INPUT_XTC%.*}"

# Optional fourth argument
if [ "$#" -ge 2 ]; then
    LIGNAME="$2"
else
    LIGNAME=""
fi

#------ ANALYSIS ------------------

# Remove PBC
echo -e "Protein \n System" | $GMX trjconv -s $FILENAME.tpr -f $FILENAME.xtc -o "$FILENAME"_noPBC.xtc -pbc mol -center >> $LOG_FILE 2>&1
# RMSD to initial structure
echo -e "Backbone \n Backbone" | $GMX rms -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -o rmsd.xvg -tu ns >> $LOG_FILE 2>&1
# RMSD to crystal structure
echo -e "Backbone \n Backbone" | $GMX rms -s em.tpr -f "$FILENAME"_noPBC.xtc -o rmsd_xtal.xvg -tu ns >> $LOG_FILE 2>&1
# RMSF
echo -e "C-alpha" | $GMX rmsf -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -o rmsf.xvg -res >> $LOG_FILE 2>&1
# Radius of gyration
echo -e "Protein" | $GMX gyrate -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -o gyrate.xvg >> $LOG_FILE 2>&1
# Hydrogen bonds
echo -e "MainChain+H \n MainChain+H" | $GMX hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_mainchain.xvg >> $LOG_FILE 2>&1
echo -e "SideChain \n SideChain" | $GMX hbond -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -tu ns -num hbnum_sidechain.xvg >> $LOG_FILE 2>&1
echo -e "Protein \n Water" | $GMX hbond -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -tu ns -num hbnum_prot_wat.xvg >> $LOG_FILE 2>&1

if [ -n "$LIGNAME" ]; then
    # Ligand-Protein hydrogen bonds
    echo -e "1 \n 13" | $GMX hbond -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -tu ns -num hbnum_prot_lig.xvg >> $LOG_FILE 2>&1
fi

# echo -e "Protein" \n "System" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center >> log 2>&1
# # RMSD to initial structure
# echo -e "Backbone" \n "Backbone" | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns >> log.log 2>&1
# # RMSD to crystal structure
# echo -e "Backbone" \n "Backbone" | $GMX rms -s em.tpr -f "$FILENAME"_noPBC.xtc -o rmsd_xtal.xvg -tu ns >> $LOG_FILE 2>&1
# echo -e "Backbone" \n "Backbone" | gmx rms -s em.tpr -f md_noPBC.xtc -o rmsd_xtal.xvg -tu ns
# # RMSF
# echo -e "C-alpha" | $GMX rmsf -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -o rmsf.xvg -res >> $LOG_FILE 2>&1
# # Radius of gyration
# echo -e "Protein" | $GMX gyrate -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -o gyrate.xvg >> $LOG_FILE 2>&1
# # Hydrogen bonds
# #echo -e "MainChain+H \n MainChain+H" | $GMX hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_mainchain.xvg
# echo -e "SideChain \n SideChain" | gmx hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_sidechain.xvg
# echo -e "Protein \n Water" | $GMX hbond -s $FILENAME.tpr -f "$FILENAME"_noPBC.xtc -tu ns -num hbnum_prot_wat.xvg >> $LOG_FILE 2>&1

# echo -e "Protein \n Water" | $GMX hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_prot_wat.xvg

# echo -e "MainChain+H \n MainChain+H" | gmx hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_mainchain.xvg