#!/bin/bash
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 input_gro npt_cpt_file log_file [ligand_name]"
    exit 1
fi

GMX='gmx'
INPUT_GRO="$1"
NPT_CPT_FILE="$2"
LOG_FILE="$3"

> $LOG_FILE 

# Optional fourth argument
if [ "$#" -ge 4 ]; then
    LIGNAME="$4"
else
    LIGNAME=""
fi

#------- EDIT MD.MDP FILE ------------

# if no ions are present, update md.mdp such that Water_and_ions group is changed to Water only
# if ligand, Protein Water becomes Protein_ligand Water
if ! grep -q "Cl-" index.ndx && ! grep -q "Na+" index.ndx; then
	md_file="md.mdp"
    original="Protein Water_and_ions"
	
	if [ -n "$LIGNAME" ]; then
		replacement="Protein_$LIGNAME Water"
		if grep "Protein_$LIGNAME" index.ndx; then
			echo "Protein_$LIGNAME already in index.ndx"
		else
			echo -e "1 | 13\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx
		fi
	else
		replacement="Protein Water"
	fi

	if grep "$original" "$md_file"; then
		sed -i "s|$original|$replacement|" "$md_file"
		echo "$replacement added successfully to tc-grps group in $md_file."
	else
		echo "tc-grps line was not found in $md_file."
	fi
else
	echo "Ions present. Keeping Water_and_ions group."
fi

# if ions are present, but ligand exists, update md.mdp accordingly
if grep -q "Cl-" index.ndx || grep -q "Na+" index.ndx; then
    if [ -n "$LIGNAME" ]; then	
        md_file="md.mdp" 

        original="Protein Water_and_ions"
        replacement="Protein_$LIGNAME Water_and_ions"
        
        if grep "Protein_$LIGNAME" index.ndx; then
            echo "Protein_$LIGNAME already in index.ndx"
        else
            echo -e "1 | 13\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx
        fi

        if grep "$original" "$md_file"; then
            sed -i "s|$original|$replacement|" "$md_file"
            echo "Protein_$LIGNAME added successfully to tc-grps group in $md_file."
        else
            echo "tc-grps line "Protein Water_and_ions" was not found in $md_file."
        fi
    fi
fi

#------- PRODUCTION MD ------------
if ! ls md.gro 1> /dev/null 2>&1; then
    $GMX grompp -f md.mdp -c $INPUT_GRO -t $NPT_CPT_FILE -p topol.top -n index.ndx -o md.tpr >> $LOG_FILE 2>&1
    echo "y" | $GMX mdrun -v -deffnm md >> $LOG_FILE 2>&1
else
    echo "'md.gro' already exists. Skipping production MD."
fi