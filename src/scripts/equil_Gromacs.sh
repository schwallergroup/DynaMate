#!/bin/bash
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 sandbox_dir input_gro [ligand_name] [ligand_file] [ligand_gro]"
    exit 1
fi

GMX='gmx'
SANDBOX_DIR="$1"
INPUT_GRO="$2"
LOG_FILE="$3"

> $LOG_FILE 
echo "Starting GROMACS Equilibration Log" >> $LOG_FILE 2>&1

# Optional fourth argument
if [ "$#" -ge 5 ]; then
    LIGNAME="$4"
	LIGFILE="$5"
	LIGGRO="$6"
else
    LIGNAME=""
	LIGFILE=""
	LIGGRO=""
fi

#------- ENERGY MINIMISATION ------------
if ! ls em.gro 1> /dev/null 2>&1; then
	$GMX grompp -f em.mdp -c $INPUT_GRO -p topol.top -o em.tpr >> $LOG_FILE 2>&1
	$GMX mdrun -v -deffnm em >> $LOG_FILE 2>&1

	if [ -f em.gro ]; then
	        echo "'em.gro' created"
		echo "11 0" | $GMX energy -f em.edr -o potential.xvg
	else
		echo "Error: Failed to create 'em.gro'"
		exit 1
	fi
else
	echo "'em.gro' already exists. Skipping energy minimisation."
fi

#----------Create posres files-----------

# Step 1: count NME residues (6 per chain)
n_nme=$(grep -c "NME" em.gro)
chains=$((n_nme / 6))
echo "Detected $chains chains based on NME residues."

if ! ls index.ndx 1> /dev/null 2>&1; then
	echo "q" | $GMX make_ndx -f em.gro -o index.ndx
fi

# Step 2: if 1 chain, create posre.itp file
if [ "$chains" -eq "1" ]; then
    if ! ls posre.itp 1> /dev/null 2>&1;then
            echo ""Protein-H"" | $GMX genrestr -f em.gro -n index.ndx -o posre.itp -fc 1000 1000 1000 
    fi
else
# Step 2: if more than 1 chain, create posre_chain(i).itp files

    # Get last residue number of each chain, by extracting residue numbers of all NME lines and selecting every 6th occurrence (last line of each NME block)
    nme_residues=($(awk '/NME/ {resnum = substr($0,1,5); gsub(/ /,"",resnum); print resnum}' em.gro))
    chain_end_residues=()
    for ((i=5; i<${#nme_residues[@]}; i+=6)); do
        chain_end_residues+=("${nme_residues[i]}")
    done
    echo "Chain end residues: ${chain_end_residues[@]}"

    # Step 3: compute chain residue ranges
    start=1
    end=0
    ranges=()
    for end_residue in "${chain_end_residues[@]}"; do
        end=$((end + end_residue))
        ranges+=("${start}-${end}")
        start=$((end + 1))
    done

    echo "Residue ranges per chain: ${ranges[@]}"

    # Step 4: create index groups for each chain
    # Start from existing index.ndx or create new
    if [ ! -f index.ndx ]; then
        $GMX make_ndx -f em.gro -o index.ndx << EOF
        q
EOF
    fi

    # Step 5: add groups per chain
    i=1
    for range in "${ranges[@]}"; do
        echo "Creating group for residues $range..."
        echo -e "ri $range\n2 & \"r_${range}\"\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx
        ((i++))
    done

    # Step 6: generate posre.itp for each chain
    i=1
    for range in "${ranges[@]}"; do
        if [ ! -f posre_chain${i}.itp ]; then
            echo "Generating position restraints for chain${i}"
            group_name="Protein-H_&_r_${range}"
            echo "$group_name" | $GMX genrestr -f em.gro -n index.ndx -o "posre_chain${i}.itp" -fc 1000 1000 1000
            
            if [ "$i" -gt "1" ]; then
                #Adjust atom indices so first = 1
                awk '
                /^\[ position_restraints \]/ { in_section=1; first_index=0; shift=0; print; next }
                /^\[/ && !/\[ position_restraints \]/ { in_section=0 }
                {
                    if (in_section && /^[0-9]/) {
                        if (first_index == 0) {
                            first_index = $1
                            if (first_index != 1) shift = first_index - 2
                        }
                        $1 = $1 - shift
                    }
                    print
                }
                ' "posre_chain${i}.itp" > "posre_chain${i}_renum.itp" && mv "posre_chain${i}_renum.itp" "posre_chain${i}.itp"
            fi

            ((i++))
        fi
    done

fi

if [ -n "$LIGNAME" ]; then
	#obabel $LIGFILE -O $LIGGRO
	if ! ls "posre_$LIGNAME.itp" 1> /dev/null 2>&1;then
			echo -e "0 & ! a H*\nq" | $GMX make_ndx -f $LIGGRO -o "index_$LIGNAME.ndx"
		echo "3" | $GMX genrestr -f $LIGGRO -n "index_$LIGNAME.ndx" -o "posre_$LIGNAME.itp" -fc 1000 1000 1000 
	fi
fi

# #-------- UPDATE TOPOLOGY FILE -----------
# TMP_FILE="topol.tmp"
	
# # Check if block already exists
# if [ -n "$LIGNAME" ]; then
# if grep -q "^; posre_$LIGNAME.itp" topol.top; then
# 	echo "Position restraint block already present in topol.top. No changes made."
# else	
# 	awk -v ligname="$LIGNAME" '
# 	BEGIN { in_dihedrals=0; inserted=0 }
# 	{
# 		if ($0 ~ /^\[ dihedrals \]/ && in_dihedrals==0) {
# 	        in_dihedrals=1
#     	}	
#     	else if (in_dihedrals==1 && $0 ~ /^$/) {
#         	# Insert blank line + restraint block
#         	print ""
#         	print "; Include Position restraint file"
#         	print "#ifdef POSRES"
#         	print "#include \"posre_" ligname ".itp\""
#         	print "#endif"
#         	inserted=1
#         	in_dihedrals=0
#     	}
#     	print
# 	}
# 	' topol.top > "$TMP_FILE"

# 	# Replace the original file with the modified one
# 	cp "topol.top" topol_old.top
# 	mv "$TMP_FILE" "topol.top"

# 	echo "Updated topol.top with position restraint block."
# fi
# fi

# Create group Water_and_ions if not exists
if grep -q "Water_and_ions" index.ndx; then
	echo "Group Water_and_ions already exists in index.ndx"	
else
	if grep -q "Cl-" index.ndx; then
		echo -e '"WAT" | "Cl-" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx
		sed -i 's/Water_Cl-/Water_and_ions/g' index.ndx
		echo "Group Water_and_ions created in index.ndx"
	fi
	if grep -q "Na+" index.ndx; then
		echo -e '"WAT" | "Na+" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx
		sed -i 's/Water_Na+/Water_and_ions/g' index.ndx
		echo "Group Water_and_ions created in index.ndx"
	fi
fi

#------ Update Water_and_ions is no ions present -----
if ! grep -q "Cl-" index.ndx && ! grep -q "Na+" index.ndx; then
	nvt_file="nvt.mdp"
	npt_file="npt.mdp"
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
	
    if grep "$original" "$nvt_file"; then
		sed -i "s|$original|$replacement|" "$nvt_file"
		echo "$replacement added successfully to tc-grps group in $nvt_file."
	else
		echo "tc-grps line was not found in $nvt_file."
	fi

	if grep "$original" "$npt_file"; then
		sed -i "s|$original|$replacement|" "$npt_file"
		echo "$replacement added successfully to tc-grps group in $npt_file."
	else
		echo "tc-grps line was not found in $npt_file."
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
#-------- UPDATE TEMP GROUPS NPT, NVT, MD.MDP FILES -----

#cp "${MDP_FILES}/nvt.mdp" .
#cp "${MDP_FILES}/npt.mdp" .
#cp "${MDP_FILES}/md.mdp" .
#echo "mdp files copied"

if [ -n "$LIGNAME" ]; then	
	nvt_file="nvt.mdp"
	npt_file="npt.mdp"
	md_file="md.mdp"

	original="Protein Water_and_ions"
	replacement="Protein_$LIGNAME Water_and_ions"
	
	if grep "Protein_$LIGNAME" index.ndx; then
		echo "Protein_$LIGNAME already in index.ndx"
	else
		echo -e "1 | 13\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx
	fi

        if grep "$original" "$nvt_file"; then
		sed -i "s|$original|$replacement|" "$nvt_file"
		echo "Protein_$LIGNAME added successfully to tc-grps group in $nvt_file."
	else
		echo "tc-grps line was not found in $nvt_file."
	fi

	if grep "$original" "$npt_file"; then
		sed -i "s|$original|$replacement|" "$npt_file"
		echo "Protein_$LIGNAME added successfully to tc-grps group in $npt_file."
	else
		echo "tc-grps line was not found in $npt_file."
	fi

	if grep "$original" "$md_file"; then
		sed -i "s|$original|$replacement|" "$md_file"
		echo "Protein_$LIGNAME added successfully to tc-grps group in $md_file."
	else
		echo "tc-grps line was not found in $md_file."
	fi
fi

#--------------- NVT --------------------

if ! ls nvt.gro 1> /dev/null 2>&1; then

	$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 2
	$GMX mdrun -v -deffnm nvt
        
	if [ -f nvt.gro ]; then
	        echo "'nvt.gro' created"
		echo -e "Temperature \n 0" | $GMX energy -f nvt.edr -o temperature.xvg
	else
		echo "Error: Failed to create 'nvt.gro'"
		exit 1
	fi
else
	echo "'nvt.gro' already exists. Skipping NVT."
fi

#--------------- NPT --------------------
if ! ls npt.gro 1> /dev/null 2>&1; then
	$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -n index.ndx -maxwarn 2
	$GMX mdrun -v -deffnm npt

	if [ -f npt.gro ]; then
	    echo "'npt.gro' created"
		echo -e "Pressure \n 0" | $GMX energy -f npt.edr -o pressure.xvg
		echo -e "Density \n 0" | $GMX energy -f npt.edr -o density.xvg
	else
		echo "Error: Failed to create 'npt.gro'"
		exit 1
	fi
else
	echo "'npt.gro' already exists. Skipping NPT."
fi