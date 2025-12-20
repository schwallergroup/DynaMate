#!/bin/bash
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 sandbox_dir input_gro log_file [ligand_name] [ligand_file] [ligand_gro]"
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
		echo "11 0" | $GMX energy -f em.edr -o potential.xvg >> $LOG_FILE 2>&1
	else
		echo "Error: Failed to create 'em.gro'" >> $LOG_FILE 2>&1
		exit 1
	fi
else
	echo "'em.gro' already exists. Skipping energy minimisation." >> $LOG_FILE 2>&1
fi

#----------Create posres files-----------

# Step 1: count NME residues (6 per chain)
n_nme=$(grep -c "NME" em.gro)
chains=$((n_nme / 6))
echo "Detected $chains chains based on NME residues." >> $LOG_FILE 2>&1

if ! ls index.ndx 1> /dev/null 2>&1; then
	echo "q" | $GMX make_ndx -f em.gro -o index.ndx >> $LOG_FILE 2>&1
fi

# Step 2: if 1 chain, create posre.itp file
if [ "$chains" -eq "1" ]; then
    if ! ls posre.itp 1> /dev/null 2>&1;then
            echo ""Protein-H"" | $GMX genrestr -f em.gro -n index.ndx -o posre.itp -fc 1000 1000 1000 >> $LOG_FILE 2>&1
    fi
else
# Step 2: if more than 1 chain, create posre_chain(i).itp files

    # Get last residue number of each chain, by extracting residue numbers of all NME lines and selecting every 6th occurrence (last line of each NME block)
    nme_residues=($(awk '/NME/ {resnum = substr($0,1,5); gsub(/ /,"",resnum); print resnum}' em.gro))
    chain_end_residues=()
    for ((i=5; i<${#nme_residues[@]}; i+=6)); do
        chain_end_residues+=("${nme_residues[i]}")
    done
    echo "Chain end residues: ${chain_end_residues[@]}" >> $LOG_FILE 2>&1

    # Step 3: compute chain residue ranges
    start=1
    end=0
    ranges=()
    for end_residue in "${chain_end_residues[@]}"; do
        end=$((end + end_residue))
        ranges+=("${start}-${end}")
        start=$((end + 1))
    done

    echo "Residue ranges per chain: ${ranges[@]}" >> $LOG_FILE 2>&1

    # Step 4: create index groups for each chain
    # Start from existing index.ndx or create new
    if [ ! -f index.ndx ]; then
        $GMX make_ndx -f em.gro -o index.ndx << EOF
        q
EOF
    fi

    # Step 5: add groups per chain
if grep -E "system1 +2" topol.top; then #special case for two identical chains named system1
	echo "You have two identical chains named system1, therefore only one position restraint file for the first chain will be created." >> $LOG_FILE 2>&1
	echo "Creating group for residues ${ranges[0]}" >> $LOG_FILE 2>&1
	echo -e "ri ${ranges[0]}\n2 & \"r_${ranges[0]}\"\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
else 
	i=1
    for range in "${ranges[@]}"; do
        echo "Creating group for residues $range..." >> $LOG_FILE 2>&1
        echo -e "ri $range\n2 & \"r_${range}\"\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
        ((i++))
    done
fi

    # Step 6: generate posre.itp for each chain
if grep -E "system1 +2" topol.top; then #special case for two identical chains named system1
	group_name="Protein-H_&_r_${ranges[0]}"
	echo "$group_name" | $GMX genrestr -f em.gro -n index.ndx -o "posre.itp" -fc 1000 1000 1000 >> $LOG_FILE 2>&1
else 
    i=1
    for range in "${ranges[@]}"; do
        if [ ! -f posre_chain${i}.itp ]; then
            echo "Generating position restraints for chain${i}" >> $LOG_FILE 2>&1
            group_name="Protein-H_&_r_${range}"
            echo "$group_name" | $GMX genrestr -f em.gro -n index.ndx -o "posre_chain${i}.itp" -fc 1000 1000 1000 >> $LOG_FILE 2>&1
            
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
fi

if [ -n "$LIGNAME" ]; then
	#obabel $LIGFILE -O $LIGGRO
	if ! ls "posre_$LIGNAME.itp" 1> /dev/null 2>&1;then
		echo "Generating position restraints for ligand $LIGNAME" >> $LOG_FILE 2>&1
		echo -e "0 & ! a H* \n q" | $GMX make_ndx -f $LIGGRO -o "index_$LIGNAME.ndx" >> $LOG_FILE 2>&1
		echo "3" | $GMX genrestr -f $LIGGRO -n "index_$LIGNAME.ndx" -o "posre_$LIGNAME.itp" -fc 1000 1000 1000 >> $LOG_FILE 2>&1
	fi
fi

# Create group Water_and_ions if not exists
if grep -q "Water_and_ions" index.ndx; then
	echo "Group Water_and_ions already exists in index.ndx"	
else
	if grep -q "Cl-" index.ndx; then
		echo -e '"WAT" | "Cl-" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
		sed -i 's/Water_Cl-/Water_and_ions/g' index.ndx
		echo "Group Water_and_ions created in index.ndx" >> $LOG_FILE 2>&1
	fi
	if grep -q "Na+" index.ndx; then
		echo -e '"WAT" | "Na+" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
		sed -i 's/Water_Na+/Water_and_ions/g' index.ndx
		echo "Group Water_and_ions created in index.ndx" >> $LOG_FILE 2>&1
	fi
fi

#------ Update Water_and_ions is no ions present -----
if ! grep -q "Cl-" index.ndx && ! grep -q "Na+" index.ndx; then
	nvt_file="nvt.mdp"
	npt_file="npt.mdp"
	original="Protein Water_and_ions"
	
	if [ -n "$LIGNAME" ]; then
		replacement="Protein_$LIGNAME Water"
		if grep "Protein_$LIGNAME" index.ndx; then
			echo "Protein_$LIGNAME already in index.ndx" >> $LOG_FILE 2>&1
		else
			echo -e "1 | 13\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
		fi
	else
		replacement="Protein Water"
	fi
	
    if grep "$original" "$nvt_file"; then
		sed -i "s|$original|$replacement|" "$nvt_file"
		echo "$replacement added successfully to tc-grps group in $nvt_file." >> $LOG_FILE 2>&1
	else
		echo "tc-grps line was not found in $nvt_file." >> $LOG_FILE 2>&1
	fi

	if grep "$original" "$npt_file"; then
		sed -i "s|$original|$replacement|" "$npt_file"
		echo "$replacement added successfully to tc-grps group in $npt_file." >> $LOG_FILE 2>&1
	else
		echo "tc-grps line was not found in $npt_file." >> $LOG_FILE 2>&1
	fi

else
	echo "Ions present. Keeping Water_and_ions group." >> $LOG_FILE 2>&1
fi
#-------- UPDATE TEMP GROUPS NPT, NVT, MD.MDP FILES if there are ions -----

if grep -q "Cl-" index.ndx || grep -q "Na+" index.ndx; then
	if [ -n "$LIGNAME" ]; then	
		nvt_file="nvt.mdp"
		npt_file="npt.mdp"

		original="Protein Water_and_ions"
		replacement="Protein_$LIGNAME Water_and_ions"
		
		if grep "Protein_$LIGNAME" index.ndx; then
			echo "Protein_$LIGNAME already in index.ndx" >> $LOG_FILE 2>&1
		else
			echo -e "1 | 13\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
		fi

		if grep "$original" "$nvt_file"; then
			sed -i "s|$original|$replacement|" "$nvt_file"
			echo "Protein_$LIGNAME added successfully to tc-grps group in $nvt_file." >> $LOG_FILE 2>&1
		else
			echo "tc-grps line "Protein Water_and_ions" was not found in $nvt_file." >> $LOG_FILE 2>&1
		fi

		if grep "$original" "$npt_file"; then
			sed -i "s|$original|$replacement|" "$npt_file"
			echo "Protein_$LIGNAME added successfully to tc-grps group in $npt_file." >> $LOG_FILE 2>&1
		else
			echo "tc-grps line "Protein Water_and_ions" was not found in $npt_file." >> $LOG_FILE 2>&1
		fi
	fi
fi
#--------------- NVT --------------------

if ! ls nvt.gro 1> /dev/null 2>&1; then

	$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 2 >> $LOG_FILE 2>&1
	$GMX mdrun -v -deffnm nvt >> $LOG_FILE 2>&1
        
	if [ -f nvt.gro ]; then
	    echo "'nvt.gro' created" >> $LOG_FILE 2>&1
		echo -e "Temperature \n 0" | $GMX energy -f nvt.edr -o temperature.xvg >> $LOG_FILE 2>&1
	else
		echo "Error: Failed to create 'nvt.gro'" >> $LOG_FILE 2>&1
		exit 1
	fi
else
	echo "'nvt.gro' already exists. Skipping NVT." >> $LOG_FILE 2>&1
fi

#--------------- NPT --------------------
if ! ls npt.gro 1> /dev/null 2>&1; then
	$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -n index.ndx -maxwarn 2 >> $LOG_FILE 2>&1
	$GMX mdrun -v -deffnm npt >> $LOG_FILE 2>&1

	if [ -f npt.gro ]; then
	    echo "'npt.gro' created" >> $LOG_FILE 2>&1
		echo -e "Pressure \n 0" | $GMX energy -f npt.edr -o pressure.xvg >> $LOG_FILE 2>&1
		echo -e "Density \n 0" | $GMX energy -f npt.edr -o density.xvg >> $LOG_FILE 2>&1
	else
		echo "Error: Failed to create 'npt.gro'" >> $LOG_FILE 2>&1
		exit 1
	fi
else
	echo "'npt.gro' already exists. Skipping NPT." >> $LOG_FILE 2>&1
fi