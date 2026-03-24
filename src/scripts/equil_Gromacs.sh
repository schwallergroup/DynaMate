#!/bin/bash
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 sandbox_dir input_gro log_file [ligand_name] [ligand_file] [ligand_gro]"
    exit 1
fi

# Detect available GROMACS executable
if command -v gmx >/dev/null 2>&1; then
    GMX='gmx'
elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX='gmx_mpi'
elif command -v gmx_d >/dev/null 2>&1; then
    GMX='gmx_d'
else
    echo "Warning: gmx not found in PATH. Please install GROMACS." >&2
    GMX='gmx'
fi
SANDBOX_DIR="$1"
INPUT_GRO="$2"
LOG_FILE="$3"

> $LOG_FILE 
echo "Starting GROMACS Equilibration Log" >> $LOG_FILE 2>&1

# Optional argument
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
if ! [ -f em.gro ]; then
	$GMX grompp -f em.mdp -c $INPUT_GRO -p topol.top -o em.tpr >> $LOG_FILE 2>&1
	$GMX mdrun -v -deffnm em >> $LOG_FILE 2>&1
fi

if [ -f em.gro ]; then
	echo "'em.gro' created"
#	echo "11 0" | $GMX energy -f em.edr -o potential.xvg >> $LOG_FILE 2>&1
else
	echo "Error: Failed to create 'em.gro'" >> $LOG_FILE 2>&1
	exit 1
fi

#----------Create posres files-----------

# ---------------------------------------------------------------------------
# 1. Parse [ molecules ] section of topol.top
#    Build two arrays:
#      MOL_TYPES   – unique system* names in first-seen order
#      CHAIN_SEQ   – expanded ordered list of molecule types, one entry per
#                   physical chain (system1 2 → system1 system1)
# ---------------------------------------------------------------------------
in_molecules=0
declare -a CHAIN_SEQ=()       # physical chain sequence
declare -a MOL_TYPES=()       # unique molecule types, ordered by first appearance
declare -A MOL_SEEN=()        # track which types we've already recorded
TOPOL_FILE="${SANDBOX_DIR%/}/topol.top"

awk_output=$(awk '
    /^\[[[:space:]]*molecules[[:space:]]*\]/ { in_mol=1; next }
    in_mol && /^\[/ { exit }
    in_mol && /^[[:space:]]*(;|$)/ { next }
    in_mol {
        mol=$1; cnt=$2
        if (mol ~ /^[Ss]ystem[0-9]*$/) {
            if (!(mol in seen)) { print ":unique:" mol; seen[mol]=1 }
            for (i=0; i<cnt; i++) print mol
        }
    }
' "$TOPOL_FILE")
 
# Check the section was found at all
if [[ -z "$awk_output" ]]; then
    # Distinguish "section missing" from "section present but no system* entries"
    if ! awk '/^\[[[:space:]]*molecules[[:space:]]*\]/{found=1} END{exit !found}' "$TOPOL_FILE"; then
        echo "ERROR: [ molecules ] section not found in $TOPOL_FILE" | tee -a "$LOG_FILE" >&2
        echo "  File: $(wc -c < "$TOPOL_FILE") bytes, $(wc -l < "$TOPOL_FILE") lines" >> "$LOG_FILE"
    else
        echo "ERROR: No protein chains (system*) found in [ molecules ] of $TOPOL_FILE" | tee -a "$LOG_FILE" >&2
        echo "  Contents of [ molecules ]:" >> "$LOG_FILE"
        awk '/^\[[[:space:]]*molecules[[:space:]]*\]/{p=1;next} p && /^\[/{exit} p' \
            "$TOPOL_FILE" | head -10 >> "$LOG_FILE"
    fi
    exit 1
fi
 
# Populate CHAIN_SEQ and MOL_TYPES from awk output
declare -a CHAIN_SEQ=()
declare -a MOL_TYPES=()
while IFS= read -r token; do
    if [[ "$token" == :unique:* ]]; then
        MOL_TYPES+=("${token#:unique:}")
    else
        CHAIN_SEQ+=("$token")
    fi
done <<< "$awk_output"
 
total_chains=${#CHAIN_SEQ[@]}
echo "Molecule type order from topol.top: ${CHAIN_SEQ[*]}" >> "$LOG_FILE"
echo "Unique protein molecule types: ${MOL_TYPES[*]}" >> "$LOG_FILE"
echo "Total physical protein chains: $total_chains" >> "$LOG_FILE"
 
# ---------------------------------------------------------------------------
# 2. Detect chain boundaries in GRO file using NME residues
#    Each chain ends with one NME capping group (6 atoms).
#    We extract the last residue number of each NME block to get the last
#    residue number of each physical chain.
# ---------------------------------------------------------------------------
n_nme=$(grep -c "NME" em.gro)
detected_chains=$((n_nme / 6))
echo "Detected $detected_chains chains in em.gro based on NME residues." >> "$LOG_FILE"
 
if [[ $detected_chains -ne $total_chains ]]; then
    echo "WARNING: topol.top implies $total_chains chains but GRO file has $detected_chains NME-capped chains." >> "$LOG_FILE"
fi
 
# Get last residue number of each chain (every 6th NME atom occurrence = last atom of that NME)
nme_chain_ends=($(awk '
{
    is_nme = ($0 ~ /NME/)
    if (in_nme && !is_nme) {
        # Transition: just left an NME block → emit the last NME atom
        print prev_resnum, prev_atomidx
    }
    if (is_nme) {
        prev_resnum = $1 + 0
        prev_atomidx = substr($0, 16, 5) + 0
    }
    in_nme = is_nme
}
END {
    # Handle case where file ends inside an NME block (unlikely but safe)
    if (in_nme) print prev_resnum, prev_atomidx
}' em.gro))
 
# nme_chain_ends: flat array of pairs "resnum atomidx" per NME residue (= per chain)
chain_end_residues=()
chain_end_atomidx=()
n_fields=${#nme_chain_ends[@]}
for ((i=0; i<n_fields; i+=2)); do
    chain_end_residues+=("${nme_chain_ends[$i]}")
    chain_end_atomidx+=("${nme_chain_ends[$((i+1))]}")
done
echo "Chain end residues:  ${chain_end_residues[*]}" >> "$LOG_FILE"
echo "Chain end atom indices: ${chain_end_atomidx[*]}" >> "$LOG_FILE"
 
# Build residue ranges per physical chain (residue numbers, 1-indexed)
start=1
end=0
declare -a CHAIN_RANGES=()
for end_residue in "${chain_end_residues[@]}"; do
    end=$((end_residue + end))
    CHAIN_RANGES+=("${start}-${end}")
    start=$((end + 1))
done
echo "Residue ranges per physical chain: ${CHAIN_RANGES[*]}" >> "$LOG_FILE"
 
# ---------------------------------------------------------------------------
# 3. Map each unique molecule type to the residue range of its FIRST
#    physical chain occurrence
# ---------------------------------------------------------------------------
declare -A TYPE_TO_RANGE=()
declare -A TYPE_TO_SHIFT=()   # atom index of last NME atom before this chain
 
for mol in "${MOL_TYPES[@]}"; do
    for ((ci=0; ci<${#CHAIN_SEQ[@]}; ci++)); do
        if [[ "${CHAIN_SEQ[$ci]}" == "$mol" ]]; then
            TYPE_TO_RANGE[$mol]="${CHAIN_RANGES[$ci]}"
            if [[ $ci -eq 0 ]]; then
                TYPE_TO_SHIFT[$mol]=0
            else
                TYPE_TO_SHIFT[$mol]="${chain_end_atomidx[$((ci-1))]}"
            fi
            echo "Molecule type $mol → first physical chain index $ci, range ${CHAIN_RANGES[$ci]}, atom shift ${TYPE_TO_SHIFT[$mol]}" >> "$LOG_FILE"
            break
        fi
    done
done

# ---------------------------------------------------------------------------
# 4. Create index.ndx if needed
# ---------------------------------------------------------------------------
if [[ ! -f index.ndx ]]; then
    echo "q" | $GMX make_ndx -f em.gro -o index.ndx >> "$LOG_FILE" 2>&1
fi
 
# ---------------------------------------------------------------------------
# 5. For each unique molecule type, generate posre_<type>.itp
#    using the residue range of its first physical chain.
#    Atom indices are renumbered to start from 1 so GROMACS can apply the
#    file to any instance of that molecule type regardless of absolute index.
# ---------------------------------------------------------------------------
for mol in "${MOL_TYPES[@]}"; do
    posre_file="posre_${mol}.itp"
    range="${TYPE_TO_RANGE[$mol]}"
    shift="${TYPE_TO_SHIFT[$mol]}"
    echo "Processing molecule type $mol with range $range and atom index shift $shift" >> "$LOG_FILE"
 
    if [[ -f "$posre_file" ]]; then
        echo "$posre_file already exists – skipping." >> "$LOG_FILE"
        continue
    fi
 
    echo "Generating $posre_file using range $range, atom index shift $shift" >> "$LOG_FILE"

    # Add residue range index group if not present
    if ! grep -Fq "Protein-H_&_r_${range}" index.ndx; then
        echo "Adding group Protein-H_&_r_${range} to index.ndx" >> "$LOG_FILE"
        echo -e "2 & ri ${range}\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> "$LOG_FILE" 2>&1
    else
        echo "Group Protein-H_&_r_${range} already exists in index.ndx" >> "$LOG_FILE"
    fi
 
    # Generate restraint file
    echo "Protein-H_&_r_${range}" | $GMX genrestr \
        -f em.gro -n index.ndx \
        -o "$posre_file" \
        -fc 1000 1000 1000 >> "$LOG_FILE" 2>&1
	
	cp "$posre_file" "$posre_file"_backup
    # Renumber atom indices to start from 1 (so the file is portable
    # across all identical chains regardless of their position in em.gro)
    awk -v shift="$shift" '
    /^\[ position_restraints \]/ { in_section=1; print; next }
    /^\[/ && !/\[ position_restraints \]/ { in_section=0 }
    {
        if (in_section && /^[[:space:]]*[0-9]/) {
            printf "%6d %4d %10d %10d %10d\n", $1-shift, $2, $3, $4, $5
        } else {
            print
        }
    }
    ' "$posre_file" > "${posre_file}.tmp" && mv "${posre_file}.tmp" "$posre_file"
 
    echo "  → $posre_file written and renumbered." >> "$LOG_FILE"
done


# ---------------------------------------------------------------------------
# 6. Ligand position restraints
# ---------------------------------------------------------------------------

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
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 2 >> $LOG_FILE 2>&1
$GMX mdrun -v -deffnm nvt >> $LOG_FILE 2>&1

if [ -f nvt.gro ]; then
	echo "'nvt.gro' created" >> $LOG_FILE 2>&1
	echo -e "Temperature \n 0" | $GMX energy -f nvt.edr -o temperature.xvg >> $LOG_FILE 2>&1
else
	echo "Error: Failed to create 'nvt.gro'" >> $LOG_FILE 2>&1
	exit 1
fi

#--------------- NPT --------------------
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