#!/bin/bash
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 sandbox_dir input_gro ionic_conc log_file [ligand_name] [ligand_file] [ligand_gro]"
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
IONIC_CONC="${3:-0.15}"  # Default to 0.15 M if not set
LOG_FILE="$4"

> $LOG_FILE 
echo "Starting GROMACS Equilibration Log" >> $LOG_FILE 2>&1

# Optional argument
if [ "$#" -ge 5 ]; then
    LIGNAME="$5"
	LIGFILE="$6"
	LIGGRO="$7"
else
    LIGNAME=""
	LIGFILE=""
	LIGGRO=""
fi

#------- FILES BACKUPS AND CLEANUPS ----
cp topol.top topol_before_equil.top && echo "Copied topol.top to topol_before_equil.top" >> $LOG_FILE 2>&1
sed -i 's/WAT/SOL/g' topol.top >> $LOG_FILE 2>&1

rm -f solv.gro solv_ions.gro em.gro nvt.gro npt.gro

# NaCl concentration (M) added on top of the neutralizing counter-ions, set by gromacs_equil
# (0.15 = physiological environment, 0 = counter-ions only)
echo "Ionic concentration: $IONIC_CONC M NaCl" >> $LOG_FILE 2>&1

#------- SOLVATE and BOX ---------------
if ! [ -f solv.gro ]; then
    gmx editconf -f $INPUT_GRO -o newbox.gro -bt cubic -c -d 1.0 >> $LOG_FILE 2>&1
    sed -i 's/WAT/SOL/g' newbox.gro >> $LOG_FILE 2>&1
    gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro >> $LOG_FILE 2>&1
    sed -i 's/OW/ O/g' solv.gro >> $LOG_FILE 2>&1
    sed -i 's/HW/ H/g' solv.gro >> $LOG_FILE 2>&1
fi

if ! [ -f solv_ions.gro ]; then
    gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr >> $LOG_FILE 2>&1
    echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname Na -nname Cl -neutral -conc $IONIC_CONC >> $LOG_FILE 2>&1
fi

#------- ENERGY MINIMISATION ------------
if ! [ -f em.gro ]; then
	$GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr >> $LOG_FILE 2>&1
	$GMX mdrun -v -deffnm em >> $LOG_FILE 2>&1
fi

if [ -f em.gro ]; then
	echo "'em.gro' created"
	echo "11 0" | $GMX energy -f em.edr -o potential.xvg >> $LOG_FILE 2>&1
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
# 2. Create index.ndx if needed
# ---------------------------------------------------------------------------
if [[ ! -f index.ndx ]]; then
    echo "q" | $GMX make_ndx -f em.gro -o index.ndx >> "$LOG_FILE" 2>&1
    cp index.ndx index_before_equil.ndx && echo "Copied index.ndx to index_before_equil.ndx" >> $LOG_FILE 2>&1
else
    cp index.ndx index_before_equil.ndx && echo "Copied index.ndx to index_before_equil.ndx" >> $LOG_FILE 2>&1
fi

# ---------------------------------------------------------------------------
# 3. For each unique molecule type, generate posre_<type>.itp directly from
#    the [ atoms ] section of its [ moleculetype ] in topol.top.
#    Atom numbers in [ atoms ] are already local to the molecule type (start
#    at 1), so no chain boundary detection or renumbering is needed and the
#    chains do not have to be capped. Heavy atoms only (same selection as
#    Protein-H: atom names whose first non-digit character is not H).
# ---------------------------------------------------------------------------
for mol in "${MOL_TYPES[@]}"; do
    posre_file="posre_${mol}.itp"

    if [[ -f "$posre_file" ]]; then
        echo "$posre_file already exists – skipping." >> "$LOG_FILE"
        continue
    fi

    echo "Generating $posre_file from [ atoms ] of $mol in $TOPOL_FILE" >> "$LOG_FILE"

    awk -v mol="$mol" '
    /^[[:space:]]*\[/ {
        section = $0
        gsub(/[][[:space:]]/, "", section)
        if (section == "moleculetype") { want_name = 1; cur = "" }
        next
    }
    /^[[:space:]]*(;|#|$)/ { next }
    want_name { cur = $1; want_name = 0; next }
    cur == mol && section == "atoms" {
        name = $5
        sub(/^[0-9]+/, "", name)
        if (name !~ /^H/) {
            printf "%6d %4d %10d %10d %10d\n", $1, 1, 1000, 1000, 1000
            n++
        }
    }
    END { exit (n > 0 ? 0 : 1) }
    ' "$TOPOL_FILE" > "${posre_file}.tmp"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: no heavy atoms found for molecule type $mol in $TOPOL_FILE – $posre_file not written." | tee -a "$LOG_FILE" >&2
        rm -f "${posre_file}.tmp"
        exit 1
    fi

    {
        echo "; position restraints for heavy atoms of $mol (generated from $TOPOL_FILE)"
        echo ""
        echo "[ position_restraints ]"
        echo ";  i funct       fcx        fcy        fcz"
        cat "${posre_file}.tmp"
    } > "$posre_file"
    rm -f "${posre_file}.tmp"

    echo "  → $posre_file written ($(grep -c '^[[:space:]]*[0-9]' "$posre_file") restrained atoms)." >> "$LOG_FILE"
done

# ---------------------------------------------------------------------------
# 6. Ligand position restraints
# ---------------------------------------------------------------------------

if [ -n "$LIGNAME" ]; then
    # ligand position retraints
	if ! ls "posre_$LIGNAME.itp" 1> /dev/null 2>&1;then
		echo "Generating position restraints for ligand $LIGNAME" >> $LOG_FILE 2>&1
		echo -e "0 & ! a H* \n q" | $GMX make_ndx -f $LIGGRO -o "index_$LIGNAME.ndx" >> $LOG_FILE 2>&1
		echo "3" | $GMX genrestr -f $LIGGRO -n "index_$LIGNAME.ndx" -o "posre_$LIGNAME.itp" -fc 1000 1000 1000 >> $LOG_FILE 2>&1
	fi

    # Create protein-ligand tcoupling group
    if grep "Protein_$LIGNAME" index.ndx; then
        echo "Protein_$LIGNAME already in index.ndx" >> $LOG_FILE 2>&1
    else
        echo -e "1 | 13\nq" | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
    fi
fi

# Create group Water_and_ions if not exists
if grep -q "Water_and_ions" index.ndx; then
	echo "Group Water_and_ions already exists in index.ndx"	
else
	if grep -q "Cl" index.ndx; then
        if grep -q "Na" index.ndx; then
            echo -e '"SOL" | "Cl" | "Na" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
            sed -i 's/SOL_Cl_Na/Water_and_ions/g' index.ndx
            echo "Group Water_and_ions created in index.ndx" >> $LOG_FILE 2>&1
        else
            echo -e '"SOL" | "Cl" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
            sed -i 's/SOL_Cl/Water_and_ions/g' index.ndx
            echo "Group Water_and_ions created in index.ndx" >> $LOG_FILE 2>&1
        fi
    elif grep -q "Na" index.ndx; then
		echo -e '"SOL" | "Na" \n q' | $GMX make_ndx -f em.gro -n index.ndx -o index.ndx >> $LOG_FILE 2>&1
		sed -i 's/SOL_Na/Water_and_ions/g' index.ndx
		echo "Group Water_and_ions created in index.ndx" >> $LOG_FILE 2>&1
	else
        echo "No ions present. Turning Water_and_ions group into Water in mdp files." >> $LOG_FILE 2>&1
        if grep "Water_and_ions" nvt.mdp; then
		    sed -i "s|Water_and_ions|Water|g" nvt.mdp
            echo "Water_and_ions replaced with Water in nvt.mdp" >> $LOG
    	else
		    echo "Water_and_ions line was not found in nvt.mdp." >> $LOG_FILE 2>&1
	    fi
        if grep "Water_and_ions" npt.mdp; then
		    sed -i "s|Water_and_ions|Water|g" npt.mdp
            echo "Water_and_ions replaced with Water in npt.mdp" >> $LOG
    	else
		    echo "Water_and_ions line was not found in npt.mdp." >> $LOG_FILE 2>&1
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
# NPT is run by gromacs_equil (gromacs_tools.py) after the NVT temperature check,
# so that NVT can be extended first if the temperature has not reached its reference.
