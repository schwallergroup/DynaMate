import subprocess
from pathlib import Path
import re
import shutil
import sys
from src import constants
from src.utils import get_class_logger
import time

logger = get_class_logger(__name__)

def gromacs_equil(sandbox_dir: str, input_gro: str, md_temp: str, ligand_name=None, ligand_files=None) -> str:
    # sometimes llm passes ligands as empty strings
    if not ligand_name:
        ligand_name = None
    if not ligand_files:
        ligand_file = None
    else:
        ligand_file=ligand_files[0]
        
    # ------------ Modify topol.top to include position restraints ------------

    input_path = Path(f"{sandbox_dir}/topol.top")
    backup_path = Path(f"{sandbox_dir}/topol_without_posre.top")

    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} not found.")

    # Make a backup
    shutil.copyfile(input_path, backup_path)
    logger.info(f"Backup created: {backup_path}")

    text = input_path.read_text(encoding="utf-8", errors="replace")

    # --- Detect all system names (system, system1, system2, etc.) ---
    system_matches = re.findall(r"^\s*(system\d*?)\b", text, flags=re.M)
    systems = sorted(set(system_matches), key=lambda x: int(re.search(r"\d*$", x).group() or 0))
    num_systems = len(systems)

    if num_systems == 0:
        logger.warning("No system entries found. No changes made.")
        sys.exit(0)
    
    with open(f"{sandbox_dir}/topol_without_posre.top", "r") as f:
        if re.search(r"system1\s+2", f.read()):
            logger.info(f"Detected 2 chain(s) because of 2 identical chains named system1 in topol.top")
        else:
            logger.info(f"Detected {num_systems} chain(s): {', '.join(systems)}")

    # --- Create posre include block ---
    def make_posre_block(posre_file):
        return f'; Include Position restraint file\n#ifdef POSRES\n#include "{posre_file}"\n#endif\n\n'

    if ligand_file is not None:
        ligand_posre_block = (
            f'; Include Position restraint file\n#ifdef POSRES\n#include "posre_{ligand_name}.itp"\n#endif\n\n'
        )

    # --- Split by [ moleculetype ] sections ---
    header_re = re.compile(r"^\[\s*moleculetype\s*\]", flags=re.I | re.M)
    headers = list(header_re.finditer(text))

    # If there are no [ moleculetype ] sections, exit
    if not headers:
        logger.warning("No [ moleculetype ] sections found. No changes made.")
        sys.exit(0)

    # Get the position of the first non-moleculetype section (e.g., [ system ])
    next_non_moleculetype = re.search(r"^\[\s*(system|molecules)\s*\]", text, flags=re.I | re.M)
    end_of_moleculetype = next_non_moleculetype.start() if next_non_moleculetype else len(text)

    # Define the boundaries for moleculetype sections
    positions = [m.start() for m in headers if m.start() < end_of_moleculetype] + [end_of_moleculetype]
    preamble = text[: positions[0]]
    segments = []
    inserted_blocks = []
    inserted_ligand = False

    # Loop over each [ moleculetype ] block
    for i in range(len(positions) - 1):
        seg = text[positions[i] : positions[i + 1]]

        # Only modify real moleculetype blocks
        match = re.search(r"^\s*(system\d*|system)\b\s+\d+", seg, flags=re.M)

        if match:
            system_name = match.group(1)
            chain_index = re.search(r"\d+$", system_name)
            chain_num = int(chain_index.group()) if chain_index else 1

            posre_file = f"posre_chain{chain_num}.itp" if num_systems > 1 else "posre.itp"
            posre_block = make_posre_block(posre_file)

            if posre_block not in seg:
                seg = seg.rstrip() + "\n\n" + posre_block
                inserted_blocks.append(system_name)

        # Check for ligand_name
        if ligand_file is not None:
            if not inserted_ligand:
                if re.search(r"^\s*" + re.escape(ligand_name) + r"\b\s+\d+", seg, flags=re.M):
                    include_lig = f'#include "posre_{ligand_name}.itp"'
                    if include_lig not in seg:
                        seg = seg.rstrip() + "\n\n" + ligand_posre_block
                    inserted_ligand = True

        segments.append(seg)

    # Reassemble modified part + untouched rest
    modified_text = preamble + "".join(segments) + text[end_of_moleculetype:]

    # Write result
    input_path.write_text(modified_text, encoding="utf-8")

    logger.info(f"Added position restraints for: {', '.join(inserted_blocks) or 'none'}")
    ## Position restraints added in topol.top

    # -------------- Create em.mdp, nvt.mdp, npt.mdp files --------------

    em_mdp_infile = open(f'{sandbox_dir}/em.mdp', 'w' )
    em_mdp_infile.write(f'''; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
''')
    em_mdp_infile.close()

    nvt_mdp_infile = open(f'{sandbox_dir}/nvt.mdp', 'w' )
    nvt_mdp_infile.write(f'''title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000     ; 2 * 5000 = 10 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = {float(md_temp)}   {float(md_temp)}                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = Berendsen                     ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT 
''')
    nvt_mdp_infile.close()

    npt_mdp_infile = open(f'{sandbox_dir}/npt.mdp', 'w' )
    npt_mdp_infile.write(f'''title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000     ; 2 * 5000 = 10 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = {float(md_temp)}   {float(md_temp)}                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = Berendsen                     ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT 
''')
    npt_mdp_infile.close()

    # -------------- Run equil_Gromacs.sh script --------------

    script = constants.SCRIPTS_DIR / "equil_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_equil.log")

    cmd = [str(script), sandbox_dir, input_gro, log_file_path]

    if ligand_file is not None:
        obabel_cmd = f"obabel {ligand_file} -O {sandbox_dir}/{ligand_name}.gro"
        obabel_result = subprocess.run(obabel_cmd, cwd=sandbox_dir, capture_output=True, text=True, shell=True)
        if obabel_result.returncode != 0:
            error_text = "\n".join(filter(None, [obabel_result.stderr, obabel_result.stdout]))
            return f"Equilibration failed with error: {error_text}" 
        ligand_gro = f"{ligand_name}.gro"
        cmd.append(ligand_name)
        cmd.append(ligand_file)
        cmd.append(ligand_gro)
        print(cmd)  

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}") 
    else:
        return (f"Equilibration ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")


def gromacs_production(sandbox_dir: str, input_gro: str, npt_cpt_file: str, md_temp: str, md_duration: str, ligand_name=None) -> str:
    """
    Run production MD with GROMACS using prod_Gromacs.sh.
    """

    # ---------- Create md.mdp file --------------
    nsteps = int(((float(md_duration)) * 1000000) / 2)  # Convert ns to number of steps (2 fs per step)
    md_mdp_infile = open(f'{sandbox_dir}/md.mdp', 'w' )
    md_mdp_infile.write(f'''title                   = Protein-ligand complex MD simulation 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = {nsteps}   ; 2 * 500,000 = 1000 ps (0.1 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000     ; save energies every 10.0 ps
nstlog                  = 5000     ; update log file every 10.0 ps
nstxout-compressed      = 5000     ; save coordinates every 10.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein Water_and_ions        ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = {float(md_temp)}   {float(md_temp)}                     ; reference temperature, one for each group, in K
; Pressure coupling 
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration 
''')
    md_mdp_infile.close()

    # ---------- Run prod_Gromacs.sh script --------------

    script = constants.SCRIPTS_DIR / "prod_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_production.log")

    cmd = [str(script), input_gro, npt_cpt_file, log_file_path]
    
    if ligand_name is not None:
        cmd.append(ligand_name)
        cmd.append(f"{sandbox_dir}/{ligand_name}.gro")

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")
    else:
        return (f"Equilibration ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")


def gromacs_analysis(sandbox_dir: str, input_xtc: str, ligand_name=None) -> str:
    """
    Run production MD with GROMACS using prod_Gromacs.sh.
    """
    script = constants.SCRIPTS_DIR / "analysis_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_analysis.log")

    cmd = [str(script), input_xtc, log_file_path]

    if ligand_name is not None:
        cmd.append(ligand_name)
        cmd.append(f"{sandbox_dir}/{ligand_name}.gro")

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}") 
    else:
        return (f"Equilibration ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")