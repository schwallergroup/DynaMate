import subprocess
from pathlib import Path
import re
import shutil
import sys
from src import constants
from src.utils import get_class_logger
from collections import defaultdict, Counter

logger = get_class_logger(__name__, log_to_file=False)


def _filter_mdrun_output(output: str) -> str:
    """Remove verbose mdrun simulation output from GROMACS log."""
    output = re.sub(
        r"Steepest Descents:.*?writing lowest energy coordinates\.\n?",
        "", output, flags=re.DOTALL
    )
    output = re.sub(
        r"starting mdrun 'Generic title'.*?Writing final coordinates\.\n?",
        "", output, flags=re.DOTALL
    )
    output = re.sub(r"^GROMACS reminds you:.*$\n?", "", output, flags=re.MULTILINE)
    output = re.sub(r"^Back Off! I just backed up.*$\n?", "", output, flags=re.MULTILINE)
    
    return output

    


def gromacs_equil(sandbox_dir: str, input_gro: str, md_temp: str, ligand_name=None, ligand_files=None) -> str:
    # sometimes llm passes ligands as empty strings or the string "None"
    if not ligand_name or str(ligand_name).strip().lower() == "none":
        ligand_name = None
    if not ligand_files:
        ligand_file = None
    else:
        ligand_file=ligand_files[0]

    # Helper functions
    def _make_posre_block(posre_file: str) -> str:
        return (
            f'; Include Position restraint file\n'
            f'#ifdef POSRES\n'
            f'#include "{posre_file}"\n'
            f'#endif\n\n'
        )

    def _parse_protein_mol_types(text: str) -> list[str]:
        """
        Parse the [ molecules ] section and return unique protein molecule type
        names (system*) in first-seen order.

        Handles both formats:
            system1    2          ← two identical chains of type system1
            system1    1          ← followed later by another system1 entry
            system2    1
        """
        in_molecules = False
        seen: dict[str, None] = {}   # ordered set (insertion-ordered dict)

        for line in text.splitlines():
            stripped = line.strip()

            # Detect section header
            if re.match(r'^\[\s*molecules\s*\]', stripped, re.I):
                in_molecules = True
                continue

            # Stop at next section header
            if in_molecules and re.match(r'^\[', stripped):
                break

            # Skip blanks and comments
            if not stripped or stripped.startswith(';'):
                continue

            if in_molecules:
                parts = stripped.split()
                if len(parts) >= 2 and re.match(r'^system\d*$', parts[0], re.I):
                    mol = parts[0]
                    if mol not in seen:
                        seen[mol] = None

        return list(seen.keys())
    
    # Modify topol.top to include position restraints
    sandbox_dir = Path(sandbox_dir)
    input_path  = sandbox_dir / "topol.top"
    backup_path = sandbox_dir / "topol_without_posre.top"
 
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} not found.")
 
    # Normalise ligand_name
    _skip_ligand = {"None", "XXX", "None_h", None, ""}
    handle_ligand = ligand_name not in _skip_ligand
 
    # Backup
    shutil.copyfile(input_path, backup_path)
    logger.info(f"Backup created: {backup_path}")
 
    text = input_path.read_text(encoding="utf-8", errors="replace")
 
    # ── 1. Discover protein molecule types from [ molecules ] ─────────────────
    protein_mol_types: list[str] = _parse_protein_mol_types(text)
    logger.info(f"Protein molecule types found in [ molecules ]: {protein_mol_types}")
 
    if not protein_mol_types:
        logger.warning("No protein chain entries (system*) found in [ molecules ]. "
                       "No position restraints inserted.")
        return
 
    # ── 2. Find the extent of [ moleculetype ] sections ──────────────────────
    # Everything from the first [ moleculetype ] up to [ system ] / [ molecules ]
    next_non_mol = re.search(
        r'^\[\s*(system|molecules)\s*\]', text, flags=re.I | re.M
    )
    end_of_moltypes = next_non_mol.start() if next_non_mol else len(text)
 
    header_re = re.compile(r'^\[\s*moleculetype\s*\]', flags=re.I | re.M)
    headers   = [m for m in header_re.finditer(text) if m.start() < end_of_moltypes]
 
    if not headers:
        logger.warning("No [ moleculetype ] sections found. No changes made.")
        return
 
    positions = [m.start() for m in headers] + [end_of_moltypes]
    preamble  = text[: positions[0]]
 
    # ── 3. Walk each [ moleculetype ] block and insert POSRES if needed ───────
    segments: list[str]    = []
    inserted_protein: set  = set()   # unique types already given a posre block
    inserted_ligand        = False
    report: list[str]      = []
 
    for i in range(len(positions) - 1):
        seg = text[positions[i] : positions[i + 1]]
 
        # Extract molecule name from the "<name>  <nrexcl>" line that follows
        # the [ moleculetype ] header.
        name_match = re.search(
            r'^\[\s*moleculetype\s*\][^\[]*?^\s*(\S+)\s+\d+',
            seg, flags=re.I | re.M | re.S
        )
        mol_name = name_match.group(1) if name_match else None
 
        # ── Protein chain ──────────────────────────────────────────────────
        if mol_name and mol_name in protein_mol_types:
            if mol_name not in inserted_protein:
                posre_file  = f"posre_{mol_name}.itp"
                posre_block = _make_posre_block(posre_file)
                if posre_block not in seg:
                    seg = seg.rstrip() + "\n\n" + posre_block
                    report.append(mol_name)
                inserted_protein.add(mol_name)
            else:
                # Duplicate moleculetype block for the same type — do NOT
                # insert another include; GROMACS applies posre per molecule
                # type, not per block.
                logger.debug(
                    f"Skipping duplicate [ moleculetype ] block for {mol_name}"
                )
 
        # ── Ligand ────────────────────────────────────────────────────────
        if handle_ligand and not inserted_ligand and mol_name == ligand_name:
            posre_block = _make_posre_block(f"posre_{ligand_name}.itp")
            if posre_block not in seg:
                seg = seg.rstrip() + "\n\n" + posre_block
                report.append(ligand_name)
            inserted_ligand = True
 
        segments.append(seg)
 
    # ── 4. Reassemble and write ───────────────────────────────────────────────
    modified_text = preamble + "".join(segments) + text[end_of_moltypes:]
    input_path.write_text(modified_text, encoding="utf-8")
 
    logger.info(
        f"Position restraint includes added for: {', '.join(report) or 'none'}"
    )
    if handle_ligand and not inserted_ligand:
        logger.warning(
            f"Ligand '{ligand_name}' was specified but no matching "
            f"[ moleculetype ] block was found — posre include NOT inserted."
        )

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

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
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
            gromacs_output = _filter_mdrun_output(log_file_path.read_text(encoding="utf-8"))
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
    if not ligand_name or str(ligand_name).strip().lower() == "none":
        ligand_name = None

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
    
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        cmd.append(ligand_name)
        cmd.append(f"{sandbox_dir}/{ligand_name}.gro")

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = _filter_mdrun_output(log_file_path.read_text(encoding="utf-8"))
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")
    else:
        return (f"Production ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")


def gromacs_analysis(sandbox_dir: str,  pdb_id: str, input_xtc: str, ligand_name=None) -> str:
    """
    Run production MD with GROMACS using prod_Gromacs.sh.
    """
    if not ligand_name or str(ligand_name).strip().lower() == "none":
        ligand_name = None

    script = constants.SCRIPTS_DIR / "analysis_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_analysis.log")

    cmd = [str(script), input_xtc, log_file_path]

<<<<<<< HEAD
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
=======
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h", ""]:
>>>>>>> b57c53e ([update] upskill and refinement fixed)
        cmd.append(ligand_name)
        cmd.append(f"{sandbox_dir}/{ligand_name}.gro")

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    # Check if ligand was not edited during the setup of simulations
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:

        logger.info("Checking if the ligand was modified during setup by comparing the original PDB of the ligand (from the input PDB) and final coordinates of the ligand (from md.gro). This is done to ensure that heavy atoms were not removed or added, which would modify the ligand.")

        ligands = defaultdict(list)
        ligands_gro = defaultdict(list)
        ligands_pdb_check = Counter()
        ligands_gro_check = Counter()

        with open(f"{sandbox_dir}/{pdb_id}.pdb", "r") as infile:
            for line in infile:
                if line.startswith("HETATM") and line[17:20].strip() == ligand_name:
                    chain = line[21].strip() or "_"
                    resnum = int(line[22:26])
                    ligands[(chain, resnum)].append(line)    

        ligand_pdb_file = f"{sandbox_dir}/{ligand_name}_check.pdb"
        with open(ligand_pdb_file, "w") as outfile:
            first_key = next(iter(ligands))
            outfile.writelines(ligands[first_key])
        logger.info(f"Extracted first occurence of ligand {ligand_name} to {ligand_pdb_file} to check if the ligand was modified during setup.")

        with open(f"{sandbox_dir}/md.gro", "r") as infile:
            ligand_gro_file = f"{sandbox_dir}/{ligand_name}_check.gro"
            for line in infile:
                if line[5:8] == ligand_name and "H" not in line[11:15]:
                    resnum = int(line[:5])
                    ligands_gro[(resnum)].append(line)
        with open(ligand_gro_file, "w") as outfile:
            first_key = next(iter(ligands_gro))
            atom_lines = ligands_gro[first_key]
            outfile.write("Ligand GRO file\n")
            outfile.write(f"{len(atom_lines)}\n")
            outfile.writelines(atom_lines)
            outfile.write("0.0 0.0 0.0\n")
            outfile.write("Ligand GRO file\n")
        logger.info(f"Extracted ligand {ligand_name} from md.gro to {ligand_gro_file} to check if the ligand was modified during setup.")
        subprocess.run(f"obabel {ligand_gro_file} -O {sandbox_dir}/{ligand_name}_gro_check.pdb", shell=True, check=False)

        #check if the 2 PDB are the same
        with open(f"{sandbox_dir}/{ligand_name}_check.pdb", "r") as pdb_file, open(f"{sandbox_dir}/{ligand_name}_gro_check.pdb", "r") as gro_pdb_file:
            pdb_lines = [line for line in pdb_file if line.startswith("HETATM") or line.startswith("ATOM")]
            gro_pdb_lines = [line for line in gro_pdb_file if line.startswith("HETATM") or line.startswith("ATOM")]

        if len(pdb_lines) != len(gro_pdb_lines):
            logger.info(f"Number of atoms differ between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro): {len(pdb_lines)} vs {len(gro_pdb_lines)}")
        else:
            logger.info(f"Number of atoms match between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro): {len(pdb_lines)}")

        for line in pdb_lines:
            if line[17:20].strip() == ligand_name:
                atom = line[75:80].strip()
                ligands_pdb_check[atom] += 1

        for line in gro_pdb_lines:
            if line[17:20].strip() == ligand_name:
                atom = line[75:80].strip()
                ligands_gro_check[atom] += 1

        # Compare the counts of each atom type
        if ligands_pdb_check == ligands_gro_check:
            logger.info(f"Atom counts (except hydrogens) match between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro).")
        else:
            logger.info(f"Atom counts (except hydrogens) differ between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro).")
            logger.info(f"The atoms that the original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro) have in common are: {ligands_pdb_check & ligands_gro_check}")
            logger.info(f"The difference between the original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro) is: {ligands_pdb_check - ligands_gro_check}")

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")
    else:
        return (f"Analysis plots produced successfully. Full GROMACS output:\n"
                f"{gromacs_output}")

