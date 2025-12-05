import subprocess
from pathlib import Path
import re
import shutil
import sys
from src import constants
from src.utils import get_class_logger
import time

logger = get_class_logger(__name__)

def gromacs_equil(sandbox_dir: str, input_gro: str, ligand_name=None, ligand_file=None) -> str:
    # sometimes llm passes ligands as empty strings
    if not ligand_name:
        ligand_name = None
    if not ligand_file:
        ligand_file = None
        
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

    # --- Loop over each [ moleculetype ] block ---
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

    script = constants.SCRIPTS_DIR / "equil_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_equil.log")

    cmd = [str(script), sandbox_dir, input_gro, log_file_path]

    print("I will convert 89G_h.pdb to 89.gro using obabel")  # Debug print
    #time.sleep(10)

    if ligand_file is not None:
        obabel_cmd = f"obabel {ligand_file} -O {sandbox_dir}/{ligand_name}.gro"
        obabel_result = subprocess.run(obabel_cmd, cwd=sandbox_dir, capture_output=True, text=True, shell=True)
        if obabel_result.returncode != 0:
            error_text = "\n".join(filter(None, [obabel_result.stderr, obabel_result.stdout]))
            return f"Equilibration failed with error: {obabel_result.stderr}" 
        ligand_gro = f"{ligand_name}.gro"
        #time.sleep(10)
        cmd.append(ligand_name)
        cmd.append(ligand_file)
        cmd.append(ligand_gro)
        print(cmd)  # Debug print

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        # Report failure and include the captured GROMACS output for debugging
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly.'}") # Note: result.stderr will be empty if we set stderr=sys.stderr, but we keep it here for safety.
    else:
        # Report success and return the captured GROMACS output
        return (f"Equilibration ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")


def gromacs_production(sandbox_dir: str, input_gro: str, npt_cpt_file: str, ligand_name=None) -> str:
    """
    Run production MD with GROMACS using prod_Gromacs.sh.
    """
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
        # Report failure and include the captured GROMACS output for debugging
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly.'}") # Note: result.stderr will be empty if we set stderr=sys.stderr, but we keep it here for safety.
    else:
        # Report success and return the captured GROMACS output
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
        # Report failure and include the captured GROMACS output for debugging
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly.'}") # Note: result.stderr will be empty if we set stderr=sys.stderr, but we keep it here for safety.
    else:
        # Report success and return the captured GROMACS output
        return (f"Equilibration ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")