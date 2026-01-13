import subprocess
from pathlib import Path
import subprocess, shlex
import re
from src.utils import get_class_logger

logger = get_class_logger(__name__)

def fix_charges(input_file, output_file=None):
            """
            Adjust charges so the total is an integer.
            Only modifies one atom's charge (last atom in coordinate/charge section).
            """
            with open(input_file) as f:
                lines = f.readlines()

            # Identify section that looks like atom definitions (before LOOP/IMPROPER)
            atom_section_end = len(lines)
            for i, line in enumerate(lines):
                if re.match(r"^\s*(LOOP|IMPROPER|DONE|STOP)\b", line):
                    atom_section_end = i
                    break

            atom_lines = []
            charges = []
            float_pattern = re.compile(r"([-+]?\d*\.\d+|\d+)(?!.*\S)")

            for i, line in enumerate(lines[:atom_section_end]):
                tokens = line.split()
                if len(tokens) >= 8:  # typical atom line length
                    last = tokens[-1]
                    try:
                        charge = float(last)
                        atom_lines.append(i)
                        charges.append(charge)
                    except ValueError:
                        pass

            if not charges:
                logger.warning("No atom charges found.")
                return

            total = sum(charges)
            target = round(total)
            delta = target - total

            if abs(delta) < 1e-6:
                logger.warning(f"Already integer total ({total:.6f}). No change.")
                return

            # Adjust last atom charge in section
            last_atom_idx = atom_lines[-1]
            old_charge = charges[-1]
            new_charge = old_charge + delta

            lines[last_atom_idx] = float_pattern.sub(f"{new_charge:.6f}", lines[last_atom_idx])

            if output_file is None:
                output_file = Path(input_file).with_name(Path(input_file).stem + "_fixed.res")

            with open(output_file, "w") as f:
                f.writelines(lines)

            logger.info(f"Adjusted total charge: {total:.6f} → {target}")
            logger.info(f"Atom line {last_atom_idx + 1}: {old_charge:.6f} → {new_charge:.6f}")
            logger.info(f"Saved to: {output_file}")

def param_ligand(sandbox_dir: str, ligand_files: str | list[str], ligand_name: str, charge_ligand: str | None = None) -> str:
    if isinstance(ligand_files, str):
        ligand_files = [ligand_files]
    else:
        ligand_files = ligand_files

    ligand_file=ligand_files[0]
    if len(ligand_files) > 1:
        logger.info(f"For now, we will only parameterize the first ligand since we will only simulate one ligand: {ligand_file}")

    logger.info(f"Parameterizing ligand file: {ligand_file}")
    ligand_stem=Path(f"{ligand_file}").stem

    # Find charge of ligand if not provided
    charges = []

    tmp_mol2file = f"{sandbox_dir}/{ligand_stem}.mol2"
    # Create temporary mol2 file using obabel
    cmd = shlex.split(f"obabel {sandbox_dir}/{ligand_file} -O {tmp_mol2file}")
    run_1 = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True)
    if run_1.returncode != 0:
        error_text = "\n".join(filter(None, [run_1.stderr, run_2.stdout]))
        return f"Ligand parameterization failed with error: {error_text}"

    # Read the mol2 file to find the charge
    in_atom_section = False

    with open(tmp_mol2file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>") and in_atom_section:
                break
            elif in_atom_section and line:
                try:
                    charge = float(line.split()[-1])
                    charges.append(charge)
                except ValueError:
                    pass  # skip malformed lines
    total_charge = sum(charges)
    charge_ligand = round(total_charge)

    # Clean up temporary mol2 file
    Path(tmp_mol2file).unlink(missing_ok=True)

    logger.info(f"Charge of ligand {ligand_stem} determined to be {charge_ligand}")

    # Create mol2 file using antechamber
    cmd = shlex.split(
        f"antechamber -i {sandbox_dir}/{ligand_file} -fi pdb -o {sandbox_dir}/{ligand_stem}.mol2 -fo mol2 -c bcc -nc {charge_ligand} -s 2"
    )
    run_2 = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True)
    if run_2.returncode != 0:
        error_text = "\n".join(filter(None, [run_2.stderr, run_2.stdout]))
        return f"Ligand parameterization failed with error: {error_text}"

    logger.info(f"Mol2 file for ligand {ligand_stem} created")

    cmd = shlex.split(f"sed -i 's/UNL/{ligand_name}/g' {sandbox_dir}/{ligand_stem}.mol2")
    run_3 = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True)
    if run_3.returncode != 0:
        error_text = "\n".join(filter(None, [run_3.stderr, run_3.stdout]))
        return f"Ligand parameterization failed with error: {error_text}"

    # Create prepi file using antechamber
    cmd = shlex.split(
        f"antechamber -i {sandbox_dir}/{ligand_stem}.mol2 -fi mol2 -o {sandbox_dir}/{ligand_stem}.prepi -fo prepi -c bcc"
    )
    run_4 = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True)
    if run_4.returncode != 0:
        error_text = "\n".join(filter(None, [run_4.stderr, run_4.stdout]))
        return f"Ligand parameterization failed with error: {error_text}"
    logger.info(f"Prepi file for ligand {ligand_file} created: {ligand_stem}.prepi")

    # Create frcmod file using parmchk2
    cmd = shlex.split(f"parmchk2 -i {sandbox_dir}/{ligand_stem}.mol2 -f mol2 -o {sandbox_dir}/{ligand_stem}.frcmod")
    run_5 = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True, check=True)
    if run_5.returncode != 0:
        error_text = "\n".join(filter(None, [run_5.stderr, run_5.stdout]))
        return f"Ligand parameterization failed with error: {error_text}"
    logger.info(f"Frcmod file for ligand {ligand_file} created: {ligand_stem}.frcmod")

    # Update charge prepi file
    fix_charges(f"{sandbox_dir}/{ligand_stem}.prepi", f"{sandbox_dir}/{ligand_stem}_fixed.prepi")

    if Path(f"{sandbox_dir}/{ligand_stem}_fixed.prepi").exists():
        prepi_file = f"{ligand_stem}_fixed.prepi"
    else:
        prepi_file = f"{ligand_stem}.prepi"

    if len(ligand_files) == 1:
        return f"Ligand parameterisation complete. Parameters saved to {sandbox_dir}/{prepi_file} and {sandbox_dir}/{ligand_stem}.frcmod"
    if len(ligand_files) > 1:
        return f"Ligand parameterisation complete for the the first ligand: {ligand_files[0]}. Parameters saved to {sandbox_dir}/{ligand_stem}.frcmod and {sandbox_dir}/{prepi_file}"
