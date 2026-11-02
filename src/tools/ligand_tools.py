import subprocess
from pathlib import Path
import subprocess, shlex
import re
from src.utils import get_class_logger
from src.tools.ligand_protonation import check_ligand_consistency, load_record

logger = get_class_logger(__name__, log_to_file=False)

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

def param_ligand(sandbox_dir: str, ligand_files: str | list[str], ligand_name: str) -> str:
    if isinstance(ligand_files, str):
        ligand_files = [ligand_files]
    else:
        ligand_files = ligand_files

    ligand_file=ligand_files[0]
    if len(ligand_files) > 1:
        logger.info(f"For now, we will only parameterize the first ligand since we will only simulate one ligand: {ligand_file}")

    logger.info(f"Parameterizing ligand file: {ligand_file}")
    ligand_stem=Path(f"{ligand_file}").stem

    # The ligand files must describe the selected protonation state (not edited outside set_ligand_protonation)
    problems = check_ligand_consistency(sandbox_dir, ligand_stem)
    if problems:
        return "Ligand parameterization failed with error: consistency check failed: " + "; ".join(problems)

    # Net charge of the selected protonation state (= sum of the formal charges in {ligand_stem}.sdf, checked above)
    charge_ligand = load_record(sandbox_dir, ligand_stem)["selected"]["net_charge"]
    logger.info(f"Charge of ligand {ligand_stem} determined to be {charge_ligand}")

    # Create mol2 file using antechamber, from the .sdf so that bond orders and formal charges are not re-perceived
    cmd = shlex.split(
        f"antechamber -i {sandbox_dir}/{ligand_stem}.sdf -fi mdl -o {sandbox_dir}/{ligand_stem}.mol2 -fo mol2 -c bcc -nc {charge_ligand} -at gaff2 -rn {ligand_name} -s 2"
    )
    run_2 = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True)
    if run_2.returncode != 0 or not Path(f"{sandbox_dir}/{ligand_stem}.mol2").exists():
        error_text = "\n".join(filter(None, [run_2.stderr, run_2.stdout]))
        logger.error(f"Ligand parameterization failed with error: {error_text}")
        return f"Ligand parameterization failed with error: {error_text}"

    logger.info(f"Mol2 file for ligand {ligand_stem} created")

    # Check: same atoms as the .sdf and AM1-BCC partial charges adding up to the net charge
    problems = check_ligand_consistency(sandbox_dir, ligand_stem, mol2_file=f"{ligand_stem}.mol2")
    if problems:
        return "Ligand parameterization failed with error: consistency check failed: " + "; ".join(problems)

    if len(ligand_files) == 1:
        return f"mol2 file created for the ligand: {ligand_files[0]} (net charge {charge_ligand:+d}), saved to {sandbox_dir}/{ligand_stem}.mol2"
    return f"mol2 file created for the first ligand: {ligand_files[0]} (net charge {charge_ligand:+d}), saved to {sandbox_dir}/{ligand_stem}.mol2"
