import subprocess
from pathlib import Path
import parmed as pmd  # type: ignore
import subprocess, shlex
import os
import re
import shutil
from src import constants
from src.utils import get_class_logger
from src.tools.ligand_protonation import check_ligand_consistency, load_record, MAX_CHARGE_CORRECTION

logger = get_class_logger(__name__, log_to_file=False)


def old_run_tleap(sandbox_dir: str, input_pdb: str, pdb_id: str) -> str:
    """
    Run tleap preparation using run_tleap.sh.
    """
    script = constants.SCRIPTS_DIR / "run_tleap.sh"
    result = subprocess.run(
        [str(script), sandbox_dir, input_pdb, pdb_id], cwd=sandbox_dir, capture_output=True, text=True
    )
    if result.returncode != 0:
        # tleap often puts errors in stdout
        error_text = "\n".join(filter(None, [result.stderr, result.stdout]))
        return f"tleap run failed with error:\n{error_text}"
    else:
        # Parmed
        try:
            prmtop_path = f"{sandbox_dir}/{pdb_id}.prmtop"
            inpcrd_path = f"{sandbox_dir}/{pdb_id}.inpcrd"

            parmed_cm = pmd.load_file(prmtop_path, inpcrd_path)
            parmed_cm.save(f"{sandbox_dir}/topol.top")
            parmed_cm.save(f"{sandbox_dir}/{pdb_id}.gro")
        except Exception as e:
            return f"ParmEd failed: {type(e).__name__}: {e}, {result}"

        return f"tleap ran successfully with output: {result.stdout}. \n New files added: {sandbox_dir}/topol.top, {sandbox_dir}/{pdb_id}.gro"


def old_run_tleap_ligand(sandbox_dir: str, input_pdb: str, pdb_id: str, ligand_files: str | list[str], ligand_name: str) -> str:
    """
    Run tleap preparation using run_tleap.sh, for a protein-ligand complex.
    """
    # make sure it's a list
    if isinstance(ligand_files, str):
        ligand_files = [ligand_files]
    else:
        ligand_files = ligand_files
    # complex.pdb
    with (
        open(f"{sandbox_dir}/{input_pdb}", "r") as pdb_infile,
        open(f"{sandbox_dir}/complex.pdb", "w") as outfile,
    ):
        for line in pdb_infile:
            if not line.startswith("END"):
                outfile.write(line)

        with open(f"{sandbox_dir}/{ligand_files[0]}", "r") as ligand_infile:
            for line in ligand_infile:
                if line.startswith("HETATM"):
                    outfile.write(line)
            logger.info(f"Only the first ligand {ligand_files[0]} was added to complex.pdb as we will only simulate one ligand.")

        outfile.write("TER\n")
        outfile.write("END\n")

    complex_pdb = f"{sandbox_dir}/complex.pdb"
    # tleap with ligand
    script = constants.SCRIPTS_DIR / "run_tleap_ligand.sh"

    ligand_stem=Path(f"{ligand_files[0]}").stem
    if Path(f"{sandbox_dir}/{ligand_stem}_fixed.prepi").exists():
        prepi_file = f"{ligand_stem}_fixed.prepi"
    else:
        prepi_file = f"{ligand_stem}.prepi"

    result = subprocess.run(
        [str(script), sandbox_dir, complex_pdb, f"{ligand_stem}.frcmod", prepi_file],
        cwd=sandbox_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        # tleap often puts errors in stdout
        error_text = "\n".join(filter(None, [result.stderr, result.stdout]))
        return f"tleap run failed with error:\n{error_text}"
    else:
        # Parmed
        try:
            prmtop_path = f"{sandbox_dir}/complex.prmtop"
            inpcrd_path = f"{sandbox_dir}/complex.inpcrd"

            parmed_cm = pmd.load_file(prmtop_path, inpcrd_path)
            parmed_cm.save(f"{sandbox_dir}/topol.top")
            parmed_cm.save(f"{sandbox_dir}/complex.gro")
        except Exception as e:
            return f"ParmEd failed: {type(e).__name__}: {e}"

        return f"tleap ran successfully with output: {result.stdout}. \n New files added: {sandbox_dir}/topol.top, {sandbox_dir}/complex.gro"

def run_tleap_ligand(sandbox_dir: str, ligand_mol2_file: str, ligand_name: str) -> str:
    """
    Run tleap preparation to generate parameters (.itp and .prm file) for a ligand.
    """
    # The mol2 file must describe the selected protonation state of {prefix}_protonation.json
    prefix = Path(ligand_mol2_file).stem
    problems = check_ligand_consistency(sandbox_dir, prefix, mol2_file=ligand_mol2_file)
    if problems:
        return "Ligand parameterization failed with error: consistency check failed: " + "; ".join(problems)
    net_charge = load_record(sandbox_dir, prefix)["selected"]["net_charge"]

    cmd = shlex.split(f"parmchk2 -i {ligand_mol2_file} -f mol2 -s gaff2 -o {ligand_name}.frcmod")
    result = subprocess.run(cmd, cwd=sandbox_dir, capture_output=True, text=True)
    if result.returncode != 0:
        error_text = "\n".join(filter(None, [result.stderr, result.stdout]))
        logger.error(f"Error running parmchk2: {error_text}")
        return f"Ligand parameterization failed with error: parmchk2 failed:\n{error_text}"
    else:
        logger.info(f"{ligand_name}.frcmod content:\n{Path(f'{sandbox_dir}/{ligand_name}.frcmod').read_text()}")

    with open(f"{sandbox_dir}/tleap_ligand.in", "w") as file:
        file.write(f"""source leaprc.gaff2
LIG = loadmol2 {ligand_mol2_file}
check LIG
loadamberparams {ligand_name}.frcmod
saveoff LIG {ligand_name}.lib
saveamberparm LIG {ligand_name}.prmtop {ligand_name}.inpcrd
quit
""")
    cmd = "tleap -f tleap_ligand.in"
    run_1 = subprocess.run(shlex.split(cmd), cwd=sandbox_dir, capture_output=True, text=True)
    if run_1.returncode != 0:
        error_text = "\n".join(filter(None, [run_1.stderr, run_1.stdout]))
        logger.error(f"Error running tleap: {error_text}")
        return f"tleap run failed with error:\n{error_text}"
    else:
        logger.info(f"tleap ran successfully with output: {run_1.stdout}. \n New files added: {sandbox_dir}/{ligand_name}.prmtop, {sandbox_dir}/{ligand_name}.inpcrd, {sandbox_dir}/{ligand_name}.lib")

    prmtop_path = f"{sandbox_dir}/{ligand_name}.prmtop"
    inpcrd_path = f"{sandbox_dir}/{ligand_name}.inpcrd"
    try:
        parmed_cm = pmd.load_file(prmtop_path, inpcrd_path)
        parmed_cm.save(f"{sandbox_dir}/{ligand_name}.top", overwrite=True)
        parmed_cm.save(f"{sandbox_dir}/{ligand_name}.gro", overwrite=True)
    except Exception as e:
        return f"ParmEd failed: {type(e).__name__}: {e}"
    os.replace(f"{sandbox_dir}/{ligand_name}.top", f"{sandbox_dir}/{ligand_name}.itp")
    logger.info(f"ParmEd converted {prmtop_path} to {sandbox_dir}/{ligand_name}.itp and {sandbox_dir}/{ligand_name}.gro")

    with open(f"{sandbox_dir}/{ligand_name}.itp", "r") as f:
        lines = f.readlines()
    lines = lines[20:-7]  # remove first 20 and last 7
    with open(f"{sandbox_dir}/{ligand_name}.itp", "w") as f:
        f.writelines(lines)
    itp_file = f"{sandbox_dir}/{ligand_name}.itp"
    prm_file = f"{sandbox_dir}/{ligand_name}.prm"
    with open(itp_file, "r") as f:
        lines = f.readlines()
    # Find the index of the first line matching "[ moleculetype ]"
    first_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "[ moleculetype ]":
            first_idx = i
            break
    if first_idx is None:
        logger.error(f"No [ moleculetype ] found in {itp_file} — {prm_file} not written.")
        raise RuntimeError(f"No [ moleculetype ] found in {itp_file}")
    else:
        prm_lines = lines[:first_idx-2]
        itp_lines = lines[first_idx:]
        with open(itp_file, "w") as f:
            f.writelines(itp_lines)
        with open(prm_file, "w") as f:
            f.writelines(prm_lines)
        logger.info(f"{prm_file} and {itp_file} written.")
    # Edit charge of itp file so that it's an integer
    try:
        fix_last_atom_charge_ligand(sandbox_dir, ligand_name, expected_charge=net_charge)
    except ValueError as e:
        return f"Ligand parameterization failed with error: {e}"

    return f"tleap ran successfully. \n New files added: {sandbox_dir}/{ligand_name}.itp, {sandbox_dir}/{ligand_name}.prm, {sandbox_dir}/{ligand_name}.gro"

def fix_last_atom_charge_ligand(sandbox_dir:str, ligand_name:str, expected_charge: int | None = None) -> None:
    """
    Reads a GROMACS .itp file, finds the last atom in the [ atoms ] section,
    and adjusts its charge so that the total charge (qtot) is an integer (expected_charge if given).
    Raises ValueError if the correction is larger than MAX_CHARGE_CORRECTION: the charges then do not
    correspond to the selected protonation state and must not be silently fixed.
    """
    shutil.copyfile(f"{sandbox_dir}/{ligand_name}.itp", f"{sandbox_dir}/{ligand_name}_backup.itp")
    logger.info(f"Backup created: {sandbox_dir}/{ligand_name}_backup.itp")

    with open(f"{sandbox_dir}/{ligand_name}.itp", 'r') as f:
        lines = f.readlines()

    # Find the [ atoms ] section and collect atom lines
    in_atoms = False
    last_atom_idx = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('['):
            if 'atoms' in stripped:
                in_atoms = True
                continue
            elif in_atoms:
                # We've hit the next section
                break
        if in_atoms and stripped and not stripped.startswith(';'):
            last_atom_idx = i

    if last_atom_idx is None:
        logger.error("Could not find any atom entries in [ atoms ] section.")

    last_line = lines[last_atom_idx]

    # Parse the qtot comment to get current total charge
    if 'qtot' not in last_line:
        logger.error("Last atom line does not contain a qtot comment.")

    qtot_str = last_line.split('qtot')[-1].strip()
    qtot = float(qtot_str)

    # Correction needed to make qtot an integer
    target = expected_charge if expected_charge is not None else round(qtot)
    correction = target - qtot
    if abs(correction) > MAX_CHARGE_CORRECTION:
        logger.error(f"Total charge of {ligand_name}.itp is {qtot:.6f}, expected {target}: correction {correction:+.6f} refused")
        raise ValueError(f"Total charge of {ligand_name}.itp is {qtot:.6f} but the selected protonation state has net charge {target}. "
                         "The ligand charges do not correspond to the selected protonation state.")

    # Extract the charge field (column index 6, 0-based)
    # Preserve formatting by working with the raw columns
    parts = last_line.split()
    old_charge = float(parts[6])
    new_charge = old_charge + correction

    # Replace charge in line, preserving spacing
    # Find the position of the charge value and replace it precisely
    old_charge_str = parts[6]
    new_charge_str = f"{new_charge:.8f}"

    # Replace only the first occurrence of the charge string after the atom index columns
    # to avoid accidentally matching the same number elsewhere
    col_start = 0
    for _ in range(6):  # skip first 6 columns
        col_start = last_line.index(parts[_], col_start) + len(parts[_])

    charge_pos = last_line.index(old_charge_str, col_start)
    new_last_line = (
        last_line[:charge_pos]
        + new_charge_str
        + last_line[charge_pos + len(old_charge_str):]
    )

    # Update qtot in the comment
    new_qtot = qtot + correction  # should now be integer
    new_last_line = new_last_line.rsplit('qtot', 1)
    new_last_line = f"qtot {round(new_qtot):.6f}".join(new_last_line)

    lines[last_atom_idx] = new_last_line

    out = f"{sandbox_dir}/{ligand_name}.itp"
    with open(out, 'w') as f:
        f.writelines(lines)

    logger.info(f"Adjusted charge of last atom by {correction:+.8f} (qtot: {qtot:.6f} → {round(new_qtot):.6f})")
    logger.info(f"Written to: {out}")

def run_tleap_protein(sandbox_dir: str, pdb_id: str, ligand_name=None) -> str:
    """
    Run tleap preparation for the protein file, regardless of whether there is a ligand or not.
    """
    logger.info("----- Protein parameterization (ff14SB) -----")
    with open(f"{sandbox_dir}/tleap_protein.in", "w") as file:
        file.write(f"""source leaprc.protein.ff14SB
source leaprc.water.tip3p   # for TIP3P
PROT = loadpdb {pdb_id}_pdb2pqr.pdb
saveamberparm PROT {pdb_id}.prmtop {pdb_id}.inpcrd
savepdb PROT {pdb_id}_tleap.pdb
quit
""")
    cmd = "tleap -f tleap_protein.in"
    run_1 = subprocess.run(shlex.split(cmd), cwd=sandbox_dir, capture_output=True, text=True)
    if run_1.returncode != 0:
        error_text = "\n".join(filter(None, [run_1.stderr, run_1.stdout]))
        logger.error(f"Error running tleap: {error_text}")
        return f"tleap run failed with error:\n{error_text}"
    else:
        logger.info(f"tleap ran successfully with output: {run_1.stdout}. \n New files added: {sandbox_dir}/{pdb_id}.prmtop, {sandbox_dir}/{pdb_id}.inpcrd, {sandbox_dir}/{pdb_id}_tleap.pdb")

    prmtop_path = f"{sandbox_dir}/{pdb_id}.prmtop"
    inpcrd_path = f"{sandbox_dir}/{pdb_id}.inpcrd"
    try:
        parmed_cm = pmd.load_file(prmtop_path, inpcrd_path)
        parmed_cm.save(f"{sandbox_dir}/{pdb_id}.top", overwrite=True)
        parmed_cm.save(f"{sandbox_dir}/{pdb_id}_tleap.gro", overwrite=True)
    except Exception as e:
        return f"ParmEd failed: {type(e).__name__}: {e}"
    os.replace(f"{sandbox_dir}/{pdb_id}.top", f"{sandbox_dir}/topol.top")
    logger.info(f"ParmEd converted {prmtop_path} to {sandbox_dir}/topol.top and {sandbox_dir}/{pdb_id}_tleap.gro")

    # Edit the topology to add the posre blocks for the protein chains

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
    
    input_path  = Path(sandbox_dir) / "topol.top"
    backup_path = Path(sandbox_dir) / "topol_without_posre.top"

    if not input_path.exists():
        logger.error(f"{input_path} not found.")
        return f"ParmEd failed: {input_path} not found after conversion of {prmtop_path}."
 
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
        return f"tleap was run for the protein but no protein chain entries (system*) were found in the [ molecules ] section of {input_path}. Inspect {input_path} and {sandbox_dir}/{pdb_id}_pdb2pqr.pdb."
 
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
        return f"tleap was run for the protein but no [ moleculetype ] sections were found in {input_path}. Inspect {input_path} and {sandbox_dir}/{pdb_id}_pdb2pqr.pdb."
 
    positions = [m.start() for m in headers] + [end_of_moltypes]
    preamble  = text[: positions[0]]
 
    # ── 3. Walk each [ moleculetype ] block and insert POSRES if needed ───────
    segments: list[str]    = []
    inserted_protein: set  = set()   # unique types already given a posre block
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
        logger.debug(f'mol_names {mol_name}')
        logger.debug(f'prot mol types {protein_mol_types}')
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
                logger.info(f"Skipping duplicate [ moleculetype ] block for {mol_name}")
 
        segments.append(seg)
 
    # ── 4. Reassemble and write ───────────────────────────────────────────────
    modified_text = preamble + "".join(segments) + text[end_of_moltypes:]
    input_path.write_text(modified_text, encoding="utf-8")
    logger.info(f"Position restraint includes added for {input_path}: {', '.join(report) or 'none'}")

    # Edit topology to include parameters for ions, water, and posre block for ligand if ligand_name is provided

    input = f"{sandbox_dir}/topol.top"
    output = f"{sandbox_dir}/topol.top"

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        ligand_prm_lines_to_insert = [
            "; Include ligand parameters\n",
            f"#include \"{ligand_name}.prm\"\n",
            "\n",
        ]

        ligand_top_lines_to_insert = [
                "; Include ligand topology\n",
            f"#include \"{ligand_name}.itp\"\n",
            "\n",
            "; Ligand position restraints\n",
            "#ifdef POSRES\n",
            f"#include \"posre_{ligand_name}.itp\"\n",
            "#endif\n",
            "\n",
        ]

        last_ligand_line =[
            f"{ligand_name}                  1\n",
        ]

    atomstypes_to_insert = [
        "OW             8  16.000000  0.00000000  A     0.31507524       0.635968\n",
        "HW             1   1.008000  0.00000000  A              0              0\n",
        "Cl            17  35.450000  0.00000000  A      0.4477657     0.14891274\n",
        "Na            11  22.990000  0.00000000  A     0.24392807     0.36584603\n",
        "\n",
    ]

    lines_to_insert = [
        "[ moleculetype ]\n",
        "; Name            nrexcl\n",
        "SOL          3\n",
        "\n",
        "[ atoms ]\n",
        ";   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB\n",
        "; residue    1 SOL rtp SOL q 0.0\n",
        "    1         OW      1    SOL      O      1 -0.83400000  16.000000   ; qtot -0.834000\n",
        "    2         HW      1    SOL     H1      2 0.41700000   1.008000   ; qtot -0.417000\n",
        "    3         HW      1    SOL     H2      3 0.41700000   1.008000   ; qtot 0.000000\n",
        "\n",
        "#ifdef FLEXIBLE\n",
        "\n",
        "[ bonds ]\n",
        ";    ai     aj funct         c0         c1         c2         c3\n",
        "      2      3     1   0.15136 462750.400000\n",
        "      1      2     1   0.09572 462750.400000\n",
        "      1      3     1   0.09572 462750.400000\n",
        "\n",
        "#else\n",
        "\n",
        "[ settles ]\n",
        "; i     funct   doh     dhh\n",
        "1     1   0.09572000   0.15136000\n",
        "\n",
        "#endif\n",
        "\n",
        "[ exclusions ]\n",
        "1  2  3\n",
        "2  1  3\n",
        "3  1  2\n",
        "\n",
        "[ moleculetype ]\n",
        "; Name            nrexcl\n",
        "Na          3\n",
        "\n",
        "[ atoms ]\n",
        ";   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB\n",
        "; residue    1 Na rtp Na q 1.0\n",
        "    1        Na      1    Na     Na       1 1.00000000  22.990000   ; qtot 1.000000\n",
        "\n",
        "[ moleculetype ]\n",
        "; Name            nrexcl\n",
        "Cl           3\n",
        "\n",
        "[ atoms ]\n",
        ";   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB\n",
        "; residue    1 Cl  rtp Cl  q -1.0\n",
        "    1        Cl       1    Cl     Cl       1 -1.00000000  35.450000   ; qtot -1.000000\n",
        "\n",
    ]

    with open(input, "r") as f:
        lines = f.readlines()

    # Find the index of the last line matching "[ system ]"
    last_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "[ system ]":
            last_idx = i

    moleculetype_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "[ moleculetype ]":
            moleculetype_idx = i
            break

    if last_idx is None or moleculetype_idx is None:
        logger.error("No [ system ] or [ moleculetype ] found in topol.top — nothing inserted.")
    else:
        if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
            new_lines = lines[:moleculetype_idx-2] + atomstypes_to_insert + ligand_prm_lines_to_insert + lines[moleculetype_idx:last_idx] + ligand_top_lines_to_insert + lines_to_insert + lines[last_idx:] + last_ligand_line
        else:
            new_lines = lines[:moleculetype_idx-2] + atomstypes_to_insert + lines[moleculetype_idx:last_idx] + lines_to_insert + lines[last_idx:]
        with open(output, "w") as f:
            f.writelines(new_lines)
        logger.info(f"{output} updated with water/ion parameters and {ligand_name} .prm/.itp/posre includes")

    # Create complex.gro if ligand
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        ligand_file = f"{sandbox_dir}/{ligand_name}.gro"
        protein_file = f"{sandbox_dir}/{pdb_id}_tleap.gro"
        complex_file = f"{sandbox_dir}/complex.gro"

        with open(ligand_file, "r") as f:
            lines_ligand = f.readlines()
        with open(protein_file, "r") as file:
            # protein lines are the lines in protein-file that don't contain HOH
            lines_protein = [line for line in file.readlines() if "HOH" not in line and "WAT" not in line and "SOL" not in line]
        with open(protein_file, "r") as file:
            lines_water = [line for line in file.readlines() if "HOH" in line or "WAT" in line or "SOL" in line]
            logger.info(f"{len(lines_water)} water atom lines found in protein file")

        atoms = int(lines_protein[1]) + int(lines_ligand[1])
        logger.info(f"Total number atoms in protein-ligand system: {atoms}")

        top_line = [
            f"Protein ligand complex {pdb_id}-{ligand_name}\n",
        ]
        atom_line = [
            f"{atoms}\n",
        ]

        new_lines = top_line + atom_line + lines_protein[2:-1] + lines_ligand[2:-1] + lines_water + list(lines_protein[-1])

        with open(complex_file, "w") as f:
            f.writelines(new_lines)
        logger.info(f"{complex_file} created")

        # Renumber the complex file since the ligand was added between the protein and water molecules
        shutil.copyfile(complex_file, f"{sandbox_dir}/complex_backup.gro")
        logger.info(f"Backup created: {sandbox_dir}/complex_backup.gro")
        with open(f'{sandbox_dir}/complex.gro', 'r') as f:
            lines = f.readlines()

        # First and last lines are header and box vectors
        header = lines[0]
        atom_count_line = lines[1]
        box_line = lines[-1]
        atom_lines = lines[2:-1]

        output_lines = []
        atom_num = 0
        res_num = 0

        for line in atom_lines:
            if len(line) < 20:
                output_lines.append(line)
                continue

            res_num = line[0:5].strip()
            res_name     = line[5:10].strip()
            atom_name    = line[10:15].strip()
            rest         = line[20:]  # coordinates onwards

            # Atom number: increment, wrap back to 1 after 99999
            atom_num += 1
            if atom_num >= 100000:
                atom_num = 1

            new_line = f'{int(res_num):5d}{res_name:<5s}{atom_name:>5s}{atom_num:5d}{rest}'
            output_lines.append(new_line)

        # Update atom count line in case it changed
        with open(f'{sandbox_dir}/complex.gro', 'w') as f:
            f.write(header)
            f.write(f'{len(output_lines):5d}\n')
            f.writelines(output_lines)
            f.write(box_line)

        logger.info(f"Done. {len(output_lines)} atoms written to complex.gro. Backup saved to complex_backup.gro")

        return f"tleap was run for the protein. Lines were added to topol.top. Protein-ligand complex created: {sandbox_dir}/complex.gro with {len(output_lines)} atoms. Backup saved to {sandbox_dir}/complex_backup.gro"

    with open(f"{sandbox_dir}/{pdb_id}_tleap.gro", "r") as f:
        n_atoms = int(f.readlines()[1])
    return f"tleap was run for the protein. Lines were added to topol.top. Protein-only system created: {sandbox_dir}/{pdb_id}_tleap.gro with {n_atoms} atoms."