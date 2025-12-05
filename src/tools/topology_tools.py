import subprocess
from pathlib import Path


def fix_topology_negative(topfile: str, sandbox_dir: str) -> str:
    """
    This script should be used if ff14sb AMBER force field is used in tleap.
    This script should be used if topfile has a net negative charge.
    This script adds the missing water and ions parameters to the topology file.
    Args:
        topfile (str): Path to the topology file to be fixed.
    Returns:
        topol.top (str): Fixed topology file.
    """
    output_file = Path(sandbox_dir) / "topol.top"

    # Read the original file
    with open(topfile, "r") as f:
        lines = f.readlines()

    # Lines to insert
    insert_lines = [
        "HW             1   1.008     0.0000      A     0.00000e+00   0.00000e+00\n"
        "OW             8   16.00     0.0000      A     3.15061e-01   6.36386e-01\n"
        "\n",
        "; Include topology for water\n",
        '#include "amber99.ff/tip3p.itp"\n',
        "\n",
        "[ moleculetype ]\n",
        "; molname       nrexcl\n",
        "NA              1\n",
        "\n",
        "[ atoms ]\n",
        "; id    at type         res nr  residu name     at name  cg nr  charge\n",
        "1       NA              1       NA              NA       1      1.00000\n",
    ]

    # Find [ atomtypes ] section
    start_index = None
    for i, line in enumerate(lines):
        if line.strip() == "[ atomtypes ]":
            start_index = i
            break

    if start_index is None:
        raise ValueError("No [ atomtypes ] section found in the file.")

    # Find last atom line (last line before blank line)
    end_index = start_index + 1
    while end_index < len(lines) and lines[end_index].strip() != "":
        end_index += 1

    # Insert lines
    lines = lines[:end_index] + insert_lines + lines[end_index:]

    # Write to output
    with open(output_file, "w") as f:
        f.writelines(lines)
    return "Successfully added missing water ions and fixed the topology. Output file is saved to topol.top"


def fix_topology_positive(topfile: str, sandbox_dir: str) -> str:
    """
    This script should be used if ff14sb AMBER force field is used in tleap.
    This script should be used if topfile has a net positive charge.
    This script adds the missing water and ions parameters to the topology file.
    Args:
        topfile (str): Path to the topology file to be fixed.
    Returns:
        topol.top (str): Fixed topology file.
    """
    output_file = sandbox_dir / "topol.top"

    # Read the original file
    with open(topfile, "r") as f:
        lines = f.readlines()

    # Lines to insert
    insert_lines = [
        "HW             1   1.008     0.0000      A     0.00000e+00   0.00000e+00\n"
        "OW             8   16.00     0.0000      A     3.15061e-01   6.36386e-01\n"
        "\n",
        "; Include topology for water\n",
        '#include "amber99.ff/tip3p.itp"\n',
        "\n",
        "[ moleculetype ]\n",
        "; molname       nrexcl\n",
        "CL              1\n",
        "\n",
        "[ atoms ]\n",
        "; id    at type         res nr  residu name     at name  cg nr  charge\n",
        "1       CL              1       CL              CL       1     -1.00000",
    ]

    # Find [ atomtypes ] section
    start_index = None
    for i, line in enumerate(lines):
        if line.strip() == "[ atomtypes ]":
            start_index = i
            break

    if start_index is None:
        raise ValueError("No [ atomtypes ] section found in the file.")

    # Find last atom line (last line before blank line)
    end_index = start_index + 1
    while end_index < len(lines) and lines[end_index].strip() != "":
        end_index += 1

    # Insert lines
    lines = lines[:end_index] + insert_lines + lines[end_index:]

    # Write to output
    with open(output_file, "w") as f:
        f.writelines(lines)
    return "Successfully added missing water ions and fixed the topology. Output file is saved to topol.top"


def analyze_Gromacs(sandbox_dir: str) -> None:
    """
    This script analyzes the trajectory from the production MD simulation using Gromacs.
    It calculates RMSD, RMSF, Radius of Gyration, and secondary structure content.
    It uses the previously generated md.xtc and topol.top files as inputs.
    The script is ran after the production_Gromacs script and requires md.gro and md.xtc files.
    Args:
        None
    Returns:
        rmsd.xvg (str): RMSD plot data.
        rmsf.xvg (str): RMSF plot data.
        gyrate.xvg (str): Radius of Gyration plot data.
    """
    rmsd_command = "echo 4 4 | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns"
    subprocess.run([rmsd_command], cwd=sandbox_dir, check=True)
    rmsf_command = "echo 4 | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg"
    subprocess.run([rmsf_command], cwd=sandbox_dir, check=True)
    rg_command = "echo 4 | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg"
    subprocess.run([rg_command], cwd=sandbox_dir, check=True)
    return None
