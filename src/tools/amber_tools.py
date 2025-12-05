import subprocess
from pathlib import Path
import parmed as pmd  # type: ignore
from src import constants
from src.utils import get_class_logger

logger = get_class_logger(__name__)


def run_tleap(sandbox_dir: str, input_pdb: str, pdb_id: str) -> str:
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


def run_tleap_ligand(sandbox_dir: str, input_pdb: str, pdb_id: str, ligand_file: str, ligand_name: str) -> str:
    """
    Run tleap preparation using run_tleap.sh, for a protein-ligand complex.
    """
    # complex.pdb
    with (
        open(f"{sandbox_dir}/{input_pdb}", "r") as pdb_infile,
        open(f"{sandbox_dir}/{ligand_file}", "r") as ligand_infile,
        open(f"{sandbox_dir}/complex.pdb", "w") as outfile,
    ):
        for line in pdb_infile:
            if not line.startswith("END"):
                outfile.write(line)

        for line in ligand_infile:
            if line.startswith("HETATM"):
                outfile.write(line)

        outfile.write("TER\n")
        outfile.write("END\n")

    complex_pdb = f"{sandbox_dir}/complex.pdb"
    # tleap with ligand
    script = constants.SCRIPTS_DIR / "run_tleap_ligand.sh"

    if Path(f"{sandbox_dir}/{ligand_name}_fixed.prepi").exists():
        prepi_file = f"{ligand_name}_fixed.prepi"
    else:
        prepi_file = f"{ligand_name}.prepi"

    result = subprocess.run(
        [str(script), sandbox_dir, complex_pdb, ligand_name, prepi_file],
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
