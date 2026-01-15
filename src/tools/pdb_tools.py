import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBList
import subprocess, shlex
import MDAnalysis as mda  # type: ignore
import numpy as np
from collections import defaultdict
import traceback
from src.utils import get_class_logger

logger = get_class_logger(__name__)


class NoLigand(Exception):
    def __init__(self, message):
        super().__init__(message)


def fix_pdb_file(input_pdb: str, fixed_pdb: str) -> str:
    """
    This is a publically available API that fixes a PDB file using PDBFixer and writes the output to a new file. This
    function should only be called when the input PDB file is known to have  issues, otherwise it may introduce new problems.
    Actions applied:
      - Removes waters
      - Finds missing residues (reports them but does not add)
      - Finds missing atoms
      - Adds missing heavy atoms (but no hydrogens)

    Args:
        input_pdb (str): The path where the faulty PDB file is located.
        fixed_pdb (str): The path where we save the fixed PDB file.

    """
    # Load a PDB file
    fixer = PDBFixer(filename=input_pdb)

    # Remove water molecules
    fixer.removeHeterogens(keepWater=False)

    # Identify problems
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    # PDBFixer cannot auto-add missing residues; just reports them.
    if fixer.missingResidues:
        logger.warning("Missing residues detected:", fixer.missingResidues)

    # Apply fixes (only adds missing atoms, no hydrogens)
    fixer.addMissingAtoms()

    # Save fixed PDB
    with open(fixed_pdb, "w") as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out)

    return f"Fixed original PDB file {input_pdb} and saved to new file {fixed_pdb}"


def fetch_and_save_pdb(sandbox_dir: str, pdb_id: str, output_pdb: str) -> str:
    """
    This is a publically available API that fetches a PDB file from the RCSB server and saves it locally.
    Then parses it with Bio.PDB.PDBParser to ensure it's valid. This function
    should only be called if you were not provided with a local PDB file.

    Args:
        pdb_id (str): The 4-character PDB ID (e.g., '1abc').
        output_pdb (str): The path to save the fetched PDB file. This file should end in ".pdb".
    """
    pdb_id = pdb_id.upper()

    # Download PDB file
    try:
        pdbl = PDBList()
        fetched_file = pdbl.retrieve_pdb_file(pdb_id, pdir=sandbox_dir, file_format="pdb")

        output_pdb = sandbox_dir / f"{pdb_id}.pdb"

        with open(fetched_file, "r") as infile, open(output_pdb, "w") as outfile:
            outfile.write(infile.read())

        return f"PDB {pdb_id} downloaded successfully to {output_pdb}"

    except Exception:
        return f"Error fetching PDB {pdb_id}: {traceback.format_exc()}"

def check_pdb_ligand(sandbox_dir: str, pdb_id: str, ligand_name: str = None) -> str:
    """
    This is a publically available API that checks if a PDB file is valid by attempting to parse it with Bio.PDB.PDBParser.

    Args:
        pdb_id (str): The 4-character PDB ID (e.g., '1abc').
        ligand_name (str): Optional: the name of the ligand if a protein-ligand complex should be simulated.
    """
    pdb_file = sandbox_dir / f"{pdb_id}.pdb"

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        # Check if ligand is present in the PDB file
        with open(pdb_file, "r") as f:
            lines = f.readlines()
            ligand_present = any(line.startswith("HETATM") and ligand_name in line for line in lines)

        if not ligand_present:
            logger.info(f"Ligand {ligand_name} not found in PDB file {pdb_file}.")
            raise NoLigand(f"Ligand {ligand_name} not found in PDB file {pdb_file}.")
            
        #Check if ligand is covalent
        with open(pdb_file, "r") as f:
            lines = f.readlines()
            ligand_covalent = any(line.startswith("LINK") and ligand_name in line for line in lines)
        
        if ligand_covalent:
            logger.info(f"Ligand {ligand_name} appears to be covalently bound in PDB file {pdb_file}. DynaMate doesn't support the parameterization of covalently bound ligands. This system cannot be processed.")
            raise NoLigand(f"Ligand {ligand_name} appears to be covalently bound in PDB file {pdb_file}. DynaMate doesn't support the parameterization of covalently bound ligands. This system cannot be processed.")

        # Count number of ligands
        with open(pdb_file, "r") as f:
            lines = f.readlines()
            ligand_count = sum(1 for line in lines if line.startswith("HET ") and ligand_name in line)
            logger.info(f"There is(are) {ligand_count} ligand(s) called {ligand_name} in {pdb_file}.")

    # Check for modified residues
    with open(pdb_file, "r") as f:
        lines = f.readlines()
        modified_residues = [line for line in lines if line.startswith("MODRES")]
        logger.info("There are ", len(modified_residues), "modified residues, which are", modified_residues, ". This should be checked and the corresponding residues modified to standard residues. If they can't be modified to standard residues, the system can't be processed.")

    return "PDB file check completed successfully."

def prepare_pdb_file_ligand(sandbox_dir: str, pdb_id: str, ligand_name: str = None) -> str:
    """
    Takes input PDB file, extract ligand to ligand_name.pdb.
    Removes HETATM, CONECT and MASTER lines from input_pdb and saves to prepared_pdb.
    Protonates ligand and pH=7 and saves to ligand_name_h.pdb.
    Args:
        input_pdb (str): The path where the input PDB file is located.
        prepared_pdb (str): The path where we save the prepared PDB file.
        ligand_name (str): The name of the ligand to extract.
        ligand_pdb (str): The path where we save the extracted ligand PDB file.
        ligand_pdb_h (str): The path where we save the protonated ligand PDB file.
    """
    # Prepare PDB
    with (
        open(f"{sandbox_dir}/{pdb_id}.pdb", "r") as infile,
        open(f"{sandbox_dir}/{pdb_id}_prepared.pdb", "w") as outfile,
    ):
        for line in infile:
            if not (line.startswith("HETATM") or line.startswith("CONECT") or line.startswith("MASTER")):
                outfile.write(line)
    logger.info(f"Prepared PDB file saved to {sandbox_dir}/{pdb_id}_prepared.pdb")

    # Extract ligand
    if (ligand_name is not None) and (ligand_name != "XXX") and (ligand_name != "None") and (ligand_name != "None_h"):
        # Count number of ligands
        resnums = set()
        with open(f"{sandbox_dir}/{pdb_id}.pdb") as f:
            for line in f:
                if line.startswith("HETATM") and ligand_name in line:
                    resnum = int(line[22:26])
                    resnums.add(resnum)
        print("resnums is: ", resnums)
        num_ligands = len(resnums)
        logger.info(f"IMPORTANT: Number of ligands {ligand_name} found: {num_ligands}")

        if num_ligands == 0:
            logger.info(f"Ligand {ligand_name} not found in PDB file {sandbox_dir}/{pdb_id}.pdb. You can either proceed without a ligand, check the ligand name provided or check the PDB file.")
            return f"Ligand {ligand_name} not found in PDB file {sandbox_dir}/{pdb_id}.pdb. You can either proceed without a ligand, check the ligand name provided or check the PDB file"

        # CASE 1: only one ligand
        if num_ligands == 1:
            ligand_pdb_file = f"{sandbox_dir}/{ligand_name}.pdb"
            ligand_pdb_files_list = [ligand_pdb_file]
            with open(f"{sandbox_dir}/{pdb_id}.pdb", "r") as infile, open(ligand_pdb_file, "w") as outfile:
                for line in infile:
                    if line.startswith("HETATM") and ligand_name in line:
                        outfile.write(line)
            logger.info(f"Extracted ligand {ligand_name} to {ligand_pdb_file}")
        
        # CASE 2: multiple ligands (split by residue number)
        else:
            ligand_pdb_files_list = []
            ligands = defaultdict(list)

            # Collect HETATM lines by residue number
            with open(f"{sandbox_dir}/{pdb_id}.pdb", "r") as infile:
                for line in infile:
                    if line.startswith("HETATM") and ligand_name in line:
                        resnum = int(line[22:26])  # residue number column
                        ligands[resnum].append(line)

            # Write one file per ligand
            for i, (resnum, atom_lines) in enumerate(ligands.items(), start=1):
                ligand_pdb_file = f"{sandbox_dir}/{ligand_name}_{i}.pdb"
                ligand_pdb_files_list.append(ligand_pdb_file)
                with open(ligand_pdb_file, "w") as outfile:
                    outfile.writelines(atom_lines)
                logger.info(f"Extracted ligand {ligand_name} residue {resnum} to {ligand_pdb_file}")

    # Protonate ligand
    list_protonated_files = []
    if (ligand_name is not None) and (ligand_name != "XXX") and (ligand_name != "None") and (ligand_name != "None_h"):
        for ligand_pdb_file in ligand_pdb_files_list: # loop over all extracted ligands
            if num_ligands == 1:
                protonated_file = f"{sandbox_dir}/{ligand_name}_h.pdb"
            else:
                index = ligand_pdb_file.split("_")[-1].split(".")[0]  # get index from filename
                protonated_file = f"{sandbox_dir}/{ligand_name}_{index}_h.pdb"
                list_protonated_files.append(f"{ligand_name}_{index}_h.pdb")
            
            cmd = shlex.split(f"obabel {ligand_pdb_file} -O {protonated_file} -p7")
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with open(protonated_file, "r") as infile:
                lines = infile.readlines()
                filtered_lines = [line for line in lines if not (line.startswith("CONECT") or line.startswith("MASTER"))]
                new_filtered_lines = [line.replace("UNL", ligand_name).replace("UNK", ligand_name) for line in filtered_lines]
            with open(protonated_file, "w") as outfile:
                outfile.writelines(new_filtered_lines)

        # Rewrite atoms names in the ligand
        ELEMENTS = {
            "H",
            "He",
            "Li",
            "Be",
            "B",
            "C",
            "N",
            "O",
            "F",
            "Ne",
            "Na",
            "Mg",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ar",
            "K",
            "Ca",
            "Sc",
            "Ti",
            "V",
            "Cr",
            "Mn",
            "Fe",
            "Co",
            "Ni",
            "Cu",
            "Zn",
            "Ga",
            "Ge",
            "As",
            "Se",
            "Br",
            "Kr",
            "Rb",
            "Sr",
            "Y",
            "Zr",
            "Nb",
            "Mo",
            "Tc",
            "Ru",
            "Rh",
            "Pd",
            "Ag",
            "Cd",
            "In",
            "Sn",
            "Sb",
            "Te",
            "I",
            "Xe",
            "Cs",
            "Ba",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
            "Hf",
            "Ta",
            "W",
            "Re",
            "Os",
            "Ir",
            "Pt",
            "Au",
            "Hg",
            "Tl",
            "Pb",
            "Bi",
            "Po",
            "At",
            "Rn",
            "Fr",
            "Ra",
            # Transition metals sometimes used in PDBs
            "U",
            "Pu",
        }

        def normalize_element(e):
            """Convert PDB element field to normalized chemical element symbol."""
            e = e.strip().capitalize()
            if len(e) == 2:
                return e[0] + e[1].lower()
            return e

        with open(protonated_file, "r") as infile:
            lines = infile.readlines()

        counters = defaultdict(int)
        new_lines = []

        for idx, line in enumerate(lines, start=1):
            if line.startswith(("ATOM", "HETATM")):
                raw_name = line[12:16].strip()
                raw_element = line[76:78]  # columns 77–78
                element = normalize_element(raw_element)

                if element not in ELEMENTS:
                    logger.error(
                        f'Unknown element "{element}" (from raw field "{raw_element.strip()}") found at line {idx}. Atom name in file: "{raw_name}". Please check the ligand PDB: unexpected element.'
                    )

                    # keep the line unchanged
                    new_lines.append(line)
                    continue

                # Known element → rename it
                counters[element] += 1
                new_name = f"{element}{counters[element]}"

                # Replace atom name in columns 13–16
                line = f"{line[:12]}{new_name:>4}{line[16:]}"

            new_lines.append(line)

        with open(protonated_file, "w") as outfile:
            outfile.writelines(new_lines)

        logger.info("Atom renaming of ligand completed.")

        if num_ligands == 1:
            return f"Successfully Prepared PDB structure with a ligand and saved the extracted protein PDB file to {sandbox_dir}/{pdb_id}_prepared.pdb and the protonated ligand PDB file to {sandbox_dir}/{ligand_name}_h.pdb. Ligand was protonated at pH=7 and atom names were cleaned (renumbered)"
        if num_ligands > 1:
            return f"Successfully Prepared PDB structure with {num_ligands} ligands and saved the extracted protein PDB file to {sandbox_dir}/{pdb_id}.pdb and the {num_ligands} protonated ligand PDB files to {sandbox_dir}/{list_protonated_files}. This list of {num_ligands} protonated files: {list_protonated_files} is IMPORTANT and should be the input parameter for future functions. The extracted pdb file was saved to {sandbox_dir}/{pdb_id}_prepared.pdb Ligands were protonated at pH=7 and atom names were cleaned (renumbered)"


    return f"Successfully Prepared PDB structure without a ligand and saved the extracted PDB file to {sandbox_dir}/{pdb_id}_prepared.pdb"


def add_caps(sandbox_dir: str, input_pdb: str, pdb_id: str) -> str:
    """
    Adds ACE and NME caps to the N- and C-termini of the protein in the input PDB file.
    Saves the modified structure to {sandbox_dir}/{pdb_id}_prepared_capped.pdb.

    Args:
        input_pdb (str): The path where the input PDB file is located.
        pdb_id (str): The PDB ID.
        sandbox_dir (str): the directory where we add and modify files.
    """
    with open(f"{sandbox_dir}/{input_pdb}") as pdbfile:
        for line in pdbfile:
            if line.startswith("HETATM") or line.startswith("CONECT") or line.startswith("MASTER"):
                pdbfile.close()
                logger.warning("Input PDB file contains HETATM, CONECT or MASTER lines. Please prepare the PDB file first to remove these lines. If the PDB file has already been prepared with the prepare_pdb_file_ligand function, use the correct parameters when calling this tool or check that it has been prepared correctly.")
                return "Error: Input PDB file contains HETATM, CONECT or MASTER lines. Please prepare the PDB file first to remove these lines. If the PDB file has already been prepared with the prepare_pdb_file_ligand function, use the correct parameters when calling this tool or check that it has been prepared correctly."    

    def create_universe(n_atoms, name, resname, positions, resids, segid):
        u_new = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=n_atoms,
            atom_resindex=np.arange(n_atoms),
            residue_segindex=np.arange(n_atoms),
            n_segments=n_atoms,
            trajectory=True,
        )  # necessary for adding coordinate

        u_new.add_TopologyAttr("name", name)
        u_new.add_TopologyAttr("resid", resids)
        u_new.add_TopologyAttr("resname", resname)
        u_new.atoms.positions = positions
        u_new.add_TopologyAttr("segid", n_atoms * [segid])
        u_new.add_TopologyAttr("chainID", n_atoms * [segid])

        return u_new

    def get_nme_pos(end_residue):
        if "OXT" in end_residue.names:
            index = np.where(end_residue.names == "OXT")[0][0]
            N_position = end_residue.positions[index]
            index_c = np.where(end_residue.names == "C")[0][0]
            carbon_position = end_residue.positions[index_c]
            vector = N_position - carbon_position
            vector /= np.sqrt(sum(vector**2))

            C_position = N_position + vector * 1.36

            return N_position, C_position

        else:
            # find midpoint of O and CA
            index_o = np.where(end_residue.names == "O")[0][0]
            index_ca = np.where(end_residue.names == "CA")[0][0]

            mid_point = (end_residue.positions[index_o] + end_residue.positions[index_ca]) / 2

            # find vector connecting mid_point and C
            index_c = np.where(end_residue.names == "C")[0][0]
            vector = end_residue.positions[index_c] - mid_point
            vector /= np.sqrt(sum(vector**2))
            N_position = end_residue.positions[index_c] + 1.36 * vector
            ##
            C_position = N_position + 1.36 * vector

        return N_position, C_position

    def get_ace_pos(end_residue):
        index_ca = np.where(end_residue.names == "CA")[0][0]
        index_n = np.where(end_residue.names == "N")[0][0]
        vector = end_residue.positions[index_n] - end_residue.positions[index_ca]
        vector /= np.sqrt(sum(vector**2))

        C1_position = end_residue.positions[index_n] + 1.36 * vector

        xa, ya, za = end_residue.positions[index_ca]
        xg, yg, zg = C1_position

        # arbritray unit vector
        # create an arbritray orientaiton for the ACE residue
        # does not really matter
        orientation = np.array([2 * np.random.rand() - 1, 2 * np.random.rand() - 1, 2 * np.random.rand() - 1])
        nx, ny, nz = orientation / np.sqrt(sum(orientation**2))

        ## The carbon and oxygen are placed on the vertices of an equilatrel triangle
        # with another vertex as the Nitrogen atom and the C as the centroid
        # The plane of the triangle is placed in an arbritrary orientation as defined before
        # The orientation does not matter
        ######################################
        x1 = xg - (xa - xg) / 2 + np.sqrt(3) * (ny * (za - zg) - nz * (ya - yg)) / 2
        y1 = yg - (ya - yg) / 2 + np.sqrt(3) * (nz * (xa - xg) - nx * (za - zg)) / 2
        z1 = zg - (za - zg) / 2 + np.sqrt(3) * (nx * (ya - yg) - ny * (xa - xg)) / 2

        ## second coordinate
        x2 = xg - (xa - xg) / 2 - np.sqrt(3) * (ny * (za - zg) - nz * (ya - yg)) / 2
        y2 = yg - (ya - yg) / 2 - np.sqrt(3) * (nz * (xa - xg) - nx * (za - zg)) / 2
        z2 = zg - (za - zg) / 2 - np.sqrt(3) * (nx * (ya - yg) - ny * (xa - xg)) / 2

        C2_position = np.array([x1, y1, z1])
        O_position = np.array([x2, y2, z2])

        ### rescale distances, the above points may be a bit far apart like 2.1 angstrom but usual bonds are 1.4 or so
        ## Therefore we shrink it
        #  C positinos

        vector = C2_position - C1_position
        vector /= np.sqrt(sum(vector**2))

        C2_position = C1_position + 1.36 * vector

        # O positions
        vector = O_position - C1_position
        vector /= np.sqrt(sum(vector**2))

        O_position = C1_position + 1.36 * vector

        return C1_position, C2_position, O_position

    # Load pdb file
    u = mda.Universe(f"{sandbox_dir}/{input_pdb}")

    # Access each fragment separately
    res_start = 0
    segment_universes = []

    for seg in u.segments:
        chain = u.select_atoms(f"segid {seg.segid}")

        # Add ACE
        resid_c = chain.residues.resids[0]
        end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")
        ace_positions = get_ace_pos(end_residue)
        ace_names = ["C", "CH3", "O"]
        resid = chain.residues.resids[0]
        kwargs = dict(
            n_atoms=len(ace_positions),
            name=ace_names,
            resname=len(ace_names) * ["ACE"],
            positions=ace_positions,
            resids=resid * np.ones(len(ace_names)),
            segid=chain.segids[0],
        )

        ace_universe = create_universe(**kwargs)

        # Add NME
        resid_c = chain.residues.resids[-1]
        end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")

        nme_positions = get_nme_pos(end_residue)
        nme_names = ["N", "C"]

        resid = chain.residues.resids[-1] + 2

        kwargs = dict(
            n_atoms=len(nme_names),
            name=nme_names,
            resname=len(nme_names) * ["NME"],
            positions=nme_positions,
            resids=resid * np.ones(len(nme_names)),
            segid=chain.segids[0],
        )

        nme_universe = create_universe(**kwargs)
        ## Merge Universe
        if "OXT" in end_residue.names:
            index = np.where(end_residue.names == "OXT")[0][0]
            OXT = end_residue[index]

            Chain = u.select_atoms(f"segid {seg.segid} and not index {OXT.index}")

        else:
            Chain = u.select_atoms(f"segid {seg.segid}")

        ### Merge ACE, Protien and NME
        u_all = mda.Merge(ace_universe.atoms, Chain, nme_universe.atoms)

        # to renumber residues
        resids_ace = [res_start + 1, res_start + 1, res_start + 1]
        resids_pro = np.arange(resids_ace[0] + 1, Chain.residues.n_residues + resids_ace[0] + 1)
        resids_nme = [resids_pro[-1] + 1, resids_pro[-1] + 1]

        u_all.atoms.residues.resids = np.concatenate(
            [resids_ace, resids_pro, resids_nme]
        )  # np.arange (1+res_start, len(u_all.atoms.residues.resids)+res_start+1)
        res_start = u_all.atoms.residues.resids[-1]
        segment_universes.append(u_all)

    ## Join all the universes
    all_uni = mda.Merge(*(seg.atoms for seg in segment_universes))
    all_uni.atoms.write(f"{sandbox_dir}/{pdb_id}_prepared_capped.pdb")

    def insert_ter(pdb_in):
        with open(pdb_in, "r") as f:
            lines = f.readlines()

        output_lines = []
        prev_resname = None

        for i, line in enumerate(lines):
            record = line[:6].strip()

            if record in ("ATOM", "HETATM"):
                resname = line[17:20].strip()

                # Insert TER if previous residue was NME and current is ACE
                if prev_resname == "NME" and resname == "ACE":
                    output_lines.append("TER\n")

                output_lines.append(line)
                prev_resname = resname

            elif record == "END":
                # Insert TER before END if the last record wasn't TER already
                if not output_lines[-1].startswith("TER"):
                    output_lines.append("TER\n")
                output_lines.append(line)

            else:
                output_lines.append(line)

        if os.path.exists(pdb_in):
            os.remove(pdb_in)

        with open(pdb_in, "w") as f:
            f.writelines(output_lines)

    insert_ter(f"{sandbox_dir}/{pdb_id}_prepared_capped.pdb")

    return f"Successfully added ACE and NME caps to the N- and C-termini of the protein and saved to {sandbox_dir}/{pdb_id}_prepared_capped.pdb"


def rename_histidines(sandbox_dir: str, input_pdb: str, pdb_id: str) -> str:
    """
    Renames histidines HIS to account for their correct protonation.
    This tool should be ran after the protein has been capped with the add_caps function for both protein-only and protein-ligand systems.
    HIS will be renamed to HID if it contains HD1 and not HE2.
    HIS will be renamed to HIE if it contains HE2 and not HD1.
    HIS will be renamed to HIP if it contains both HD1 and HE2.

    Args:
        input_pdb (str): The path where the input PDB file is located.
        pdb_id (str): pdb id of the protein.
    """

    # Read all lines
    with open(f"{sandbox_dir}/{input_pdb}", "r") as f:
        pdb_lines = f.readlines()

    # Collect atom names per residue (keyed by resname, chain, resnum)
    residues = defaultdict(list)
    for i, line in enumerate(pdb_lines):
        if line.startswith(("ATOM", "HETATM")):
            res_name = line[17:20]
            chain_id = line[21]
            res_num = line[22:26]
            key = (chain_id, res_num)
            residues[key].append((i, line))

    # Determine new residue names for HIS residues
    his_renames = {}
    for key, atom_list in residues.items():
        # get residue name (all atoms share it)
        res_name = atom_list[0][1][17:20].strip()
        if res_name == "HIS":
            atom_names = {l[12:16].strip() for _, l in atom_list}
            if "HD1" in atom_names and "HE2" not in atom_names:
                new_name = "HID"
            elif "HD1" in atom_names and "HE2" in atom_names:
                new_name = "HIP"
            elif "HD1" not in atom_names and "HE2" in atom_names:
                new_name = "HIE"
            else:
                new_name = "HIS"
            his_renames[key] = new_name

    # Modify the lines in place
    for key, new_name in his_renames.items():
        for i, line in residues[key]:
            pdb_lines[i] = line[:17] + new_name.ljust(3) + line[20:]

    # Write output file
    output_path = os.path.join(sandbox_dir, f"{pdb_id}_prepared_capped_his.pdb")
    with open(output_path, "w") as f:
        f.writelines(pdb_lines)

    return f"Successfully renamed histidines HIS to account for their correct protonation in the PDB files and saved to {output_path}"
