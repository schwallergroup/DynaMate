import os
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBList
import subprocess, shlex
import MDAnalysis as mda  # type: ignore
import numpy as np
from collections import defaultdict
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
import traceback
from pathlib import Path
from src.utils import get_class_logger
from src.tools.ligand_protonation import protonate_ligand, summarize_record, LigandProtonationError

logger = get_class_logger(__name__, log_to_file=False)


class NoLigand(Exception):
    def __init__(self, message):
        super().__init__(message)


def fix_pdb_file(input_pdb: str, fixed_pdb: str) -> str:
    """
    This is a publically available API that fixes a PDB file using PDBFixer and writes the output to a new file. This
    function should only be called when the input PDB file is known to have  issues, otherwise it may introduce new problems.
    Actions applied:
      - Removes heterogens (ligands, ions), keeps waters
      - Finds missing residues (from the SEQRES records)
      - Finds missing atoms
      - Adds missing residues (including missing terminal residues) and missing heavy atoms (but no hydrogens)

    Args:
        input_pdb (str): The path where the faulty PDB file is located.
        fixed_pdb (str): The path where we save the fixed PDB file.

    """
    # Load a PDB file
    fixer = PDBFixer(filename=input_pdb)

    # Remove heterogens (ligands, ions) but keep water molecules
    fixer.removeHeterogens(keepWater=True)

    # Identify problems
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    if fixer.missingResidues:
        logger.warning(f"Missing residues detected, they are added by PDBFixer: {fixer.missingResidues}")

    # Apply fixes (adds missing residues and heavy atoms, no hydrogens)
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

def prepare_pdb_file_ligand(sandbox_dir: str, pdb_id: str, ph: float, ligand_name: str = None) -> str:
    """
    Takes input PDB file. Fixes it with PDBFixer, removes HETATM, CONECT and MASTER lines and saves the protein to {pdb_id}_prepared.pdb.
    The protein is protonated later, after capping, by protonate_protein.
    If there is one ligand, extract it from the original PDB file to ligand_name.pdb. If there are multiple ligands, extract each to a separate file named ligand_name_{index}.pdb.
    Protonates ligand at the given pH and saves to ligand_name_h.pdb (or ligand_name_{index}_h.pdb if there are multiple ligands).
    Args:
        sandbox_dir (str): Directory where the input PDB file is located and where output files will be saved.
        pdb_id (str): The 4-character PDB ID (e.g., '1abc').
        ph (float): pH at which the protonation state is chosen.
        ligand_name (str): Optional: the name of the ligand if a protein-ligand complex should be simulated.
    """
    # Fix PDB file if necessary
    fix_pdb_file(f"{sandbox_dir}/{pdb_id}.pdb", f"{sandbox_dir}/{pdb_id}_fixed.pdb")

    # Prepare PDB
    with (
        open(f"{sandbox_dir}/{pdb_id}_fixed.pdb", "r") as infile,
        open(f"{sandbox_dir}/{pdb_id}_prepared.pdb", "w") as outfile,
    ):
        for line in infile:
            if not (line.startswith("HETATM") or line.startswith("CONECT") or line.startswith("MASTER")):
                outfile.write(line)
    logger.info(f"Prepared PDB file saved to {sandbox_dir}/{pdb_id}_prepared.pdb")

    # Extract ligand from the original PDB file (fix_pdb_file removes all heterogens)
    if (ligand_name is not None) and (ligand_name != "XXX") and (ligand_name != "None") and (ligand_name != "None_h"):
        # Count number of ligands
        ligand_keys = set()
        with open(f"{sandbox_dir}/{pdb_id}.pdb") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand_name and line[16] in (" ", "A"):
                    chain = line[21]
                    resnum = int(line[22:26])
                    ligand_keys.add((chain, resnum))
        num_ligands = len(ligand_keys)
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
                    if line.startswith("HETATM") and line[17:20].strip() == ligand_name and line[16] in (" ", "A"):
                        outfile.write(line)
            logger.info(f"Extracted ligand {ligand_name} to {ligand_pdb_file}")
        
        # CASE 2: multiple ligands (split by residue number)
        else:
            ligand_pdb_files_list = []
            ligands = defaultdict(list)

            # Collect HETATM lines by residue number
            with open(f"{sandbox_dir}/{pdb_id}.pdb", "r") as infile:
                for line in infile:
                    if line.startswith("HETATM") and line[17:20].strip() == ligand_name and line[16] in (" ", "A"):
                        chain = line[21].strip() or "_"
                        resnum = int(line[22:26])
                        ligands[(chain, resnum)].append(line)

            # Write one file per ligand
            for i, ((chain, resnum), atom_lines) in enumerate(ligands.items(), start=1):
                ligand_pdb_file = f"{sandbox_dir}/{ligand_name}_{i}.pdb"
                ligand_pdb_files_list.append(ligand_pdb_file)
                with open(ligand_pdb_file, "w") as outfile:
                    outfile.writelines(atom_lines)
                logger.info(f"Extracted ligand {ligand_name} residue {resnum} to {ligand_pdb_file}")

    # Protonate ligand
    list_protonated_files = []
    records = []
    if (ligand_name is not None) and (ligand_name != "XXX") and (ligand_name != "None") and (ligand_name != "None_h"):
        for ligand_pdb_file in ligand_pdb_files_list: # loop over all extracted ligands
            if num_ligands == 1:
                protonated_prefix = f"{ligand_name}_h"
            else:
                index = Path(ligand_pdb_file).stem.split("_")[-1]  # get index from filename
                protonated_prefix = f"{ligand_name}_{index}_h"

            list_protonated_files.append(f"{protonated_prefix}.pdb")
            # protonate_ligand works with file names relative to sandbox_dir
            try:
                records.append(protonate_ligand(sandbox_dir, Path(ligand_pdb_file).name, protonated_prefix, ligand_name, ph=ph))
            except (LigandProtonationError, RuntimeError) as e:
                logger.error(f"Ligand protonation failed with error: {e}")
                return f"Ligand protonation failed with error: {e}"

        # Only the first ligand is simulated: report its protonation (candidates, agreement between predictors)
        protonation_report = summarize_record(records[0])

        if num_ligands == 1:
            message = f"Successfully Prepared PDB structure with a ligand and saved the extracted protein PDB file to {sandbox_dir}/{pdb_id}_prepared.pdb and the protonated ligand to {sandbox_dir}/{ligand_name}_h.sdf and {sandbox_dir}/{ligand_name}_h.pdb. Ligand was protonated at pH={ph} and atom names were cleaned (renumbered). The protonation decision is logged in {ligand_name}_h_protonation.json.\n{protonation_report}"
        else:
            message = f"Successfully Prepared PDB structure with {num_ligands} ligands and saved the extracted protein PDB file to {sandbox_dir}/{pdb_id}_prepared.pdb and the {num_ligands} protonated ligand PDB files to {sandbox_dir}/{list_protonated_files}. This list of {num_ligands} protonated files: {list_protonated_files} is IMPORTANT and should be the input parameter for future functions. Ligands were protonated at pH={ph} (.sdf and .pdb files) and atom names were cleaned (renumbered). Only the first ligand is simulated, its protonation decision is logged in {Path(list_protonated_files[0]).stem}_protonation.json.\n{protonation_report}"
        logger.info(message)
        return message


    return f"Successfully Prepared PDB structure without a ligand and saved the extracted PDB file to {sandbox_dir}/{pdb_id}_prepared.pdb. The protein is not protonated yet: run add_caps and then protonate_protein."

# Residue names that differ from the SEQRES name after preparation/protonation
SEQRES_EQUIVALENT_RESNAMES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
                              "CYX": "CYS", "CYM": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "MSE": "MET"}


def find_truncated_termini(pdb_file: str, chains: dict[str, list[str]], n_compare: int = 3) -> dict[str, tuple[bool, bool]]:
    """
    Finds which chain ends of a structure are truncated, by comparing the residues of each chain with the SEQRES
    records of the original PDB file. A chain end is a real terminus if its first (last) n_compare residues are the
    first (last) residues of the SEQRES sequence, e.g. when PDBFixer has rebuilt the missing terminal residues.
    Chains without SEQRES records (e.g. AlphaFold or modelled structures) are treated as truncated at both ends.

    Args:
        pdb_file (str): Path of the original PDB file (with SEQRES records).
        chains (dict): {chain_id: residue names of the chain in the structure to cap}.
        n_compare (int): Number of terminal residues compared with the SEQRES sequence.
    Returns:
        {chain_id: (n_terminus_truncated, c_terminus_truncated)}
    """
    seqres = {sequence.chainId: sequence.residues for sequence in PDBFixer(filename=pdb_file).sequences}

    truncated = {}
    for chain_id, resnames in chains.items():
        if chain_id not in seqres:
            logger.warning(f"Chain {chain_id}: no SEQRES record in {pdb_file}, the chain ends are treated as truncated and capped with ACE/NME")
            truncated[chain_id] = (True, True)
            continue
        observed = [SEQRES_EQUIVALENT_RESNAMES.get(name, name) for name in resnames]
        reference = [SEQRES_EQUIVALENT_RESNAMES.get(name, name) for name in seqres[chain_id]]
        n = min(n_compare, len(observed), len(reference))
        n_truncated = observed[:n] != reference[:n]
        c_truncated = observed[-n:] != reference[-n:]
        if n_truncated:
            logger.info(f"Chain {chain_id}: N-terminus {'-'.join(resnames[:n])} differs from the SEQRES start {'-'.join(seqres[chain_id][:n])}, it will be capped with ACE")
        if c_truncated:
            logger.info(f"Chain {chain_id}: C-terminus {'-'.join(resnames[-n:])} differs from the SEQRES end {'-'.join(seqres[chain_id][-n:])}, it will be capped with NME")
        truncated[chain_id] = (n_truncated, c_truncated)

    return truncated


def insert_ter_after_chain_ends(pdb_file: str) -> int:
    """
    Makes sure there is one TER record after the end of each chain, so that tleap treats each chain end as a
    terminus (charged or capped) and does not bond consecutive chains.
    A chain end is detected when the chain ID changes, after an NME cap, or after a residue containing an OXT atom.
    Existing TER records are kept, duplicated ones are removed. The file is overwritten.

    Args:
        pdb_file (str): Path of the PDB file to edit.
    Returns:
        Number of TER records in the written file.
    """
    with open(pdb_file, "r") as f:
        lines = f.readlines()

    output_lines = []
    prev_chain = None
    prev_residue = None
    prev_is_chain_end = False  # previous residue is an NME cap or contains an OXT atom
    has_end_record = False

    def add_ter():
        if output_lines and not output_lines[-1].startswith("TER"):
            output_lines.append("TER\n")

    for line in lines:
        record = line[:6].strip()

        if record in ("ATOM", "HETATM"):
            chain_id = line[21]
            residue = line[17:27]  # residue name, chain ID, residue number, insertion code

            # New residue: insert TER if the previous residue ended a chain
            if prev_residue is not None and residue != prev_residue:
                if chain_id != prev_chain or prev_is_chain_end:
                    add_ter()
                prev_is_chain_end = False

            if line[17:20].strip() == "NME" or line[12:16].strip() == "OXT":
                prev_is_chain_end = True

            output_lines.append(line)
            prev_chain, prev_residue = chain_id, residue

        elif record == "TER":
            add_ter()

        elif record in ("END", "ENDMDL"):
            add_ter()
            output_lines.append(line)
            has_end_record = True

        else:
            output_lines.append(line)

    if not has_end_record:
        add_ter()

    with open(pdb_file, "w") as f:
        f.writelines(output_lines)

    n_ter = sum(1 for line in output_lines if line.startswith("TER"))
    logger.info(f"{n_ter} TER records written in {pdb_file}")
    return n_ter


def add_crystal_waters(water_source_pdb: str, pdb_file: str) -> int:
    """
    Appends the crystal waters (HOH oxygens) of water_source_pdb to the end of pdb_file, after the protein chains and
    followed by a TER record. Waters are renumbered after the last residue (and atom) of pdb_file. Hydrogens are not
    needed: tleap adds them from the TIP3P template. The file is overwritten.

    Args:
        water_source_pdb (str): PDB file containing the crystal waters (e.g. {pdb_id}_fixed.pdb, alternate locations resolved by PDBFixer).
        pdb_file (str): PDB file to which the waters are added (e.g. {pdb_id}_pdb2pqr.pdb).
    Returns:
        Number of waters added.
    """
    with open(water_source_pdb, "r") as f:
        water_lines = [line for line in f if line.startswith(("ATOM", "HETATM")) and line[17:20] == "HOH" and line[12:16].strip() == "O"]
    if not water_lines:
        return 0

    with open(pdb_file, "r") as f:
        lines = f.readlines()

    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    serial = int(atom_lines[-1][6:11]) if atom_lines else 0
    resid = int(atom_lines[-1][22:26]) if atom_lines else 0

    new_lines = []
    for line in water_lines:
        serial += 1
        resid += 1
        new_lines.append(f"HETATM{serial % 100000:5d}{line[11:22]}{resid % 10000:4d}{line[26:].rstrip()}\n")
    new_lines.append("TER\n")

    # Insert the waters before END (or at the end of the file if there is no END record)
    end_idx = next((i for i, line in enumerate(lines) if line.startswith(("END", "ENDMDL"))), len(lines))
    if end_idx > 0 and not lines[end_idx - 1].endswith("\n"):
        lines[end_idx - 1] += "\n"
    lines = lines[:end_idx] + new_lines + lines[end_idx:]

    with open(pdb_file, "w") as f:
        f.writelines(lines)

    logger.info(f"{len(water_lines)} crystal waters from {water_source_pdb} added to {pdb_file}")
    return len(water_lines)


def add_caps(sandbox_dir: str, input_pdb: str, pdb_id: str) -> str:
    """
    Adds ACE and NME caps to the N- and C-termini of the protein chains in the input PDB file, only where the chain
    end is truncated (chain end differs from the SEQRES records of the original {pdb_id}.pdb).
    Real chain ends are left uncapped so that tleap treats them as charged termini.
    Chains without SEQRES records (e.g. AlphaFold structures) are capped at both ends.
    Saves the modified structure to {sandbox_dir}/{pdb_id}_prepared_capped.pdb (also written if no cap is added).

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

    # Find which chain ends are truncated by comparing the chains with the SEQRES records of the original PDB file
    chain_resnames = {}
    for seg in u.segments:
        chain = u.select_atoms(f"segid {seg.segid}")
        chain_id = chain.chainIDs[0] if hasattr(chain, "chainIDs") else seg.segid
        chain_resnames[chain_id] = list(chain.residues.resnames)
    truncated = find_truncated_termini(f"{sandbox_dir}/{pdb_id}.pdb", chain_resnames)

    # Access each fragment separately
    res_start = 0
    segment_universes = []
    report = []

    for seg in u.segments:
        chain = u.select_atoms(f"segid {seg.segid}")
        chain_id = chain.chainIDs[0] if hasattr(chain, "chainIDs") else seg.segid
        cap_n, cap_c = truncated[chain_id]
        report.append(f"chain {chain_id}: N-terminus {'capped with ACE' if cap_n else 'left uncapped'}, C-terminus {'capped with NME' if cap_c else 'left uncapped'}")

        # Add ACE
        if cap_n:
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

        # Add NME (methyl carbon named C as in ff14SB)
        resid_c = chain.residues.resids[-1]
        end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")
        if cap_c:
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

        ## Merge Universe: the OXT of a capped C-terminus is replaced by the NME nitrogen
        if cap_c and "OXT" in end_residue.names:
            index = np.where(end_residue.names == "OXT")[0][0]
            OXT = end_residue[index]

            Chain = u.select_atoms(f"segid {seg.segid} and not index {OXT.index}")

        else:
            Chain = u.select_atoms(f"segid {seg.segid}")

        ### Merge ACE, Protein and NME
        parts = ([ace_universe.atoms] if cap_n else []) + [Chain] + ([nme_universe.atoms] if cap_c else [])
        u_all = mda.Merge(*parts)

        # to renumber residues
        first_resid = res_start + 1
        resids_ace = [first_resid] * 3 if cap_n else []
        first_pro = first_resid + 1 if cap_n else first_resid
        resids_pro = np.arange(first_pro, first_pro + Chain.residues.n_residues)
        resids_nme = [resids_pro[-1] + 1] * 2 if cap_c else []

        u_all.atoms.residues.resids = np.concatenate(
            [resids_ace, resids_pro, resids_nme]
        )
        res_start = u_all.atoms.residues.resids[-1]
        segment_universes.append(u_all)

    ## Join all the universes
    all_uni = mda.Merge(*(seg.atoms for seg in segment_universes))
    all_uni.atoms.write(f"{sandbox_dir}/{pdb_id}_prepared_capped.pdb")

    insert_ter_after_chain_ends(f"{sandbox_dir}/{pdb_id}_prepared_capped.pdb")

    logger.info("Capping: " + "; ".join(report))
    return f"Capping done ({'; '.join(report)}). Only truncated chain ends (or chains without SEQRES records) were capped with ACE/NME, real chain ends were left uncapped so that they are charged by tleap. Saved to {sandbox_dir}/{pdb_id}_prepared_capped.pdb"


def protonate_protein(sandbox_dir: str, input_pdb: str, pdb_id: str, ph: float) -> str:
    """
    Protonates the protein at the given pH with PDB2PQR (pKa prediction with PROPKA, AMBER residue and atom names,
    e.g. HID/HIE/HIP for histidines). Should be run after add_caps.
    Writes {input_pdb stem}.pqr and {pdb_id}_pdb2pqr.pdb, the input of run_tleap_protein.

    Args:
        sandbox_dir (str): Directory where the input PDB file is located and where output files will be saved.
        input_pdb (str): The capped protein PDB file, relative to sandbox_dir (e.g. {pdb_id}_prepared_capped.pdb).
        pdb_id (str): The PDB ID.
        ph (float): pH at which the protonation states are chosen.
    """
    pqr_file = f"{Path(input_pdb).stem}.pqr"
    pdb_output = f"{pdb_id}_pdb2pqr.pdb"
    cmd = f"pdb2pqr30 -ff=AMBER --titration-state-method=propka --with-ph={ph} --ffout=AMBER --keep-chain --pdb-output={pdb_output} {input_pdb} {pqr_file}"
    run_1 = subprocess.run(shlex.split(cmd), cwd=sandbox_dir, capture_output=True, text=True)
    if run_1.returncode != 0 or not Path(f"{sandbox_dir}/{pdb_output}").exists():
        error_text = "\n".join(filter(None, [run_1.stderr, run_1.stdout]))
        logger.error(f"PDB2PQR failed with error: {error_text}")
        return f"PDB2PQR failed with error: {error_text}"

    logger.info(f"PDB2PQR output saved to {sandbox_dir}/{pqr_file} and {sandbox_dir}/{pdb_output}")

    # Make sure each chain end is followed by a TER record for tleap
    n_ter = insert_ter_after_chain_ends(f"{sandbox_dir}/{pdb_output}")

    # Add back the crystal waters (removed by prepare_pdb_file_ligand), they are not protonated by PDB2PQR: tleap builds their hydrogens
    n_waters = 0
    fixed_pdb = f"{sandbox_dir}/{pdb_id}_fixed.pdb"
    if Path(fixed_pdb).exists():
        n_waters = add_crystal_waters(fixed_pdb, f"{sandbox_dir}/{pdb_output}")
    else:
        logger.warning(f"{fixed_pdb} not found: no crystal waters added to {pdb_output}")

    return f"Protein protonated at pH={ph} with PDB2PQR (PROPKA). Saved to {sandbox_dir}/{pdb_output} (input for run_tleap_protein, {n_ter} TER records, one per chain end) and {sandbox_dir}/{pqr_file}. {n_waters} crystal waters from {pdb_id}_fixed.pdb were added at the end of {pdb_output}."


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
