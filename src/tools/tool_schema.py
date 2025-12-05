from pydantic import BaseModel
from typing import Any, Dict


class Tool(BaseModel):
    name: str
    description: str
    parameters: Dict[str, Any]


def create_tool_schema_prep(sandbox_dir):
    tools = [
        Tool(
            name="find_input",
            description="Find the uploaded file from the user. This always searches the sandbox directory automatically.",
            parameters={
                "type": "object",
                "properties": {},
                "required": [],
            },
        ),
        Tool(
            name="read_file",
            description="Read the contents of a file at the specified path",
            parameters={
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "The path to the file to read (this is the sandbox directory located at {sandbox_dir})",
                    },
                },
                "required": ["path"],
            },
        ),
        Tool(
            name="fetch_and_save_pdb",
            description="Fetch a PDB file from the RCSB server using the PDB ID and save it locally in the sandbox directory",
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": f"The directory path where input files are located and output files are produced (this is the sandbox directory): {sandbox_dir}",
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": "The 4-character PDB ID (e.g., '1abc')",
                    },
                    "output_pdb": {
                        "type": "string",
                        "description": "The path to save the fetched PDB file. This should be named with all capital letters from the four letter pdb_id",
                    },
                },
                "required": ["sandbox_dir", "pdb_id", "output_pdb"],
            },
        ),
        Tool(
            name="search_papers",
            description=(
                """Useful to answer system specific questions that may be found in literature. Ask a specific question as the input."""
            ),
            parameters={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": ("What to query the papers for."),
                    },
                },
                "required": ["query"],
            },
        ),
    ]

    return [
        {
            "type": "function",
            "function": {
                "name": tool.name,
                "description": tool.description,
                "parameters": tool.parameters,
            },
        }
        for tool in tools
    ]


def create_tool_schema_md(sandbox_dir, ligand_name, pdb_id):
    tools = [
        Tool(
            name="prepare_pdb_file_ligand",
            description=(
                "For protein-ligand complexes and protein-only comlexes. This tool will prepare the PDB files:"
                "One pdb file containing the protein only, without solvent and ligand."
                "If a ligand is present, the pdb file containing the ligand only, without solvent and protein."
                f"The ligand file is named {ligand_name}.pdb and is then protonated and saved to {ligand_name}_h.pdb. You should check the protonation of the ligand in the literature to make sure it has been properly protonated at the appropriate pH."
                f"The output protein PDB file is named {pdb_id}_prepared.pdb."
                f"The output ligand PDB file is named {ligand_name}_h.pdb."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": f"The directory path where input files are located and output files are produced (this is the sandbox directory): {sandbox_dir}",
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (
                            f"The PDB ID of the structure to be prepared: {pdb_id}. The PDB file {pdb_id}.pdb must exist in the sandbox_dir."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three-letter residue name of the ligand to extract from the PDB file, in capital letters, called {ligand_name}."
                        ),
                    },
                },
                "required": ["path", "pdb_id", "ligand_name"],
            },
        ),
        Tool(
            name="add_caps",
            description=(
                f"Add ACE and NME caps to the N- and C-termini of a protein PDB file. "
                "This script is used if the pdb file was imported from the PDB database and is missing caps. "
                "This script should be run after the PDB file has been prepared from prepare_pdb_file_ligand."
                f"The input PDB file should be named {pdb_id}_prepared.pdb and the output file will be named {pdb_id}_prepared_capped.pdb."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": f"The directory path where input files are located and output files are produced (this is the sandbox directory): {sandbox_dir}",
                    },
                    "input_pdb": {
                        "type": "string",
                        "description": (
                            f"Input pdb file to be capped (without the path). This file must exist in the sandbox_dir ({sandbox_dir})."
                            f" The input file should be named {pdb_id}_prepared.pdb if it has previously been prepared with the prepare_pdb_file_ligand tool."
                        ),
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (
                            f"The PDB ID of the structure to be prepared: {pdb_id}. The PDB file {pdb_id}_prepared.pdb must exist in the sandbox_dir ({sandbox_dir})."
                        ),
                    },
                },
                "required": ["sandbox_dir", "input_pdb", "pdb_id"],
            },
        ),
        Tool(
            name="rename_histidines",
            description=(
                "Renames histidines HIS to account for their correct protonation."
                "This tool should be ran after the protein has been capped with the add_caps function for both protein-only and protein-ligand systems."
                f"The input PDB file should be named {pdb_id}_prepared_capped.pdb and the output file will be named {pdb_id}_prepared_capped_his.pdb."
                "HIS will be renamed to HID if it contains HD1 and not HE2."
                "HIS will be renamed to HIE if it contains HE2 and not HD1."
                "HIS will be renamed to HIP if it contains both HD1 and HE2."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": f"The directory path where input files are located and output files are produced (this is the sandbox directory): {sandbox_dir}",
                    },
                    "input_pdb": {
                        "type": "string",
                        "description": (
                            f"Input pdb file for which HIS will be modified if needed (without the path). This file must exist in the sandbox_dir ({sandbox_dir})."
                            f" The input file should be named {pdb_id}_prepared_capped.pdb or {pdb_id}_prepared.pdb if it has previously been prepared with the prepare_pdb_file_ligand tool and potentially capped with the add_caps tool."
                        ),
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (
                            f"The PDB ID of the structure to be prepared: {pdb_id}. The PDB file {pdb_id}_prepared.pdb must exist in the sandbox_dir ({sandbox_dir})."
                        ),
                    },
                },
                "required": ["sandbox_dir", "input_pdb", "pdb_id"],
            },
        ),
        Tool(
            name="run_tleap",
            description=(
                "Prepare the molecular system using Amber’s tleap utility. "
                "The script should ONLY be used for protein-only systems. Ligands are NOT included and should be processed with the run_tleap_ligand tool. "
                "This script should be run after the PDB file has been prepared and capped if necessary. "
                "This is the FIRST step in the molecular dynamics (MD) pipeline when the force field ff14sb is chosen and the system is a protein alone (without a ligand). "
                "The tleap process parameterizes the biomolecule using the Amber force field ff14sb, "
                f"generates Amber topology ({pdb_id}.prmtop) and coordinate ({pdb_id}.inpcrd) files. "
                f"ParmEd is then used to convert these files to GROMACS format (topol.top, {pdb_id}.gro). "
                "It executes the shell script 'run_tleap.sh' in the working directory specified by sandbox_dir ({sandbox_dir}). "
                "All input and output files are read from and written to this directory. "
                "If tleap fails, inspect the PDB for missing atoms or nonstandard residues."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            "Absolute path to the working directory where tleap should be executed. All output files (.prmtop, .inpcrd) will be saved here."
                        ),
                    },
                    "input_pdb": {
                        "type": "string",
                        "description": (
                            "Input PDB structure file to be processed by tleap. "
                            "This file must exist in the provided sandbox_dir."
                            f"If the PDB file was previously prepared and capped, input_pdb will be {pdb_id}_prepared_capped.pdb."
                        ),
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (f"The PDB ID of the structure to be prepared: {pdb_id}. "),
                    },
                },
                "required": ["sandbox_dir", "input_pdb", "pdb_id"],
            },
        ),
        Tool(
            name="param_ligand",
            description=(
                "Parameterize a small molecule ligand using Amber's antechamber and parmchk2 utilities. "
                "This tool generates the necessary files for including a ligand in molecular dynamics simulations. "
                "It creates a mol2 file with assigned charges, a prepi file for the ligand, and a frcmod file containing any missing force field parameters. "
                f"The input ligand file should be in PDB format, protonated ({ligand_name}_h.pdb) and located in the specified sandbox_dir ({sandbox_dir}). "
                "If the ligand's net charge is not provided, it will be calculated from the structure. "
                f"The output files ({ligand_name}.mol2, {ligand_name}.prepi, {ligand_name}.frcmod) will be saved in the same sandbox_dir ({sandbox_dir})."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": f"The directory path where input files are located and output files are produced (this is the sandbox directory): {sandbox_dir}",
                    },
                    "ligand_file": {
                        "type": "string",
                        "description": (
                            f"Input PDB file of the protonated ligand to be parameterized (without full path). "
                            f"The ligand file is named {ligand_name}_h.pdb if previously protonated with the prepare_pdb_file_ligand tool."
                            "This file must exist in the provided sandbox_dir ({sandbox_dir})."
                        ),
                    },
                    # "charge_ligand": {
                    #     "type": ["integer", "null"],
                    #     "description": (
                    #         "The net charge of the ligand. If not provided, the charge will be calculated from the structure."
                    #     ),
                    # },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three character residue name of the ligand to extract from the PDB file, in capital letters, called {ligand_name}."
                        ),
                    },
                },
                "required": ["sandbox_dir", "ligand_file", "ligand_name"],
            },
        ),
        Tool(
            name="run_tleap_ligand",
            description=(
                "Prepare the molecular system using Amber’s tleap utility."
                "This script should ONLY be used for a protein-ligand complex, so if the user specifies that ligand should be included in the simulation. "
                "This is the FIRST step in the molecular dynamics (MD) pipeline when the force field ff14sb is chosen. "
                "First a file complex.gro is created by combining the protein PDB file and the ligand PDB file. "
                f"To process the ligand, previously the ligand must have been parameterized with the param_ligand tool, generating the files {ligand_name}.mol2, {ligand_name}.prepi and {ligand_name}.frcmod. "
                "The tleap process parameterizes the protein-ligand complex using the Amber force field ff14sb, "
                f"generates Amber topology ({pdb_id}.prmtop) and coordinate ({pdb_id}.inpcrd) files. "
                f"ParmEd is then used to convert these files to GROMACS format (topol.top, {pdb_id}.gro). "
                "It executes the shell script 'run_tleap_ligand.sh' in the working directory specified by 'sandbox_dir' ({sandbox_dir}). "
                "All input and output files are read from and written to this directory. "
                "If tleap fails, inspect the PDB for missing atoms or nonstandard residues."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "input_pdb": {
                        "type": "string",
                        "description": (
                            f"Input PDB of the protein (without ligand) to be merged with ligand file (don't provide full path). "
                            f"If the protein has been prepared and capped, it will be called {pdb_id}_prepared_capped.pdb. "
                            "This file must exist in the provided sandbox_dir."
                        ),
                    },
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            "Absolute path to the working directory where tleap should be executed. "
                            "All output files (.prmtop, .inpcrd) will be saved here."
                        ),
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (f"The PDB ID of the structure to be prepared: {pdb_id}. "),
                    },
                    "ligand_file": {
                        "type": "string",
                        "description": (
                            f"Input PDB file of the ligand to be merged with the protein (without full path). "
                            f"The ligand file is named {ligand_name}_h.pdb if previously protonated with the prepare_pdb_file_ligand tool."
                            "This file must exist in the provided sandbox_dir."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three-letter residue name of the ligand to extract from the PDB file, in capital letters, called {ligand_name}."
                        ),
                    },
                },
                "required": [
                    "sandbox_dir",
                    "input_pdb",
                    "pdb_id",
                    "ligand_file",
                    "ligand_name",
                ],
            },
        ),
        Tool(
            name="gromacs_equil",
            description=(
                """
                Run GROMACS equilibration and optionally the production MD phase.
                This tool should be executed AFTER solvation and ion addition are complete.

                It executes 'equil_Gromacs.sh' inside the working directory provided in
                'sandbox_dir'. The script performs energy minimization, NVT and NPT
                equilibration, and possibly production MD.

                Parameter (.mdp) files can be found in the sandbox_dir. All intermediate and output files —
                including md.tpr, md.xtc, md.edr, and md.log — are generated in the same
                'sandbox_dir'.
                """
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            "Absolute path to the working directory where GROMACS will perform equilibration and production MD. "
                            "Outputs are written here."
                        ),
                    },
                    "input_gro": {
                        "type": "string",
                        "description": (
                            "File of the structure to be equilibrated. "
                            "If a protein-ligand complex is simulated, the file will be complex.gro. "
                            f"If a protein alone is simulated, the file will be {pdb_id}.gro"
                            "This file must be located within sandbox_dir."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            "Three-character residue name of the ligand (capital letters or numbers). If no ligand was provided by user, do not input this argument."
                        ),
                    },
                    "ligand_file": {
                        "type": "string",
                        "description": (
                            "Optional argument for the input PDB file of the ligand to be merged with the protein (without full path). "
                            "The ligand file is typically named '<ligand_name>_h.pdb' if previously protonated. "
                            "If no ligand was provided by user, do not input this argument."
                        ),
                    },
                },
                "required": ["sandbox_dir", "input_gro"],
            },
        ),
        Tool(
            name="gromacs_production",
            description=(
                """
                    Run GROMACS production MD phase.
                    This tool should be executed AFTER equilibration with the gromacs_equil tool is COMPLETE.
                    The script performs production MD using the provided checkpoint file from NPT equilibration, npt.cpt.
                    All intermediate and output files — including md.tpr, md.xtc, md.edr, and md.log are generated in the same 'sandbox_dir'.
                    If the user requests another production run length than the default 0.1 ns, the md.mdp file must be edited accordingly before running this tool.
                    The output trajectory file is md.xtc.
                    """
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            "Absolute path to the working directory where GROMACS will perform production MD. Outputs are written here."
                        ),
                    },
                    "input_gro": {
                        "type": "string",
                        "description": (
                            "Equilibrated structure file (e.g., npt.gro) to be used as input for production MD (without full file path). "
                            "Must be located within sandbox_dir."
                        ),
                    },
                    "npt_cpt_file": {
                        "type": "string",
                        "description": (
                            "Checkpoint file (e.g., npt.cpt) from NPT equilibration (without full file path). "
                            "Must be located within sandbox_dir."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            "Three-character residue name of the ligand (capital letters or numbers). If no ligand was provided by user, do not input this argument."
                        ),
                    }
                },
                "required": ["sandbox_dir", "input_gro", "npt_cpt_file"],
            },
        ),
        Tool(
            name="gromacs_analysis",
            description=(
                """
                    Analyse GROMACS production MD phase.
                    This tool should be executed AFTER production with the gromacs_production tool is COMPLETE.
                    The input is the production trajectory md.xtc.
                    The output trajectory file is md.xtc.
                    PBC is removed to generate md_noPBC.xtc which is used for analysis. 
                    RMSD analysis is performed to create rmsd.xvg (with respect to the initial structure of the trajectory) and rmsd_xtal.xvg (with respect to the crystal structure). 
                    RMSF analysis is performed to create rmsf.xvg.
                    Radius of gyration analysis is performed to create gyrate.xvg.
                    Hydrogen bond analysis is performed to create hbnum_mainchain.xvg (number of hydrogen bonds in the protein backbone), hbnum_sidechain.xvg (number of hydrogen bonds in the protein side chains), and hbnum_prot_wat.xvg (total number of hydrogen bonds between the protein and water molecules).
                    If a ligand is present, hbnum_prot_lig.xvg (number of hydrogen bonds between the protein and the ligand) is also created.
                    """
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            "Absolute path to the working directory where GROMACS will perform production MD. Outputs are written here."
                        ),
                    },
                    "input_xtc": {
                        "type": "string",
                        "description": (
                            "Trajectory file from the production run. It should be called md.xtc (without full file path). "
                            "Must be located within sandbox_dir."
                        ),
                    },
                    "ligand_name": {
                        "type": ["string", "null"],
                        "description": (
                            f"The three-letter residue name of the ligand in capital letters {ligand_name}. "
                            "If no ligand is present, this can be left null."
                        ),
                    },
                },
                "required": ["sandbox_dir", "input_xtc"],
            },
        ),
        Tool(
            name="run_gmxMMPBSA",
            description=(
                """
                    Run MMPBSA analysis on the GROMACS production MD trajectory using gmx_MMPBSA for protein-ligand systems only.
                    This tool should be executed AFTER production with the gromacs_production tool is COMPLETE, and the gromacs_analysis tool successfully created the md_noPBC.xtc trajectory file.
                    nsteps and nstxout values can be found in the md.mdp file used for the production run, located in sandbox_dir.
                    The temperature used during the MD simulation is also required for the MMPBSA calculation.
                    A new directory called gmx_MMPBSA is created and all MMPBSA output files are saved there, including the final binding energy summary file called FINAL_RESULTS_MMPBSA.dat.
                    The final binding energy can be found in FINAL_RESULTS_MMPBSA.dat at the last ΔTOTAL line.
                    """
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            "Absolute path to the working directory where GROMACS will perform production MD. Outputs are written here."
                        ),
                    },
                    "pdb_id": {
                        "type": ["string"],
                        "description": (f"The PDB ID of the structure to be prepared: {pdb_id}. "),
                    },
                    "nsteps": {
                        "type": ["string"],
                        "description": ("Number of MD steps performed during the production run, found in the md.mdp file located in sandbox_dir."),
                    },
                    "nstxout_compressed": {
                        "type": ["string"],
                        "description": ("Number of MD steps that elapse between writing position coordinates using lossy compression, found in the md.mdp file located in sandbox_dir."),
                    },
                    "temp": {
                        "type": ["string"],
                        "description": ("Temperature used during the MD simulation, found in the md.mdp file located in sandbox_dir. This value is an integer with base 10."),
                    },
                },
                "required": ["sandbox_dir", "pdb_id", "nsteps", "nstxout_compressed", "temp"],
            },
        ),
        Tool(
            name="find_input",
            description="Find the uploaded file from the user. This always searches the sandbox directory automatically.",
            parameters={
                "type": "object",
                "properties": {},
                "required": [],
            },
        ),
        Tool(
            name="read_file",
            description=(
                f"""Read and return the contents of a chosen text file from the {sandbox_dir}. Useful for inspecting
                     logs, configuration files, or generated input/output data. Use this ONLY for diagnostic 
                     purposes. This does NOT run simulations or modify files."""
            ),
            parameters={
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "Absolute path to the file to be read. Must exist inside the sandbox directory."
                        ),
                    },
                },
                "required": ["path"],
            },
        ),
        Tool(
            name="list_files",
            description=(
                """List all files and subdirectories in the sandbox directory. 
                    Intended for verifying the presence of expected outputs (e.g., topol.top, md.xtc) 
                    after each simulation stage. Use this ONLY for diagnostic purposes. This does NOT run 
                    simulations or modify files."""
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": ("Absolute path to the directory to list."),
                    },
                },
                "required": ["sandbox_dir"],
            },
        ),
        Tool(
            name="edit_file",
            description=(
                """Edit or create a text file inside the working directory. 
                    Use this tool to modify .mdp files, correct input scripts, or patch configuration files. 
                    The 'sandbox_dir' argument should be an absolute file path within sandbox_dir."""
            ),
            parameters={
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": ("Absolute path to the file to edit or create. Must be within sandbox_dir."),
                    },
                    "old_text": {
                        "type": "string",
                        "description": ("Exact text to replace. Leave blank if creating a new file."),
                    },
                    "new_text": {
                        "type": "string",
                        "description": ("New content or replacement text to write to the file."),
                    },
                },
                "required": ["path", "old_text", "new_text"],
            },
        ),
        Tool(
            name="search_papers",
            description=(
                """Useful to answer system specific questions that may be found in literature. Ask a specific question as the input.
                Helpful things to search for in papers include:
                1. protonation of the ligand
                2. protonation of the active site residues and check that the protonation is correct with our current protein
                3. if the protein works better at a specific pH (influences protonation of the protein)
                4. if the protein works better at a specific temperature
                5. if the protein was crystallised under specific conditions that would required some modifications to the pdb file

                For example, you could ask: What protonation states are typical for protein simulations at pH 7?
                """
            ),
            parameters={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": ("What to query the papers for."),
                    },
                },
                "required": ["query"],
            },
        ),
    ]

    return [
        {
            "type": "function",
            "function": {
                "name": tool.name,
                "description": tool.description,
                "parameters": tool.parameters,
            },
        }
        for tool in tools
    ]
