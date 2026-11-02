from pydantic import BaseModel
from typing import Any, Dict


class Tool(BaseModel):
    name: str
    description: str
    parameters: Dict[str, Any]


def create_tool_schema_prep(sandbox_dir):
    tools = [
        Tool(
            name="grep_file",
            description=(
                "Search a file for lines matching a regex pattern. Returns matching lines with their line numbers. "
                "Useful for inspecting specific parts of large files (e.g. finding HETATM lines in a PDB, checking specific residues, "
                "finding errors in log files) without reading the entire file."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": f"Absolute path to the file to search. Must be inside the sandbox directory ({sandbox_dir}).",
                    },
                    "pattern": {
                        "type": "string",
                        "description": "A regex pattern to match against each line of the file (e.g. '^HETATM', 'ERROR|FATAL', 'ILE A   2').",
                    },
                },
                "required": ["path", "pattern"],
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


def create_tool_schema_md(sandbox_dir, ligand_name, pdb_id, ph):
    tools = [
        Tool(
            name="prepare_pdb_file_ligand",
            description=(
                "For protein-ligand complexes and protein-only comlexes. This tool will prepare the PDB files:"
                "One pdb file containing the protein only, without solvent and ligand. The protein is NOT protonated by this tool: it is protonated later with protonate_protein, after add_caps. "
                "If one or more ligands are present, one pdb file will be created per ligand, without solvent and protein."
                f"If there is one ligand, the ligand file is named {ligand_name}.pdb and is then protonated and saved to {ligand_name}_h.sdf and {ligand_name}_h.pdb. "
                f"If there are multiple ligands (i ligands), the ligand files are named {ligand_name}_i.pdb and are then protonated and saved to {ligand_name}_i_h.sdf and {ligand_name}_i_h.pdb. "
                "Ligand protonation: bond orders are taken from the PDB Chemical Component Dictionary (or perceived by Open Babel), the state at the pH predicted by Dimorphite-DL is applied "
                "and compared with the state predicted by pkasolver. All candidate states, the agreement between the two predictors and the selected state are logged in "
                f"{ligand_name}_h_protonation.json and reported in the output of this tool. If the output reports a disagreement or an ambiguous state, look for literature evidence "
                "and, if needed, change the state with set_ligand_protonation. Never edit the ligand files directly. "
                f"The output protein PDB file is named {pdb_id}_prepared.pdb."
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
                    "ph": {
                        "type": "number",
                        "description": (
                            f"The pH at which the protonation state of the ligand is chosen: {ph}."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three-letter residue name of the ligand to extract from the PDB file, in capital letters, called {ligand_name}."
                        ),
                    },
                },
                "required": ["path", "pdb_id", "ph", "ligand_name"],
            },
        ),
        Tool(
            name="add_caps",
            description=(
                "Add ACE and NME caps to the N- and C-termini of the protein chains, ONLY where the chain end is truncated "
                f"(the first or last residues of the chain differ from the SEQRES records of {pdb_id}.pdb; terminal residues rebuilt by PDBFixer count as real ends). "
                "Real chain ends are left uncapped so that tleap treats them as charged termini. Chains without SEQRES records (e.g. AlphaFold structures) are capped at both ends. "
                "This tool should be run after prepare_pdb_file_ligand and before protonate_protein, for both protein-only and protein-ligand systems. "
                f"The input PDB file should be named {pdb_id}_prepared.pdb and the output file will be named {pdb_id}_prepared_capped.pdb (written even if no cap is added)."
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
            name="protonate_protein",
            description=(
                "Protonate the protein at the given pH with PDB2PQR (pKa prediction with PROPKA, AMBER residue and atom names). "
                "Histidines are renamed to HID, HIE or HIP according to their protonation state, so rename_histidines is not needed. "
                "This tool must be run after add_caps and before run_tleap_protein, for both protein-only and protein-ligand systems. "
                f"The input PDB file should be named {pdb_id}_prepared_capped.pdb. The output files are {pdb_id}_pdb2pqr.pdb (the input of run_tleap_protein) "
                f"and {pdb_id}_prepared_capped.pqr. The crystal waters of {pdb_id}_fixed.pdb are added back at the end of {pdb_id}_pdb2pqr.pdb (oxygens only, tleap adds the hydrogens). "
                f"You should check in the literature that the protonation states of key residues are appropriate at pH {ph}."
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
                            f"Capped protein PDB file to be protonated (without the path), named {pdb_id}_prepared_capped.pdb. "
                            f"This file must exist in the sandbox_dir ({sandbox_dir})."
                        ),
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (f"The PDB ID of the structure to be prepared: {pdb_id}."),
                    },
                    "ph": {
                        "type": "number",
                        "description": (f"The pH at which the protonation states of the protein are chosen: {ph}."),
                    },
                },
                "required": ["sandbox_dir", "input_pdb", "pdb_id", "ph"],
            },
        ),
        Tool(
            name="rename_histidines",
            description=(
                "Renames histidines HIS to account for their correct protonation."
                "This tool is NOT part of the standard workflow: protonate_protein already names histidines HID, HIE or HIP. Only use it if a PDB file still contains HIS residues with hydrogens. "
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
            name="set_ligand_protonation",
            description=(
                "Change the protonation state (charge state or tautomer) of the ligand. This is the ONLY way to change the ligand protonation: "
                "the ligand files cannot be edited directly. Use it when prepare_pdb_file_ligand reports a disagreement between Dimorphite-DL and pkasolver "
                "or an ambiguous state and you have decided which state is correct, or when the literature supports another state. "
                f"Choose one of the candidate states listed in {ligand_name}_h_protonation.json (candidate_id), or give a SMILES of the full ligand with all its "
                "hydrogens and charges implicit in the SMILES (a literature reference is then required if the state is not one of the candidates). "
                "The state must have the same heavy atoms and connectivity as the crystal pose; only hydrogens, charges and bond orders can change. "
                "The tool rebuilds the ligand .sdf/.pdb, logs your decision and rationale, and moves files built from the previous state to a stale_protonation_* folder. "
                "Afterwards param_ligand, run_tleap_ligand, run_tleap_protein and the GROMACS steps must be run again."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "ligand_file": {
                        "type": "string",
                        "description": (
                            f"The protonated ligand file (without path): {ligand_name}_h.pdb, or {ligand_name}_1_h.pdb if there are multiple ligands."
                        ),
                    },
                    "candidate_id": {
                        "type": "string",
                        "description": "Id of one of the candidate states (e.g. 'C2'), as listed in the output of prepare_pdb_file_ligand. Give either candidate_id or smiles.",
                    },
                    "smiles": {
                        "type": "string",
                        "description": "SMILES of the ligand in the new protonation state. Give either candidate_id or smiles.",
                    },
                    "rationale": {
                        "type": "string",
                        "description": "Why this state is chosen (e.g. which predictor is trusted and why, experimental pKa, literature evidence).",
                    },
                    "reference": {
                        "type": "string",
                        "description": "Literature reference (DOI or citation) supporting the state. Required when the SMILES is not one of the candidate states.",
                    },
                },
                "required": ["ligand_file", "rationale"],
            },
        ),
        Tool(
            name="param_ligand",
            description=(
                "If the protein-ligand system contains more than one ligand, only the first one will be kept and parameterized. Nevertheless, the full list of ligands should be parsed as ligand_files."
                "Assign partial charges to a small molecule ligand using Amber's antechamber utility (AM1-BCC charges, GAFF2 atom types). "
                "This tool must be run before run_tleap_ligand. "
                f"The input ligand file is the protonated ligand ({ligand_name}_h.pdb) located in the specified sandbox_dir ({sandbox_dir}); antechamber reads the matching {ligand_name}_h.sdf "
                "so that bond orders and formal charges are kept. The net charge is the one of the selected protonation state. "
                "Before and after antechamber, the tool checks that the ligand files are consistent with the selected protonation state and were not modified outside the tools; "
                "if the check fails, the ligand must be protonated again with prepare_pdb_file_ligand or set_ligand_protonation. "
                f"The output mol2 file has the same name as the input ligand file with a .mol2 extension (e.g. {ligand_name}_h.mol2) and is saved in the same sandbox_dir ({sandbox_dir})."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": f"The directory path where input files are located and output files are produced (this is the sandbox directory): {sandbox_dir}",
                    },
                    "ligand_files": {
                        "type": "array",
                        "items": {
                            "type": "string"
                        },
                        "description": (
                            "List of input PDB file(s) of the protonated ligand(s) to be parameterized "
                            "(without full path). If there is only one ligand, provide a list with one element, "
                            f"e.g. ['{ligand_name}_h.pdb']. If there are multiple ligands, provide a list such as "
                            f"['{ligand_name}_1_h.pdb', '{ligand_name}_2_h.pdb'] where the ligands are indexed starting at 1. "
                            f"All files must exist in the sandbox directory ({sandbox_dir})."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three character residue name of the ligand to extract from the PDB file, in capital letters, called {ligand_name}."
                        ),
                    },
                },
                "required": ["sandbox_dir", "ligand_files", "ligand_name"],
            },
        ),
        Tool(
            name="run_tleap_ligand",
            description=(
                "Parameterize the ligand with the GAFF2 force field using Amber's parmchk2 and tleap utilities. "
                "This tool should ONLY be used for a protein-ligand complex, so if the user specifies that a ligand should be included in the simulation. "
                "It must be run AFTER param_ligand (which creates the ligand mol2 file) and BEFORE run_tleap_protein. "
                f"parmchk2 generates the missing force field parameters ({ligand_name}.frcmod), then tleap generates {ligand_name}.prmtop, {ligand_name}.inpcrd and {ligand_name}.lib. "
                f"ParmEd is then used to convert these files to GROMACS format: the ligand topology {ligand_name}.itp, the ligand atom types {ligand_name}.prm and the ligand coordinates {ligand_name}.gro. "
                "The charge of the last ligand atom is adjusted so that the total ligand charge is an integer. "
                f"All input and output files are read from and written to the sandbox_dir ({sandbox_dir}). "
                "If tleap fails, inspect the mol2 file for wrong atom types or missing atoms."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            f"Absolute path to the working directory where tleap should be executed: {sandbox_dir}. "
                            "All output files will be saved here."
                        ),
                    },
                    "ligand_mol2_file": {
                        "type": "string",
                        "description": (
                            "The ligand mol2 file created by the param_ligand tool (don't provide full path). "
                            f"If there is one ligand, it is called {ligand_name}_h.mol2. If there are multiple ligands, only the first one is used: {ligand_name}_1_h.mol2. "
                            "This file must exist in the provided sandbox_dir."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three-letter residue name of the ligand, in capital letters, called {ligand_name}."
                        ),
                    },
                },
                "required": ["sandbox_dir", "ligand_mol2_file", "ligand_name"],
            },
        ),
        Tool(
            name="run_tleap_protein",
            description=(
                "Parameterize the protein with the Amber force field ff14SB using Amber's tleap utility. "
                "This tool must ALWAYS be run, for both protein-only systems and protein-ligand complexes. "
                "For a protein-ligand complex, it must be run AFTER run_tleap_ligand, because it uses the ligand files "
                f"{ligand_name}.itp, {ligand_name}.prm and {ligand_name}.gro created by run_tleap_ligand. "
                f"It must be run after protonate_protein: tleap loads the protonated protein {pdb_id}_pdb2pqr.pdb. "
                f"tleap generates the Amber topology ({pdb_id}.prmtop), coordinates ({pdb_id}.inpcrd) and {pdb_id}_tleap.pdb. "
                f"ParmEd is then used to convert these files to GROMACS format (topol.top, {pdb_id}_tleap.gro). "
                "The crystal waters are renamed from WAT to SOL in topol.top and the .gro file, to match the GROMACS waters. "
                "topol.top is then edited to add the position restraint includes for the protein chains, the water and ion parameters "
                "and, if a ligand is given, the ligand parameters, topology and position restraint includes. "
                f"For a protein-only system the output structure is {pdb_id}_tleap.gro. "
                "For a protein-ligand complex, the protein and ligand coordinates are combined into complex.gro. "
                f"All input and output files are read from and written to the sandbox_dir ({sandbox_dir}). "
                "If tleap fails, inspect the input file for missing atoms or nonstandard residues."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "sandbox_dir": {
                        "type": "string",
                        "description": (
                            f"Absolute path to the working directory where tleap should be executed: {sandbox_dir}. "
                            "All output files will be saved here."
                        ),
                    },
                    "pdb_id": {
                        "type": "string",
                        "description": (f"The PDB ID of the structure to be prepared: {pdb_id}. "),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            f"The three-letter residue name of the ligand, in capital letters, called {ligand_name}. "
                            "Only provide this argument for a protein-ligand complex. If no ligand was provided by the user, do not input this argument."
                        ),
                    },
                },
                "required": ["sandbox_dir", "pdb_id"],
            },
        ),
        Tool(
            name="gromacs_equil",
            description=(
                """
                Run GROMACS equilibration and optionally the production MD phase.
                This tool should be executed AFTER solvation and ion addition are complete.

                It executes 'equil_Gromacs.sh' inside the working directory provided in
                'sandbox_dir' (solvation, ions, energy minimization, NVT), then runs NPT.
                NVT and NPT last 100 ps each and are checked automatically:
                - after NVT, the mean temperature over the second half must be within 5 K of md_temp
                  (temperature.xvg), otherwise NVT is extended by 100 ps;
                - after NPT, the reference pressure (1 bar) must be within one standard deviation of
                  the mean pressure (pressure.xvg) and the density must have stabilized between
                  900 and 1100 kg/m3 (density.xvg), otherwise NPT is extended by 100 ps.
                The results are reported in the output and saved to equilibration_checks.json.
                If the output starts with WARNING, the checks still failed after the extensions:
                tell the user. In test mode (--test) both stages are kept at 10 ps without extension.
                The system is neutralized with Na+/Cl- counter-ions; the NaCl concentration of the ionic
                environment chosen in the plan (physiological 0.15 M, or neutral: counter-ions only) is applied automatically.

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
                            f"If a protein alone is simulated, the file will be {pdb_id}_tleap.gro. "
                            "This file must be located within sandbox_dir."
                        ),
                    },
                    "md_temp": {
                        "type": "string",
                        "description": (
                            "Temperature in Kelvin for the MD simulation."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            "Three-character residue name of the ligand (capital letters or numbers). If no ligand was provided by user, do not input this argument."
                        ),
                    },
                    "ligand_files": {
                        "type": "array",
                        "items": {
                            "type": "string"
                        },
                        "description": (
                             "List of input PDB file(s) of the protonated ligand(s) to be parameterized "
                            "(without full path). If there is only one ligand, provide a list with one element, "
                            f"e.g. ['{ligand_name}_h.pdb']. If there are multiple ligands, provide a list such as "
                            f"['{ligand_name}_h_1.pdb', '{ligand_name}_h_2.pdb']. "
                            f"All files must exist in the sandbox directory ({sandbox_dir})."
                        ),
                    },
                },
                "required": ["sandbox_dir", "input_gro", "md_temp"],
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
                    "md_temp": {
                        "type": "string",
                        "description": (
                            "Temperature in Kelvin for the MD simulation."
                        ),
                    },
                    "md_duration": {
                        "type": "string",
                        "description": (
                            "Duration of the MD simulation in nanoseconds."
                        ),
                    },
                    "ligand_name": {
                        "type": "string",
                        "description": (
                            "Three-character residue name of the ligand (capital letters or numbers). If no ligand was provided by user, do not input this argument."
                        ),
                    }
                },
                "required": ["sandbox_dir", "input_gro", "npt_cpt_file", "md_temp", "md_duration"],
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
                    "pdb_id": {
                        "type": "string",
                        "description": (f"The PDB ID of the structure used: {pdb_id}. "),
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
                "required": ["sandbox_dir", "pdb_id", "input_xtc"] + (["ligand_name"] if ligand_name else []),
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
                    The salt concentration of the GB (saltcon) and PB (istrng) solvent models is set automatically to the NaCl
                    concentration of the ionic environment chosen in the plan (0.15 M physiological, 0 M neutral).
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
                    "md_temp": {
                        "type": ["string"],
                        "description": ("Temperature used during the MD simulation, found in the md.mdp file located in sandbox_dir. This value is an integer with base 10."),
                    },
                },
                "required": ["sandbox_dir", "pdb_id", "nsteps", "nstxout_compressed", "md_temp"],
            },
        ),
        Tool(
            name="grep_file",
            description=(
                "Search a file for lines matching a regex pattern. Returns matching lines with their line numbers. "
                "Useful for inspecting specific parts of large files (e.g. finding HETATM lines in a PDB, checking specific residues, "
                "finding errors in log files) without reading the entire file."
            ),
            parameters={
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": f"Absolute path to the file to search. Must be inside the sandbox directory ({sandbox_dir}).",
                    },
                    "pattern": {
                        "type": "string",
                        "description": "A regex pattern to match against each line of the file (e.g. '^HETATM', 'ERROR|FATAL', 'ILE A   2').",
                    },
                },
                "required": ["path", "pattern"],
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
