from src.tools.amber_tools import run_tleap_protein, run_tleap_ligand
from src.tools.gromacs_tools import gromacs_equil, gromacs_production, gromacs_analysis
from src.tools.ligand_tools import param_ligand
from src.tools.ligand_protonation import set_ligand_protonation
from src.tools.pdb_tools import prepare_pdb_file_ligand, add_caps, protonate_protein, rename_histidines, fetch_and_save_pdb
from src.tools.coding_tools import read_file, edit_file, list_files, find_input, grep_file
from src.tools.RAG_tools import search_papers
from src.tools.MMPBSA import run_gmxMMPBSA
from src import constants


def truncate_file_output(full_content: str) -> str:
    """Truncate long file outputs to prevent LLM token overload."""
    if len(full_content) <= 2 * constants.MAX_CHARACTERS_TO_LOG:
        return full_content
    return f"{full_content[:constants.MAX_CHARACTERS_TO_LOG]}... truncated ... {full_content[-constants.MAX_CHARACTERS_TO_LOG:]}" if full_content else ""


TOOL_MAP = {
    # Basic File Operations
    "read_file": lambda _, i: read_file(i['path']),
    "grep_file": lambda _, i: grep_file(i['path'], i['pattern']),
    "list_files": lambda s, _: list_files(s.sandbox_dir),
    "find_input": lambda s, _: find_input(s.sandbox_dir),
    "edit_file": lambda s, i: edit_file(i["path"], i["old_text"], i["new_text"], ligand_name=s.ligand_name),
    # Protein prep
    "fetch_and_save_pdb": lambda s, i: fetch_and_save_pdb(s.sandbox_dir, i["pdb_id"], i["output_pdb"]),
    "prepare_pdb_file_ligand": lambda s, i: prepare_pdb_file_ligand(s.sandbox_dir, i["pdb_id"], i["ph"], i["ligand_name"]),
    "add_caps": lambda s, i: add_caps(s.sandbox_dir, i["input_pdb"], i["pdb_id"]),
    "protonate_protein": lambda s, i: protonate_protein(s.sandbox_dir, i["input_pdb"], i["pdb_id"], i["ph"]),
    "rename_histidines": lambda s, i: rename_histidines(s.sandbox_dir, i["input_pdb"], i["pdb_id"]),
    # Ligand handling
    "set_ligand_protonation": lambda s, i: set_ligand_protonation(
        s.sandbox_dir, i["ligand_file"], i["rationale"], candidate_id=i.get("candidate_id"), smiles=i.get("smiles"), reference=i.get("reference")
    ),
    "param_ligand": lambda s, i: param_ligand(s.sandbox_dir, i["ligand_files"] if isinstance(i["ligand_files"], list) else [i["ligand_files"]], i["ligand_name"]),
    # AMBER-related
    "run_tleap_ligand": lambda s, i: run_tleap_ligand(s.sandbox_dir, i["ligand_mol2_file"], i["ligand_name"]),
    "run_tleap_protein": lambda s, i: run_tleap_protein(s.sandbox_dir, i["pdb_id"], ligand_name=i.get("ligand_name")),
    # GROMACS-related
    "gromacs_equil": lambda s, i: gromacs_equil(
        s.sandbox_dir, i["input_gro"], i["md_temp"], ligand_name=i.get("ligand_name"), ligand_files=i.get("ligand_files"),
        test=getattr(s, "test_mode", False), ionic_conc=getattr(s, "ionic_conc", 0.15),
    ),
    "gromacs_production": lambda s, i: gromacs_production(
        s.sandbox_dir, i["input_gro"], i["npt_cpt_file"], i["md_temp"], i["md_duration"], ligand_name=i.get("ligand_name")
    ),
    "gromacs_analysis": lambda s, i: gromacs_analysis(s.sandbox_dir, i["pdb_id"], i["input_xtc"], ligand_name=i.get("ligand_name")),
    # MMPBSA-related
    "run_gmxMMPBSA": lambda s, i: run_gmxMMPBSA(
        s.sandbox_dir, i["pdb_id"], i["nsteps"], i["nstxout_compressed"], i["md_temp"],
        ionic_conc=getattr(s, "ionic_conc", 0.15),
    ),
    # # RAG tools
    "search_papers": lambda _, i: search_papers(i["query"]),
}
