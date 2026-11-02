import os
import re
from pathlib import Path


def find_input(directory: Path) -> str:
    try:
        items = []
        for file in sorted(os.listdir(directory)):
            item_path = os.path.join(directory, file)
            if os.path.isdir(item_path):
                pass
            else:
                items.append(file)

        if not items:
            return f"Empty directory: {directory}. You forgot to upload your PDB file"

        pdb_files = list(directory.glob("*.pdb"))
        pdb_files_str = "\n".join(str(f) for f in pdb_files)

        return f"User uploads to {directory}:\n{pdb_files_str}"

    except Exception as e:
        return f"Error listing files: {str(e)}"


def read_file(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
        return f"File contents of {path}:\n{content}"
    except FileNotFoundError:
        return f"File not found: {path}"
    except Exception as e:
        return f"Error reading file: {str(e)}"


def list_files(path: str) -> str:
    try:
        if not os.path.exists(path):
            return f"Path not found: {path}"

        items = []
        for item in sorted(os.listdir(path)):
            item_path = os.path.join(path, item)
            if os.path.isdir(item_path):
                items.append(f"[DIR]  {item}/")
            else:
                items.append(f"[FILE] {item}")

        if not items:
            return f"Empty directory: {path}"

        return f"Contents of {path}:\n" + "\n".join(items)
    except Exception as e:
        return f"Error listing files: {str(e)}"


def grep_file(path: str, pattern: str) -> str:
    """Search a file for lines matching a regex pattern and return them with line numbers."""
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        regex = re.compile(pattern)
        matches = []
        for i, line in enumerate(lines, start=1):
            if regex.search(line):
                matches.append(f"{i}: {line.rstrip()}")

        if not matches:
            return f"No matches for pattern '{pattern}' in {path}"

        return f"Matches for '{pattern}' in {path} ({len(matches)} hits):\n" + "\n".join(matches)
    except FileNotFoundError:
        return f"File not found: {path}"
    except re.error as e:
        return f"Invalid regex pattern '{pattern}': {e}"
    except Exception as e:
        return f"Error searching file: {str(e)}"


# Ligand structure, parameter and provenance files: only changed through prepare_pdb_file_ligand / set_ligand_protonation
PROTECTED_LIGAND_SUFFIXES = {".pdb", ".sdf", ".mol2", ".frcmod", ".lib", ".prmtop", ".inpcrd", ".itp", ".prm", ".gro", ".json", ".smi"}


def edit_file(path: str, old_text: str, new_text: str, ligand_name: str | None = None) -> str:
    name = Path(path).name
    if ligand_name and name.upper().startswith(ligand_name.upper()) and Path(path).suffix in PROTECTED_LIGAND_SUFFIXES:
        return (
            f"Editing {name} is not allowed: ligand structure and parameter files can only be changed through the tools. "
            "To change the protonation state of the ligand, use set_ligand_protonation (with a candidate_id, or a SMILES "
            "and a literature reference) and a rationale, then run param_ligand, run_tleap_ligand and run_tleap_protein again."
        )
    try:
        if os.path.exists(path) and old_text:
            with open(path, "r", encoding="utf-8") as f:
                content = f.read()

            if old_text not in content:
                return f"Text not found in file: {old_text}"

            content = content.replace(old_text, new_text)

            with open(path, "w", encoding="utf-8") as f:
                f.write(content)

            return f"Successfully edited {path}"
        else:
            # Only create directory if path contains subdirectories
            dir_name = os.path.dirname(path)
            if dir_name:
                os.makedirs(dir_name, exist_ok=True)

            with open(path, "w", encoding="utf-8") as f:
                f.write(new_text)

            return f"Successfully created {path}"
    except Exception as e:
        return f"Error editing file: {str(e)}"
