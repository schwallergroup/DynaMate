import hashlib
import json
import shlex
import shutil
import subprocess
import urllib.request
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import dimorphite_dl
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from src.utils import get_class_logger

logger = get_class_logger(__name__, log_to_file=False)

# Pseudo-codes that are not real CCD entries for a specific molecule
GENERIC_LIGAND_CODES = {"UNL", "UNK", "LIG", "MOL", "XXX", "DRG", "INH"}
# Max number of Dimorphite-DL alternative states kept as candidates
MAX_NEARBY_STATES = 10
# Max correction (e) accepted when making the ligand total charge an integer
MAX_CHARGE_CORRECTION = 0.01


class LigandProtonationError(ValueError):
    pass


# ---------------------------------------------------------------- helpers

def record_path(sandbox_dir: str, output_prefix: str) -> Path:
    """Provenance record of the ligand protonation, e.g. LIG_h_protonation.json."""
    return Path(sandbox_dir) / f"{output_prefix}_protonation.json"


def load_record(sandbox_dir: str, output_prefix: str) -> dict | None:
    path = record_path(sandbox_dir, output_prefix)
    if not path.exists():
        return None
    return json.loads(path.read_text())


def save_record(sandbox_dir: str, output_prefix: str, record: dict) -> None:
    record_path(sandbox_dir, output_prefix).write_text(json.dumps(record, indent=2))


def file_hash(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_smiles(smiles: str) -> str | None:
    """Canonical SMILES without stereochemistry: protonation states are compared on connectivity, charges and hydrogens only."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def smiles_charge(smiles: str) -> int:
    return Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))


def neutralize_smiles(smiles: str) -> str:
    """Neutral form of the molecule (where possible), so that the predictors start from the same reference state."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    return Chem.MolToSmiles(rdMolStandardize.Uncharger().uncharge(mol))


# ---------------------------------------------------------------- bond orders

def get_ccd_smiles(ligand_name: str) -> str | None:
    """
    SMILES of the ligand from the PDB Chemical Component Dictionary (RCSB), used as bond-order template.
    Returns None for generic codes or if the entry cannot be fetched.
    """
    if ligand_name.upper() in GENERIC_LIGAND_CODES:
        return None
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_name.upper()}"
    try:
        with urllib.request.urlopen(url, timeout=15) as response:
            data = json.loads(response.read().decode())
    except Exception as e:
        logger.warning(f"Could not fetch CCD entry {ligand_name} from {url}: {e}")
        return None

    descriptor = data.get("rcsb_chem_comp_descriptor", {})
    smiles = descriptor.get("smiles_stereo") or descriptor.get("smiles")
    if not smiles:
        for entry in data.get("pdbx_chem_comp_descriptor", []):
            if entry.get("type", "").upper().startswith("SMILES"):
                smiles = entry.get("descriptor")
                break
    if smiles:
        logger.info(f"CCD SMILES of {ligand_name}: {smiles}")
    return smiles


def perceive_smiles_openbabel(sandbox_dir: str, ligand_pdb_file: str, output_prefix: str) -> str:
    """Bond orders perceived by Open Babel from the 3D pose (fallback when no CCD template is available)."""
    smi_file = f"{output_prefix}_input.smi"
    cmd = f"obabel {ligand_pdb_file} -O {smi_file}"
    run_1 = subprocess.run(shlex.split(cmd), cwd=sandbox_dir, capture_output=True, text=True)
    if run_1.returncode != 0:
        error_text = "\n".join(filter(None, [run_1.stderr, run_1.stdout]))
        logger.error(f"Error perceiving SMILES from Open Babel: {error_text}")
        raise RuntimeError(f"Open Babel could not perceive a SMILES from {ligand_pdb_file}: {error_text}")
    smiles = Path(f"{sandbox_dir}/{smi_file}").read_text().split()[0]
    logger.info(f"SMILES perceived by Open Babel from {ligand_pdb_file}: {smiles}")
    return smiles


# ---------------------------------------------------------------- protonation predictors

def dimorphite_states(smiles: str, ph: float) -> tuple[str, list[str]]:
    """Most probable state at pH (precision=0) and states within 1 pKa unit of pH (precision=1)."""
    states = dimorphite_dl.protonate_smiles(smiles, ph_min=ph, ph_max=ph, precision=0.0)
    if not states:
        logger.error(f"Dimorphite-DL returned no protonation state for {smiles}")
        raise RuntimeError(f"Dimorphite-DL returned no protonation state for {smiles}")
    if len(states) > 1:
        logger.warning(f"Dimorphite-DL returned several states at pH {ph}: {states}. Using the first one.")
    nearby_states = dimorphite_dl.protonate_smiles(smiles, ph_min=ph, ph_max=ph, precision=1.0)
    logger.info(f"Dimorphite-DL: state at pH {ph}: {states[0]}; states within 1 pKa unit: {nearby_states}")
    return states[0], list(nearby_states)[:MAX_NEARBY_STATES]


def pkasolver_state(smiles: str, ph: float) -> dict:
    """
    State at pH predicted by pkasolver. pkasolver predicts the pKa of each sequential (de)protonation step;
    sorted by pKa, the deprotonated form of step i is the protonated form of step i+1. The state at pH is the
    deprotonated form of the last step with pKa < pH, or the protonated form of the first step if all pKa > pH.
    Returns {"available", "smiles", "pka_values", "error"}; pkasolver is optional and never stops the pipeline.
    """
    try:
        from pkasolver.query import calculate_microstate_pka_values
    except Exception as e:
        logger.warning(f"pkasolver is not available, only Dimorphite-DL is used: {e}")
        return {"available": False, "smiles": None, "pka_values": [], "error": f"pkasolver not available: {e}"}

    try:
        steps = calculate_microstate_pka_values(Chem.MolFromSmiles(smiles))
        steps = sorted(steps, key=lambda step: step.pka)
        pka_values = [round(float(step.pka), 2) for step in steps]
        if not steps:
            state_smiles = smiles
        else:
            below = [step for step in steps if step.pka < ph]
            state_mol = below[-1].deprotonated_mol if below else steps[0].protonated_mol
            state_smiles = Chem.MolToSmiles(Chem.RemoveHs(state_mol))
        logger.info(f"pkasolver: pKa values {pka_values}, state at pH {ph}: {state_smiles}")
        return {"available": True, "smiles": state_smiles, "pka_values": pka_values, "error": None}
    except Exception as e:
        logger.warning(f"pkasolver failed for {smiles}: {type(e).__name__}: {e}")
        return {"available": True, "smiles": None, "pka_values": [], "error": f"{type(e).__name__}: {e}"}


# ---------------------------------------------------------------- build ligand files from a state

def build_ligand_from_smiles(sandbox_dir: str, ligand_pdb_file: str, prot_smiles: str, output_prefix: str, ligand_name: str) -> Chem.Mol:
    """
    Puts the protonation state prot_smiles on the crystal pose: the heavy atoms and connectivity of the pose are kept,
    bond orders, formal charges and hydrogens come from prot_smiles (RDKit AssignBondOrdersFromTemplate),
    hydrogens get coordinates (Chem.AddHs(addCoords=True)).
    Writes {output_prefix}.sdf (bond orders + formal charges) and {output_prefix}.pdb (atom names C1, C2, ..., H1, ...).
    Raises LigandProtonationError if the state does not match the pose or is not chemically valid.
    """
    pose = Chem.MolFromPDBFile(f"{sandbox_dir}/{ligand_pdb_file}", removeHs=True)
    if pose is None:
        raise LigandProtonationError(f"RDKit could not read {ligand_pdb_file}")
    template = Chem.MolFromSmiles(prot_smiles)
    if template is None:
        raise LigandProtonationError(f"Invalid SMILES: {prot_smiles}")
    # AssignBondOrdersFromTemplate also accepts a template matching only part of the pose: check the whole molecule
    if template.GetNumAtoms() != pose.GetNumAtoms():
        raise LigandProtonationError(f"SMILES {prot_smiles} has {template.GetNumAtoms()} heavy atoms but {ligand_pdb_file} has {pose.GetNumAtoms()}")
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, pose)
    except ValueError as e:
        raise LigandProtonationError(f"SMILES {prot_smiles} does not match the heavy atoms of {ligand_pdb_file}: {e}")
    if Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False) != Chem.MolToSmiles(template, isomericSmiles=False):
        raise LigandProtonationError(f"SMILES {prot_smiles} does not match the connectivity of {ligand_pdb_file}")
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        raise LigandProtonationError(f"State {prot_smiles} is not chemically valid on the pose: {e}")
    if Descriptors.NumRadicalElectrons(mol) != 0:
        raise LigandProtonationError(f"State {prot_smiles} has radical electrons: hydrogens are missing")

    mol = Chem.AddHs(mol, addCoords=True)
    Chem.AssignStereochemistryFrom3D(mol)
    mol.SetProp("_Name", ligand_name)

    # Atom names element+counter (C1, C2, ..., H1, ...), right format for antechamber/GROMACS
    counters = defaultdict(int)
    for atom in mol.GetAtoms():
        element = atom.GetSymbol()
        counters[element] += 1
        name = f"{element}{counters[element]}"
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f" {name:<3}" if len(element) == 1 and len(name) < 4 else f"{name:<4}")
        info.SetResidueName(ligand_name)
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)

    Chem.MolToMolFile(mol, f"{sandbox_dir}/{output_prefix}.sdf")
    Chem.MolToPDBFile(mol, f"{sandbox_dir}/{output_prefix}.pdb")
    logger.info(f"Protonated ligand ({mol.GetNumAtoms()} atoms, net formal charge {Chem.GetFormalCharge(mol)}) "
                f"written to {sandbox_dir}/{output_prefix}.sdf and {sandbox_dir}/{output_prefix}.pdb")
    return mol


def _selection(candidate: dict, source: str, rationale: str, reference: str | None) -> dict:
    return {
        "id": candidate["id"],
        "smiles": candidate["smiles"],
        "net_charge": candidate["net_charge"],
        "source": source,
        "matches_predictor": bool(candidate["sources"]),
        "rationale": rationale,
        "reference": reference,
        "time": datetime.now().isoformat(timespec="seconds"),
    }


def _store_file_hashes(sandbox_dir: str, output_prefix: str, record: dict) -> None:
    record["files"] = {
        f"{output_prefix}.sdf": file_hash(Path(sandbox_dir) / f"{output_prefix}.sdf"),
        f"{output_prefix}.pdb": file_hash(Path(sandbox_dir) / f"{output_prefix}.pdb"),
    }


def summarize_record(record: dict) -> str:
    """Short text for the agent: chosen state, agreement between predictors, alternatives."""
    selected = record["selected"]
    lines = [
        f"Ligand {record['ligand_name']} protonation at pH {record['ph']} (bond orders from {record['bond_order_source']}): "
        f"selected state {selected['id']} {selected['smiles']} (net charge {selected['net_charge']:+d}, source: {selected['source']}).",
        f"Comparison Dimorphite-DL vs pkasolver: {record['comparison'].upper()}."
        + (" Ambiguous: several states are within 1 pKa unit of the pH." if record["ambiguous"] else ""),
        "Candidate states: " + "; ".join(
            f"{c['id']} {c['smiles']} (charge {c['net_charge']:+d}, from {', '.join(c['sources']) or 'agent'})" for c in record["candidates"]
        ),
    ]
    if record["pkasolver"]["error"]:
        lines.append(f"pkasolver: {record['pkasolver']['error']}")
    if record["comparison"] in ("minor_disagreement", "major_disagreement") or record["ambiguous"]:
        lines.append(
            "ACTION: the predictors disagree or the state is ambiguous. Look for literature evidence (search_papers) on the protonation "
            "state of this ligand, then either keep the selected state or call set_ligand_protonation with a candidate_id (or a SMILES "
            "and a reference) and a rationale. Do NOT edit the ligand files directly."
        )
    return "\n".join(lines)


# ---------------------------------------------------------------- main entry points

def protonate_ligand(sandbox_dir: str, ligand_pdb_file: str, output_prefix: str, ligand_name: str, ph: float = 7.0) -> dict:
    """
    Protonates a ligand at the given pH while keeping its heavy-atom coordinates.
      1. Bond orders: from the PDB Chemical Component Dictionary entry of ligand_name, otherwise perceived by Open Babel from the pose.
      2. Protonation state: Dimorphite-DL state at `ph` (applied by default) compared with the pkasolver state at `ph`.
      3. All distinct states (Dimorphite-DL state and states within 1 pKa unit, pkasolver state) are stored as candidates
         in {output_prefix}_protonation.json, with the agreement between predictors.
      4. The selected state is put on the pose (build_ligand_from_smiles): writes {output_prefix}.sdf and {output_prefix}.pdb.

    Args:
        sandbox_dir (str): Directory where the ligand PDB file is located and where output files will be saved.
        ligand_pdb_file (str): Extracted ligand pose (heavy atoms, hydrogens are ignored), relative to sandbox_dir.
        output_prefix (str): Name of the output files without extension, relative to sandbox_dir (e.g. UNL_h).
        ligand_name (str): Residue name written in the output PDB.
        ph (float): pH at which the protonation state is chosen.
    Returns:
        The provenance record (dict).
    """
    # 1. Bond orders: CCD template if it matches the pose, otherwise Open Babel
    input_smiles, bond_order_source = None, None
    ccd_smiles = get_ccd_smiles(ligand_name)
    if ccd_smiles:
        try:
            build_ligand_from_smiles(sandbox_dir, ligand_pdb_file, ccd_smiles, output_prefix, ligand_name)
            input_smiles, bond_order_source = ccd_smiles, "ccd"
        except LigandProtonationError as e:
            logger.warning(f"CCD SMILES of {ligand_name} does not match the pose in {ligand_pdb_file} ({e}), using Open Babel bond orders")
    if input_smiles is None:
        input_smiles, bond_order_source = perceive_smiles_openbabel(sandbox_dir, ligand_pdb_file, output_prefix), "openbabel"

    # 2. Predictors, both starting from the neutral form
    neutral_smiles = neutralize_smiles(input_smiles)
    dimorphite_smiles, nearby_states = dimorphite_states(neutral_smiles, ph)
    pkasolver_result = pkasolver_state(neutral_smiles, ph)

    # 3. Candidates (deduplicated on canonical SMILES)
    candidates: list[dict] = []

    def add_candidate(smiles: str, source: str) -> str | None:
        canonical = canonical_smiles(smiles)
        if canonical is None:
            return None
        for candidate in candidates:
            if candidate["canonical"] == canonical:
                if source not in candidate["sources"]:
                    candidate["sources"].append(source)
                return candidate["id"]
        candidate_id = f"C{len(candidates) + 1}"
        candidates.append({"id": candidate_id, "smiles": smiles, "canonical": canonical,
                           "net_charge": smiles_charge(smiles), "sources": [source]})
        return candidate_id

    dimorphite_id = add_candidate(dimorphite_smiles, "dimorphite")
    pkasolver_id = add_candidate(pkasolver_result["smiles"], "pkasolver") if pkasolver_result["smiles"] else None
    nearby_ids = [add_candidate(s, "dimorphite_nearby") for s in nearby_states]

    dimorphite_candidate = next(c for c in candidates if c["id"] == dimorphite_id)
    if pkasolver_id is None:
        comparison = "single_predictor"
    elif pkasolver_id == dimorphite_id:
        comparison = "agree"
    else:
        pkasolver_candidate = next(c for c in candidates if c["id"] == pkasolver_id)
        same_charge = pkasolver_candidate["net_charge"] == dimorphite_candidate["net_charge"]
        comparison = "minor_disagreement" if same_charge else "major_disagreement"

    ambiguous = len({i for i in nearby_ids if i}) > 1 or any(abs(p - ph) < 1.0 for p in pkasolver_result["pka_values"])

    # 4. Apply the Dimorphite-DL state by default (raises LigandProtonationError if it does not fit the pose)
    build_ligand_from_smiles(sandbox_dir, ligand_pdb_file, dimorphite_smiles, output_prefix, ligand_name)

    record = {
        "ligand_name": ligand_name,
        "pose_file": ligand_pdb_file,
        "output_prefix": output_prefix,
        "ph": ph,
        "bond_order_source": bond_order_source,
        "input_smiles": input_smiles,
        "neutral_smiles": neutral_smiles,
        "candidates": candidates,
        "dimorphite": {"state": dimorphite_id, "nearby_states": [i for i in nearby_ids if i]},
        "pkasolver": {**pkasolver_result, "state": pkasolver_id},
        "comparison": comparison,
        "ambiguous": ambiguous,
        "selected": _selection(dimorphite_candidate, "dimorphite (default)", "Default: most probable Dimorphite-DL state at the simulation pH", None),
        "history": [],
        "manual_override": False,
    }
    _store_file_hashes(sandbox_dir, output_prefix, record)
    save_record(sandbox_dir, output_prefix, record)
    logger.info(summarize_record(record))
    return record


def set_ligand_protonation(sandbox_dir: str, ligand_file: str, rationale: str, candidate_id: str | None = None,
                           smiles: str | None = None, reference: str | None = None) -> str:
    """
    The only way for the agent to change the protonation state of a ligand. The new state is either one of the candidate
    states of {prefix}_protonation.json (candidate_id) or a SMILES (e.g. from the literature, reference required if it is
    not a candidate). The state is checked against the crystal pose, the ligand .sdf/.pdb are rebuilt, the decision is
    logged in the provenance record and all files derived from the previous state are moved to a stale_protonation_* folder.
    """
    prefix = Path(ligand_file).stem
    record = load_record(sandbox_dir, prefix)
    if record is None:
        return f"Ligand protonation change refused: no protonation record {record_path(sandbox_dir, prefix).name} found. Run prepare_pdb_file_ligand first."
    if not rationale or not rationale.strip():
        return "Ligand protonation change refused: a rationale is required."

    candidates = record["candidates"]
    if candidate_id:
        candidate = next((c for c in candidates if c["id"] == candidate_id), None)
        if candidate is None:
            return f"Ligand protonation change refused: unknown candidate_id {candidate_id}. Candidates: {[c['id'] for c in candidates]}"
    elif smiles:
        canonical = canonical_smiles(smiles)
        if canonical is None:
            return f"Ligand protonation change refused: invalid SMILES {smiles}"
        candidate = next((c for c in candidates if c["canonical"] == canonical), None)
        if candidate is None:
            if not reference or not reference.strip():
                return ("Ligand protonation change refused: this state was not predicted by Dimorphite-DL or pkasolver, "
                        "a literature reference (DOI or citation) is required.")
            candidate = {"id": f"C{len(candidates) + 1}", "smiles": smiles, "canonical": canonical,
                         "net_charge": smiles_charge(smiles), "sources": []}
    else:
        return "Ligand protonation change refused: provide a candidate_id or a SMILES."

    if candidate["id"] == record["selected"]["id"]:
        return f"State {candidate['id']} is already the selected state, nothing changed."

    try:
        build_ligand_from_smiles(sandbox_dir, record["pose_file"], candidate["smiles"], prefix, record["ligand_name"])
    except LigandProtonationError as e:
        return f"Ligand protonation change refused: {e}"

    if candidate not in candidates:
        candidates.append(candidate)
    source = "literature" if not candidate["sources"] else ", ".join(candidate["sources"])
    record["history"].append(record["selected"])
    record["selected"] = _selection(candidate, f"agent: {source}", rationale.strip(), reference)
    record["manual_override"] = record["manual_override"] or not candidate["sources"]
    _store_file_hashes(sandbox_dir, prefix, record)
    moved = _move_stale_files(sandbox_dir, prefix, record["ligand_name"])
    save_record(sandbox_dir, prefix, record)

    logger.info(f"Ligand protonation changed to {candidate['id']} {candidate['smiles']} by the agent: {rationale} (reference: {reference})")
    return (
        f"Ligand protonation changed to {candidate['id']} {candidate['smiles']} (net charge {candidate['net_charge']:+d}). "
        f"{prefix}.sdf and {prefix}.pdb were rebuilt and the decision was logged in {record_path(sandbox_dir, prefix).name}. "
        + (f"Files derived from the previous state were moved to {moved}. " if moved else "")
        + "You must now run param_ligand, run_tleap_ligand and run_tleap_protein again, then the GROMACS steps."
    )


def _move_stale_files(sandbox_dir: str, prefix: str, ligand_name: str) -> str | None:
    """Moves files built from the previous ligand state (parameters, topology, MD outputs) to a stale_protonation_* folder."""
    sandbox = Path(sandbox_dir)
    patterns = [
        f"{prefix}.mol2", f"{ligand_name}.frcmod", f"{ligand_name}.lib", f"{ligand_name}.prmtop", f"{ligand_name}.inpcrd",
        f"{ligand_name}.itp", f"{ligand_name}_backup.itp", f"{ligand_name}.prm", f"{ligand_name}.gro",
        f"posre_{ligand_name}.itp", f"index_{ligand_name}.ndx",
        "complex.gro", "complex_backup.gro", "topol.top", "topol_without_posre.top",
        "em.*", "nvt.*", "npt.*", "md.*", "*.xvg",
    ]
    stale = sorted({p for pattern in patterns for p in sandbox.glob(pattern) if p.is_file() and p.suffix != ".mdp"})
    if not stale:
        return None
    stale_dir = sandbox / f"stale_protonation_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    stale_dir.mkdir()
    for path in stale:
        shutil.move(str(path), stale_dir / path.name)
    logger.info(f"Moved {len(stale)} files built from the previous ligand state to {stale_dir}")
    return str(stale_dir)


def check_ligand_consistency(sandbox_dir: str, prefix: str, mol2_file: str | None = None) -> list[str]:
    """
    Checks that the ligand files all describe the selected protonation state of {prefix}_protonation.json:
      - {prefix}.sdf and {prefix}.pdb were not modified outside protonate_ligand/set_ligand_protonation (file hashes)
      - net charge of the .sdf equals the selected state, same elements in the same order in .sdf, .pdb and .mol2
      - sum of the partial charges of the .mol2 equals the net charge
    Returns the list of problems (empty if consistent).
    """
    sandbox = Path(sandbox_dir)
    record = load_record(sandbox_dir, prefix)
    if record is None:
        return [f"No protonation record {record_path(sandbox_dir, prefix).name}: the ligand was not protonated with prepare_pdb_file_ligand."]

    problems = []
    for name, expected_hash in record.get("files", {}).items():
        path = sandbox / name
        if not path.exists():
            problems.append(f"{name} is missing")
        elif file_hash(path) != expected_hash:
            problems.append(f"{name} was modified outside prepare_pdb_file_ligand/set_ligand_protonation")

    sdf = Chem.MolFromMolFile(str(sandbox / f"{prefix}.sdf"), removeHs=False)
    if sdf is None:
        return problems + [f"RDKit could not read {prefix}.sdf"]
    net_charge = record["selected"]["net_charge"]
    if Chem.GetFormalCharge(sdf) != net_charge:
        problems.append(f"{prefix}.sdf has net charge {Chem.GetFormalCharge(sdf)}, selected state has {net_charge}")
    sdf_elements = [a.GetSymbol() for a in sdf.GetAtoms()]

    pdb = Chem.MolFromPDBFile(str(sandbox / f"{prefix}.pdb"), removeHs=False, sanitize=False)
    if pdb is None:
        problems.append(f"RDKit could not read {prefix}.pdb")
    elif [a.GetSymbol() for a in pdb.GetAtoms()] != sdf_elements:
        problems.append(f"{prefix}.pdb and {prefix}.sdf do not contain the same atoms ({pdb.GetNumAtoms()} vs {len(sdf_elements)} atoms)")

    if mol2_file is not None:
        mol2_path = sandbox / mol2_file
        if not mol2_path.exists():
            problems.append(f"{mol2_file} is missing")
        else:
            atom_block = mol2_path.read_text().split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0]
            atom_lines = [l.split() for l in atom_block.strip().splitlines() if l.strip()]
            # GAFF/GAFF2 hydrogen atom types all start with "h" (ha, hc, hn, ho, h1...), heavy-atom types never do
            mol2_n_h = sum(1 for fields in atom_lines if fields[5].lower().startswith("h"))
            sdf_n_h = sdf_elements.count("H")
            if len(atom_lines) != len(sdf_elements):
                problems.append(f"{mol2_file} has {len(atom_lines)} atoms, {prefix}.sdf has {len(sdf_elements)}")
            if mol2_n_h != sdf_n_h:
                problems.append(f"{mol2_file} has {mol2_n_h} hydrogens, {prefix}.sdf has {sdf_n_h}")
            total_charge = sum(float(fields[8]) for fields in atom_lines)
            if abs(total_charge - net_charge) > MAX_CHARGE_CORRECTION:
                problems.append(f"Sum of partial charges in {mol2_file} is {total_charge:.4f}, selected state has net charge {net_charge}")

    for problem in problems:
        logger.error(f"Ligand consistency check: {problem}")
    return problems
