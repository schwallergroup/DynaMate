from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean
from typing import Optional


def parse_system_name(system_name: str) -> tuple[str, Optional[str]]:
    """Return (pdb_id, ligand_name) from a system name like '1FDH_GLY' or '1FDH_None'."""
    parts = system_name.split("_")
    pdb_id = parts[0]
    ligand_name = parts[1] if len(parts) > 1 and parts[1].upper() != "NONE" else None
    return pdb_id, ligand_name


def build_necessary_steps_protein() -> dict[str, list[str]]:
    """Expected output files per pipeline step for a protein-only system."""
    return {
        "fetch_and_save_pdb": ["{pdb_id}.pdb"],
        "prepare_pdb_file_ligand": ["{pdb_id}_prepared.pdb"],
        "add_caps": ["{pdb_id}_prepared_capped.pdb"],
        "rename_histidines": ["{pdb_id}_prepared_capped_his.pdb"],
        "run_tleap": [
            "{pdb_id}.prmtop",
            "{pdb_id}.inpcrd",
            "{pdb_id}_tleap.pdb",
            "topol.top",
            "{pdb_id}.gro",
        ],
        "gromacs_equil": [
            "topol_without_posre.top",
            "em.gro",
            "nvt.gro",
            "nvt.xtc",
            "npt.gro",
            "npt.xtc",
            "temperature.xvg",
            "pressure.xvg",
            "density.xvg",
            "potential.xvg",
        ],
        "gromacs_production": ["md.gro", "md.xtc"],
        "gromacs_analysis": [
            "rmsd.xvg",
            "rmsd_xtal.xvg",
            "rmsf.xvg",
            "gyrate.xvg",
            "hbnum_prot_wat.xvg",
            "hbnum_sidechain.xvg",
        ],
    }


def build_necessary_steps_ligand() -> dict[str, list[str]]:
    """Expected output files per pipeline step for a protein-ligand system."""
    return {
        "fetch_and_save_pdb": ["{pdb_id}.pdb"],
        "prepare_pdb_file_ligand": [
            "{pdb_id}_prepared.pdb",
            "{ligand_name}.pdb",
            "{ligand_name}_h.pdb",
        ],
        "add_caps": ["{pdb_id}_prepared_capped.pdb"],
        "rename_histidines": ["{pdb_id}_prepared_capped_his.pdb"],
        "param_ligand": ["{ligand_name}.prepi", "{ligand_name}.frcmod"],
        "run_tleap_ligand": [
            "complex.pdb",
            "complex.prmtop",
            "complex.inpcrd",
            "complex_tleap.pdb",
            "topol.top",
            "complex.gro",
        ],
        "gromacs_equil": [
            "topol_without_posre.top",
            "{ligand_name}.gro",
            "em.gro",
            "nvt.gro",
            "nvt.xtc",
            "npt.gro",
            "npt.xtc",
            "temperature.xvg",
            "pressure.xvg",
            "density.xvg",
            "potential.xvg",
        ],
        "gromacs_production": ["md.gro", "md.xtc"],
        "gromacs_analysis": [
            "rmsd.xvg",
            "rmsd_xtal.xvg",
            "rmsf.xvg",
            "gyrate.xvg",
            "hbnum_prot_lig.xvg",
            "hbnum_prot_wat.xvg",
            "hbnum_sidechain.xvg",
        ],
    }


def list_nonempty_files(directory: str | Path) -> list[str]:
    """Return basenames of all non-empty files in *directory* (non-recursive)."""
    directory = Path(directory)
    if not directory.exists():
        return []
    return [
        f
        for f in os.listdir(directory)
        if os.path.isfile(directory / f) and os.path.getsize(directory / f) > 0
    ]


def _template_satisfied(template: str, pdb_id: str, ligand_name: str, nonempty_files: list[str]) -> bool:
    """Check whether a single file template is satisfied by any file in nonempty_files.

    For ligand templates (containing {ligand_name}), matching is flexible:
    the filename just needs to contain the ligand name and share the same extension.
    This handles hydrogenation suffixes (_h), multiple chains (_1, _1_h), etc.

    For protein/fixed templates, exact match after placeholder substitution is used.
    """
    ext = os.path.splitext(template)[1].lower()

    if "{ligand_name}" in template:
        return any(
            os.path.splitext(f)[1].lower() == ext and ligand_name.lower() in f.lower()
            for f in nonempty_files
        )

    expanded = template.replace("{pdb_id}", pdb_id)
    return expanded in nonempty_files


def score_steps(
    steps: dict[str, list[str]],
    pdb_id: str,
    ligand_name: Optional[str],
    nonempty_files: list[str],
) -> dict[str, float]:
    """Compute per-step completion ratios.

    Each step score = (templates satisfied) / (total templates).
    Ligand file templates use flexible matching (contains ligand name + correct extension).
    Protein/fixed file templates use exact matching after {pdb_id} substitution.
    """
    ligand_name = ligand_name or ""
    step_scores: dict[str, float] = {}

    for step, templates in steps.items():
        satisfied = sum(
            1 for t in templates
            if _template_satisfied(t, pdb_id, ligand_name, nonempty_files)
        )
        step_scores[step] = satisfied / len(templates) if templates else 0.0

    return step_scores


@dataclass
class MDEvalResult:
    system_name: str
    pdb_id: str
    ligand_name: Optional[str]
    is_protein_ligand: bool
    step_scores: dict[str, float] = field(default_factory=dict)
    overall_score: float = 0.0
    pipeline_successful: bool = False
    nonempty_files: list[str] = field(default_factory=list)


def evaluate_md_run(
    sandbox_dir: str | Path,
    system_name: str,
) -> MDEvalResult:
    """Evaluate a single MD run by checking non-empty output files.

    Args:
        sandbox_dir: Directory where the agent wrote output files.
        system_name: System identifier in 'PDBID_LIGAND' or 'PDBID_None' format.
                     The protein-only vs protein-ligand step definitions are
                     selected automatically from this name.

    Returns:
        MDEvalResult with per-step scores and overall metrics.
    """
    pdb_id, ligand_name = parse_system_name(system_name)
    steps = build_necessary_steps_ligand() if ligand_name else build_necessary_steps_protein()

    nonempty = list_nonempty_files(sandbox_dir)
    step_scores = score_steps(steps, pdb_id, ligand_name, nonempty)
    overall = mean(step_scores.values()) if step_scores else 0.0
    successful = bool(step_scores) and all(v == 1.0 for v in step_scores.values())

    return MDEvalResult(
        system_name=system_name,
        pdb_id=pdb_id,
        ligand_name=ligand_name,
        is_protein_ligand=ligand_name is not None,
        step_scores=step_scores,
        overall_score=overall,
        pipeline_successful=successful,
        nonempty_files=nonempty,
    )


def evaluate_md_runs(
    systems: list[tuple[str | Path, str]],
) -> list[MDEvalResult]:
    """Evaluate multiple MD runs, each toggling protein/ligand independently.

    Args:
        systems: List of (sandbox_dir, system_name) pairs.

    Returns:
        List of MDEvalResult, one per system.
    """
    return [evaluate_md_run(sandbox_dir, system_name) for sandbox_dir, system_name in systems]
