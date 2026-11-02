import os
import subprocess
from pathlib import Path
import re
import shutil
import sys
import json
import numpy as np
from src import constants
from src.utils import get_class_logger
from collections import defaultdict, Counter

logger = get_class_logger(__name__, log_to_file=False)


def _filter_mdrun_output(output: str) -> str:
    """Remove verbose mdrun simulation output from GROMACS log."""
    output = re.sub(
        r"Steepest Descents:.*?writing lowest energy coordinates\.\n?",
        "", output, flags=re.DOTALL
    )
    output = re.sub(
        r"starting mdrun 'Generic title'.*?Writing final coordinates\.\n?",
        "", output, flags=re.DOTALL
    )
    output = re.sub(r"^GROMACS reminds you:.*$\n?", "", output, flags=re.MULTILINE)
    output = re.sub(r"^Back Off! I just backed up.*$\n?", "", output, flags=re.MULTILINE)
    
    return output


# -------------- Equilibration checks --------------

def _gmx_executable() -> str:
    for exe in ("gmx", "gmx_mpi", "gmx_d"):
        if shutil.which(exe):
            return exe
    return "gmx"


def _run_gmx(sandbox_dir: str, args: list[str], log_file: Path, stdin: str | None = None) -> int:
    """Runs a GROMACS command in sandbox_dir, appending its output to log_file. Returns the return code."""
    with open(log_file, "a") as log:
        result = subprocess.run([_gmx_executable(), *args], cwd=sandbox_dir, input=stdin, stdout=log, stderr=subprocess.STDOUT, text=True)
    return result.returncode


def read_xvg(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Time (ps) and first data column of a GROMACS .xvg file."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith(("#", "@")):
            rows.append([float(x) for x in line.split()[:2]])
    data = np.array(rows).reshape(-1, 2)
    return data[:, 0], data[:, 1]


def _second_half(path: str | Path, quantity: str) -> dict:
    """Values of the second half of the run (the first half is still relaxing), or an error if the file is missing or too short."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {"produced": False, "passed": False, "detail": f"{path.name} was not produced: the {quantity} cannot be checked."}
    time, values = read_xvg(path)
    if len(time) < 4:
        return {"produced": True, "passed": False, "detail": f"{path.name} has only {len(time)} points: the {quantity} cannot be checked."}
    mask = time >= time[0] + (time[-1] - time[0]) / 2
    return {"produced": True, "time": time[mask], "values": values[mask], "window_ps": [float(time[mask][0]), float(time[-1])]}


def check_nvt_temperature(xvg: str | Path, ref_t: float) -> dict:
    """NVT passes if the mean temperature over the second half is within NVT_TEMP_TOLERANCE_K of the reference."""
    half = _second_half(xvg, "temperature")
    if "values" not in half:
        return half
    mean, std = float(np.mean(half["values"])), float(np.std(half["values"]))
    passed = abs(mean - ref_t) <= constants.NVT_TEMP_TOLERANCE_K
    return {"produced": True, "passed": passed, "mean": round(mean, 2), "std": round(std, 2), "window_ps": half["window_ps"],
            "detail": f"T = {mean:.1f} ± {std:.1f} K over {half['window_ps'][0]:.0f}-{half['window_ps'][1]:.0f} ps, reference {ref_t:.1f} K "
                      f"(tolerance ±{constants.NVT_TEMP_TOLERANCE_K} K): {'passed' if passed else 'not reached'}"}


def check_npt_pressure(xvg: str | Path, ref_p: float = constants.REF_PRESSURE_BAR) -> dict:
    """NPT pressure passes if the reference pressure is within one standard deviation of the mean over the second half."""
    half = _second_half(xvg, "pressure")
    if "values" not in half:
        return half
    mean, std = float(np.mean(half["values"])), float(np.std(half["values"]))
    passed = abs(mean - ref_p) <= std
    return {"produced": True, "passed": passed, "mean": round(mean, 2), "std": round(std, 2), "window_ps": half["window_ps"],
            "detail": f"P = {mean:.1f} ± {std:.1f} bar over {half['window_ps'][0]:.0f}-{half['window_ps'][1]:.0f} ps, reference {ref_p:.1f} bar "
                      f"(reference within one standard deviation): {'passed' if passed else 'not reached'}"}


def check_npt_density(xvg: str | Path) -> dict:
    """
    NPT density passes if it has stabilized: mean over the second half within DENSITY_RANGE_KG_M3 (protein + ions + water,
    so no exact match with pure water is expected) and linear drift over that window below DENSITY_MAX_DRIFT_KG_M3.
    """
    half = _second_half(xvg, "density")
    if "values" not in half:
        return half
    mean, std = float(np.mean(half["values"])), float(np.std(half["values"]))
    slope = float(np.polyfit(half["time"], half["values"], 1)[0])
    drift = slope * (half["window_ps"][1] - half["window_ps"][0])
    low, high = constants.DENSITY_RANGE_KG_M3
    in_range = low <= mean <= high
    stable = abs(drift) <= constants.DENSITY_MAX_DRIFT_KG_M3
    passed = in_range and stable
    reasons = [] if passed else [r for r, ok in ((f"outside {low:.0f}-{high:.0f} kg/m3", in_range), ("still drifting", stable)) if not ok]
    return {"produced": True, "passed": passed, "mean": round(mean, 2), "std": round(std, 2), "drift": round(drift, 2), "window_ps": half["window_ps"],
            "detail": f"density = {mean:.1f} ± {std:.1f} kg/m3, drift {drift:+.1f} kg/m3 over {half['window_ps'][0]:.0f}-{half['window_ps'][1]:.0f} ps "
                      f"(range {low:.0f}-{high:.0f}, max drift ±{constants.DENSITY_MAX_DRIFT_KG_M3:.0f}): {'stabilized' if passed else ', '.join(reasons)}"}


def _extend_stage(sandbox_dir: str, stage: str, log_file: Path) -> bool:
    """Extends a finished stage (nvt or npt) by EQUIL_EXTENSION_PS, continuing from its checkpoint (appends to its .edr/.xtc/.log)."""
    logger.info(f"Extending {stage} by {constants.EQUIL_EXTENSION_PS} ps")
    if _run_gmx(sandbox_dir, ["convert-tpr", "-s", f"{stage}.tpr", "-extend", str(constants.EQUIL_EXTENSION_PS), "-o", f"{stage}.tpr"], log_file) != 0:
        return False
    return _run_gmx(sandbox_dir, ["mdrun", "-v", "-deffnm", stage, "-cpi", f"{stage}.cpt"], log_file) == 0


def _write_energy(sandbox_dir: str, edr: str, term: str, xvg: str, log_file: Path) -> None:
    _run_gmx(sandbox_dir, ["energy", "-f", edr, "-o", xvg], log_file, stdin=f"{term}\n0\n")


def gromacs_equil(sandbox_dir: str, input_gro: str, md_temp: str, ligand_name=None, ligand_files=None, test: bool = False, ionic_conc: float = 0.15) -> str:
    # sometimes llm passes ligands as empty strings or the string "None"
    if not ligand_name or str(ligand_name).strip().lower() == "none":
        ligand_name = None
    if not ligand_files:
        ligand_file = None
    else:
        ligand_file=ligand_files[0]

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        protein_tcoupling = f"Protein_{ligand_name}"
    else:
        protein_tcoupling = "Protein"

    # NVT and NPT length: 100 ps, or 10 ps in test mode (no extension)
    equil_nvt_ps = constants.EQUIL_TEST_TIME_PS if test else constants.EQUIL_TIME_NVT_PS
    equil_npt_ps = constants.EQUIL_TEST_TIME_PS if test else constants.EQUIL_TIME_NPT_PS
    equil_nvt_nsteps = int(equil_nvt_ps / 0.002)  # 2 fs per step
    equil_npt_nsteps = int(equil_npt_ps / 0.002)  # 2 fs per step

    # -------------- Create ions.mdp, em.mdp, nvt.mdp, npt.mdp files --------------

    ions_mdp_infile = open(f'{sandbox_dir}/ions.mdp', 'w')
    ions_mdp_infile.write(f'''; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = cutoff	; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
''')
    ions_mdp_infile.close()

    em_mdp_infile = open(f'{sandbox_dir}/em.mdp', 'w' )
    em_mdp_infile.write(f'''; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
''')
    em_mdp_infile.close()

    nvt_mdp_infile = open(f'{sandbox_dir}/nvt.mdp', 'w' )
    nvt_mdp_infile.write(f'''title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = {equil_nvt_nsteps}     ; 2 fs * {equil_nvt_nsteps} = {equil_nvt_ps} ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = no       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = {protein_tcoupling} Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = {float(md_temp)}   {float(md_temp)}                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = no
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = {float(md_temp)}       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
''')
    nvt_mdp_infile.close()

    npt_mdp_infile = open(f'{sandbox_dir}/npt.mdp', 'w' )
    npt_mdp_infile.write(f'''title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = {equil_npt_nsteps}     ; 2 fs * {equil_npt_nsteps} = {equil_npt_ps} ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = {protein_tcoupling} Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = {float(md_temp)}   {float(md_temp)}                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = C-rescale                     ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT 
''')
    npt_mdp_infile.close()

    # if there are backups copy them to the original file
    
    if os.path.exists(f"{sandbox_dir}/topol_before_equil.top"):
        shutil.copy(f"{sandbox_dir}/topol_before_equil.top", f"{sandbox_dir}/topol.top")
    if os.path.exists(f"{sandbox_dir}/index_before_equil.ndx"):
        shutil.copy(f"{sandbox_dir}/index_before_equil.ndx", f"{sandbox_dir}/index.ndx")

    # -------------- Run equil_Gromacs.sh script --------------

    script = constants.SCRIPTS_DIR / "equil_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_equil.log")

    cmd = [str(script), sandbox_dir, input_gro, str(ionic_conc), log_file_path]

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        obabel_cmd = f"obabel {ligand_file} -O {sandbox_dir}/{ligand_name}.gro"
        obabel_result = subprocess.run(obabel_cmd, cwd=sandbox_dir, capture_output=True, text=True, shell=True)
        if obabel_result.returncode != 0:
            error_text = "\n".join(filter(None, [obabel_result.stderr, obabel_result.stdout]))
            return f"Equilibration failed with error: {error_text}" 
        ligand_gro = f"{ligand_name}.gro"
        cmd.append(ligand_name)
        cmd.append(ligand_file)
        cmd.append(ligand_gro)
        print(cmd)  

    # equil_Gromacs.sh: solvation, ions (counter-ions + ionic_conc M NaCl), energy minimization and NVT
    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    def read_log() -> str:
        if not log_file_path.exists():
            return ""
        try:
            return _filter_mdrun_output(log_file_path.read_text(encoding="utf-8"))
        except Exception as e:
            return f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{read_log()}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")

    report = {"test": test, "ionic_conc_m": float(ionic_conc), "equil_nvt_ps": equil_nvt_ps, "equil_npt_ps": equil_npt_ps, "md_temp": float(md_temp)}

    # -------------- NVT check: temperature reached, otherwise extend NVT --------------
    nvt = check_nvt_temperature(f"{sandbox_dir}/temperature.xvg", float(md_temp))
    nvt_extensions = 0
    while nvt["produced"] and not nvt["passed"] and not test and nvt_extensions < constants.MAX_NVT_EXTENSIONS:
        if not _extend_stage(sandbox_dir, "nvt", log_file_path):
            return f"Equilibration failed: the NVT extension did not run.\n--- Full GROMACS Log ---\n{read_log()}"
        nvt_extensions += 1
        _write_energy(sandbox_dir, "nvt.edr", "Temperature", "temperature.xvg", log_file_path)
        nvt = check_nvt_temperature(f"{sandbox_dir}/temperature.xvg", float(md_temp))
    if not nvt["produced"]:
        return f"Equilibration failed: {nvt['detail']}\n--- Full GROMACS Log ---\n{read_log()}"
    report["nvt"] = {**nvt, "extensions": nvt_extensions, "length_ps": equil_nvt_ps + nvt_extensions * constants.EQUIL_EXTENSION_PS}
    logger.info(f"NVT check: {nvt['detail']} ({nvt_extensions} extension(s))")

    # -------------- NPT --------------
    npt_cmds = [
        ["grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-t", "nvt.cpt", "-r", "nvt.gro", "-p", "topol.top", "-o", "npt.tpr", "-n", "index.ndx", "-maxwarn", "2"],
        ["mdrun", "-v", "-deffnm", "npt"],
    ]
    for args in npt_cmds:
        if _run_gmx(sandbox_dir, args, log_file_path) != 0 or (args[0] == "mdrun" and not Path(f"{sandbox_dir}/npt.gro").exists()):
            return f"Equilibration failed: NPT step 'gmx {args[0]}' did not complete.\n--- Full GROMACS Log ---\n{read_log()}"

    # -------------- NPT checks: pressure around 1 bar and stabilized density, otherwise extend NPT --------------
    def check_npt():
        _write_energy(sandbox_dir, "npt.edr", "Pressure", "pressure.xvg", log_file_path)
        _write_energy(sandbox_dir, "npt.edr", "Density", "density.xvg", log_file_path)
        return check_npt_pressure(f"{sandbox_dir}/pressure.xvg"), check_npt_density(f"{sandbox_dir}/density.xvg")

    pressure, density = check_npt()
    npt_extensions = 0
    while (pressure["produced"] and density["produced"] and not (pressure["passed"] and density["passed"])
           and not test and npt_extensions < constants.MAX_NPT_EXTENSIONS):
        if not _extend_stage(sandbox_dir, "npt", log_file_path):
            return f"Equilibration failed: the NPT extension did not run.\n--- Full GROMACS Log ---\n{read_log()}"
        npt_extensions += 1
        pressure, density = check_npt()
    if not (pressure["produced"] and density["produced"]):
        return f"Equilibration failed: {pressure['detail'] if not pressure['produced'] else density['detail']}\n--- Full GROMACS Log ---\n{read_log()}"
    report["npt"] = {"pressure": pressure, "density": density, "extensions": npt_extensions,
                     "length_ps": equil_npt_ps + npt_extensions * constants.EQUIL_EXTENSION_PS}
    logger.info(f"NPT checks: {pressure['detail']}; {density['detail']} ({npt_extensions} extension(s))")

    passed = nvt["passed"] and pressure["passed"] and density["passed"]
    report["passed"] = passed
    Path(f"{sandbox_dir}/equilibration_checks.json").write_text(json.dumps(report, indent=2))

    summary = (
        f"Ions: " + (f" {float(ionic_conc)} M NaCl added to have a neutral system.\n" if float(ionic_conc) > 0 else ", Na or Cl ions added to have a neutral system.\n") +
        f"NVT ({report['nvt']['length_ps']} ps, {nvt_extensions} extension(s)): {nvt['detail']}\n"
        f"NPT ({report['npt']['length_ps']} ps, {npt_extensions} extension(s)): {pressure['detail']}; {density['detail']}\n"
        f"Details saved to {sandbox_dir}/equilibration_checks.json."
    )
    if test:
        header = (f"Equilibration ran in TEST mode: NVT and NPT kept at {equil_nvt_ps} ps without extension. "
                  "The checks below are reported for information only and do not demonstrate equilibration.")
    elif passed:
        header = "Equilibration ran successfully and passed the temperature, pressure and density checks."
    else:
        header = ("WARNING: equilibration ran but did NOT pass all checks after the maximum number of extensions. "
                  "Report this to the user before relying on the production run; a longer equilibration may be needed.")
    return f"{header}\n{summary}\nFull GROMACS output:\n{read_log()}"


def gromacs_production(sandbox_dir: str, input_gro: str, npt_cpt_file: str, md_temp: str, md_duration: str, ligand_name=None) -> str:
    """
    Run production MD with GROMACS using prod_Gromacs.sh.
    """
    if not ligand_name or str(ligand_name).strip().lower() == "none":
        ligand_name = None

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        protein_tcoupling = f"Protein_{ligand_name}"
    else:
        protein_tcoupling = "Protein"

    # ---------- Create md.mdp file --------------
    nsteps = int(((float(md_duration)) * 1000000) / 2)  # Convert ns to number of steps (2 fs per step)
    duration_ps = int(nsteps)*2/1000
    duration_ns = int(duration_ps)/1000
    md_mdp_infile = open(f'{sandbox_dir}/md.mdp', 'w' )
    md_mdp_infile.write(f'''title                   = Protein-ligand complex MD simulation 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = {nsteps}   ; 2 * {nsteps} = {duration_ps} ps ({duration_ns} ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000     ; save energies every 10.0 ps
nstlog                  = 5000     ; update log file every 10.0 ps
nstxout-compressed      = 5000     ; save coordinates every 10.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = {protein_tcoupling} Water_and_ions        ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = {float(md_temp)}   {float(md_temp)}                     ; reference temperature, one for each group, in K
; Pressure coupling 
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no 
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration 
''')
    md_mdp_infile.close()

    # ---------- Run prod_Gromacs.sh script --------------

    script = constants.SCRIPTS_DIR / "prod_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_production.log")

    cmd = [str(script), input_gro, npt_cpt_file, log_file_path]
    
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:
        cmd.append(ligand_name)
        cmd.append(f"{sandbox_dir}/{ligand_name}.gro")

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = _filter_mdrun_output(log_file_path.read_text(encoding="utf-8"))
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")
    else:
        return (f"Production ran successfully. Full GROMACS output:\n"
                f"{gromacs_output}")


def gromacs_analysis(sandbox_dir: str,  pdb_id: str, input_xtc: str, ligand_name=None) -> str:
    """
    Run production MD with GROMACS using prod_Gromacs.sh.
    """
    if not ligand_name or str(ligand_name).strip().lower() == "none":
        ligand_name = None

    script = constants.SCRIPTS_DIR / "analysis_Gromacs.sh"
    log_file_path = Path(f"{sandbox_dir}/gromacs_analysis.log")

    cmd = [str(script), input_xtc, log_file_path]

    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h", ""]:
        cmd.append(ligand_name)
        cmd.append(f"{sandbox_dir}/{ligand_name}.gro")

    result = subprocess.run(cmd, cwd=sandbox_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    gromacs_output = ""

    if log_file_path.exists():
        try:
            gromacs_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            gromacs_output = f"Could not read GROMACS log file: {e}"

    # Check if ligand was not edited during the setup of simulations
    if ligand_name is not None and ligand_name not in ["XXX", "None", "None_h"]:

        logger.info("Checking if the ligand was modified during setup by comparing the original PDB of the ligand (from the input PDB) and final coordinates of the ligand (from md.gro). This is done to ensure that heavy atoms were not removed or added, which would modify the ligand.")

        ligands = defaultdict(list)
        ligands_gro = defaultdict(list)
        ligands_pdb_check = Counter()
        ligands_gro_check = Counter()

        with open(f"{sandbox_dir}/{pdb_id}.pdb", "r") as infile:
            for line in infile:
                if line.startswith("HETATM") and line[17:20].strip() == ligand_name:
                    chain = line[21].strip() or "_"
                    resnum = int(line[22:26])
                    ligands[(chain, resnum)].append(line)    

        ligand_pdb_file = f"{sandbox_dir}/{ligand_name}_check.pdb"
        with open(ligand_pdb_file, "w") as outfile:
            first_key = next(iter(ligands))
            outfile.writelines(ligands[first_key])
        logger.info(f"Extracted first occurence of ligand {ligand_name} to {ligand_pdb_file} to check if the ligand was modified during setup.")

        with open(f"{sandbox_dir}/md.gro", "r") as infile:
            ligand_gro_file = f"{sandbox_dir}/{ligand_name}_check.gro"
            for line in infile:
                if line[5:8] == ligand_name and "H" not in line[11:15]:
                    resnum = int(line[:5])
                    ligands_gro[(resnum)].append(line)
        with open(ligand_gro_file, "w") as outfile:
            first_key = next(iter(ligands_gro))
            atom_lines = ligands_gro[first_key]
            outfile.write("Ligand GRO file\n")
            outfile.write(f"{len(atom_lines)}\n")
            outfile.writelines(atom_lines)
            outfile.write("0.0 0.0 0.0\n")
            outfile.write("Ligand GRO file\n")
        logger.info(f"Extracted ligand {ligand_name} from md.gro to {ligand_gro_file} to check if the ligand was modified during setup.")
        subprocess.run(f"obabel {ligand_gro_file} -O {sandbox_dir}/{ligand_name}_gro_check.pdb", shell=True, check=False)

        #check if the 2 PDB are the same
        with open(f"{sandbox_dir}/{ligand_name}_check.pdb", "r") as pdb_file, open(f"{sandbox_dir}/{ligand_name}_gro_check.pdb", "r") as gro_pdb_file:
            pdb_lines = [line for line in pdb_file if line.startswith("HETATM") or line.startswith("ATOM")]
            gro_pdb_lines = [line for line in gro_pdb_file if line.startswith("HETATM") or line.startswith("ATOM")]

        if len(pdb_lines) != len(gro_pdb_lines):
            logger.info(f"Number of atoms differ between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro): {len(pdb_lines)} vs {len(gro_pdb_lines)}")
        else:
            logger.info(f"Number of atoms match between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro): {len(pdb_lines)}")

        for line in pdb_lines:
            if line[17:20].strip() == ligand_name:
                atom = line[75:80].strip()
                ligands_pdb_check[atom] += 1

        for line in gro_pdb_lines:
            if line[17:20].strip() == ligand_name:
                atom = line[75:80].strip()
                ligands_gro_check[atom] += 1

        # Compare the counts of each atom type
        if ligands_pdb_check == ligands_gro_check:
            logger.info(f"Atom counts (except hydrogens) match between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro).")
        else:
            logger.info(f"Atom counts (except hydrogens) differ between original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro).")
            logger.info(f"The atoms that the original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro) have in common are: {ligands_pdb_check & ligands_gro_check}")
            logger.info(f"The difference between the original PDB of the ligand (from {pdb_id}.pdb) and final coordinates of the ligand (from md.gro) is: {ligands_pdb_check - ligands_gro_check}")

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full GROMACS Log ---\n"
                f"{gromacs_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")
    else:
        return (f"Analysis plots produced successfully. Full GROMACS output:\n"
                f"{gromacs_output}")

