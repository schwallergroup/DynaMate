import subprocess
from pathlib import Path
import subprocess, shlex
import parmed as pmd # type: ignore
import sys
import re
import os
from src import constants

def run_gmxMMPBSA(sandbox_dir: str, pdb_id: str, nsteps:str, nstxout_compressed:str, md_temp:str, ionic_conc: float = 0.15) -> str:
    nframes=int(nsteps)/int(nstxout_compressed)
    os.makedirs(f"{sandbox_dir}/gmx_MMPBSA", exist_ok=True)
    MMPBSA_dir=f"{sandbox_dir}/gmx_MMPBSA"
    start_frame = int(int(nframes)/4)
    mmpbsa_infile = open(f"{MMPBSA_dir}/mmpbsa.in", 'w' )
    mmpbsa_infile.write(f'''&general
sys_name={pdb_id}
startframe={start_frame}
endframe={int(float(nframes))}
interval=5
temperature={int(float(md_temp))}
verbose=2
/
&gb
igb=5
saltcon={ionic_conc}
/
&pb
istrng={ionic_conc}
fillratio=4.0
inp=1
/
''')
    mmpbsa_infile.close()

    tpr_file=f"{sandbox_dir}/md.tpr"
    xtc_file=f"{sandbox_dir}/md_noPBC.xtc"
    index_file=f"{sandbox_dir}/index.ndx"
    topol_file=f"{sandbox_dir}/topol.top"
    MMPBSA_ENV_DIR=os.environ.get("MMPBSA_ENV_DIR", str(constants.MMPBSA_ENV_DIR))

    cmd = [
        MMPBSA_ENV_DIR,
        "-O",
        "-i", "mmpbsa.in",
        "-cs", tpr_file,
        "-ct", xtc_file,
        "-ci", index_file,
        "-cg", "1", "13",
        "-cp", topol_file,
        "-o", "FINAL_RESULTS_MMPBSA.dat",
        "-eo", "FINAL_RESULTS_MMPBSA.csv",
        "-nogui"
    ]

    result = subprocess.run(cmd, cwd=MMPBSA_dir, stdout=sys.stdout, stderr=sys.stderr, text=True)

    MMPBSA_output = ""

    log_file_path = Path(MMPBSA_dir) / "gmx_MMPBSA.log"

    if log_file_path.exists():
        try:
            MMPBSA_output = log_file_path.read_text(encoding="utf-8")
        except Exception as e:
            MMPBSA_output = f"Could not read gmx_MMPBSA log file: {e}"

    if result.returncode != 0:
        return (f"Equilibration script failed with return code {result.returncode}.\n"
                f"--- Full gmx_MMPBSA Log ---\n"
                f"{MMPBSA_output}\n"
                f"--- Shell Script Stderr ---\n"
                f"{result.stderr or 'None captured directly'}")
    else:
        return (f"MMPBSA complete! Files created: {MMPBSA_dir}/FINAL_RESULTS_MMPBSA.dat and {MMPBSA_dir}/FINAL_RESULTS_MMPBSA.csv."
                f"Full gmx_MMPBSA output:\n"
                f"{MMPBSA_output}")