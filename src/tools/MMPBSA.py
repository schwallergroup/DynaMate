import subprocess
from pathlib import Path
import subprocess, shlex
from tkinter import constants
import parmed as pmd # type: ignore
import sys
import re
import os

def run_gmxMMPBSA(sandbox_dir: str, pdb_id: str, nsteps:str, nstxout_compressed:str, temp=str) -> str:
    nframes=int(nsteps)/int(nstxout_compressed)
    os.makedirs(f"{sandbox_dir}/gmx_MMPBSA", exist_ok=True)
    os.chdir(f"{sandbox_dir}/gmx_MMPBSA")
    mmpbsa_infile = open('mmpbsa.in', 'w' )
    mmpbsa_infile.write(f'''&general
sys_name={pdb_id}
startframe=1
endframe={int(float(nframes))}
interval=5
temperature={int(float(temp))}
verbose=2
/
&pb
  ipb                  = 2                                              # Dielectric model for PB
  inp                  = 1                                              # Nonpolar solvation method
  sander_apbs          = 0                                              # Use sander.APBS?
  indi                 = 1.0                                            # Internal dielectric constant
  exdi                 = 80.0                                           # External dielectric constant
  emem                 = 4.0                                            # Membrane dielectric constant
  smoothopt            = 1                                              # Set up dielectric values for finite-difference grid edges that are located across the solute/solvent dielectric boundary
  istrng               = 0.0                                            # Ionic strength (M)
  radiopt              = 1                                              # Use optimized radii?
  prbrad               = 1.4                                            # Probe radius
  iprob                = 2.0                                            # Mobile ion probe radius (Angstroms) for ion accessible surface used to define the Stern layer
  sasopt               = 0                                              # Molecular surface in PB implict model
  arcres               = 0.25                                           # The resolution (Å) to compute solvent accessible arcs
  memopt               = 0                                              # Use PB optimization for membrane
  poretype             = 1                                              # Use exclusion region for channel proteins
  npbopt               = 0                                              # Use NonLinear PB solver?
  solvopt              = 1                                              # Select iterative solver
  accept               = 0.001                                          # Sets the iteration convergence criterion (relative to the initial residue)
  linit                = 1000                                           # Number of SCF iterations
  fillratio            = 4.0                                            # Ratio between the longest dimension of the rectangular finite-difference grid and that of the solute
  scale                = 2.0                                            # 1/scale = grid spacing for the finite difference solver (default = 1/2 Å)
  nbuffer              = 0.0                                            # Sets how far away (in grid units) the boundary of the finite difference grid is away from the solute surface
  nfocus               = 2                                              # Electrostatic focusing calculation
  fscale               = 8                                              # Set the ratio between the coarse and fine grid spacings in an electrostatic focussing calculation
  npbgrid              = 1                                              # Sets how often the finite-difference grid is regenerated
  bcopt                = 5                                              # Boundary condition option
  eneopt               = 2                                              # Compute electrostatic energy and forces
  frcopt               = 0                                              # Output for computing electrostatic forces
  scalec               = 0                                              # Option to compute reaction field energy and forces
  cutfd                = 5.0                                            # Cutoff for finite-difference interactions
  cutnb                = 0.0                                            # Cutoff for nonbonded interations
  nsnba                = 1                                              # Sets how often atom-based pairlist is generated
  decompopt            = 2                                              # Option to select different decomposition schemes when INP = 2
  use_rmin             = 1                                              # The option to set up van der Waals radii
  sprob                = 0.557                                          # Solvent probe radius for SASA used to compute the dispersion term
  vprob                = 1.3                                            # Solvent probe radius for molecular volume (the volume enclosed by SASA)
  rhow_effect          = 1.129                                          # Effective water density used in the non-polar dispersion term calculation
  use_sav              = 1                                              # Use molecular volume (the volume enclosed by SASA) for cavity term calculation
  cavity_surften       = 0.0378                                         # Surface tension
  cavity_offset        = -0.5692                                        # Offset for nonpolar solvation calc
  maxsph               = 400                                            # Approximate number of dots to represent the maximum atomic solvent accessible surface
  maxarcdot            = 1500                                           # Number of dots used to store arc dots per atom
  npbverb              = 0                                             # Option to turn on verbose mode
/
''')
    mmpbsa_infile.close()

    #run_gmxMMPBSA("6JJ3","10000000","5000","300")

    tpr_file=f"{sandbox_dir}/md.tpr"
    xtc_file=f"{sandbox_dir}/md_noPBC.xtc"
    index_file=f"{sandbox_dir}/index.ndx"
    topol_file=f"{sandbox_dir}/topol.top"

    GMXMMPBSA_PATH = "/home/hackathon/miniforge3/envs/gmxMMPBSA/bin/gmx_MMPBSA"

    cmd = [
        GMXMMPBSA_PATH,
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

    subprocess.run(cmd, check=True)
    return f"MMPBSA complete! Files created: {sandbox_dir}/gmx_MMPBSA/FINAL_RESULTS_MMPBSA.dat and {sandbox_dir}/gmx_MMPBSA/FINAL_RESULTS_MMPBSA.csv"
