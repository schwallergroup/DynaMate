import os
import glob
from pathlib import Path
from dotenv import load_dotenv, set_key
import shutil
import tyro
from dataclasses import dataclass

from src.agents import MDAgent, PrepAgent
from src import utils
from src import constants

if constants.ENV_FILE.exists():
    load_dotenv(dotenv_path=constants.ENV_FILE)


def _ensure_api_key(env_var: str, prompt_name: str) -> str | None:
    """Ensures a required API key is set, prompting the user if necessary."""
    key = os.environ.get(env_var)

    if not key:
        print(f"--- Missing API Key: {prompt_name} ---")
        key = input(f"Please enter your {prompt_name} API key: ").strip()
        if key:
            os.environ[env_var] = key
            set_key(str(constants.ENV_FILE), env_var, key)

    return key


def update_mdp_file(src_path, dest_path, md_temp, md_duration):
    def extract_comment(part):
        comment = ""

        if ";" in part:
            _, comment = part.split(";", 1)
            comment = ";" + comment

        return comment

    with open(src_path, "r") as f:
        lines = f.readlines()

    dt = None

    for line in lines:
        if line.strip().startswith("dt"):
            dt = float(line.split("=")[1].split(";")[0].strip())
            break
    if dt is None:
        raise ValueError(f"dt not found in {src_path}")

    # Calculate nsteps
    # md_duration is in ns, dt is in ps, so convert: 1 ns = 1000 ps
    nsteps = int((md_duration * 1000) / dt)

    # Update the lines
    new_lines = []

    for line in lines:
        if line.strip().startswith("ref_t"):
            # replace the value with md_temp
            parts = line.split("=")
            comment = extract_comment(parts[1])
            new_line = f"{parts[0]}= {md_temp} {md_temp} {comment}\n"

        elif line.strip().startswith("nsteps"):
            parts = line.split("=")
            comment = extract_comment(parts[1])

            new_line = f"{parts[0]}= {nsteps} {comment}\n"
        else:
            new_line = line

        new_lines.append(new_line)

    with open(dest_path, "w") as f:
        f.writelines(new_lines)


def process_mdp_files(mdp_dir, sandbox_dir, md_temp, md_duration):
    if not mdp_dir.exists():
        raise FileNotFoundError(f"MDP files directory not found: {mdp_dir}")

    mdp_files = os.listdir(mdp_dir)

    for file in mdp_files:
        full_file_name = os.path.join(mdp_dir, file)
        dest_file_name = os.path.join(sandbox_dir, file)

        if not os.path.isfile(full_file_name):
            continue

        if file in ("md.mdp", "npt.mdp", "nvt.mdp"):
            update_mdp_file(full_file_name, dest_file_name, md_temp, md_duration)
        else:
            shutil.copy(full_file_name, sandbox_dir)


@dataclass
class CommandLineArgs:
    """
    Tyro automatically generates a command line interface from this class.
    """

    pdb_id: str
    "PDB ID."

    model: str
    "Model name to use for the MD pipeline."

    ligand: str | None = None
    "Ligand ID (optional; defaults to no ligand)."

    model_supports_system_messages: bool = True


def main(config: CommandLineArgs):
    root_logger = utils.get_class_logger("Main")

    # create a run directory inside of sandbox
    run_name = f"run_{utils.time_now()}"
    sandbox_dir = constants.DATA_DIR / run_name
    sandbox_dir.mkdir(parents=True, exist_ok=True)

    root_logger.info("autoMD - your assistant for running molecular dynamics")
    root_logger.info("================")
    root_logger.info(
        "I can read, fetch, and prepare PDB files to run MD simulations. Is there a particular system I can help you with today?"
    )

    try:
        _ensure_api_key("OPENROUTER_API_KEY", "OPENROUTER_API_KEY")
    except ValueError as e:
        root_logger.error(str(e))
        return

    root_logger.info("\n=== 1. Starting PrepAgent (Planning & Parameter Determination) ===")
    prep_agent = PrepAgent(
        model_name=config.model,
        temperature=constants.TEMPERATURE,
        sandbox_dir=sandbox_dir,
        pdb_id=config.pdb_id,
        ligand_name=config.ligand,
        model_supports_system_messages=config.model_supports_system_messages,
    )
    # Start the prep agent to understand request from user and generate a plan
    prep_agent.setup_tools()
    pdb_file_path, ligand_name, plan, llm_cost = prep_agent.run()
    root_logger.info("PrepAgent completed. Plan generated.")

    # Copy required mdp_files into sandbox for the agent and update ref_t and nsteps
    md_temp, md_duration = plan["parameters"]["temperature_k"], plan["parameters"]["duration_ns"]

    md_duration=0.01

    root_logger.info(f"Applying parameters: Temp={md_temp}K, Duration={md_duration}ns")
    process_mdp_files(constants.MDP_FILES, sandbox_dir, md_temp, md_duration)

    root_logger.info("\n=== 2. Starting MDAgent (Execution & Tool Loop) ===")

    # Copy PDB into run directory
    md_agent = MDAgent(
        model_name=config.model,
        temperature=constants.TEMPERATURE,
        sandbox_dir=sandbox_dir,
        structure_path=pdb_file_path,
        pdb_id=Path(pdb_file_path).stem,
        ligand_name=ligand_name,
        model_supports_system_messages=config.model_supports_system_messages,
        plan=plan,
    )
    md_agent.setup_tools()

    # Run the MD pipeline (handles user input and retries internally)
    result = md_agent.run()

    # Print result and summary
    if not result:
        root_logger.error("=== MD Pipeline failed or incomplete ===")
    else:
        root_logger.info("=== MD Pipeline completed successfully ===")



if __name__ == "__main__":
    config = tyro.cli(CommandLineArgs)
    main(config)
