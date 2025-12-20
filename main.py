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

    temp: float | None = None
    "Simulation temperature in Kelvin."

    duration: float | None = None
    "Simulation length in nanoseconds."

    model_supports_system_messages: bool = True


def main(config: CommandLineArgs):
    root_logger = utils.get_class_logger("Main")

    # create a run directory inside of sandbox
    run_name = f"run_{utils.time_now()}"
    sandbox_dir = constants.DATA_DIR / run_name
    sandbox_dir.mkdir(parents=True, exist_ok=True)

    root_logger.info("DynaMate - your assistant for running molecular dynamics")

    try:
        _ensure_api_key("OPENROUTER_API_KEY", "OPENROUTER_API_KEY")
    except ValueError as e:
        root_logger.error(str(e))
        return

    root_logger.info("\n=== Starting PrepAgent (Planning & Parameter Determination) ===")
    prep_agent = PrepAgent(
        model_name=config.model,
        temperature=constants.TEMPERATURE,
        sandbox_dir=sandbox_dir,
        pdb_id=config.pdb_id,
        ligand_name=config.ligand,
        md_temp=config.temp,
        md_duration=config.duration,
        model_supports_system_messages=config.model_supports_system_messages,
    )
    prep_agent.setup_tools()
    pdb_file_path, ligand_name, plan, llm_cost = prep_agent.run()
    root_logger.info("PrepAgent completed. Plan generated.")

    md_temp, md_duration = plan["parameters"]["temperature_k"], plan["parameters"]["duration_ns"]

    root_logger.info(f"Applying parameters: Temp={md_temp}K, Duration={md_duration}ns")

    root_logger.info("\n=== Starting MDAgent (Execution & Tool Loop) ===")

    md_agent = MDAgent(
        model_name=config.model,
        temperature=constants.TEMPERATURE,
        sandbox_dir=sandbox_dir,
        structure_path=pdb_file_path,
        pdb_id=Path(pdb_file_path).stem,
        ligand_name=ligand_name,
        md_temp=md_temp,
        md_duration=md_duration,
        model_supports_system_messages=config.model_supports_system_messages,
        plan=plan,
    )
    md_agent.setup_tools()

    # Run the MD pipeline (handles user input and retries internally)
    result = md_agent.run()

    if not result:
        root_logger.error("=== MD Pipeline failed or incomplete ===")
    else:
        root_logger.info("=== MD Pipeline completed successfully ===")


if __name__ == "__main__":
    config = tyro.cli(CommandLineArgs)
    main(config)
