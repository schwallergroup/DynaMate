from pathlib import Path
from typing import Dict, Any
import json
import litellm

from src.agents.agent import BaseAgent, ToolOutputError
from src.prompts import MD_SYSTEM_PROMPT
from src.tools import tool_schema
from src.constants import MAX_ITERATIONS

litellm.drop_params = True 


class MDAgent(BaseAgent):
    ESSENTIAL_STEPS = {
        "prepare_pdb_file_ligand",
        "add_caps",
        "rename_histidines",
        "run_tleap_ligand",
        "run_gromacs_equil",
        "gromacs_production",
        "gromacs_analysis",
    }

    def __init__(
        self,
        model_name,
        temperature,
        sandbox_dir,
        structure_path,
        pdb_id=None,
        ligand_name=None,
        md_temp=None,
        md_duration=None,
        model_supports_system_messages=True,
        plan: Dict[str, Any] = None,
        post_pipeline_chat: bool = True,
    ):
        super().__init__(
            model_name,
            temperature,
            sandbox_dir,
            pdb_id,
            ligand_name,
            md_temp,
            md_duration,
            model_supports_system_messages
        )

        self.structure_path = Path(structure_path)
        self.plan = plan
        self.post_pipeline_chat = post_pipeline_chat

        self.completed_steps = []
        self.completed_summary = ""
        self.EXPECTED_FILES_PROTEIN_TEMPLATE = [
            "{pdb_id}.pdb",
            "{pdb_id}_prepared.pdb",
            "{pdb_id}_prepared_capped.pdb",
            "{pdb_id}_prepared_capped_his.pdb",
            "{pdb_id}.prmtop",
            "{pdb_id}.inpcrd",
            "{pdb_id}_tleap.pdb",
            "topol.top",
            "{pdb_id}.gro",
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
            "md.gro",
            "md.xtc",
            "rmsd.xvg",
            "rmsd_xtal.xvg",
            "rmsf.xvg",
            "gyrate.xvg",
            "hbnum_prot_wat.xvg",
            "hbnum_sidechain.xvg",
        ]
        self.EXPECTED_FILES_TEMPLATE = [
            "{pdb_id}.pdb",
            "{pdb_id}_prepared.pdb",
            ["{ligand_name}.pdb", "{ligand_name}_1.pdb"],
            ["{ligand_name}_h.pdb", "{ligand_name}_1_h.pdb"],
            "{pdb_id}_prepared_capped.pdb",
            "{pdb_id}_prepared_capped_his.pdb",
            ["{ligand_name}_h.prepi", "{ligand_name}_1_h.prepi"],
            ["{ligand_name}.frcmod", "{ligand_name}_h.frcmod", "{ligand_name}_1_h.frcmod"],
            "complex.pdb",
            "complex.prmtop",
            "complex.inpcrd",
            "complex_tleap.pdb",
            "topol.top",
            "complex.gro",
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
            "md.gro",
            "md.xtc",
            "rmsd.xvg",
            "rmsd_xtal.xvg",
            "rmsf.xvg",
            "gyrate.xvg",
            "hbnum_prot_lig.xvg",
            "hbnum_prot_wat.xvg",
            "hbnum_sidechain.xvg",
        ]

        self.logger.info(f"MDAgent initialized.")

    def _additional_check_for_errors_tool_output(self, tool_name, tool_call):
        if (tool_name in ("gromacs_production", "gromacs_equil", "gromacs_analysis")) and (
            " failed with return code " in tool_call
        ):
            self.logger.error(f"{tool_call}")
            return False

        if tool_name in ("run_tleap", "run_tleap_ligand", "param_ligand"):
            if "tleap run failed with error:" in tool_call:
                self.logger.error(f"{tool_call}")
                return False
            
            if "Ligand parameterization failed with error:" in tool_call:
                self.logger.error(f"{tool_call}")
                return False

            if "ParmEd failed:" in tool_call:
                self.logger.error(f"{tool_call}")
                return False

        return True

    def _reset_pipeline(self):
        from src.agents.agent import _TrackedMessageList
        self.completed_steps = []
        self.messages = _TrackedMessageList()
        self.completed_summary = ""

    def _validate_and_setup(self) -> bool:
        if not self.tool_schemas:
            raise NameError("Tool schema is not defined, call setup_tools() first.")

        self._reset_pipeline()

        if not self.structure_path.exists():
            self.logger.error(f"PDB file not found at {self.structure_path}")
            return False

        self.logger.info("=== Starting Molecular Dynamics Workflow ===")

        return True

    def _setup_system_prompt(self) -> None:
        plan_text = json.dumps(self.plan, indent=2) if self.plan else "No plan provided."
        system_prompt_text = MD_SYSTEM_PROMPT.format(
            pdb_path=self.structure_path,
            sandbox_dir=self.sandbox_dir,
        )

        system_context = (
            f"{system_prompt_text}\n\n"
            f"Here is the structured preparation plan generated by a specialized scientist:\n"
            f"{plan_text}\n\n"
            f"Use this plan to guide the molecular dynamics workflow. "
            f"Ensure the same parameters (temperature, duration) are respected.\n\n"
        )

        system_prompt = {
            "role": "system",
            "content": system_context,
        }

        self.messages.append(system_prompt)
        self.logger.info(f"Completed steps summary:\n{self.completed_summary}")

    def _process_tool_results(self, name, exec_result, remaining_steps):
        if exec_result["ok"]:
            self.completed_steps.append(name)
            self.completed_summary += f"{name} succeeded;\n"

            if name in remaining_steps:
                remaining_steps.remove(name)
        else:
            self.completed_summary += f"{name} failed;\n"

    def _process_tool_results_bfe(self, name, exec_result):
        if exec_result["ok"]:
            self.completed_steps.append(name)
            self.completed_summary += f"{name} succeeded;\n"
        else:
            self.completed_summary += f"{name} failed;\n"

    def _run_agent(self, remaining_steps: list[str]):
        # Inject the task prompt once before the loop
        remaining_steps_string = "\n".join(remaining_steps)
        self.messages.append({
            "role": "user",
            "content": (
                "Execute the molecular dynamics pipeline. The required steps are:\n"
                f"{remaining_steps_string}\n\n"
                "Work through each step in order using the available tools. "
                "Do not ask the user anything."
            ),
        })

        # Phase 1: run until the pipeline is successful
        iteration = 0
        while not self._pipeline_successful():
            iteration += 1

            if iteration == MAX_ITERATIONS:
                print(f"\n[Warning] Pipeline has reached {MAX_ITERATIONS} iterations without completing.")
                # answer = input("Would you like to quit and get a summary of what went wrong? (yes/no): ").strip().lower()
                answer = "yes"
                if answer in ("yes", "y"):
                    self.logger.info(f"User chose to quit after {iteration} iterations.")
                    return

            try:
                response = self._call_llm(self.messages)
                tool_calls = response.tool_calls

                if tool_calls:
                    self.logger.info(f"Iteration {iteration}: {len(tool_calls)} tool call(s)")
                    self.messages.append(response)
                    for tool_call in tool_calls:
                        exec_result = self._process_tool_call(tool_call)
                        self._process_tool_results(tool_call.function.name, exec_result, remaining_steps)
                else:
                    self.messages.append({"role": "assistant", "content": response.content})
                    if response.content:
                        print(f"\nAgent: {response.content}")
                    else:
                        print("\nAgent: I need your input to continue.")
                    # user_input = input("You: ").strip()
                    user_input = "No input provided, please make your own decision to help the pipeline continue."
                    self.messages.append({"role": "user", "content": user_input})

            except Exception as e:
                self.logger.error(str(e))
                self.messages.append({"role": "user", "content": f"An error occurred: {e}"})

        if not self.post_pipeline_chat:
            return

        # Phase 2: pipeline complete — open conversation until user presses Enter to exit
        self.logger.info("Pipeline successful. Entering post-pipeline conversation.")
        self.messages.append({
            "role": "user",
            "content": "The pipeline has completed successfully. Let the user know and ask if there is anything else you can help with (e.g. further analysis, binding free energy calculation).",
        })

        while True:
            try:
                response = self._call_llm(self.messages)
                tool_calls = response.tool_calls

                if tool_calls:
                    self.logger.info(f"Post-pipeline: {len(tool_calls)} tool call(s)")
                    self.messages.append(response)
                    for tool_call in tool_calls:
                        exec_result = self._process_tool_call(tool_call)
                        self._process_tool_results_bfe(tool_call.function.name, exec_result)
                else:
                    self.messages.append({"role": "assistant", "content": response.content})
                    if response.content:
                        print(f"\nAgent: {response.content}")
                    # user_input = input("You (press Enter to exit): ").strip()
                    user_input = None
                    if not user_input:
                        self.logger.info("User exited post-pipeline conversation.")
                        break
                    self.messages.append({"role": "user", "content": user_input})

            except Exception as e:
                self.logger.error(str(e))
                raise


    def _get_initial_steps(self) -> list[str]:
        if self.plan and "plan" in self.plan:
            return [step["step"] for step in self.plan["plan"]]

        return list(self.ESSENTIAL_STEPS)
    
    def _format_expected_files(self, templates):
        values = {
            "pdb_id": self.pdb_id.upper(),
            "ligand_name": self.ligand_name.upper() if self.ligand_name else "",
        }

        formatted = []

        for item in templates:
            if isinstance(item, (list, tuple)):
                formatted.append(
                    [s.format(**values) for s in item]
                )
            else:
                formatted.append(
                    item.format(**values)
                )

        return formatted

    def _resolve_file(self, f):
        if isinstance(f, (list, tuple)):
            for candidate in f:
                path = self.sandbox_dir / candidate
                if path.exists():
                    return candidate
            return f[0]
        return f

    def _pipeline_successful(self) -> bool:
        """Check whether the full MD pipeline completed successfully."""
        expected_files = (
            self._format_expected_files(self.EXPECTED_FILES_TEMPLATE)
            if self.ligand_name
            else self._format_expected_files(self.EXPECTED_FILES_PROTEIN_TEMPLATE)
        )

        missing_or_empty = []

        for f in expected_files:
            resolved = self._resolve_file(f)
            path = self.sandbox_dir / resolved

            if not path.exists() or path.stat().st_size == 0:
                missing_or_empty.append(resolved)

        if missing_or_empty:
            return False

        return True

    def _generate_and_log_summary(self, success: bool):
        if success:
            self.logger.info("=== MD Pipeline completed successfully ===")
            summary_prompt = "All tools succeeded. Summarize what was accomplished and key output files that can be used for human evaluation of the production run. Do not write this summary to any files, just respond with language. If errors occured during the pipeline, explain what the error was and how you fixed it."
        else:
            self.logger.error("=== MD Pipeline failed or incomplete ===")
            summary_prompt = "The MD pipeline failed or is incomplete. Summarize what went wrong and what should be investigated further. Do not write this summary to any files, just respond with language. If errors occured during the pipeline, explain what the error was and how you tried to fix it."

        self.messages.append(
            {
                "role": "user",
                "content": summary_prompt,
            }
        )

        try:
            summary_response = self._call_llm(self.messages)

            if summary_response.tool_calls:
                self.messages.append(summary_response)
                for tool_call in summary_response.tool_calls:
                    self._process_tool_call(tool_call)
            else:
                assistant_message = {"role": "assistant", "content": summary_response.content}
                self.messages.append(assistant_message)
                self.logger.info(f"Final summary for task:\n{summary_response.content}")

        except Exception as e:
            self.logger.error(str(e))
            raise

    def setup_tools(self):
        self.tool_schemas = tool_schema.create_tool_schema_md(self.sandbox_dir, self.ligand_name, self.pdb_id)

    def run(self):
        self._reset_pipeline()

        if not self._validate_and_setup():
            return False

        self._setup_system_prompt()

        remaining_steps = self._get_initial_steps()
        self._run_agent(remaining_steps)

        # _run_agent only returns after the user exits — pipeline is guaranteed successful here
        success = self._pipeline_successful()

        self._generate_and_log_summary(success)
        self._create_logs()
        self._final_log(self.llm_cost)

        return success
