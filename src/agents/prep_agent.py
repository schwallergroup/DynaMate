import re
import json
import sys
from typing import Dict, Any, List
from pydantic import BaseModel
import litellm
from src.tools import tool_schema
from src.agents.agent import BaseAgent
from src.prompts import PREP_SYSTEM_PROMPT

litellm.drop_params = True


class Tool(BaseModel):
    name: str
    description: str
    parameters: Dict[str, Any]


class PrepAgent(BaseAgent):
    def __init__(
        self,
        model_name,
        temperature,
        sandbox_dir,
        pdb_id=None,
        ligand_name=None,
        md_temp=None,
        md_duration=None,
        model_supports_system_messages=True,
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

        self.pdb_file_path = None
        self.user_temp = None
        self.agent_plan = ""

        self.logger.info(f"PrepAgent initialized.")

    def setup_tools(self):
        self.tool_schemas = tool_schema.create_tool_schema_prep(self.sandbox_dir)

    def _additional_check_for_errors_tool_output(self, tool_name, tool_call) -> bool:
        # No additional domain-specific checks needed for PrepAgent's tools based on output string
        return True

    def _reset_pipeline(self):
        self.agent_plan = ""

    def _setup_system_prompt(self) -> None:
        system_prompt_text = PREP_SYSTEM_PROMPT.format(
            sandbox_dir=self.sandbox_dir,
        )

        system_prompt = {
            "role": "system",
            "content": system_prompt_text,
        }

        self.messages.append(system_prompt)

    def _ask_for_system(self):
        """Have the LLM ask the user for the PDB ID and optional ligand, parse via LLM, then confirm."""
        # Step 1: LLM asks the user
        self.messages.append({
            "role": "user",
            "content": "Ask the user what molecular system they would like to simulate (PDB ID or file upload) and whether they have a ligand to include (3-letter code).",
        })
        response = self._call_llm(self.messages)
        self.messages.append({"role": "assistant", "content": response.content})
        print(f"\nAgent: {response.content}")

        while True:
            user_input = input("You: ").strip()
            self.messages.append({"role": "user", "content": user_input})

            # Step 2: LLM extracts PDB ID and ligand as JSON
            parse_messages = [
                {
                    "role": "user",
                    "content": (
                        f"Extract the PDB ID (4-character alphanumeric code) and ligand ID "
                        f"(3-character alphanumeric code, if present) from this user response:\n\n"
                        f"\"{user_input}\"\n\n"
                        f"Reply with only valid JSON in this exact format: "
                        f'{{\"pdb_id\": \"XXXX\", \"ligand\": \"XXX\"}} or '
                        f'{{\"pdb_id\": \"XXXX\", \"ligand\": null}} if no ligand was mentioned.'
                    ),
                }
            ]
            parse_response = self._call_llm(parse_messages)

            try:
                raw = parse_response.content.strip()
                raw = re.sub(r'^```(?:json)?\s*|\s*```$', '', raw, flags=re.MULTILINE).strip()
                extracted = json.loads(raw)
                pdb_id = (extracted.get("pdb_id") or "").strip().upper() or None
                ligand = (extracted.get("ligand") or "").strip().upper() or None
            except (json.JSONDecodeError, AttributeError):
                pdb_id, ligand = None, None

            # Step 3: Confirm with user
            confirm_parts = [f"PDB ID: {pdb_id or 'not found'}"]
            confirm_parts.append(f"Ligand: {ligand}" if ligand else "Ligand: none")
            print(f"\nAgent: I understood the following — {', '.join(confirm_parts)}. Is that correct? (yes/no)")

            confirmation = input("You: ").strip().lower()
            if confirmation in ("yes", "y"):
                break

            print("Agent: No problem, please provide the PDB ID and ligand again.")

        self.pdb_id = pdb_id
        if self.ligand_name is None:
            self.ligand_name = ligand

    def _get_pdb_file_path(self, prompt):
        # Inject the task prompt once before the loop
        self.messages.append({"role": "user", "content": prompt})

        for _ in range(5):
            response = self._call_llm(self.messages)
            self.logger.info(f"Response: {response.content}")

            tool_calls = response.tool_calls
            if tool_calls:
                self.logger.info(f"Length of tool calls: {len(tool_calls)}")
                self.messages.append(response)
                for tool_call in tool_calls:
                    self._process_tool_call(tool_call)
            else:
                self.messages.append({"role": "assistant", "content": response.content})
                if response.content:
                    print(f"\nAgent: {response.content}")
                user_input = "No input provided, please make your own decision to help the pipeline continue."
                # user_input = input("You: ").strip()
                if user_input:
                    self.messages.append({"role": "user", "content": user_input})

            pdb_files = list(self.sandbox_dir.glob("*.pdb"))
            if pdb_files:
                return str(pdb_files[0])

        return None

    def _find_ligand(self):
        while True:
            user_input = self.ligand_name

            if not user_input:
                self.logger.info("Defining system without a ligand")
                return

            lig_name = re.search(r"^[A-Z0-9]{3}$", user_input)
            if not lig_name:
                self.logger.error(
                    f"'{user_input}' is not a valid 3-character ligand code. Continuing without a ligand."
                )
                self.ligand_name = None
                return

            self.ligand_name = lig_name.group()
            self.logger.info(f"User requested ligand: {self.ligand_name}")

            lig_found = False
            with open(self.pdb_file_path, "r") as infile:
                for line in infile:
                    if line.startswith("HETATM") and self.ligand_name in line:
                        lig_found = True
                        break

            if lig_found:
                return

            self.logger.error(
                f"Ligand '{self.ligand_name}' not found in {self.pdb_file_path}."
            )
            self.ligand_name = None

    def _find_simulation_temperature(self):
        temperature = self.md_temp

        while temperature is None:
            user_input = "Please pick a suitable temperature (in Kelvins) for running the molecular dynamics simulation. Also provide a rational for why you selected this temperature."
            response = self._prompt_llm(user_input)
            tool_calls = response.tool_calls

            if tool_calls:
                self.logger.info(f"Length of tool calls: {len(tool_calls)}")
                self.messages.append(response)

                for tool_call in tool_calls:
                    self._process_tool_call(tool_call)

            else:
                assistant_message = {"role": "assistant", "content": response.content}
                self.messages.append(assistant_message)
                match = re.search(r"(\d+\.?\d*)\s+K", response.content)
                temperature = float(match.group(1)) if match else 310.0
                self.logger.info(f"Rational for picking simulation temperature {temperature}: {response.content}")

        return temperature

    def _calculate_duration(self):
        duration = self.md_duration

        while duration is None:
            user_input = "Please pick a suitable simulation duration (in nanoseconds) for running a short molecular dynamics of this system. Also provide a rational for why you selected this duration. Note that dynamics simulations of 10 ns take 1 hour to complete, so keep the experiment brief (less than 1 ns duration)."
            response = self._prompt_llm(user_input)
            tool_calls = response.tool_calls

            if tool_calls:
                self.logger.info(f"Length of tool calls: {len(tool_calls)}")
                self.messages.append(response)
                for tool_call in tool_calls:
                    self._process_tool_call(tool_call)
            else:
                assistant_message = {"role": "assistant", "content": response.content}
                self.messages.append(assistant_message)
                matches = re.findall(r"(\d+\.?\d*)\s+ns", response.content)
                duration = float(matches[-1]) if matches else 0.01
                self.logger.info(f"Rational for picking simulation duration {duration}: {response.content}")

        return duration

    def _generate_plan(self, temperature, duration):
        if self.ligand_name:
            steps = [
                {
                    "step": "prepare_pdb_file_ligand",
                    "description": "Clean and preprocess PDB file for protein-ligand system.",
                },
                {"step": "add_caps", "description": "Add N- and C-terminal capping groups."},
                {"step": "rename_histidines", "description": "Rename HIS to HIE, HIP or HID."},
                {"step": "param_ligand", "description": "Generate ligand parameters using antechamber or acpype."},
                {"step": "run_tleap_ligand", "description": "Build system topology and solvate complex using tleap."},
                {"step": "gromacs_equil", "description": "Perform energy minimization and equilibration."},
                {"step": "gromacs_production", "description": "Run production MD simulation."},
                {"step": "gromacs_analysis", "description": "Analyse production MD simulation."},
            ]
        else:
            steps = [
                {
                    "step": "prepare_pdb_file_ligand",
                    "description": "Clean and preprocess PDB file for protein-only system.",
                },
                {"step": "add_caps", "description": "Add N- and C-terminal capping groups."},
                {"step": "rename_histidines", "description": "Rename HIS to HIE, HIP or HID."},
                {"step": "run_tleap", "description": "Build system topology and solvate complex using tleap."},
                {"step": "gromacs_equil", "description": "Perform energy minimization and equilibration."},
                {"step": "gromacs_production", "description": "Run production MD simulation."},
                {"step": "gromacs_analysis", "description": "Analyse production MD simulation."},
            ]

        plan = {
            "sandbox_dir": str(self.sandbox_dir),
            "pdb_file_path": self.pdb_file_path,
            "ligand": self.ligand_name,
            "plan": steps,
            "parameters": {"temperature_k": float(temperature), "duration_ns": float(duration)},
        }

        return plan

    def run(self):
        self._reset_pipeline()

        if not self.tool_schemas:
            raise NameError("Tool schema is not defined, call setup_tools() first.")

        self._setup_system_prompt()

        interactive = self.pdb_id is None
        if interactive:
            self._ask_for_system()

        if interactive:
            prompt = "Please proceed to fetch and prepare the PDB file for the system we just discussed."
        else:
            prompt = f"Fetch and prepare the PDB file for {self.pdb_id}."

        self.logger.info(f"Starting prep for PDB: {self.pdb_id}, Ligand: {self.ligand_name or 'none'}")
        self.pdb_file_path = self._get_pdb_file_path(prompt)
        self.logger.info(f"Structure ready: {self.pdb_file_path}")

        self._find_ligand()

        if self.md_temp is None:
            self.md_temp = self._find_simulation_temperature()
        self.logger.info(f"Using simulation temperature: {self.md_temp} K")

        if self.md_duration is None:
            self.md_duration = self._calculate_duration()
        self.logger.info(f"Using simulation duration: {self.md_duration} ns")

        # Build plan steps depending on ligand
        plan = self._generate_plan(self.md_temp, self.md_duration)
        self.agent_plan = json.dumps(plan, indent=2)
        self.logger.info(f"Generated plan: {self.agent_plan}")

        self._create_logs()
        return self.pdb_file_path, self.ligand_name, plan, self.llm_cost
