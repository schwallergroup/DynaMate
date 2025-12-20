import os
from pathlib import Path
from typing import Dict, Any, List
import json
from litellm import completion
from abc import ABC, abstractmethod
import traceback
from src import constants

from src.tools.map import TOOL_MAP
from src import utils, constants
import tiktoken

ENC = tiktoken.get_encoding("cl100k_base")


class ToolOutputError(Exception):
    """Custom exception raised when a tool returns a known failure string."""

    pass


class BaseAgent(ABC):
    def __init__(
        self,
        model_name: str,
        temperature: float,
        sandbox_dir: str,
        pdb_id: str | None = None,
        ligand_name: str | None = None,
        md_temp: float | None = None,
        md_duration: float | None = None,
        model_supports_system_messages: bool = True,
    ):
        self.model_name = model_name
        self.temperature = temperature
        self.sandbox_dir = Path(sandbox_dir)
        self.pdb_id = pdb_id
        self.ligand_name = ligand_name
        self.md_temp = md_temp
        self.md_duration = md_duration

        self.model_supports_system_messages = model_supports_system_messages

        self.tool_schemas = None
        self.messages: List[Dict[str, Any]] = []
        self.llm_cost = 0

        self.logger = utils.get_class_logger(self.__class__.__name__)

    @abstractmethod
    def setup_tools(self):
        pass

    @abstractmethod
    def _additional_check_for_errors_tool_output(self, tool_name, tool_call) -> bool:
        pass

    import json

    def _count_tokens(self, text):
        if not isinstance(text, str):
            text = json.dumps(text)
        return len(ENC.encode(text))
    
    def _find_recent_block(self, messages):
        """
        Return the smallest suffix of messages that forms a logically valid block.
        """
        if not messages:
            return []

        block = []
        i = len(messages) - 1

        while i >= 0:
            m = messages[i]
            block.insert(0, m)

            role = m.get("role")

            # Stop conditions for different last-message types:
            if role == "tool":
                i -= 1
                while i >= 0 and messages[i].get("role") != "assistant":
                    block.insert(0, messages[i])
                    i -= 1
                if i >= 0:
                    block.insert(0, messages[i])
                return block
            if role == "assistant" and m.get("tool_calls"):
                return block
            if role == "assistant":
                if i - 1 >= 0 and messages[i-1].get("role") == "user":
                    block.insert(0, messages[i-1])
                return block

            if role == "user":
                if i - 1 >= 0 and messages[i-1].get("role") == "assistant":
                    block.insert(0, messages[i-1])
                return block

            i -= 1

        return block


    def _summarize_history(self, messages):
        recent_messages = self._find_recent_block(messages)
        history_to_summarize = messages[:-len(recent_messages)]
        
        # Turn history into text
        history_text = "\n".join(
            f"{m['role']}: {m['content']}"
            for m in history_to_summarize
        )

        # Summarize
        summary_response = completion(
            model=self.model_name,
            temperature=0.1,
            messages=[
                {"role": "system", "content": "Summarize the conversation concisely but fully."},
                {"role": "user", "content": history_text}
            ],
            max_tokens=constants.SUMMARY_OUTPUT_TOKENS,
        )

        summary_text = summary_response.choices[0].message["content"]

        # NEW conversation = summary + recent logical block
        return [
            {"role": "assistant", "content": f"[Conversation Summary]\n{summary_text}"},
            *recent_messages
        ]


    def _validate_tool_path(self, tool_input) -> None:
        # if "path" in tool_input and not utils.is_path_child_dir(tool_input["path"], self.sandbox_dir):
        #     raise PermissionError(f"Access outside sandbox not allowed: {tool_input['path']}, {self.sandbox_dir}")

        # if "sandbox_dir" in tool_input and not utils.is_path_child_dir(tool_input["sandbox_dir"], self.sandbox_dir):
            # raise PermissionError(f"Access outside sandbox not allowed: {tool_input['sandbox_dir']}, {self.sandbox_dir}")
        pass

    def _safe_execute_tool(self, tool_name: str, tool_input: Dict[str, Any]) -> dict:
        """Executes a tool and catches any errors to pass back to the LLM."""

        self.logger.info(f"Executing tool: {tool_name} with input: {tool_input}")

        if not isinstance(tool_input, dict):
            try:
                tool_input = json.loads(tool_input)
            except Exception:
                error_msg = f"Invalid tool input format (not a dict): {tool_input}"
                self.logger.error(error_msg)
                return {"ok": False, "output": error_msg}

        tool_output = None

        self._validate_tool_path(tool_input)
        func = TOOL_MAP.get(tool_name)

        if not func:
            raise ValueError(f"Unknown tool: {tool_name}")

        tool_output = func(self, tool_input)
        passed = self._additional_check_for_errors_tool_output(tool_name, tool_output)

        return {"ok": passed, "output": tool_output}


    def _format_tool_usage_ouput(self, id_, tool_name, arguments, output):
        return {
            "tool_call_id": id_,
            "role": "tool",
            "name": tool_name,
            "arguments": arguments,
            "content": output,
        }

    def _create_logs(self):
        def to_dict_safe(msg):
            if isinstance(msg, dict):
                return msg
            elif hasattr(msg, "model_dump"):
                return msg.model_dump()
            else:
                return str(msg)

        safe_messages = [to_dict_safe(m) for m in self.messages]
        created_files = [f for f in os.listdir(self.sandbox_dir) if os.path.isfile(os.path.join(self.sandbox_dir, f))]

        # Build run log
        run_data = {
            "timestamp": utils.time_now(),
            "model": self.model_name,
            "temperature": self.temperature,
            "tools": self.tool_schemas,
            "messages": safe_messages,
            "files_created": created_files,
        }

        utils.append_jsonl(run_data, constants.JSON_LOG_FILE)

    def _call_llm(self, messages):
        # periodically summarize
        total_tokens = 0
        for m in messages:
            total_tokens += self._count_tokens(m.get("content", ""))

        if total_tokens > constants.MAX_CONTEXT_TOKENS and len(messages) > 3:
            messages = self._summarize_history(messages)
            self.messages = messages

        response = completion(
            model=self.model_name,
            temperature=self.temperature,
            supports_system_message=self.model_supports_system_messages,
            messages=messages,
            web_search_options={"search_context_size": "medium"},
            tools=self.tool_schemas,
            tool_choice="auto",
        )
        self.llm_cost += response._hidden_params["response_cost"]

        message = response.choices[0].message

        return message

    def _prompt_llm(self, prompt):
        self.messages.append({"role": "user", "content": prompt})
        return self._call_llm(self.messages)

    def _final_log(self, llm_cost):
        final_logs = {
            "timestamp_final": utils.time_now(),
            "total_completion_cost": llm_cost,
        }

        utils.append_jsonl(final_logs, constants.JSON_LOG_FILE)

    def _process_tool_call(self, tool_call):
        function_name = tool_call.function.name
        raw_args = tool_call.function.arguments

        if not raw_args:
            self.logger.warning(f"Tool '{function_name}' returned empty arguments. Using empty dict.")
            function_args = {}
        else:
            try:
                function_args = json.loads(raw_args)
            except json.JSONDecodeError:
                self.logger.warning(f"Tool '{function_name}' returned invalid JSON: {raw_args}. Using empty dict.")
                function_args = {}

        exec_result = self._safe_execute_tool(function_name, function_args)
        tool_output = utils.truncate_string(exec_result["output"])
        exec_result["output"] = tool_output
        
        self.logger.info(f"Tool result: {tool_output}.")

        tool_message = self._format_tool_usage_ouput(tool_call.id, function_name, function_args, tool_output)
        self.messages.append(tool_message)

        return exec_result

    @abstractmethod
    def _reset_pipeline(self) -> None:
        pass

    @abstractmethod
    def run(self):
        pass