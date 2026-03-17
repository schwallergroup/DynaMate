import json
from pathlib import Path
from src import constants


def _format_message(msg: dict) -> str:
    """Format a single message dict into readable trace text."""
    role = msg.get("role", "unknown")
    content = msg.get("content") or ""

    if role == "system":
        return f"[SYSTEM]\n{content}\n"

    if role == "user":
        return f"[USER]\n{content}\n"

    if role == "assistant":
        lines = []
        if content:
            lines.append(f"[ASSISTANT]\n{content}")
        tool_calls = msg.get("tool_calls") or []
        for tc in tool_calls:
            if isinstance(tc, dict):
                fn = tc.get("function", {})
                name = fn.get("name", "unknown_tool")
                args = fn.get("arguments", "{}")
            else:
                try:
                    fn = tc.get("function", {}) if isinstance(tc, dict) else {}
                    name = fn.get("name", "unknown_tool")
                    args = fn.get("arguments", "{}")
                except Exception:
                    name, args = "unknown_tool", "{}"
            lines.append(f"[ASSISTANT TOOL CALL: {name}]\n{args}")
        return "\n".join(lines) + "\n"

    if role == "tool":
        name = msg.get("name", "tool")
        return f"[TOOL RESULT: {name}]\n{content}\n"

    return f"[{role.upper()}]\n{content}\n"


def _load_runs() -> list[dict]:
    """Load all full run records (those with 'messages') from the JSONL log."""
    log_file = constants.JSON_LOG_FILE
    if not log_file.exists():
        return []
    runs = []
    with open(log_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
                if "messages" in record:
                    runs.append(record)
            except json.JSONDecodeError:
                continue
    return runs


def _is_error_result(content: str) -> bool:
    low = content.lower()
    return any(kw in low for kw in constants.ERROR_KEYWORDS)


def _extract_error_snippet(content: str, max_len: int = 150) -> str:
    """Return the first line that looks like an error message."""
    for line in content.splitlines():
        if line.strip() and _is_error_result(line):
            return line.strip()[:max_len]
    return content.strip()[:max_len]


def extract_error_examples(run_index: int = -1) -> list[str]:
    """
    Extract error-recovery patterns from an agent run as concise example
    strings suitable for upskill's generate_skill(examples=...).

    Each string captures what tool failed, the error, and what the agent
    did to recover — the most informative signal for skill generation in
    an agentic MD pipeline.

    Args:
        run_index: Which run to read (0 = first, -1 = last).

    Returns:
        List of strings like:
          "run_tleap_ligand failed with 'Could not find atom type: f1'
           → agent called param_ligand to regenerate force field parameters"
    """
    runs = _load_runs()
    if not runs:
        return []

    messages = runs[run_index].get("messages", [])
    examples: list[str] = []

    i = 0
    while i < len(messages):
        msg = messages[i]

        if msg.get("role") == "tool":
            content = msg.get("content", "")
            if _is_error_result(content):
                failed_tool = msg.get("name", "unknown_tool")
                error_snippet = _extract_error_snippet(content)

                # Look ahead up to 3 messages for the recovery action.
                recovery_tools: list[str] = []
                reasoning: str = ""
                for j in range(i + 1, min(i + 4, len(messages))):
                    next_msg = messages[j]
                    if next_msg.get("role") != "assistant":
                        continue
                    # Capture first sentence of reasoning text.
                    if not reasoning and next_msg.get("content"):
                        first_sentence = next_msg["content"].split(".")[0].strip()
                        if first_sentence and len(first_sentence) < 200:
                            reasoning = first_sentence
                    # Capture fix tool calls.
                    for tc in next_msg.get("tool_calls") or []:
                        fix_name = ""
                        if isinstance(tc, dict):
                            fix_name = tc.get("function", {}).get("name", "")
                        if fix_name and fix_name != failed_tool:
                            recovery_tools.append(fix_name)

                if recovery_tools or reasoning:
                    recovery = reasoning or f"called {', '.join(recovery_tools)}"
                    if recovery_tools:
                        recovery += f" (tools: {', '.join(recovery_tools)})"
                    example = (
                        f"{failed_tool} failed with '{error_snippet}' "
                        f"→ {recovery}"
                    )
                    examples.append(example)

        i += 1

    return examples


def export_run(run_index: int = -1, output_path: Path | None = None) -> Path:
    """
    Export a specific run from the JSONL log as a plain-text trace file.

    Args:
        run_index: Which run to export (0 = first, -1 = last). Selects among
                   records that contain a 'messages' key (skips cost-only records).
        output_path: Where to write the trace file. Defaults to
                     agent_logs/agent-trace.txt.

    Returns:
        Path to the written trace file.
    """
    runs = _load_runs()
    if not runs:
        raise ValueError("No completed run records found in agent log.")

    run = runs[run_index]
    messages = run.get("messages", [])
    model = run.get("model", "unknown")
    timestamp = run.get("timestamp", "unknown")
    files_created = run.get("files_created", [])

    lines = [
        f"# DynaMate Agent Trace",
        f"# Model: {model}",
        f"# Timestamp: {timestamp}",
        f"# Files created: {', '.join(files_created) if files_created else 'none'}",
        "",
    ]

    for msg in messages:
        lines.append(_format_message(msg))
        lines.append("---")

    trace_text = "\n".join(lines)

    if output_path is None:
        output_path = constants.AGENT_LOGS / "agent-trace.txt"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(trace_text, encoding="utf-8")

    return output_path


# Convenience alias used by main.py
def export_last_run(output_path: Path | None = None) -> Path:
    return export_run(run_index=-1, output_path=output_path)
