from pathlib import Path

MAX_CHARACTERS_TO_LOG = 5000
SUMMARY_OUTPUT_TOKENS = 6000
MAX_CONTEXT_TOKENS = 32000
PAPER_DIR = Path(__file__).resolve().parent.parent / "my_papers"
MODEL_NAME = "openrouter/openai/gpt-4.1-2025-04-14"
TEMPERATURE = 0.1

SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
MDP_FILES = SCRIPTS_DIR / "mdp_files"

ENV_FILE = Path(__file__).resolve().parent.parent / ".env"
DATA_DIR = Path(__file__).resolve().parent.parent / "sandbox"

AGENT_LOGS = Path(__file__).resolve().parent.parent / "agent_logs"
JSON_LOG_FILE = AGENT_LOGS / "agent_runs.jsonl"